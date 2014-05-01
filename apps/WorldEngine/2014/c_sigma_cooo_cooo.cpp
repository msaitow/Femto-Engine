                                                                                
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
#include <sci/icmr/Femto/elems/c_sigma_cooo_cooo.h>                                  
                                                                                
using std::cout;                                                                
using std::endl;                                                                
                                                                                
#define FLOPCOUNT                                                               
                                                                                
using orz::mr::Femto::Fc0;                                                      
using orz::mr::Femto::h1_int;                                                   
using orz::mr::Femto::h6_int;                                                   
                                                                                
//        :::::::::: ::::::::::   :::   ::: ::::::::::: :::::::: 
//       :+:        :+:         :+:+: :+:+:    :+:    :+:    :+: 
//      +:+        +:+        +:+ +:+:+ +:+   +:+    +:+    +:+  
//     :#::+::#   +#++:++#   +#+  +:+  +#+   +#+    +#+    +:+   
//    +#+        +#+        +#+       +#+   +#+    +#+    +#+    
//   #+#        #+#        #+#       #+#   #+#    #+#    #+#     
//  ###        ########## ###       ###   ###     ########       

//                                   Generated date : Sun Apr 20 10:26:14 2014

                                                                                
// ***************************************************************************  
// orz::mr::mrci                                                                
// ***************************************************************************  
									    /*!        
   @brief CT input                                                              
									     */        
                                                                                
orz::mr::BareAmpPack orz::mr::Femto::sigma_cooo_cooo(const orz::mr::Input &ctinp,                                    
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
  // |-- [    0] --| W0(k,j,a2,a1) += (    1.00000000) D2(k,j,a0,a2) Fc1(a1,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -0.50000000) T2(w,a2,a1,i) W0(k,j,a2,a1) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W0aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no0_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO0_X0_TYPE0_NOERI)
      (sj, ij, W0aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    for(int si = 0;si < nir;++si){ 
    for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
      T2b = T2.get_amp2(ii);
      FC_FUNC(g_if_sigma_cooo_cooo_no0_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO0_X1_TYPE0_NOERI)
        (si, ii, sj, ij, T2b.cptr(), W0aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- No. [    1] -- 
  // |-- [    0] --| W1(w,k,j,a1) += (    1.00000000) T2(w,a2,a1,a0) D2(k,j,a0,a2) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) Fc1(i,a1) W1(w,k,j,a1) 
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W1caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
    for(int sa0 = 0;sa0 < nir;++sa0){ 
    for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
      T2b = T2.get_amp2(ia0);
      FC_FUNC(g_if_sigma_cooo_cooo_no1_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO1_X0_TYPE0_NOERI)
        (sa0, ia0, sj, ij, T2b.cptr(), W1caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
    } // End ia0
    } // End sa0
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no1_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO1_X1_TYPE0_NOERI)
      (sj, ij, W1caa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W2(k,a2,i,a1) += (    1.00000000) D2(k,a2,a0,i) Fc1(a1,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -0.50000000) T2(w,a2,a1,j) W2(k,a2,i,a1) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W2aaaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, 0));
  FC_FUNC(g_if_sigma_cooo_cooo_no2_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO2_X0_TYPE0_NOERI)
    (W2aaaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no2_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO2_X1_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W2aaaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W3(w,k,i,a1) += (    1.00000000) T2(w,a2,a1,a0) D2(k,i,a0,a2) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) Fc1(j,a1) W3(w,k,i,a1) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W3caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_cooo_cooo_no3_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO3_X0_TYPE0_NOERI)
      (sa0, ia0, T2b.cptr(), W3caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no3_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO3_X1_TYPE0_NOERI)
      (sj, ij, W3caaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W4(w,k,a0,i) += (    1.00000000) T2(w,a1,a0,i) D1(k,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) Fc1(j,a0) W4(w,k,a0,i) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    orz::DTensor W4caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, si));
    FC_FUNC(g_if_sigma_cooo_cooo_no4_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO4_X0_TYPE0_NOERI)
      (si, ii, T2b.cptr(), W4caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_cooo_no4_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO4_X1_TYPE0_NOERI)
        (si, ii, sj, ij, W4caa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- No. [    5] -- 
  // |-- [    0] --| W5(w,k,a0,j) += (    1.00000000) T2(w,a1,a0,j) D1(k,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) Fc1(i,a0) W5(w,k,a0,j) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W5caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no5_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO5_X0_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W5caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no5_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO5_X1_TYPE0_NOERI)
      (sj, ij, W5caa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W6(k,a1,j,a3,i,a2) += (    1.00000000) D3(k,j,a1,a3,a0,i) Fc1(a2,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -0.50000000) T2(w,a3,a2,a1) W6(k,a1,j,a3,i,a2) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      orz::DTensor W6aaaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sa1^sj));
      FC_FUNC(g_if_sigma_cooo_cooo_no6_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO6_X0_TYPE0_NOERI)
        (sa1, ia1, sj, ij, W6aaaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_cooo_no6_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO6_X1_TYPE0_NOERI)
        (sa1, ia1, sj, ij, T2b.cptr(), W6aaaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W7(k,j,a2,a1) += (    1.00000000) D2(k,j,a0,a2) Fc1(a1,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) T2(w,a2,i,a1) W7(k,j,a2,a1) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      orz::DTensor W7aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sa1));
      FC_FUNC(g_if_sigma_cooo_cooo_no7_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO7_X0_TYPE0_NOERI)
        (sa1, ia1, sj, ij, W7aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_cooo_no7_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO7_X1_TYPE0_NOERI)
        (sa1, ia1, sj, ij, T2b.cptr(), W7aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| W8(w,k,j,a0) += (    1.00000000) T2(w,a2,a1,a0) D2(k,j,a1,a2) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) Fc1(i,a0) W8(w,k,j,a0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      orz::DTensor W8ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sa0));
      FC_FUNC(g_if_sigma_cooo_cooo_no8_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO8_X0_TYPE0_NOERI)
        (sa0, ia0, sj, ij, T2b.cptr(), W8ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
      FC_FUNC(g_if_sigma_cooo_cooo_no8_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO8_X1_TYPE0_NOERI)
        (sa0, ia0, sj, ij, W8ca_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- No. [    9] -- 
  // |-- [    0] --| W9(k,i,a2,a1) += (    1.00000000) D2(k,i,a0,a2) Fc1(a1,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -0.50000000) T2(a2,w,a1,j) W9(k,i,a2,a1) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W9aaaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, 0));
  FC_FUNC(g_if_sigma_cooo_cooo_no9_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO9_X0_TYPE0_NOERI)
    (W9aaaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no9_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO9_X1_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W9aaaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| W10(w,k,i,a0) += (    1.00000000) T2(w,a2,a1,a0) D2(k,a2,a1,i) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) Fc1(j,a0) W10(w,k,i,a0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W10caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa0));
    FC_FUNC(g_if_sigma_cooo_cooo_no10_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO10_X0_TYPE0_NOERI)
      (sa0, ia0, T2b.cptr(), W10caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_cooo_no10_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO10_X1_TYPE0_NOERI)
        (sa0, ia0, sj, ij, W10caa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- No. [   11] -- 
  // |-- [    0] --| W11(w,k,i,a0) += (    1.00000000) T2(w,a1,i,a0) D1(k,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) Fc1(j,a0) W11(w,k,i,a0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W11caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa0));
    FC_FUNC(g_if_sigma_cooo_cooo_no11_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO11_X0_TYPE0_NOERI)
      (sa0, ia0, T2b.cptr(), W11caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_cooo_no11_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO11_X1_TYPE0_NOERI)
        (sa0, ia0, sj, ij, W11caa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- No. [   12] -- 
  // |-- [    0] --| W12(w,k,j,a0) += (    1.00000000) T2(a1,w,a0,j) D1(k,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) Fc1(i,a0) W12(w,k,j,a0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W12caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no12_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO12_X0_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W12caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no12_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO12_X1_TYPE0_NOERI)
      (sj, ij, W12caa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| W13(k,a2,j,i,a3,a1) += (    1.00000000) D3(k,j,a2,i,a0,a3) Fc1(a1,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -0.50000000) T2(w,a3,a2,a1) W13(k,a2,j,i,a3,a1) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      orz::DTensor W13aaaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sj^sa1));
      FC_FUNC(g_if_sigma_cooo_cooo_no13_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO13_X0_TYPE0_NOERI)
        (sa1, ia1, sj, ij, W13aaaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_cooo_no13_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO13_X1_TYPE0_NOERI)
        (sa1, ia1, sj, ij, T2b.cptr(), W13aaaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   14] -- 
  // |-- [    0] --| W14(c0,k,j,i) += (    1.00000000) T2(c0,a1,a0,i) D2(k,j,a0,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) Fc1(w,c0) W14(c0,k,j,i) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      orz::DTensor W14ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^si));
      FC_FUNC(g_if_sigma_cooo_cooo_no14_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO14_X0_TYPE0_NOERI)
        (si, ii, sj, ij, T2b.cptr(), W14ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_cooo_no14_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO14_X1_TYPE0_NOERI)
        (si, ii, sj, ij, W14ca_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- No. [   15] -- 
  // |-- [    0] --| W15(c0,k,j,i) += (    1.00000000) T2(c0,a1,i,a0) D2(k,j,a0,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (   -2.00000000) Fc1(w,c0) W15(c0,k,j,i) 
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W15caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
    for(int sa0 = 0;sa0 < nir;++sa0){ 
    for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
      T2b = T2.get_amp2(ia0);
      FC_FUNC(g_if_sigma_cooo_cooo_no15_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO15_X0_TYPE0_NOERI)
        (sa0, ia0, sj, ij, T2b.cptr(), W15caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
    } // End ia0
    } // End sa0
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no15_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO15_X1_TYPE0_NOERI)
      (sj, ij, W15caa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   16] -- 
  // |-- [    0] --| W16(c0,k,i,j) += (    1.00000000) T2(c0,a1,a0,j) D2(k,a1,a0,i) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) Fc1(w,c0) W16(c0,k,i,j) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W16caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no16_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO16_X0_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W16caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no16_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO16_X1_TYPE0_NOERI)
      (sj, ij, W16caa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   17] -- 
  // |-- [    0] --| W17(c0,k,i,j) += (    1.00000000) T2(a1,c0,a0,j) D2(k,i,a0,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) Fc1(w,c0) W17(c0,k,i,j) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W17caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no17_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO17_X0_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W17caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no17_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO17_X1_TYPE0_NOERI)
      (sj, ij, W17caa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   18] -- 
  // |-- [    0] --| W18(c0,k,j,i) += (    1.00000000) T2(a0,c0,i,j) D1(k,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) Fc1(w,c0) W18(c0,k,j,i) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W18caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no18_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO18_X0_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W18caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no18_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO18_X1_TYPE0_NOERI)
      (sj, ij, W18caa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   19] -- 
  // |-- [    0] --| W19(c0,k,i,j) += (    1.00000000) T2(c0,a0,i,j) D1(k,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -2.00000000) Fc1(w,c0) W19(c0,k,i,j) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W19caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no19_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO19_X0_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W19caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no19_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO19_X1_TYPE0_NOERI)
      (sj, ij, W19caa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   20] -- 
  // |-- [    0] --| W20(c0,k,j,i) += (    1.00000000) T2(c0,a2,a1,a0) D3(k,j,a1,i,a0,a2) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) Fc1(w,c0) W20(c0,k,j,i) 
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W20caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
    for(int sa0 = 0;sa0 < nir;++sa0){ 
    for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
      T2b = T2.get_amp2(ia0);
      FC_FUNC(g_if_sigma_cooo_cooo_no20_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO20_X0_TYPE0_NOERI)
        (sa0, ia0, sj, ij, T2b.cptr(), W20caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
    } // End ia0
    } // End sa0
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no20_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO20_X1_TYPE0_NOERI)
      (sj, ij, W20caa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   21] -- 
  // |-- [    0] --| W21(k,a1,j,a2) += (    1.00000000) D2(k,j,a1,a0) Fc1(a2,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) T2(w,a2,i,a1) W21(k,a1,j,a2) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      orz::DTensor W21aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa1^sj));
      FC_FUNC(g_if_sigma_cooo_cooo_no21_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO21_X0_TYPE0_NOERI)
        (sa1, ia1, sj, ij, W21aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_cooo_no21_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO21_X1_TYPE0_NOERI)
        (sa1, ia1, sj, ij, T2b.cptr(), W21aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   22] -- 
  // |-- [    0] --| W22(k,a1,j,a2) += (    1.00000000) D2(k,j,a1,a0) Fc1(a2,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.50000000) T2(w,a2,a1,i) W22(k,a1,j,a2) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W22aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no22_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO22_X0_TYPE0_NOERI)
      (sj, ij, W22aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    for(int si = 0;si < nir;++si){ 
    for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
      T2b = T2.get_amp2(ii);
      FC_FUNC(g_if_sigma_cooo_cooo_no22_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO22_X1_TYPE0_NOERI)
        (si, ii, sj, ij, T2b.cptr(), W22aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- No. [   23] -- 
  // |-- [    0] --| W23(k,a1,i,a2) += (    1.00000000) D2(k,i,a1,a0) Fc1(a2,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.50000000) T2(a2,w,a1,j) W23(k,a1,i,a2) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W23aaaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, 0));
  FC_FUNC(g_if_sigma_cooo_cooo_no23_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO23_X0_TYPE0_NOERI)
    (W23aaaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no23_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO23_X1_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W23aaaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   24] -- 
  // |-- [    0] --| W24(k,a1,i,a2) += (    1.00000000) D2(k,a0,a1,i) Fc1(a2,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.50000000) T2(w,a2,a1,j) W24(k,a1,i,a2) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W24aaaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, 0));
  FC_FUNC(g_if_sigma_cooo_cooo_no24_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO24_X0_TYPE0_NOERI)
    (W24aaaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no24_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO24_X1_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W24aaaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   25] -- 
  // |-- [    0] --| W25(k,a1) += (    1.00000000) D1(k,a0) Fc1(a1,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) T2(w,a1,i,j) W25(k,a1) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W25aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, 0));
  FC_FUNC(g_if_sigma_cooo_cooo_no25_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO25_X0_TYPE0_NOERI)
    (W25aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no25_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO25_X1_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W25aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   26] -- 
  // |-- [    0] --| W26(k,a1) += (    1.00000000) D1(k,a0) Fc1(a1,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.50000000) T2(a1,w,i,j) W26(k,a1) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W26aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, 0));
  FC_FUNC(g_if_sigma_cooo_cooo_no26_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO26_X0_TYPE0_NOERI)
    (W26aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no26_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO26_X1_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W26aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   27] -- 
  // |-- [    0] --| W27(k,a2,a1,j,i,a3) += (    1.00000000) D3(k,j,a2,i,a1,a0) Fc1(a3,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.50000000) T2(w,a3,a2,a1) W27(k,a2,a1,j,i,a3) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      orz::DTensor W27aaaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sa1^sj));
      FC_FUNC(g_if_sigma_cooo_cooo_no27_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO27_X0_TYPE0_NOERI)
        (sa1, ia1, sj, ij, W27aaaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_cooo_no27_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO27_X1_TYPE0_NOERI)
        (sa1, ia1, sj, ij, T2b.cptr(), W27aaaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   28] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -1.00000000) E0 T2(w,a0,a2,a1) D3(k,j,a2,i,a1,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_cooo_no28_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO28_X0_TYPE0_NOERI)
        (sa1, ia1, sj, ij, &E0, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   29] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    2.00000000) E0 T2(w,a0,i,a1) D2(k,j,a1,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_cooo_no29_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO29_X0_TYPE0_NOERI)
        (sa1, ia1, sj, ij, &E0, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   30] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -1.00000000) E0 T2(w,a0,a1,i) D2(k,j,a1,a0) 
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
      FC_FUNC(g_if_sigma_cooo_cooo_no30_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO30_X0_TYPE0_NOERI)
        (si, ii, sj, ij, &E0, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- No. [   31] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -1.00000000) E0 T2(a0,w,a1,j) D2(k,i,a1,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no31_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO31_X0_TYPE0_NOERI)
      (sj, ij, &E0, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   32] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -1.00000000) E0 T2(w,a0,a1,j) D2(k,a0,a1,i) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no32_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO32_X0_TYPE0_NOERI)
      (sj, ij, &E0, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   33] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    2.00000000) E0 T2(w,a0,i,j) D1(k,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no33_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO33_X0_TYPE0_NOERI)
      (sj, ij, &E0, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   34] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -1.00000000) E0 T2(a0,w,i,j) D1(k,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no34_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO34_X0_TYPE0_NOERI)
      (sj, ij, &E0, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   35] -- 
  // |-- [    0] --| W403(w,k,a0,j) += (    1.00000000) T2(a2,w,a1,j) D2(k,a0,a1,a2) 
  // |-- [    1] --| S2(w,k,i,j) += (   -0.50000000) Fc1(i,a0) W403(w,k,a0,j) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W403caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no35_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO35_X0_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W403caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no35_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO35_X1_TYPE0_NOERI)
      (sj, ij, W403caa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   36] -- 
  // |-- [    0] --| W404(w,k,a0,j) += (    1.00000000) T2(w,a2,a1,j) D2(k,a2,a1,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -0.50000000) Fc1(i,a0) W404(w,k,a0,j) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W404caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no36_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO36_X0_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W404caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no36_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO36_X1_TYPE0_NOERI)
      (sj, ij, W404caa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   37] -- 
  // |-- [    0] --| W405(w,k,j,a0) += (    1.00000000) T2(w,a3,a2,a1) D3(k,j,a2,a0,a1,a3) 
  // |-- [    1] --| S2(w,k,i,j) += (   -0.50000000) Fc1(i,a0) W405(w,k,j,a0) 
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W405caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
    for(int sa1 = 0;sa1 < nir;++sa1){ 
    for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
      T2b = T2.get_amp2(ia1);
      FC_FUNC(g_if_sigma_cooo_cooo_no37_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO37_X0_TYPE0_NOERI)
        (sa1, ia1, sj, ij, T2b.cptr(), W405caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
    } // End ia1
    } // End sa1
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no37_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO37_X1_TYPE0_NOERI)
      (sj, ij, W405caa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   38] -- 
  // |-- [    0] --| W406(w,k,a0,i) += (    1.00000000) T2(w,a2,i,a1) D2(k,a0,a1,a2) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) Fc1(j,a0) W406(w,k,a0,i) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W406caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no38_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO38_X0_TYPE0_NOERI)
      (sa1, ia1, T2b.cptr(), W406caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no38_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO38_X1_TYPE0_NOERI)
      (sj, ij, W406caaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   39] -- 
  // |-- [    0] --| W407(w,k,a0,i) += (    1.00000000) T2(w,a2,a1,i) D2(k,a0,a1,a2) 
  // |-- [    1] --| S2(w,k,i,j) += (   -0.50000000) Fc1(j,a0) W407(w,k,a0,i) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    orz::DTensor W407caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, si));
    FC_FUNC(g_if_sigma_cooo_cooo_no39_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO39_X0_TYPE0_NOERI)
      (si, ii, T2b.cptr(), W407caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_cooo_no39_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO39_X1_TYPE0_NOERI)
        (si, ii, sj, ij, W407caa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- No. [   40] -- 
  // |-- [    0] --| W408(w,k,a0,i) += (    1.00000000) T2(w,a3,a2,a1) D3(k,a0,a2,i,a1,a3) 
  // |-- [    1] --| S2(w,k,i,j) += (   -0.50000000) Fc1(j,a0) W408(w,k,a0,i) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W408caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no40_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO40_X0_TYPE0_NOERI)
      (sa1, ia1, T2b.cptr(), W408caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no40_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO40_X1_TYPE0_NOERI)
      (sj, ij, W408caaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   41] -- 
  // |-- [    0] --| W409(w,j,a0,i) += (    1.00000000) T2(w,a2,i,a1) D2(j,a0,a2,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) Fc1(k,a0) W409(w,j,a0,i) 
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W409caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
    for(int sa1 = 0;sa1 < nir;++sa1){ 
    for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
      T2b = T2.get_amp2(ia1);
      FC_FUNC(g_if_sigma_cooo_cooo_no41_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO41_X0_TYPE0_NOERI)
        (sa1, ia1, sj, ij, T2b.cptr(), W409caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
    } // End ia1
    } // End sa1
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no41_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO41_X1_TYPE0_NOERI)
      (sj, ij, W409caa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   42] -- 
  // |-- [    0] --| W410(w,i,a0,j) += (    1.00000000) T2(a2,w,a1,j) D2(i,a0,a2,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.50000000) Fc1(k,a0) W410(w,i,a0,j) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W410caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no42_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO42_X0_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W410caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no42_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO42_X1_TYPE0_NOERI)
      (sj, ij, W410caa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   43] -- 
  // |-- [    0] --| W411(w,j,a0,i) += (    1.00000000) T2(w,a2,a1,i) D2(j,a0,a2,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.50000000) Fc1(k,a0) W411(w,j,a0,i) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      orz::DTensor W411ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^si));
      FC_FUNC(g_if_sigma_cooo_cooo_no43_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO43_X0_TYPE0_NOERI)
        (si, ii, sj, ij, T2b.cptr(), W411ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_cooo_no43_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO43_X1_TYPE0_NOERI)
        (si, ii, sj, ij, W411ca_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- No. [   44] -- 
  // |-- [    0] --| W412(w,i,a0,j) += (    1.00000000) T2(w,a2,a1,j) D2(i,a1,a2,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.50000000) Fc1(k,a0) W412(w,i,a0,j) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W412caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no44_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO44_X0_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W412caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no44_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO44_X1_TYPE0_NOERI)
      (sj, ij, W412caa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   45] -- 
  // |-- [    0] --| W413(w,a0,i,j) += (    1.00000000) T2(w,a1,i,j) D1(a1,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) Fc1(k,a0) W413(w,a0,i,j) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W413caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no45_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO45_X0_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W413caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no45_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO45_X1_TYPE0_NOERI)
      (sj, ij, W413caa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   46] -- 
  // |-- [    0] --| W414(w,a0,j,i) += (    1.00000000) T2(a1,w,i,j) D1(a1,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.50000000) Fc1(k,a0) W414(w,a0,j,i) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W414caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no46_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO46_X0_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W414caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no46_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO46_X1_TYPE0_NOERI)
      (sj, ij, W414caa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   47] -- 
  // |-- [    0] --| W415(w,j,i,a0) += (    1.00000000) T2(w,a3,a2,a1) D3(j,a0,i,a2,a3,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.50000000) Fc1(k,a0) W415(w,j,i,a0) 
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W415caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
    for(int sa1 = 0;sa1 < nir;++sa1){ 
    for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
      T2b = T2.get_amp2(ia1);
      FC_FUNC(g_if_sigma_cooo_cooo_no47_x0_type0_noeri,G_IF_SIGMA_COOO_COOO_NO47_X0_TYPE0_NOERI)
        (sa1, ia1, sj, ij, T2b.cptr(), W415caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
    } // End ia1
    } // End sa1
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no47_x1_type0_noeri,G_IF_SIGMA_COOO_COOO_NO47_X1_TYPE0_NOERI)
      (sj, ij, W415caa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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

//-@type(2).declaration(begin)
  // --  Title : sigma_cooo_cooo
  //  >> Intermediates for the external contractions are defined here << 
  orz::DTensor W135ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W137ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W143ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W147caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W151caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W161caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W163ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W167caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W169caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W181caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W183ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W185caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W189caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W191caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W193caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W195caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W197ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W201caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W205caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W215caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W217ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W221caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W223caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W235caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W237ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W239caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W243caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W245caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W247caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W249caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W275ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W277ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W283ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W285caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W289caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W293caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W295caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W297caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W299caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W303ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W305caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W311caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W313caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W315caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W317caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W321caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W323ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W325caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W329caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W331caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W333caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W335caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W337ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W339caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W343caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W347caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W349caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W351caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W353caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W357ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W359caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W365caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W367caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W369caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W371caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W375caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W377ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W379caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W383caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W385caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W387caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W389caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W252ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W274ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W302ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W308ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W356ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W362ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
//-@type(2).declaration(end)

  //*-- Entering to take the type 0 contractions --*//
//-@type(0).contraction(begin)

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W135(c0,a2) += (    1.00000000) T2(c0,a0,a2,a1) D1(a0,a1) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no0_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO0_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W135ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W137(c0,a1) += (    1.00000000) T2(c0,a0,a2,a1) D1(a0,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no1_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO1_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W137ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W143(c0,a3) += (    1.00000000) T2(c0,a0,a2,a1) D2(a0,a1,a3,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no2_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO2_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W143ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W147(c0,a3,k,a2) += (    1.00000000) T2(c0,a0,a2,a1) D2(a0,a1,a3,k) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no3_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO3_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W147caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W151(c0,a3,k,a1) += (    1.00000000) T2(c0,a0,a2,a1) D2(a0,k,a3,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no4_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO4_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W151caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W161(c0,i,a4,a2) += (    1.00000000) T2(c0,a0,a2,a1) D2(a0,a1,i,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no5_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO5_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W161caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W163(c0,i) += (    1.00000000) T2(c0,a0,a2,a1) D2(a0,a1,i,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no6_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO6_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W163ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W167(c0,i,a4,a1) += (    1.00000000) T2(c0,a0,a2,a1) D2(a0,a2,i,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no7_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO7_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W167caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| W169(c0,a3,j,a0) += (    1.00000000) T2(c0,a0,a2,a1) D2(a3,a2,j,a1) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no8_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO8_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W169caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W181(c0,j,a4,a2) += (    1.00000000) T2(c0,a0,a2,a1) D2(a0,a1,j,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no9_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO9_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W181caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| W183(c0,j) += (    1.00000000) T2(c0,a0,a2,a1) D2(a0,a1,j,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no10_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO10_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W183ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| W185(c0,i,a3,a0) += (    1.00000000) T2(c0,a0,a2,a1) D2(i,a1,a3,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no11_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO11_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W185caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| W189(c0,j,a4,a1) += (    1.00000000) T2(c0,a0,a2,a1) D2(a0,a4,j,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no12_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO12_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W189caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| W191(c0,j,a4,a1) += (    1.00000000) T2(c0,a0,a2,a1) D2(a0,a4,j,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no13_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO13_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W191caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   14] -- 
  // |-- [    0] --| W193(c0,j,a4,a1) += (    1.00000000) T2(c0,a0,a2,a1) D2(a0,a2,j,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no14_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO14_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W193caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   15] -- 
  // |-- [    0] --| W195(c0,j,a4,a1) += (    1.00000000) T2(c0,a0,a2,a1) D2(a0,a2,j,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no15_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO15_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W195caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   16] -- 
  // |-- [    0] --| W197(c0,a3) += (    1.00000000) T2(c0,a0,a2,a1) C2(a0,a1,a3,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no16_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO16_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W197ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   17] -- 
  // |-- [    0] --| W201(c0,a3,k,a2) += (    1.00000000) T2(c0,a0,a2,a1) C2(a0,a1,a3,k) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no17_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO17_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W201caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   18] -- 
  // |-- [    0] --| W205(c0,a3,k,a1) += (    1.00000000) T2(c0,a0,a2,a1) C2(a0,k,a3,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no18_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO18_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W205caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   19] -- 
  // |-- [    0] --| W215(c0,i,a4,a2) += (    1.00000000) T2(c0,a0,a2,a1) C2(a0,a1,i,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no19_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO19_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W215caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   20] -- 
  // |-- [    0] --| W217(c0,i) += (    1.00000000) T2(c0,a0,a2,a1) C2(a0,a1,i,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no20_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO20_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W217ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   21] -- 
  // |-- [    0] --| W221(c0,i,a4,a1) += (    1.00000000) T2(c0,a0,a2,a1) C2(a0,a2,i,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no21_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO21_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W221caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   22] -- 
  // |-- [    0] --| W223(c0,a3,j,a0) += (    1.00000000) T2(c0,a0,a2,a1) C2(a1,j,a2,a3) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no22_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO22_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W223caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   23] -- 
  // |-- [    0] --| W235(c0,j,a4,a2) += (    1.00000000) T2(c0,a0,a2,a1) C2(a0,a1,j,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no23_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO23_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W235caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   24] -- 
  // |-- [    0] --| W237(c0,j) += (    1.00000000) T2(c0,a0,a2,a1) C2(a0,a1,j,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no24_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO24_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W237ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   25] -- 
  // |-- [    0] --| W239(c0,i,a3,a0) += (    1.00000000) T2(c0,a0,a2,a1) C2(a1,i,a2,a3) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no25_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO25_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W239caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   26] -- 
  // |-- [    0] --| W243(c0,j,a4,a1) += (    1.00000000) T2(c0,a0,a2,a1) C2(a0,a4,j,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no26_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO26_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W243caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   27] -- 
  // |-- [    0] --| W245(c0,j,a4,a1) += (    1.00000000) T2(c0,a0,a2,a1) C2(a0,a4,j,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no27_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO27_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W245caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   28] -- 
  // |-- [    0] --| W247(c0,j,a4,a1) += (    1.00000000) T2(c0,a0,a2,a1) C2(a0,a2,j,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no28_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO28_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W247caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   29] -- 
  // |-- [    0] --| W249(c0,j,a4,a1) += (    1.00000000) T2(c0,a0,a2,a1) C2(a0,a2,j,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no29_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO29_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W249caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   30] -- 
  // |-- [    0] --| W275(c0,a2) += (    1.00000000) T2(c0,a0,a2,a1) D1(a0,a1) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no30_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO30_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W275ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   31] -- 
  // |-- [    0] --| W277(c0,a1) += (    1.00000000) T2(c0,a0,a2,a1) D1(a0,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no31_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO31_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W277ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   32] -- 
  // |-- [    0] --| W283(c0,i) += (    1.00000000) T2(c0,a0,a2,a1) D2(a0,a1,i,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no32_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO32_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W283ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   33] -- 
  // |-- [    0] --| W285(c0,i,a4,a2) += (    1.00000000) T2(c0,a0,a2,a1) D2(a0,a1,i,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no33_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO33_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W285caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   34] -- 
  // |-- [    0] --| W289(c0,i,a4,a1) += (    1.00000000) T2(c0,a0,a2,a1) D2(a0,a4,i,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no34_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO34_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W289caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   35] -- 
  // |-- [    0] --| W293(c0,a3,j,a0) += (    1.00000000) T2(c0,a0,a2,a1) D2(a3,a1,j,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no35_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO35_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W293caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   36] -- 
  // |-- [    0] --| W295(c0,a3,j,a0) += (    1.00000000) T2(c0,a0,a2,a1) D2(a3,a2,j,a1) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no36_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO36_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W295caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   37] -- 
  // |-- [    0] --| W297(c0,a3,j,a0) += (    1.00000000) T2(c0,a0,a2,a1) D2(a3,a1,j,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no37_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO37_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W297caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   38] -- 
  // |-- [    0] --| W299(c0,a3,j,a0) += (    1.00000000) T2(c0,a0,a2,a1) D2(a3,a2,j,a1) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no38_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO38_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W299caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   39] -- 
  // |-- [    0] --| W303(c0,a3) += (    1.00000000) T2(c0,a0,a2,a1) D2(a0,a1,a3,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no39_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO39_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W303ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   40] -- 
  // |-- [    0] --| W305(c0,a3,k,a2) += (    1.00000000) T2(c0,a0,a2,a1) D2(a0,a1,a3,k) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no40_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO40_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W305caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   41] -- 
  // |-- [    0] --| W311(c0,a3,k,a1) += (    1.00000000) T2(c0,a0,a2,a1) D2(a0,a2,a3,k) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no41_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO41_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W311caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   42] -- 
  // |-- [    0] --| W313(c0,a3,k,a1) += (    1.00000000) T2(c0,a0,a2,a1) D2(a0,a2,a3,k) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no42_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO42_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W313caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   43] -- 
  // |-- [    0] --| W315(c0,a3,k,a1) += (    1.00000000) T2(c0,a0,a2,a1) D2(a0,k,a3,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no43_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO43_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W315caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   44] -- 
  // |-- [    0] --| W317(c0,a3,k,a1) += (    1.00000000) T2(c0,a0,a2,a1) D2(a0,k,a3,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no44_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO44_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W317caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   45] -- 
  // |-- [    0] --| W321(c0,j,a4,a2) += (    1.00000000) T2(c0,a0,a2,a1) D2(a0,a1,j,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no45_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO45_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W321caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   46] -- 
  // |-- [    0] --| W323(c0,j) += (    1.00000000) T2(c0,a0,a2,a1) D2(a0,a1,j,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no46_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO46_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W323ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   47] -- 
  // |-- [    0] --| W325(c0,a3,i,a0) += (    1.00000000) T2(c0,a0,a2,a1) D2(a3,a1,i,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no47_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO47_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W325caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   48] -- 
  // |-- [    0] --| W329(c0,j,a4,a1) += (    1.00000000) T2(c0,a0,a2,a1) D2(a0,a4,j,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no48_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO48_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W329caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   49] -- 
  // |-- [    0] --| W331(c0,j,a4,a1) += (    1.00000000) T2(c0,a0,a2,a1) D2(a0,a4,j,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no49_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO49_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W331caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   50] -- 
  // |-- [    0] --| W333(c0,j,a4,a1) += (    1.00000000) T2(c0,a0,a2,a1) D2(a0,a2,j,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no50_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO50_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W333caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   51] -- 
  // |-- [    0] --| W335(c0,j,a4,a1) += (    1.00000000) T2(c0,a0,a2,a1) D2(a0,a2,j,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no51_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO51_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W335caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   52] -- 
  // |-- [    0] --| W337(c0,i) += (    1.00000000) T2(c0,a0,a2,a1) C2(a0,a1,i,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no52_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO52_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W337ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   53] -- 
  // |-- [    0] --| W339(c0,i,a4,a2) += (    1.00000000) T2(c0,a0,a2,a1) C2(a0,a1,i,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no53_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO53_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W339caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   54] -- 
  // |-- [    0] --| W343(c0,i,a4,a1) += (    1.00000000) T2(c0,a0,a2,a1) C2(a0,a4,i,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no54_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO54_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W343caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   55] -- 
  // |-- [    0] --| W347(c0,a3,j,a0) += (    1.00000000) T2(c0,a0,a2,a1) C2(a1,a3,a2,j) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no55_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO55_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W347caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   56] -- 
  // |-- [    0] --| W349(c0,a3,j,a0) += (    1.00000000) T2(c0,a0,a2,a1) C2(a1,j,a2,a3) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no56_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO56_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W349caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   57] -- 
  // |-- [    0] --| W351(c0,a3,j,a0) += (    1.00000000) T2(c0,a0,a2,a1) C2(a1,a3,a2,j) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no57_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO57_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W351caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   58] -- 
  // |-- [    0] --| W353(c0,a3,j,a0) += (    1.00000000) T2(c0,a0,a2,a1) C2(a1,j,a2,a3) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no58_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO58_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W353caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   59] -- 
  // |-- [    0] --| W357(c0,a3) += (    1.00000000) T2(c0,a0,a2,a1) C2(a0,a1,a3,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no59_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO59_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W357ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   60] -- 
  // |-- [    0] --| W359(c0,a3,k,a2) += (    1.00000000) T2(c0,a0,a2,a1) C2(a0,a1,a3,k) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no60_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO60_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W359caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   61] -- 
  // |-- [    0] --| W365(c0,a3,k,a1) += (    1.00000000) T2(c0,a0,a2,a1) C2(a0,a2,a3,k) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no61_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO61_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W365caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   62] -- 
  // |-- [    0] --| W367(c0,a3,k,a1) += (    1.00000000) T2(c0,a0,a2,a1) C2(a0,a2,a3,k) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no62_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO62_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W367caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   63] -- 
  // |-- [    0] --| W369(c0,a3,k,a1) += (    1.00000000) T2(c0,a0,a2,a1) C2(a0,k,a3,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no63_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO63_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W369caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   64] -- 
  // |-- [    0] --| W371(c0,a3,k,a1) += (    1.00000000) T2(c0,a0,a2,a1) C2(a0,k,a3,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no64_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO64_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W371caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   65] -- 
  // |-- [    0] --| W375(c0,j,a4,a2) += (    1.00000000) T2(c0,a0,a2,a1) C2(a0,a1,j,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no65_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO65_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W375caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   66] -- 
  // |-- [    0] --| W377(c0,j) += (    1.00000000) T2(c0,a0,a2,a1) C2(a0,a1,j,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no66_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO66_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W377ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   67] -- 
  // |-- [    0] --| W379(c0,a3,i,a0) += (    1.00000000) T2(c0,a0,a2,a1) C2(a1,a3,a2,i) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no67_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO67_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W379caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   68] -- 
  // |-- [    0] --| W383(c0,j,a4,a1) += (    1.00000000) T2(c0,a0,a2,a1) C2(a0,a4,j,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no68_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO68_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W383caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   69] -- 
  // |-- [    0] --| W385(c0,j,a4,a1) += (    1.00000000) T2(c0,a0,a2,a1) C2(a0,a4,j,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no69_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO69_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W385caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   70] -- 
  // |-- [    0] --| W387(c0,j,a4,a1) += (    1.00000000) T2(c0,a0,a2,a1) C2(a0,a2,j,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no70_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO70_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W387caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   71] -- 
  // |-- [    0] --| W389(c0,j,a4,a1) += (    1.00000000) T2(c0,a0,a2,a1) C2(a0,a2,j,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no71_x0_type0_eri_c,G_IF_SIGMA_COOO_COOO_NO71_X0_TYPE0_ERI_C)
      (sa1, ia1, T2b.cptr(), W389caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

//-@type(0).contraction(end)

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
  // -- Title : sigma_cooo_cooo

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W37(w,a0,j,a2) += (    1.00000000) V2(w,a2,c0,a1) T2(c0,a0,a1,j) 
  // |-- [    1] --| S2(w,k,i,j) += (   -2.00000000) D2(k,a0,a2,i) W37(w,a0,j,a2) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W37aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no0_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO0_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), V2_sym.cptr(), W37aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no0_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO0_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W37aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W38(w,a0,j,a2) += (    1.00000000) V2(w,c0,a2,a1) T2(c0,a0,a1,j) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) D2(k,a0,a2,i) W38(w,a0,j,a2) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W38aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no1_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO1_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), V2_sym.cptr(), W38aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no1_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO1_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W38aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W42(w,a0,j,i) += (    1.00000000) V2(w,i,c0,a1) T2(c0,a0,a1,j) 
  // |-- [    1] --| S2(w,k,i,j) += (    4.00000000) D1(k,a0) W42(w,a0,j,i) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W42aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no2_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO2_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), V2_sym.cptr(), W42aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no2_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO2_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W42aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W50(w,a0,j,i) += (    1.00000000) V2(w,c0,i,a1) T2(c0,a0,a1,j) 
  // |-- [    1] --| S2(w,k,i,j) += (   -2.00000000) D1(k,a0) W50(w,a0,j,i) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W50aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no3_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO3_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), V2_sym.cptr(), W50aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no3_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO3_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W50aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W63(w,a0,j,a2) += (    1.00000000) V2(w,a2,c0,a1) T2(a0,c0,a1,j) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) D2(k,a0,a2,i) W63(w,a0,j,a2) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W63aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no4_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO4_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), V2_sym.cptr(), W63aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no4_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO4_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W63aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W64(w,a0,j,a2) += (    1.00000000) V2(w,c0,a2,a1) T2(a0,c0,a1,j) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) D2(k,i,a2,a0) W64(w,a0,j,a2) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W64aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no5_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO5_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), V2_sym.cptr(), W64aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no5_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO5_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W64aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W68(w,a0,j,i) += (    1.00000000) V2(w,i,c0,a1) T2(a0,c0,a1,j) 
  // |-- [    1] --| S2(w,k,i,j) += (   -2.00000000) D1(k,a0) W68(w,a0,j,i) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W68aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no6_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO6_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), V2_sym.cptr(), W68aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no6_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO6_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W68aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W76(w,a0,j,i) += (    1.00000000) V2(w,c0,i,a1) T2(a0,c0,a1,j) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) D1(k,a0) W76(w,a0,j,i) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W76aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no7_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO7_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), V2_sym.cptr(), W76aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no7_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO7_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W76aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| W87(w,a0,a1,i,a3,a2) += (    1.00000000) V2(w,c0,a3,a2) T2(c0,a0,a1,i) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) D3(k,j,a3,a2,a1,a0) W87(w,a0,a1,i,a3,a2) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    orz::DTensor W87aaaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sw^si));
    FC_FUNC(g_if_sigma_cooo_cooo_no8_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO8_X0_TYPE1_ERI_C)
      (si, ii, sw, iw, T2b.cptr(), V2_sym.cptr(), W87aaaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_cooo_no8_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO8_X1_TYPE1_ERI_C)
        (si, ii, sj, ij, sw, iw, W87aaaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W88(w,a0,a1,i,a3,a2) += (    1.00000000) V2(w,a3,c0,a2) T2(c0,a0,a1,i) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) D3(k,j,a3,a0,a1,a2) W88(w,a0,a1,i,a3,a2) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    orz::DTensor W88aaaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sw^si));
    FC_FUNC(g_if_sigma_cooo_cooo_no9_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO9_X0_TYPE1_ERI_C)
      (si, ii, sw, iw, T2b.cptr(), V2_sym.cptr(), W88aaaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_cooo_no9_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO9_X1_TYPE1_ERI_C)
        (si, ii, sj, ij, sw, iw, W88aaaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| W89(c0,k,j,a3) += (    1.00000000) T2(c0,a0,a2,a1) D3(k,j,a2,a3,a1,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) V2(c0,w,i,a3) W89(c0,k,j,a3) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <3, 2> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W89aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sj));
    for(int sa1 = 0;sa1 < nir;++sa1){ 
    for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
      T2b = T2.get_amp2(ia1);
      FC_FUNC(g_if_sigma_cooo_cooo_no10_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO10_X0_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, T2b.cptr(), W89aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
    } // End ia1
    } // End sa1
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no10_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO10_X1_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, V2_sym.cptr(), W89aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| W90(c0,k,j,a3) += (    1.00000000) T2(c0,a0,a2,a1) D3(k,j,a2,a3,a1,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -2.00000000) V2(c0,a3,w,i) W90(c0,k,j,a3) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <3, 2> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W90aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sj));
    for(int sa1 = 0;sa1 < nir;++sa1){ 
    for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
      T2b = T2.get_amp2(ia1);
      FC_FUNC(g_if_sigma_cooo_cooo_no11_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO11_X0_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, T2b.cptr(), W90aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
    } // End ia1
    } // End sa1
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no11_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO11_X1_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, V2_sym.cptr(), W90aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| W91(w,a0,i,a1,a3,a2) += (    1.00000000) V2(w,c0,a3,a2) T2(c0,a0,i,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (   -2.00000000) D3(k,j,a3,a2,a1,a0) W91(w,a0,i,a1,a3,a2) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W91aaaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sw^sa1));
    FC_FUNC(g_if_sigma_cooo_cooo_no12_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO12_X0_TYPE1_ERI_C)
      (sa1, ia1, sw, iw, T2b.cptr(), V2_sym.cptr(), W91aaaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_cooo_no12_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO12_X1_TYPE1_ERI_C)
        (sa1, ia1, sj, ij, sw, iw, W91aaaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| W92(w,a0,i,a1,a3,a2) += (    1.00000000) V2(w,a3,c0,a2) T2(c0,a0,i,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) D3(k,j,a3,a2,a1,a0) W92(w,a0,i,a1,a3,a2) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W92aaaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sw^sa1));
    FC_FUNC(g_if_sigma_cooo_cooo_no13_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO13_X0_TYPE1_ERI_C)
      (sa1, ia1, sw, iw, T2b.cptr(), V2_sym.cptr(), W92aaaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_cooo_no13_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO13_X1_TYPE1_ERI_C)
        (sa1, ia1, sj, ij, sw, iw, W92aaaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   14] -- 
  // |-- [    0] --| W93(w,a0,a1,j,a3,a2) += (    1.00000000) V2(w,c0,a3,a2) T2(c0,a0,a1,j) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) D3(k,a0,a3,a2,a1,i) W93(w,a0,a1,j,a3,a2) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W93aaaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sw^sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no14_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO14_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), V2_sym.cptr(), W93aaaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no14_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO14_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W93aaaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   15] -- 
  // |-- [    0] --| W94(w,a0,a1,j,a3,a2) += (    1.00000000) V2(w,a3,c0,a2) T2(c0,a0,a1,j) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) D3(k,a0,a3,i,a1,a2) W94(w,a0,a1,j,a3,a2) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W94aaaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sw^sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no15_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO15_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), V2_sym.cptr(), W94aaaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no15_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO15_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W94aaaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   16] -- 
  // |-- [    0] --| W97(w,a0,j,a1,a3,a2) += (    1.00000000) V2(w,c0,a3,a2) T2(a0,c0,a1,j) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) D3(k,i,a3,a2,a1,a0) W97(w,a0,j,a1,a3,a2) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W97aaaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sw^sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no16_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO16_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), V2_sym.cptr(), W97aaaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no16_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO16_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W97aaaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   17] -- 
  // |-- [    0] --| W98(w,a0,j,a1,a3,a2) += (    1.00000000) V2(w,a3,c0,a2) T2(a0,c0,a1,j) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) D3(k,a2,a3,i,a1,a0) W98(w,a0,j,a1,a3,a2) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W98aaaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sw^sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no17_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO17_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), V2_sym.cptr(), W98aaaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no17_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO17_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W98aaaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   18] -- 
  // |-- [    0] --| W101(w,c0,k,a0) += (    1.00000000) V2(w,c0,a2,a1) D2(k,a0,a2,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) T2(a0,c0,i,j) W101(w,c0,k,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W101caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_cooo_no18_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO18_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W101caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no18_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO18_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), W101caa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   19] -- 
  // |-- [    0] --| W102(w,c0,k,a0) += (    1.00000000) V2(w,a2,c0,a1) D2(k,a1,a2,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) T2(a0,c0,i,j) W102(w,c0,k,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W102caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_cooo_no19_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO19_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W102caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no19_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO19_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), W102caa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   20] -- 
  // |-- [    0] --| W103(c0,k,a2,j) += (    1.00000000) T2(c0,a0,a1,j) D2(k,a0,a1,a2) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) V2(c0,w,i,a2) W103(c0,k,a2,j) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W103aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no20_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO20_X0_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, T2b.cptr(), W103aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no20_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO20_X1_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, V2_sym.cptr(), W103aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   21] -- 
  // |-- [    0] --| W104(c0,k,a2,j) += (    1.00000000) T2(c0,a0,a1,j) D2(k,a0,a1,a2) 
  // |-- [    1] --| S2(w,k,i,j) += (   -2.00000000) V2(c0,a2,w,i) W104(c0,k,a2,j) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W104aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no21_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO21_X0_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, T2b.cptr(), W104aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no21_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO21_X1_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, V2_sym.cptr(), W104aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   22] -- 
  // |-- [    0] --| W105(c0,k,a2,j) += (    1.00000000) T2(a0,c0,a1,j) D2(k,a2,a1,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) V2(c0,w,i,a2) W105(c0,k,a2,j) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W105aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no22_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO22_X0_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, T2b.cptr(), W105aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no22_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO22_X1_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, V2_sym.cptr(), W105aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   23] -- 
  // |-- [    0] --| W106(c0,k,a2,j) += (    1.00000000) T2(a0,c0,a1,j) D2(k,a2,a1,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -2.00000000) V2(c0,a2,w,i) W106(c0,k,a2,j) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W106aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no23_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO23_X0_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, T2b.cptr(), W106aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no23_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO23_X1_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, V2_sym.cptr(), W106aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   24] -- 
  // |-- [    0] --| W107(w,c0,k,a0) += (    1.00000000) V2(w,c0,a2,a1) D2(k,a0,a2,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (   -2.00000000) T2(c0,a0,i,j) W107(w,c0,k,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W107caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_cooo_no24_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO24_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W107caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no24_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO24_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), W107caa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   25] -- 
  // |-- [    0] --| W108(w,c0,k,a0) += (    1.00000000) V2(w,a2,c0,a1) D2(k,a0,a2,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) T2(c0,a0,i,j) W108(w,c0,k,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W108caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_cooo_no25_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO25_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W108caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no25_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO25_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), W108caa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   26] -- 
  // |-- [    0] --| W113(w,a0,a2,a1,a4,a3) += (    1.00000000) V2(w,a4,c0,a3) T2(c0,a0,a2,a1) 
  // |-- [    1] --| W114(w,j,a0,a2,a4,a3) += (    1.00000000) D1(j,a1) W113(w,a0,a2,a1,a4,a3) 
  // |-- [    2] --| S2(w,k,i,j) += (   -0.50000000) D3(a0,k,a3,a2,i,a4) W114(w,j,a0,a2,a4,a3) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 1 <3, 2> 
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 2
  // Type1 anormality found for external indices of the contraction >> 0 <<
  orz::DTensor W113aaaaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaaaa(symblockinfo, sw));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no26_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO26_X0_TYPE1_ERI_C)
      (sa1, ia1, sw, iw, T2b.cptr(), V2_sym.cptr(), W113aaaaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W114aaaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sw^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no26_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO26_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W113aaaaa_sigma_cooo_cooo.cptr(), W114aaaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no26_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO26_X2_TYPE1_ERI_C)
      (sj, ij, sw, iw, W114aaaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   27] -- 
  // |-- [    0] --| W115(w,a0,a2,a1,a4,a3) += (    1.00000000) V2(w,a4,c0,a3) T2(c0,a0,a2,a1) 
  // |-- [    1] --| W116(w,j,a0,a1,a4,a3) += (    1.00000000) D1(j,a2) W115(w,a0,a2,a1,a4,a3) 
  // |-- [    2] --| S2(w,k,i,j) += (   -0.50000000) D3(a0,a1,a3,k,i,a4) W116(w,j,a0,a1,a4,a3) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 0
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W115aaaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sw^sa1));
    FC_FUNC(g_if_sigma_cooo_cooo_no27_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO27_X0_TYPE1_ERI_C)
      (sa1, ia1, sw, iw, T2b.cptr(), V2_sym.cptr(), W115aaaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      orz::DTensor W116aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw^sj^sa1));
      FC_FUNC(g_if_sigma_cooo_cooo_no27_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO27_X1_TYPE1_ERI_C)
        (sa1, ia1, sj, ij, sw, iw, W115aaaa_sigma_cooo_cooo.cptr(), W116aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_cooo_no27_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO27_X2_TYPE1_ERI_C)
        (sa1, ia1, sj, ij, sw, iw, W116aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   28] -- 
  // |-- [    0] --| W118(c0,a3,i,k) += (    1.00000000) T2(c0,a0,a2,a1) D3(a0,a1,a3,a2,i,k) 
  // |-- [    1] --| W117(w,c0,j,a3) += (    1.00000000) V2(c0,a3,w,a4) D1(j,a4) 
  // |-- [    2] --| S2(w,k,i,j) += (   -0.50000000) W117(w,c0,j,a3) W118(c0,a3,i,k) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W118aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sc0));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no28_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO28_X0_TYPE1_ERI_C)
      (sa1, ia1, sc0, ic0, T2b.cptr(), W118aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W117ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sc0^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no28_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO28_X1_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, V2_sym.cptr(), W117ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no28_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO28_X2_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, W117ca_sigma_cooo_cooo.cptr(), W118aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   29] -- 
  // |-- [    0] --| W119(w,c0,i,a3) += (    1.00000000) V2(c0,a3,w,a4) D1(i,a4) 
  // |-- [    1] --| W120(c0,a3,j,k) += (    1.00000000) T2(c0,a0,a2,a1) D3(a0,a1,a3,a2,j,k) 
  // |-- [    2] --| S2(w,k,i,j) += (    1.00000000) W119(w,c0,i,a3) W120(c0,a3,j,k) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// ==> makeContractions: [1] Special care is taken for the contraction 1 <3, 2> 
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 1
  orz::DTensor W119caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sc0));
  FC_FUNC(g_if_sigma_cooo_cooo_no29_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO29_X0_TYPE1_ERI_C)
    (sc0, ic0, V2_sym.cptr(), W119caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W120aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sj));
    for(int sa1 = 0;sa1 < nir;++sa1){ 
    for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
      T2b = T2.get_amp2(ia1);
      FC_FUNC(g_if_sigma_cooo_cooo_no29_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO29_X1_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, T2b.cptr(), W120aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
    } // End ia1
    } // End sa1
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no29_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO29_X2_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, W119caa_sigma_cooo_cooo.cptr(), W120aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   30] -- 
  // |-- [    0] --| W121(w,a0,a2,a1,a4,a3) += (    1.00000000) V2(w,a4,c0,a3) T2(c0,a0,a2,a1) 
  // |-- [    1] --| W122(w,i,a0,a2,a4,a3) += (    1.00000000) D1(i,a1) W121(w,a0,a2,a1,a4,a3) 
  // |-- [    2] --| S2(w,k,i,j) += (   -0.50000000) D3(a0,a4,a3,a2,j,k) W122(w,i,a0,a2,a4,a3) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 1 <2, 1> 
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 2
  orz::DTensor W122aaaaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaaaa(symblockinfo, sw));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W121aaaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sw^sa1));
    FC_FUNC(g_if_sigma_cooo_cooo_no30_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO30_X0_TYPE1_ERI_C)
      (sa1, ia1, sw, iw, T2b.cptr(), V2_sym.cptr(), W121aaaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no30_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO30_X1_TYPE1_ERI_C)
      (sa1, ia1, sw, iw, W121aaaa_sigma_cooo_cooo.cptr(), W122aaaaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no30_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO30_X2_TYPE1_ERI_C)
      (sj, ij, sw, iw, W122aaaaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   31] -- 
  // |-- [    0] --| W123(w,a0,a2,a1,a4,a3) += (    1.00000000) V2(w,a4,c0,a3) T2(c0,a0,a2,a1) 
  // |-- [    1] --| W124(w,i,a0,a1,a4,a3) += (    1.00000000) D1(i,a2) W123(w,a0,a2,a1,a4,a3) 
  // |-- [    2] --| S2(w,k,i,j) += (   -0.50000000) D3(a0,a1,a3,a4,j,k) W124(w,i,a0,a1,a4,a3) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 0
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W123aaaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sw^sa1));
    orz::DTensor W124aaaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sw^sa1));
    FC_FUNC(g_if_sigma_cooo_cooo_no31_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO31_X0_TYPE1_ERI_C)
      (sa1, ia1, sw, iw, T2b.cptr(), V2_sym.cptr(), W123aaaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no31_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO31_X1_TYPE1_ERI_C)
      (sa1, ia1, sw, iw, W123aaaa_sigma_cooo_cooo.cptr(), W124aaaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_cooo_no31_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO31_X2_TYPE1_ERI_C)
        (sa1, ia1, sj, ij, sw, iw, W124aaaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   32] -- 
  // |-- [    0] --| W127(w,c0,a2,a4) += (    1.00000000) V2(w,a4,c0,a3) D1(a3,a2) 
  // |-- [    1] --| W128(w,a0,a1,a4) += (    1.00000000) T2(c0,a0,a2,a1) W127(w,c0,a2,a4) 
  // |-- [    2] --| S2(w,k,i,j) += (    1.00000000) D3(a0,a1,i,a4,j,k) W128(w,a0,a1,a4) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 0
  orz::DTensor W127caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_cooo_no32_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO32_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W127caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W128aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sa1));
    FC_FUNC(g_if_sigma_cooo_cooo_no32_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO32_X1_TYPE1_ERI_C)
      (sa1, ia1, sw, iw, T2b.cptr(), W127caa_sigma_cooo_cooo.cptr(), W128aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_cooo_no32_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO32_X2_TYPE1_ERI_C)
        (sa1, ia1, sj, ij, sw, iw, W128aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   33] -- 
  // |-- [    0] --| W129(w,c0,a1,a4) += (    1.00000000) V2(w,a4,c0,a3) D1(a3,a1) 
  // |-- [    1] --| W130(w,a0,a2,a4) += (    1.00000000) T2(c0,a0,a2,a1) W129(w,c0,a1,a4) 
  // |-- [    2] --| S2(w,k,i,j) += (   -0.50000000) D3(a0,a2,i,a4,j,k) W130(w,a0,a2,a4) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 1 <2, 1> 
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 2
  orz::DTensor W130aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W129ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sw^sa1));
    FC_FUNC(g_if_sigma_cooo_cooo_no33_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO33_X0_TYPE1_ERI_C)
      (sa1, ia1, sw, iw, V2_sym.cptr(), W129ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no33_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO33_X1_TYPE1_ERI_C)
      (sa1, ia1, sw, iw, T2b.cptr(), W129ca_sigma_cooo_cooo.cptr(), W130aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no33_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO33_X2_TYPE1_ERI_C)
      (sj, ij, sw, iw, W130aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   34] -- 
  // |-- [    0] --| W131(w,c0) += (    1.00000000) V2(c0,a3,w,a4) D1(a3,a4) 
  // |-- [    1] --| W132(c0,i,j,k) += (    1.00000000) T2(c0,a0,a2,a1) D3(a0,a1,i,a2,j,k) 
  // |-- [    2] --| S2(w,k,i,j) += (   -0.50000000) W131(w,c0) W132(c0,i,j,k) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// ==> makeContractions: [1] Special care is taken for the contraction 1 <3, 2> 
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 1
  orz::DTensor W131c_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, sc0));
  FC_FUNC(g_if_sigma_cooo_cooo_no34_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO34_X0_TYPE1_ERI_C)
    (sc0, ic0, V2_sym.cptr(), W131c_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W132aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sj));
    for(int sa1 = 0;sa1 < nir;++sa1){ 
    for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
      T2b = T2.get_amp2(ia1);
      FC_FUNC(g_if_sigma_cooo_cooo_no34_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO34_X1_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, T2b.cptr(), W132aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
    } // End ia1
    } // End sa1
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no34_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO34_X2_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, W131c_sigma_cooo_cooo.cptr(), W132aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   35] -- 
  // |-- [    0] --| W136(w,a4,a3,a2) += (    1.00000000) V2(w,a4,c0,a3) W135(c0,a2) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) D3(a3,a2,i,a4,j,k) W136(w,a4,a3,a2) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W136aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_cooo_no35_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO35_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W135ca_sigma_cooo_cooo.cptr(), W136aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no35_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO35_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W136aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   36] -- 
  // |-- [    0] --| W138(w,a4,a3,a1) += (    1.00000000) V2(w,a4,c0,a3) W137(c0,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (   -0.50000000) D3(a3,a1,i,a4,j,k) W138(w,a4,a3,a1) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W138aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_cooo_no36_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO36_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W137ca_sigma_cooo_cooo.cptr(), W138aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no36_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO36_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W138aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   37] -- 
  // |-- [    0] --| W139(w,c0,a0,a3) += (    1.00000000) V2(w,a4,c0,a3) D1(a0,a4) 
  // |-- [    1] --| W140(w,a2,a1,a3) += (    1.00000000) T2(c0,a0,a2,a1) W139(w,c0,a0,a3) 
  // |-- [    2] --| S2(w,k,i,j) += (   -0.50000000) D3(a3,a2,i,a1,j,k) W140(w,a2,a1,a3) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 0
  orz::DTensor W139caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_cooo_no37_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO37_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W139caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W140aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sa1));
    FC_FUNC(g_if_sigma_cooo_cooo_no37_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO37_X1_TYPE1_ERI_C)
      (sa1, ia1, sw, iw, T2b.cptr(), W139caa_sigma_cooo_cooo.cptr(), W140aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_cooo_no37_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO37_X2_TYPE1_ERI_C)
        (sa1, ia1, sj, ij, sw, iw, W140aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   38] -- 
  // |-- [    0] --| W141(w,a0,a2,a1,a4,a3) += (    1.00000000) V2(w,a4,c0,a3) T2(c0,a0,a2,a1) 
  // |-- [    1] --| W142(w,k,a2,a1,a4,a3) += (    1.00000000) D1(a0,k) W141(w,a0,a2,a1,a4,a3) 
  // |-- [    2] --| S2(w,k,i,j) += (   -0.50000000) D3(a3,a2,i,a4,j,a1) W142(w,k,a2,a1,a4,a3) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 0
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W141aaaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sw^sa1));
    orz::DTensor W142aaaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sw^sa1));
    FC_FUNC(g_if_sigma_cooo_cooo_no38_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO38_X0_TYPE1_ERI_C)
      (sa1, ia1, sw, iw, T2b.cptr(), V2_sym.cptr(), W141aaaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no38_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO38_X1_TYPE1_ERI_C)
      (sa1, ia1, sw, iw, W141aaaa_sigma_cooo_cooo.cptr(), W142aaaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_cooo_no38_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO38_X2_TYPE1_ERI_C)
        (sa1, ia1, sj, ij, sw, iw, W142aaaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   39] -- 
  // |-- [    0] --| W144(w,a4) += (    1.00000000) V2(w,a4,c0,a3) W143(c0,a3) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) D2(i,a4,j,k) W144(w,a4) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W144a_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_cooo_no39_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO39_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W143ca_sigma_cooo_cooo.cptr(), W144a_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no39_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO39_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W144a_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   40] -- 
  // |-- [    0] --| W148(w,a4,k,a2) += (    1.00000000) V2(w,a4,c0,a3) W147(c0,a3,k,a2) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.50000000) D2(i,a4,j,a2) W148(w,a4,k,a2) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W148aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_cooo_no40_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO40_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W147caaa_sigma_cooo_cooo.cptr(), W148aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no40_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO40_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W148aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   41] -- 
  // |-- [    0] --| W152(w,a4,k,a1) += (    1.00000000) V2(w,a4,c0,a3) W151(c0,a3,k,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.50000000) D2(i,a4,j,a1) W152(w,a4,k,a1) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W152aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_cooo_no41_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO41_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W151caaa_sigma_cooo_cooo.cptr(), W152aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no41_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO41_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W152aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   42] -- 
  // |-- [    0] --| W153(w,c0,a0,k) += (    1.00000000) V2(c0,a3,w,a4) D2(a0,a4,a3,k) 
  // |-- [    1] --| W154(c0,i,j,a0) += (    1.00000000) T2(c0,a0,a2,a1) D2(i,a1,j,a2) 
  // |-- [    2] --| S2(w,k,i,j) += (   -0.33333333) W153(w,c0,a0,k) W154(c0,i,j,a0) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// ==> makeContractions: [1] Special care is taken for the contraction 1 <3, 2> 
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 1
  orz::DTensor W153caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sc0));
  FC_FUNC(g_if_sigma_cooo_cooo_no42_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO42_X0_TYPE1_ERI_C)
    (sc0, ic0, V2_sym.cptr(), W153caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W154aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sj));
    for(int sa1 = 0;sa1 < nir;++sa1){ 
    for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
      T2b = T2.get_amp2(ia1);
      FC_FUNC(g_if_sigma_cooo_cooo_no42_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO42_X1_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, T2b.cptr(), W154aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
    } // End ia1
    } // End sa1
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no42_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO42_X2_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, W153caa_sigma_cooo_cooo.cptr(), W154aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   43] -- 
  // |-- [    0] --| W155(w,c0,a0,k) += (    1.00000000) V2(c0,a3,w,a4) D2(a0,a4,a3,k) 
  // |-- [    1] --| W156(c0,i,j,a0) += (    1.00000000) T2(c0,a0,a2,a1) D2(i,a2,j,a1) 
  // |-- [    2] --| S2(w,k,i,j) += (   -0.16666667) W155(w,c0,a0,k) W156(c0,i,j,a0) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// ==> makeContractions: [1] Special care is taken for the contraction 1 <3, 2> 
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 1
  orz::DTensor W155caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sc0));
  FC_FUNC(g_if_sigma_cooo_cooo_no43_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO43_X0_TYPE1_ERI_C)
    (sc0, ic0, V2_sym.cptr(), W155caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W156aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sj));
    for(int sa1 = 0;sa1 < nir;++sa1){ 
    for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
      T2b = T2.get_amp2(ia1);
      FC_FUNC(g_if_sigma_cooo_cooo_no43_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO43_X1_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, T2b.cptr(), W156aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
    } // End ia1
    } // End sa1
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no43_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO43_X2_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, W155caa_sigma_cooo_cooo.cptr(), W156aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   44] -- 
  // |-- [    0] --| W157(w,c0,a0,k) += (    1.00000000) V2(c0,a3,w,a4) D2(a0,k,a3,a4) 
  // |-- [    1] --| W158(c0,i,j,a0) += (    1.00000000) T2(c0,a0,a2,a1) D2(i,a1,j,a2) 
  // |-- [    2] --| S2(w,k,i,j) += (   -0.16666667) W157(w,c0,a0,k) W158(c0,i,j,a0) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// ==> makeContractions: [1] Special care is taken for the contraction 1 <3, 2> 
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 1
  orz::DTensor W157caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sc0));
  FC_FUNC(g_if_sigma_cooo_cooo_no44_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO44_X0_TYPE1_ERI_C)
    (sc0, ic0, V2_sym.cptr(), W157caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W158aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sj));
    for(int sa1 = 0;sa1 < nir;++sa1){ 
    for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
      T2b = T2.get_amp2(ia1);
      FC_FUNC(g_if_sigma_cooo_cooo_no44_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO44_X1_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, T2b.cptr(), W158aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
    } // End ia1
    } // End sa1
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no44_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO44_X2_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, W157caa_sigma_cooo_cooo.cptr(), W158aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   45] -- 
  // |-- [    0] --| W159(w,c0,a0,k) += (    1.00000000) V2(c0,a3,w,a4) D2(a0,k,a3,a4) 
  // |-- [    1] --| W160(c0,i,j,a0) += (    1.00000000) T2(c0,a0,a2,a1) D2(i,a2,j,a1) 
  // |-- [    2] --| S2(w,k,i,j) += (   -0.33333333) W159(w,c0,a0,k) W160(c0,i,j,a0) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// ==> makeContractions: [1] Special care is taken for the contraction 1 <3, 2> 
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 1
  orz::DTensor W159caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sc0));
  FC_FUNC(g_if_sigma_cooo_cooo_no45_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO45_X0_TYPE1_ERI_C)
    (sc0, ic0, V2_sym.cptr(), W159caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W160aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sj));
    for(int sa1 = 0;sa1 < nir;++sa1){ 
    for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
      T2b = T2.get_amp2(ia1);
      FC_FUNC(g_if_sigma_cooo_cooo_no45_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO45_X1_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, T2b.cptr(), W160aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
    } // End ia1
    } // End sa1
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no45_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO45_X2_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, W159caa_sigma_cooo_cooo.cptr(), W160aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   46] -- 
  // |-- [    0] --| W162(w,a3,i,a2) += (    1.00000000) V2(w,a4,c0,a3) W161(c0,i,a4,a2) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) D2(a3,a2,j,k) W162(w,a3,i,a2) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W162aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_cooo_no46_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO46_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W161caaa_sigma_cooo_cooo.cptr(), W162aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no46_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO46_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W162aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   47] -- 
  // |-- [    0] --| W164(w,a4,a3,i) += (    1.00000000) V2(w,a4,c0,a3) W163(c0,i) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.50000000) D2(a3,a4,j,k) W164(w,a4,a3,i) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W164aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_cooo_no47_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO47_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W163ca_sigma_cooo_cooo.cptr(), W164aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no47_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO47_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W164aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   48] -- 
  // |-- [    0] --| W166(c0,i,k,a2) += (    1.00000000) T2(c0,a0,a2,a1) D2(a0,a1,i,k) 
  // |-- [    1] --| W165(w,c0,j,a2) += (    1.00000000) V2(c0,a3,w,a4) D2(a3,a2,j,a4) 
  // |-- [    2] --| S2(w,k,i,j) += (    0.50000000) W165(w,c0,j,a2) W166(c0,i,k,a2) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W166aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sc0));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no48_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO48_X0_TYPE1_ERI_C)
      (sa1, ia1, sc0, ic0, T2b.cptr(), W166aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W165ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sc0^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no48_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO48_X1_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, V2_sym.cptr(), W165ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no48_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO48_X2_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, W165ca_sigma_cooo_cooo.cptr(), W166aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   49] -- 
  // |-- [    0] --| W168(w,a3,i,a1) += (    1.00000000) V2(w,a4,c0,a3) W167(c0,i,a4,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.50000000) D2(a3,a1,j,k) W168(w,a3,i,a1) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W168aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_cooo_no49_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO49_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W167caaa_sigma_cooo_cooo.cptr(), W168aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no49_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO49_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W168aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   50] -- 
  // |-- [    0] --| W170(w,a4,j,a0) += (    1.00000000) V2(w,a4,c0,a3) W169(c0,a3,j,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.50000000) D2(a0,k,i,a4) W170(w,a4,j,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W170aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no50_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO50_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, V2_sym.cptr(), W169caaa_sigma_cooo_cooo.cptr(), W170aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no50_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO50_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W170aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   51] -- 
  // |-- [    0] --| W172(c0,i,k,a1) += (    1.00000000) T2(c0,a0,a2,a1) D2(a0,a2,i,k) 
  // |-- [    1] --| W171(w,c0,j,a1) += (    1.00000000) V2(c0,a3,w,a4) D2(a3,a1,j,a4) 
  // |-- [    2] --| S2(w,k,i,j) += (   -0.33333333) W171(w,c0,j,a1) W172(c0,i,k,a1) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 0
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W172aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sa1));
    FC_FUNC(g_if_sigma_cooo_cooo_no51_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO51_X0_TYPE1_ERI_C)
      (sa1, ia1, sc0, ic0, T2b.cptr(), W172aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      orz::DTensor W171c_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, sc0^sj^sa1));
      FC_FUNC(g_if_sigma_cooo_cooo_no51_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO51_X1_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, V2_sym.cptr(), W171c_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_cooo_no51_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO51_X2_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, W171c_sigma_cooo_cooo.cptr(), W172aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   52] -- 
  // |-- [    0] --| W174(c0,i,k,a1) += (    1.00000000) T2(c0,a0,a2,a1) D2(a0,a2,i,k) 
  // |-- [    1] --| W173(w,c0,j,a1) += (    1.00000000) V2(c0,a3,w,a4) D2(a3,a4,j,a1) 
  // |-- [    2] --| S2(w,k,i,j) += (   -0.16666667) W173(w,c0,j,a1) W174(c0,i,k,a1) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W174aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sa1));
    FC_FUNC(g_if_sigma_cooo_cooo_no52_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO52_X0_TYPE1_ERI_C)
      (sa1, ia1, sc0, ic0, T2b.cptr(), W174aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      orz::DTensor W173c_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, sc0^sj^sa1));
      FC_FUNC(g_if_sigma_cooo_cooo_no52_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO52_X1_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, V2_sym.cptr(), W173c_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
      FC_FUNC(g_if_sigma_cooo_cooo_no52_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO52_X2_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, W173c_sigma_cooo_cooo.cptr(), W174aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   53] -- 
  // |-- [    0] --| W176(c0,i,k,a1) += (    1.00000000) T2(c0,a0,a2,a1) D2(a0,k,i,a2) 
  // |-- [    1] --| W175(w,c0,j,a1) += (    1.00000000) V2(c0,a3,w,a4) D2(a3,a1,j,a4) 
  // |-- [    2] --| S2(w,k,i,j) += (   -0.16666667) W175(w,c0,j,a1) W176(c0,i,k,a1) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 0
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W176aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sa1));
    FC_FUNC(g_if_sigma_cooo_cooo_no53_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO53_X0_TYPE1_ERI_C)
      (sa1, ia1, sc0, ic0, T2b.cptr(), W176aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      orz::DTensor W175c_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, sc0^sj^sa1));
      FC_FUNC(g_if_sigma_cooo_cooo_no53_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO53_X1_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, V2_sym.cptr(), W175c_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_cooo_no53_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO53_X2_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, W175c_sigma_cooo_cooo.cptr(), W176aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   54] -- 
  // |-- [    0] --| W178(c0,i,k,a1) += (    1.00000000) T2(c0,a0,a2,a1) D2(a0,k,i,a2) 
  // |-- [    1] --| W177(w,c0,j,a1) += (    1.00000000) V2(c0,a3,w,a4) D2(a3,a4,j,a1) 
  // |-- [    2] --| S2(w,k,i,j) += (   -0.33333333) W177(w,c0,j,a1) W178(c0,i,k,a1) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W178aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sa1));
    FC_FUNC(g_if_sigma_cooo_cooo_no54_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO54_X0_TYPE1_ERI_C)
      (sa1, ia1, sc0, ic0, T2b.cptr(), W178aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      orz::DTensor W177c_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, sc0^sj^sa1));
      FC_FUNC(g_if_sigma_cooo_cooo_no54_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO54_X1_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, V2_sym.cptr(), W177c_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
      FC_FUNC(g_if_sigma_cooo_cooo_no54_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO54_X2_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, W177c_sigma_cooo_cooo.cptr(), W178aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   55] -- 
  // |-- [    0] --| W179(w,c0,i,a2) += (    1.00000000) V2(c0,a3,w,a4) D2(i,a4,a3,a2) 
  // |-- [    1] --| W180(c0,j,k,a2) += (    1.00000000) T2(c0,a0,a2,a1) D2(a0,a1,j,k) 
  // |-- [    2] --| S2(w,k,i,j) += (   -1.00000000) W179(w,c0,i,a2) W180(c0,j,k,a2) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// ==> makeContractions: [1] Special care is taken for the contraction 1 <3, 2> 
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 1
  orz::DTensor W179caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sc0));
  FC_FUNC(g_if_sigma_cooo_cooo_no55_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO55_X0_TYPE1_ERI_C)
    (sc0, ic0, V2_sym.cptr(), W179caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W180aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sj));
    for(int sa1 = 0;sa1 < nir;++sa1){ 
    for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
      T2b = T2.get_amp2(ia1);
      FC_FUNC(g_if_sigma_cooo_cooo_no55_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO55_X1_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, T2b.cptr(), W180aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
    } // End ia1
    } // End sa1
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no55_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO55_X2_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, W179caa_sigma_cooo_cooo.cptr(), W180aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   56] -- 
  // |-- [    0] --| W182(w,a3,j,a2) += (    1.00000000) V2(w,a4,c0,a3) W181(c0,j,a4,a2) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.50000000) D2(i,k,a3,a2) W182(w,a3,j,a2) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W182aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no56_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO56_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, V2_sym.cptr(), W181caaa_sigma_cooo_cooo.cptr(), W182aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no56_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO56_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W182aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   57] -- 
  // |-- [    0] --| W184(w,a4,a3,j) += (    1.00000000) V2(w,a4,c0,a3) W183(c0,j) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.50000000) D2(i,a4,a3,k) W184(w,a4,a3,j) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W184aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no57_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO57_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, V2_sym.cptr(), W183ca_sigma_cooo_cooo.cptr(), W184aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no57_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO57_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W184aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   58] -- 
  // |-- [    0] --| W186(w,a4,i,a0) += (    1.00000000) V2(w,a4,c0,a3) W185(c0,i,a3,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.50000000) D2(a0,a4,j,k) W186(w,a4,i,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W186aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_cooo_no58_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO58_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W185caaa_sigma_cooo_cooo.cptr(), W186aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no58_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO58_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W186aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   59] -- 
  // |-- [    0] --| W187(w,c0,i,a1) += (    1.00000000) V2(c0,a3,w,a4) D2(i,a4,a3,a1) 
  // |-- [    1] --| W188(c0,j,k,a1) += (    1.00000000) T2(c0,a0,a2,a1) D2(a0,a2,j,k) 
  // |-- [    2] --| S2(w,k,i,j) += (    0.50000000) W187(w,c0,i,a1) W188(c0,j,k,a1) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W187ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sc0^sa1));
    FC_FUNC(g_if_sigma_cooo_cooo_no59_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO59_X0_TYPE1_ERI_C)
      (sa1, ia1, sc0, ic0, V2_sym.cptr(), W187ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      orz::DTensor W188a_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xa(symblockinfo, sc0^sj^sa1));
      FC_FUNC(g_if_sigma_cooo_cooo_no59_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO59_X1_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, T2b.cptr(), W188a_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
      FC_FUNC(g_if_sigma_cooo_cooo_no59_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO59_X2_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, W187ca_sigma_cooo_cooo.cptr(), W188a_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   60] -- 
  // |-- [    0] --| W190(w,a3,j,a1) += (    1.00000000) V2(w,a4,c0,a3) W189(c0,j,a4,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (   -0.33333333) D2(i,a1,a3,k) W190(w,a3,j,a1) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W190aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no60_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO60_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, V2_sym.cptr(), W189caaa_sigma_cooo_cooo.cptr(), W190aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no60_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO60_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W190aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   61] -- 
  // |-- [    0] --| W192(w,a3,j,a1) += (    1.00000000) V2(w,a4,c0,a3) W191(c0,j,a4,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (   -0.16666667) D2(i,k,a3,a1) W192(w,a3,j,a1) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W192aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no61_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO61_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, V2_sym.cptr(), W191caaa_sigma_cooo_cooo.cptr(), W192aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no61_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO61_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W192aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   62] -- 
  // |-- [    0] --| W194(w,a3,j,a1) += (    1.00000000) V2(w,a4,c0,a3) W193(c0,j,a4,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (   -0.16666667) D2(i,a1,a3,k) W194(w,a3,j,a1) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W194aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no62_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO62_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, V2_sym.cptr(), W193caaa_sigma_cooo_cooo.cptr(), W194aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no62_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO62_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W194aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   63] -- 
  // |-- [    0] --| W196(w,a3,j,a1) += (    1.00000000) V2(w,a4,c0,a3) W195(c0,j,a4,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (   -0.33333333) D2(i,k,a3,a1) W196(w,a3,j,a1) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W196aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no63_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO63_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, V2_sym.cptr(), W195caaa_sigma_cooo_cooo.cptr(), W196aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no63_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO63_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W196aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   64] -- 
  // |-- [    0] --| W198(w,a4) += (    1.00000000) V2(w,a4,c0,a3) W197(c0,a3) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) C2(a4,i,k,j) W198(w,a4) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W198a_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_cooo_no64_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO64_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W197ca_sigma_cooo_cooo.cptr(), W198a_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no64_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO64_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W198a_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   65] -- 
  // |-- [    0] --| W202(w,a4,k,a2) += (    1.00000000) V2(w,a4,c0,a3) W201(c0,a3,k,a2) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) C2(a2,j,a4,i) W202(w,a4,k,a2) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W202aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_cooo_no65_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO65_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W201caaa_sigma_cooo_cooo.cptr(), W202aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no65_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO65_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W202aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   66] -- 
  // |-- [    0] --| W206(w,a4,k,a1) += (    1.00000000) V2(w,a4,c0,a3) W205(c0,a3,k,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) C2(a1,j,a4,i) W206(w,a4,k,a1) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W206aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_cooo_no66_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO66_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W205caaa_sigma_cooo_cooo.cptr(), W206aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no66_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO66_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W206aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   67] -- 
  // |-- [    0] --| W207(w,c0,a0,k) += (    1.00000000) V2(c0,a3,w,a4) C2(a0,a4,a3,k) 
  // |-- [    1] --| W208(c0,i,j,a0) += (    1.00000000) T2(c0,a0,a2,a1) C2(a1,i,a2,j) 
  // |-- [    2] --| S2(w,k,i,j) += (    0.66666667) W207(w,c0,a0,k) W208(c0,i,j,a0) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// ==> makeContractions: [1] Special care is taken for the contraction 1 <3, 2> 
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 1
  orz::DTensor W207caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sc0));
  FC_FUNC(g_if_sigma_cooo_cooo_no67_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO67_X0_TYPE1_ERI_C)
    (sc0, ic0, V2_sym.cptr(), W207caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W208aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sj));
    for(int sa1 = 0;sa1 < nir;++sa1){ 
    for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
      T2b = T2.get_amp2(ia1);
      FC_FUNC(g_if_sigma_cooo_cooo_no67_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO67_X1_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, T2b.cptr(), W208aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
    } // End ia1
    } // End sa1
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no67_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO67_X2_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, W207caa_sigma_cooo_cooo.cptr(), W208aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   68] -- 
  // |-- [    0] --| W209(w,c0,a0,k) += (    1.00000000) V2(c0,a3,w,a4) C2(a0,a4,a3,k) 
  // |-- [    1] --| W210(c0,i,j,a0) += (    1.00000000) T2(c0,a0,a2,a1) C2(a1,j,a2,i) 
  // |-- [    2] --| S2(w,k,i,j) += (    0.33333333) W209(w,c0,a0,k) W210(c0,i,j,a0) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// ==> makeContractions: [1] Special care is taken for the contraction 1 <3, 2> 
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 1
  orz::DTensor W209caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sc0));
  FC_FUNC(g_if_sigma_cooo_cooo_no68_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO68_X0_TYPE1_ERI_C)
    (sc0, ic0, V2_sym.cptr(), W209caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W210aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sj));
    for(int sa1 = 0;sa1 < nir;++sa1){ 
    for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
      T2b = T2.get_amp2(ia1);
      FC_FUNC(g_if_sigma_cooo_cooo_no68_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO68_X1_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, T2b.cptr(), W210aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
    } // End ia1
    } // End sa1
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no68_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO68_X2_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, W209caa_sigma_cooo_cooo.cptr(), W210aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   69] -- 
  // |-- [    0] --| W211(w,c0,a0,k) += (    1.00000000) V2(c0,a3,w,a4) C2(a0,k,a3,a4) 
  // |-- [    1] --| W212(c0,i,j,a0) += (    1.00000000) T2(c0,a0,a2,a1) C2(a1,i,a2,j) 
  // |-- [    2] --| S2(w,k,i,j) += (    0.33333333) W211(w,c0,a0,k) W212(c0,i,j,a0) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// ==> makeContractions: [1] Special care is taken for the contraction 1 <3, 2> 
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 1
  orz::DTensor W211caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sc0));
  FC_FUNC(g_if_sigma_cooo_cooo_no69_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO69_X0_TYPE1_ERI_C)
    (sc0, ic0, V2_sym.cptr(), W211caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W212aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sj));
    for(int sa1 = 0;sa1 < nir;++sa1){ 
    for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
      T2b = T2.get_amp2(ia1);
      FC_FUNC(g_if_sigma_cooo_cooo_no69_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO69_X1_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, T2b.cptr(), W212aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
    } // End ia1
    } // End sa1
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no69_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO69_X2_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, W211caa_sigma_cooo_cooo.cptr(), W212aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   70] -- 
  // |-- [    0] --| W213(w,c0,a0,k) += (    1.00000000) V2(c0,a3,w,a4) C2(a0,k,a3,a4) 
  // |-- [    1] --| W214(c0,i,j,a0) += (    1.00000000) T2(c0,a0,a2,a1) C2(a1,j,a2,i) 
  // |-- [    2] --| S2(w,k,i,j) += (    0.66666667) W213(w,c0,a0,k) W214(c0,i,j,a0) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// ==> makeContractions: [1] Special care is taken for the contraction 1 <3, 2> 
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 1
  orz::DTensor W213caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sc0));
  FC_FUNC(g_if_sigma_cooo_cooo_no70_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO70_X0_TYPE1_ERI_C)
    (sc0, ic0, V2_sym.cptr(), W213caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W214aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sj));
    for(int sa1 = 0;sa1 < nir;++sa1){ 
    for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
      T2b = T2.get_amp2(ia1);
      FC_FUNC(g_if_sigma_cooo_cooo_no70_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO70_X1_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, T2b.cptr(), W214aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
    } // End ia1
    } // End sa1
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no70_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO70_X2_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, W213caa_sigma_cooo_cooo.cptr(), W214aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   71] -- 
  // |-- [    0] --| W216(w,a3,i,a2) += (    1.00000000) V2(w,a4,c0,a3) W215(c0,i,a4,a2) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) C2(a2,a3,k,j) W216(w,a3,i,a2) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W216aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_cooo_no71_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO71_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W215caaa_sigma_cooo_cooo.cptr(), W216aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no71_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO71_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W216aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   72] -- 
  // |-- [    0] --| W218(w,a4,a3,i) += (    1.00000000) V2(w,a4,c0,a3) W217(c0,i) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) C2(a3,a4,j,k) W218(w,a4,a3,i) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W218aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_cooo_no72_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO72_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W217ca_sigma_cooo_cooo.cptr(), W218aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no72_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO72_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W218aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   73] -- 
  // |-- [    0] --| W220(c0,i,k,a2) += (    1.00000000) T2(c0,a0,a2,a1) C2(a0,a1,i,k) 
  // |-- [    1] --| W219(w,c0,j,a2) += (    1.00000000) V2(c0,a3,w,a4) C2(a2,a3,a4,j) 
  // |-- [    2] --| S2(w,k,i,j) += (   -1.00000000) W219(w,c0,j,a2) W220(c0,i,k,a2) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W220aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sc0));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no73_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO73_X0_TYPE1_ERI_C)
      (sa1, ia1, sc0, ic0, T2b.cptr(), W220aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W219ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sc0^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no73_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO73_X1_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, V2_sym.cptr(), W219ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no73_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO73_X2_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, W219ca_sigma_cooo_cooo.cptr(), W220aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   74] -- 
  // |-- [    0] --| W222(w,a3,i,a1) += (    1.00000000) V2(w,a4,c0,a3) W221(c0,i,a4,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) C2(a1,a3,k,j) W222(w,a3,i,a1) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W222aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_cooo_no74_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO74_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W221caaa_sigma_cooo_cooo.cptr(), W222aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no74_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO74_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W222aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   75] -- 
  // |-- [    0] --| W224(w,a4,j,a0) += (    1.00000000) V2(w,a4,c0,a3) W223(c0,a3,j,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) C2(a0,k,i,a4) W224(w,a4,j,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W224aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no75_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO75_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, V2_sym.cptr(), W223caaa_sigma_cooo_cooo.cptr(), W224aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no75_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO75_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W224aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   76] -- 
  // |-- [    0] --| W226(c0,i,k,a1) += (    1.00000000) T2(c0,a0,a2,a1) C2(a0,a2,i,k) 
  // |-- [    1] --| W225(w,c0,j,a1) += (    1.00000000) V2(c0,a3,w,a4) C2(a1,a3,a4,j) 
  // |-- [    2] --| S2(w,k,i,j) += (    0.66666667) W225(w,c0,j,a1) W226(c0,i,k,a1) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 0
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W226aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sa1));
    FC_FUNC(g_if_sigma_cooo_cooo_no76_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO76_X0_TYPE1_ERI_C)
      (sa1, ia1, sc0, ic0, T2b.cptr(), W226aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      orz::DTensor W225c_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, sc0^sj^sa1));
      FC_FUNC(g_if_sigma_cooo_cooo_no76_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO76_X1_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, V2_sym.cptr(), W225c_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_cooo_no76_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO76_X2_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, W225c_sigma_cooo_cooo.cptr(), W226aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   77] -- 
  // |-- [    0] --| W228(c0,i,k,a1) += (    1.00000000) T2(c0,a0,a2,a1) C2(a0,a2,i,k) 
  // |-- [    1] --| W227(w,c0,j,a1) += (    1.00000000) V2(c0,a3,w,a4) C2(a1,j,a4,a3) 
  // |-- [    2] --| S2(w,k,i,j) += (    0.33333333) W227(w,c0,j,a1) W228(c0,i,k,a1) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 0
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W228aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sa1));
    FC_FUNC(g_if_sigma_cooo_cooo_no77_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO77_X0_TYPE1_ERI_C)
      (sa1, ia1, sc0, ic0, T2b.cptr(), W228aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      orz::DTensor W227c_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, sc0^sj^sa1));
      FC_FUNC(g_if_sigma_cooo_cooo_no77_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO77_X1_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, V2_sym.cptr(), W227c_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_cooo_no77_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO77_X2_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, W227c_sigma_cooo_cooo.cptr(), W228aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   78] -- 
  // |-- [    0] --| W230(c0,i,k,a1) += (    1.00000000) T2(c0,a0,a2,a1) C2(a0,k,i,a2) 
  // |-- [    1] --| W229(w,c0,j,a1) += (    1.00000000) V2(c0,a3,w,a4) C2(a1,a3,a4,j) 
  // |-- [    2] --| S2(w,k,i,j) += (    0.33333333) W229(w,c0,j,a1) W230(c0,i,k,a1) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 0
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W230aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sa1));
    FC_FUNC(g_if_sigma_cooo_cooo_no78_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO78_X0_TYPE1_ERI_C)
      (sa1, ia1, sc0, ic0, T2b.cptr(), W230aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      orz::DTensor W229c_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, sc0^sj^sa1));
      FC_FUNC(g_if_sigma_cooo_cooo_no78_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO78_X1_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, V2_sym.cptr(), W229c_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_cooo_no78_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO78_X2_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, W229c_sigma_cooo_cooo.cptr(), W230aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   79] -- 
  // |-- [    0] --| W232(c0,i,k,a1) += (    1.00000000) T2(c0,a0,a2,a1) C2(a0,k,i,a2) 
  // |-- [    1] --| W231(w,c0,j,a1) += (    1.00000000) V2(c0,a3,w,a4) C2(a1,j,a4,a3) 
  // |-- [    2] --| S2(w,k,i,j) += (    0.66666667) W231(w,c0,j,a1) W232(c0,i,k,a1) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 0
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W232aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sa1));
    FC_FUNC(g_if_sigma_cooo_cooo_no79_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO79_X0_TYPE1_ERI_C)
      (sa1, ia1, sc0, ic0, T2b.cptr(), W232aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      orz::DTensor W231c_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, sc0^sj^sa1));
      FC_FUNC(g_if_sigma_cooo_cooo_no79_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO79_X1_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, V2_sym.cptr(), W231c_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_cooo_no79_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO79_X2_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, W231c_sigma_cooo_cooo.cptr(), W232aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   80] -- 
  // |-- [    0] --| W233(w,c0,i,a2) += (    1.00000000) V2(c0,a3,w,a4) C2(a2,a3,a4,i) 
  // |-- [    1] --| W234(c0,j,k,a2) += (    1.00000000) T2(c0,a0,a2,a1) C2(a0,a1,j,k) 
  // |-- [    2] --| S2(w,k,i,j) += (    2.00000000) W233(w,c0,i,a2) W234(c0,j,k,a2) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// ==> makeContractions: [1] Special care is taken for the contraction 1 <3, 2> 
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 1
  orz::DTensor W233caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sc0));
  FC_FUNC(g_if_sigma_cooo_cooo_no80_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO80_X0_TYPE1_ERI_C)
    (sc0, ic0, V2_sym.cptr(), W233caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W234aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sj));
    for(int sa1 = 0;sa1 < nir;++sa1){ 
    for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
      T2b = T2.get_amp2(ia1);
      FC_FUNC(g_if_sigma_cooo_cooo_no80_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO80_X1_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, T2b.cptr(), W234aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
    } // End ia1
    } // End sa1
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no80_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO80_X2_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, W233caa_sigma_cooo_cooo.cptr(), W234aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   81] -- 
  // |-- [    0] --| W236(w,a3,j,a2) += (    1.00000000) V2(w,a4,c0,a3) W235(c0,j,a4,a2) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) C2(a2,a3,k,i) W236(w,a3,j,a2) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W236aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no81_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO81_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, V2_sym.cptr(), W235caaa_sigma_cooo_cooo.cptr(), W236aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no81_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO81_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W236aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   82] -- 
  // |-- [    0] --| W238(w,a4,a3,j) += (    1.00000000) V2(w,a4,c0,a3) W237(c0,j) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) C2(a3,k,i,a4) W238(w,a4,a3,j) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W238aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no82_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO82_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, V2_sym.cptr(), W237ca_sigma_cooo_cooo.cptr(), W238aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no82_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO82_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W238aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   83] -- 
  // |-- [    0] --| W240(w,a4,i,a0) += (    1.00000000) V2(w,a4,c0,a3) W239(c0,i,a3,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) C2(a0,a4,j,k) W240(w,a4,i,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W240aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_cooo_no83_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO83_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W239caaa_sigma_cooo_cooo.cptr(), W240aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no83_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO83_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W240aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   84] -- 
  // |-- [    0] --| W241(w,c0,i,a1) += (    1.00000000) V2(c0,a3,w,a4) C2(a1,a3,a4,i) 
  // |-- [    1] --| W242(c0,j,k,a1) += (    1.00000000) T2(c0,a0,a2,a1) C2(a0,a2,j,k) 
  // |-- [    2] --| S2(w,k,i,j) += (   -1.00000000) W241(w,c0,i,a1) W242(c0,j,k,a1) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W241ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sc0^sa1));
    FC_FUNC(g_if_sigma_cooo_cooo_no84_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO84_X0_TYPE1_ERI_C)
      (sa1, ia1, sc0, ic0, V2_sym.cptr(), W241ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      orz::DTensor W242a_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xa(symblockinfo, sc0^sj^sa1));
      FC_FUNC(g_if_sigma_cooo_cooo_no84_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO84_X1_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, T2b.cptr(), W242a_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
      FC_FUNC(g_if_sigma_cooo_cooo_no84_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO84_X2_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, W241ca_sigma_cooo_cooo.cptr(), W242a_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   85] -- 
  // |-- [    0] --| W244(w,a3,j,a1) += (    1.00000000) V2(w,a4,c0,a3) W243(c0,j,a4,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.66666667) C2(a1,i,k,a3) W244(w,a3,j,a1) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W244aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no85_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO85_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, V2_sym.cptr(), W243caaa_sigma_cooo_cooo.cptr(), W244aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no85_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO85_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W244aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   86] -- 
  // |-- [    0] --| W246(w,a3,j,a1) += (    1.00000000) V2(w,a4,c0,a3) W245(c0,j,a4,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.33333333) C2(a1,a3,k,i) W246(w,a3,j,a1) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W246aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no86_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO86_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, V2_sym.cptr(), W245caaa_sigma_cooo_cooo.cptr(), W246aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no86_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO86_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W246aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   87] -- 
  // |-- [    0] --| W248(w,a3,j,a1) += (    1.00000000) V2(w,a4,c0,a3) W247(c0,j,a4,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.33333333) C2(a1,i,k,a3) W248(w,a3,j,a1) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W248aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no87_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO87_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, V2_sym.cptr(), W247caaa_sigma_cooo_cooo.cptr(), W248aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no87_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO87_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W248aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   88] -- 
  // |-- [    0] --| W250(w,a3,j,a1) += (    1.00000000) V2(w,a4,c0,a3) W249(c0,j,a4,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.66666667) C2(a1,a3,k,i) W250(w,a3,j,a1) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W250aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no88_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO88_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, V2_sym.cptr(), W249caaa_sigma_cooo_cooo.cptr(), W250aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no88_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO88_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W250aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   89] -- 
  // |-- [    0] --| W251(c0,i,a3,a4) += (    1.00000000) T2(c0,a0,a2,a1) D3(a0,a1,i,a2,a3,a4) 
  // |-- [    1] --| W252(w,i) += (    1.00000000) V2(c0,w,a4,a3) W251(c0,i,a3,a4) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
  // Pref: 2
  orz::DTensor W251aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sc0));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no89_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO89_X0_TYPE1_ERI_C)
      (sa1, ia1, sc0, ic0, T2b.cptr(), W251aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  FC_FUNC(g_if_sigma_cooo_cooo_no89_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO89_X1_TYPE1_ERI_C)
    (sc0, ic0, V2_sym.cptr(), W251aaa_sigma_cooo_cooo.cptr(), W252ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   90] -- 
  // |-- [    0] --| W253(w,a0,a2,a1,a4,a3) += (    1.00000000) V2(w,c0,a4,a3) T2(c0,a0,a2,a1) 
  // |-- [    1] --| W254(w,j,a0,a2,a4,a3) += (    1.00000000) D1(j,a1) W253(w,a0,a2,a1,a4,a3) 
  // |-- [    2] --| S2(w,k,i,j) += (   -0.50000000) D3(a0,k,i,a2,a3,a4) W254(w,j,a0,a2,a4,a3) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 1 <3, 2> 
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 2
  // Type1 anormality found for external indices of the contraction >> 0 <<
  orz::DTensor W253aaaaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaaaa(symblockinfo, sw));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no90_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO90_X0_TYPE1_ERI_C)
      (sa1, ia1, sw, iw, T2b.cptr(), V2_sym.cptr(), W253aaaaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W254aaaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sw^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no90_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO90_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W253aaaaa_sigma_cooo_cooo.cptr(), W254aaaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no90_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO90_X2_TYPE1_ERI_C)
      (sj, ij, sw, iw, W254aaaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   91] -- 
  // |-- [    0] --| W255(w,a0,a2,a1,a4,a3) += (    1.00000000) V2(w,c0,a4,a3) T2(c0,a0,a2,a1) 
  // |-- [    1] --| W256(w,j,a0,a1,a4,a3) += (    1.00000000) D1(j,a2) W255(w,a0,a2,a1,a4,a3) 
  // |-- [    2] --| S2(w,k,i,j) += (   -0.50000000) D3(a0,a1,i,k,a3,a4) W256(w,j,a0,a1,a4,a3) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 0
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W255aaaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sw^sa1));
    FC_FUNC(g_if_sigma_cooo_cooo_no91_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO91_X0_TYPE1_ERI_C)
      (sa1, ia1, sw, iw, T2b.cptr(), V2_sym.cptr(), W255aaaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      orz::DTensor W256aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw^sj^sa1));
      FC_FUNC(g_if_sigma_cooo_cooo_no91_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO91_X1_TYPE1_ERI_C)
        (sa1, ia1, sj, ij, sw, iw, W255aaaa_sigma_cooo_cooo.cptr(), W256aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_cooo_no91_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO91_X2_TYPE1_ERI_C)
        (sa1, ia1, sj, ij, sw, iw, W256aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   92] -- 
  // |-- [    0] --| W258(c0,i,a3,k) += (    1.00000000) T2(c0,a0,a2,a1) D3(a0,a1,i,a2,a3,k) 
  // |-- [    1] --| W257(w,c0,j,a3) += (    1.00000000) V2(c0,w,a4,a3) D1(j,a4) 
  // |-- [    2] --| S2(w,k,i,j) += (   -0.50000000) W257(w,c0,j,a3) W258(c0,i,a3,k) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W258aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sc0));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no92_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO92_X0_TYPE1_ERI_C)
      (sa1, ia1, sc0, ic0, T2b.cptr(), W258aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W257ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sc0^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no92_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO92_X1_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, V2_sym.cptr(), W257ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no92_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO92_X2_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, W257ca_sigma_cooo_cooo.cptr(), W258aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   93] -- 
  // |-- [    0] --| W259(w,c0) += (    1.00000000) V2(c0,w,a4,a3) D1(a3,a4) 
  // |-- [    1] --| W260(c0,i,j,k) += (    1.00000000) T2(c0,a0,a2,a1) D3(a0,a1,i,a2,j,k) 
  // |-- [    2] --| S2(w,k,i,j) += (    1.00000000) W259(w,c0) W260(c0,i,j,k) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// ==> makeContractions: [1] Special care is taken for the contraction 1 <3, 2> 
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 1
  orz::DTensor W259c_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, sc0));
  FC_FUNC(g_if_sigma_cooo_cooo_no93_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO93_X0_TYPE1_ERI_C)
    (sc0, ic0, V2_sym.cptr(), W259c_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W260aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sj));
    for(int sa1 = 0;sa1 < nir;++sa1){ 
    for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
      T2b = T2.get_amp2(ia1);
      FC_FUNC(g_if_sigma_cooo_cooo_no93_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO93_X1_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, T2b.cptr(), W260aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
    } // End ia1
    } // End sa1
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no93_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO93_X2_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, W259c_sigma_cooo_cooo.cptr(), W260aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   94] -- 
  // |-- [    0] --| W261(w,c0,a1,a4) += (    1.00000000) V2(w,c0,a4,a3) D1(a3,a1) 
  // |-- [    1] --| W262(w,a0,a2,a4) += (    1.00000000) T2(c0,a0,a2,a1) W261(w,c0,a1,a4) 
  // |-- [    2] --| S2(w,k,i,j) += (   -0.50000000) D3(a0,a4,i,a2,j,k) W262(w,a0,a2,a4) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 1 <2, 1> 
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 2
  orz::DTensor W262aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W261ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sw^sa1));
    FC_FUNC(g_if_sigma_cooo_cooo_no94_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO94_X0_TYPE1_ERI_C)
      (sa1, ia1, sw, iw, V2_sym.cptr(), W261ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no94_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO94_X1_TYPE1_ERI_C)
      (sa1, ia1, sw, iw, T2b.cptr(), W261ca_sigma_cooo_cooo.cptr(), W262aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no94_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO94_X2_TYPE1_ERI_C)
      (sj, ij, sw, iw, W262aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   95] -- 
  // |-- [    0] --| W263(w,c0,a2,a4) += (    1.00000000) V2(w,c0,a4,a3) D1(a3,a2) 
  // |-- [    1] --| W264(w,a0,a1,a4) += (    1.00000000) T2(c0,a0,a2,a1) W263(w,c0,a2,a4) 
  // |-- [    2] --| S2(w,k,i,j) += (   -0.50000000) D3(a0,a1,i,a4,j,k) W264(w,a0,a1,a4) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 0
  orz::DTensor W263caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_cooo_no95_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO95_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W263caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W264aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sa1));
    FC_FUNC(g_if_sigma_cooo_cooo_no95_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO95_X1_TYPE1_ERI_C)
      (sa1, ia1, sw, iw, T2b.cptr(), W263caa_sigma_cooo_cooo.cptr(), W264aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_cooo_no95_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO95_X2_TYPE1_ERI_C)
        (sa1, ia1, sj, ij, sw, iw, W264aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   96] -- 
  // |-- [    0] --| W265(w,c0,k,a4) += (    1.00000000) V2(c0,w,a4,a3) D1(a3,k) 
  // |-- [    1] --| W266(c0,i,j,a4) += (    1.00000000) T2(c0,a0,a2,a1) D3(a0,a1,i,a2,j,a4) 
  // |-- [    2] --| S2(w,k,i,j) += (   -0.50000000) W265(w,c0,k,a4) W266(c0,i,j,a4) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// ==> makeContractions: [1] Special care is taken for the contraction 1 <3, 2> 
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 1
  orz::DTensor W265caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sc0));
  FC_FUNC(g_if_sigma_cooo_cooo_no96_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO96_X0_TYPE1_ERI_C)
    (sc0, ic0, V2_sym.cptr(), W265caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W266aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sj));
    for(int sa1 = 0;sa1 < nir;++sa1){ 
    for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
      T2b = T2.get_amp2(ia1);
      FC_FUNC(g_if_sigma_cooo_cooo_no96_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO96_X1_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, T2b.cptr(), W266aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
    } // End ia1
    } // End sa1
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no96_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO96_X2_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, W265caa_sigma_cooo_cooo.cptr(), W266aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   97] -- 
  // |-- [    0] --| W267(w,a0,a2,a1,a4,a3) += (    1.00000000) V2(w,c0,a4,a3) T2(c0,a0,a2,a1) 
  // |-- [    1] --| W268(w,i,a0,a1,a4,a3) += (    1.00000000) D1(i,a2) W267(w,a0,a2,a1,a4,a3) 
  // |-- [    2] --| S2(w,k,i,j) += (    1.00000000) D3(a0,a1,a3,a4,j,k) W268(w,i,a0,a1,a4,a3) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 0
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W267aaaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sw^sa1));
    orz::DTensor W268aaaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sw^sa1));
    FC_FUNC(g_if_sigma_cooo_cooo_no97_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO97_X0_TYPE1_ERI_C)
      (sa1, ia1, sw, iw, T2b.cptr(), V2_sym.cptr(), W267aaaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no97_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO97_X1_TYPE1_ERI_C)
      (sa1, ia1, sw, iw, W267aaaa_sigma_cooo_cooo.cptr(), W268aaaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_cooo_no97_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO97_X2_TYPE1_ERI_C)
        (sa1, ia1, sj, ij, sw, iw, W268aaaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   98] -- 
  // |-- [    0] --| W269(w,a0,a2,a1,a4,a3) += (    1.00000000) V2(w,c0,a4,a3) T2(c0,a0,a2,a1) 
  // |-- [    1] --| W270(w,i,a0,a2,a4,a3) += (    1.00000000) D1(i,a1) W269(w,a0,a2,a1,a4,a3) 
  // |-- [    2] --| S2(w,k,i,j) += (   -0.50000000) D3(a0,a2,a3,a4,j,k) W270(w,i,a0,a2,a4,a3) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 1 <2, 1> 
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 2
  orz::DTensor W270aaaaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaaaa(symblockinfo, sw));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W269aaaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sw^sa1));
    FC_FUNC(g_if_sigma_cooo_cooo_no98_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO98_X0_TYPE1_ERI_C)
      (sa1, ia1, sw, iw, T2b.cptr(), V2_sym.cptr(), W269aaaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no98_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO98_X1_TYPE1_ERI_C)
      (sa1, ia1, sw, iw, W269aaaa_sigma_cooo_cooo.cptr(), W270aaaaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no98_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO98_X2_TYPE1_ERI_C)
      (sj, ij, sw, iw, W270aaaaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   99] -- 
  // |-- [    0] --| W271(w,c0,i,a3) += (    1.00000000) V2(c0,w,a4,a3) D1(i,a4) 
  // |-- [    1] --| W272(c0,a3,j,k) += (    1.00000000) T2(c0,a0,a2,a1) D3(a0,a1,a3,a2,j,k) 
  // |-- [    2] --| S2(w,k,i,j) += (   -0.50000000) W271(w,c0,i,a3) W272(c0,a3,j,k) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// ==> makeContractions: [1] Special care is taken for the contraction 1 <3, 2> 
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 1
  orz::DTensor W271caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sc0));
  FC_FUNC(g_if_sigma_cooo_cooo_no99_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO99_X0_TYPE1_ERI_C)
    (sc0, ic0, V2_sym.cptr(), W271caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W272aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sj));
    for(int sa1 = 0;sa1 < nir;++sa1){ 
    for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
      T2b = T2.get_amp2(ia1);
      FC_FUNC(g_if_sigma_cooo_cooo_no99_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO99_X1_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, T2b.cptr(), W272aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
    } // End ia1
    } // End sa1
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no99_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO99_X2_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, W271caa_sigma_cooo_cooo.cptr(), W272aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  100] -- 
  // |-- [    0] --| W273(c0,a3,j,a4) += (    1.00000000) T2(c0,a0,a2,a1) D3(a0,a1,a3,a4,j,a2) 
  // |-- [    1] --| W274(w,j) += (    1.00000000) V2(c0,w,a4,a3) W273(c0,a3,j,a4) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
  // Pref: 2
  orz::DTensor W273aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sc0));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no100_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO100_X0_TYPE1_ERI_C)
      (sa1, ia1, sc0, ic0, T2b.cptr(), W273aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  FC_FUNC(g_if_sigma_cooo_cooo_no100_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO100_X1_TYPE1_ERI_C)
    (sc0, ic0, V2_sym.cptr(), W273aaa_sigma_cooo_cooo.cptr(), W274ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  101] -- 
  // |-- [    0] --| W276(w,a4,a3,a2) += (    1.00000000) V2(w,c0,a4,a3) W275(c0,a2) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) D3(i,a2,a3,a4,j,k) W276(w,a4,a3,a2) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W276aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_cooo_no101_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO101_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W275ca_sigma_cooo_cooo.cptr(), W276aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no101_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO101_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W276aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  102] -- 
  // |-- [    0] --| W278(w,a4,a3,a1) += (    1.00000000) V2(w,c0,a4,a3) W277(c0,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (   -0.50000000) D3(i,a1,a3,a4,j,k) W278(w,a4,a3,a1) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W278aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_cooo_no102_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO102_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W277ca_sigma_cooo_cooo.cptr(), W278aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no102_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO102_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W278aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  103] -- 
  // |-- [    0] --| W279(w,c0,a0,a3) += (    1.00000000) V2(w,c0,a4,a3) D1(a0,a4) 
  // |-- [    1] --| W280(w,a2,a1,a3) += (    1.00000000) T2(c0,a0,a2,a1) W279(w,c0,a0,a3) 
  // |-- [    2] --| S2(w,k,i,j) += (   -0.50000000) D3(i,a2,a3,a1,j,k) W280(w,a2,a1,a3) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 0
  orz::DTensor W279caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_cooo_no103_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO103_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W279caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W280aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sa1));
    FC_FUNC(g_if_sigma_cooo_cooo_no103_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO103_X1_TYPE1_ERI_C)
      (sa1, ia1, sw, iw, T2b.cptr(), W279caa_sigma_cooo_cooo.cptr(), W280aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_cooo_no103_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO103_X2_TYPE1_ERI_C)
        (sa1, ia1, sj, ij, sw, iw, W280aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  104] -- 
  // |-- [    0] --| W281(w,a0,a2,a1,a4,a3) += (    1.00000000) V2(w,c0,a4,a3) T2(c0,a0,a2,a1) 
  // |-- [    1] --| W282(w,k,a2,a1,a4,a3) += (    1.00000000) D1(a0,k) W281(w,a0,a2,a1,a4,a3) 
  // |-- [    2] --| S2(w,k,i,j) += (   -0.50000000) D3(i,a2,a3,a4,j,a1) W282(w,k,a2,a1,a4,a3) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 0
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W281aaaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sw^sa1));
    orz::DTensor W282aaaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sw^sa1));
    FC_FUNC(g_if_sigma_cooo_cooo_no104_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO104_X0_TYPE1_ERI_C)
      (sa1, ia1, sw, iw, T2b.cptr(), V2_sym.cptr(), W281aaaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no104_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO104_X1_TYPE1_ERI_C)
      (sa1, ia1, sw, iw, W281aaaa_sigma_cooo_cooo.cptr(), W282aaaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_cooo_no104_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO104_X2_TYPE1_ERI_C)
        (sa1, ia1, sj, ij, sw, iw, W282aaaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  105] -- 
  // |-- [    0] --| W284(w,a4,a3,i) += (    1.00000000) V2(w,c0,a4,a3) W283(c0,i) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) D2(a3,a4,j,k) W284(w,a4,a3,i) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W284aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_cooo_no105_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO105_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W283ca_sigma_cooo_cooo.cptr(), W284aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no105_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO105_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W284aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  106] -- 
  // |-- [    0] --| W286(w,a3,i,a2) += (    1.00000000) V2(w,c0,a4,a3) W285(c0,i,a4,a2) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.50000000) D2(a3,a2,j,k) W286(w,a3,i,a2) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W286aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_cooo_no106_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO106_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W285caaa_sigma_cooo_cooo.cptr(), W286aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no106_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO106_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W286aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  107] -- 
  // |-- [    0] --| W288(c0,i,k,a2) += (    1.00000000) T2(c0,a0,a2,a1) D2(a0,a1,i,k) 
  // |-- [    1] --| W287(w,c0,j,a2) += (    1.00000000) V2(c0,w,a4,a3) D2(a3,a4,j,a2) 
  // |-- [    2] --| S2(w,k,i,j) += (    0.50000000) W287(w,c0,j,a2) W288(c0,i,k,a2) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W288aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sc0));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no107_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO107_X0_TYPE1_ERI_C)
      (sa1, ia1, sc0, ic0, T2b.cptr(), W288aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W287ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sc0^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no107_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO107_X1_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, V2_sym.cptr(), W287ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no107_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO107_X2_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, W287ca_sigma_cooo_cooo.cptr(), W288aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  108] -- 
  // |-- [    0] --| W290(w,a3,i,a1) += (    1.00000000) V2(w,c0,a4,a3) W289(c0,i,a4,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.50000000) D2(a3,a1,j,k) W290(w,a3,i,a1) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W290aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_cooo_no108_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO108_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W289caaa_sigma_cooo_cooo.cptr(), W290aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no108_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO108_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W290aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  109] -- 
  // |-- [    0] --| W292(c0,i,k,a1) += (    1.00000000) T2(c0,a0,a2,a1) D2(a0,k,i,a2) 
  // |-- [    1] --| W291(w,c0,j,a1) += (    1.00000000) V2(c0,w,a4,a3) D2(a3,a4,j,a1) 
  // |-- [    2] --| S2(w,k,i,j) += (    0.50000000) W291(w,c0,j,a1) W292(c0,i,k,a1) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W292aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sa1));
    FC_FUNC(g_if_sigma_cooo_cooo_no109_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO109_X0_TYPE1_ERI_C)
      (sa1, ia1, sc0, ic0, T2b.cptr(), W292aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      orz::DTensor W291c_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, sc0^sj^sa1));
      FC_FUNC(g_if_sigma_cooo_cooo_no109_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO109_X1_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, V2_sym.cptr(), W291c_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
      FC_FUNC(g_if_sigma_cooo_cooo_no109_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO109_X2_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, W291c_sigma_cooo_cooo.cptr(), W292aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  110] -- 
  // |-- [    0] --| W294(w,a4,j,a0) += (    1.00000000) V2(w,c0,a4,a3) W293(c0,a3,j,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -0.33333333) D2(a0,a4,i,k) W294(w,a4,j,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W294aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no110_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO110_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, V2_sym.cptr(), W293caaa_sigma_cooo_cooo.cptr(), W294aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no110_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO110_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W294aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  111] -- 
  // |-- [    0] --| W296(w,a4,j,a0) += (    1.00000000) V2(w,c0,a4,a3) W295(c0,a3,j,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -0.16666667) D2(a0,a4,i,k) W296(w,a4,j,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W296aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no111_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO111_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, V2_sym.cptr(), W295caaa_sigma_cooo_cooo.cptr(), W296aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no111_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO111_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W296aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  112] -- 
  // |-- [    0] --| W298(w,a4,j,a0) += (    1.00000000) V2(w,c0,a4,a3) W297(c0,a3,j,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -0.16666667) D2(a0,k,i,a4) W298(w,a4,j,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W298aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no112_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO112_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, V2_sym.cptr(), W297caaa_sigma_cooo_cooo.cptr(), W298aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no112_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO112_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W298aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  113] -- 
  // |-- [    0] --| W300(w,a4,j,a0) += (    1.00000000) V2(w,c0,a4,a3) W299(c0,a3,j,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -0.33333333) D2(a0,k,i,a4) W300(w,a4,j,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W300aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no113_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO113_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, V2_sym.cptr(), W299caaa_sigma_cooo_cooo.cptr(), W300aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no113_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO113_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W300aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  114] -- 
  // |-- [    0] --| W301(c0,a3,a4,a2) += (    1.00000000) T2(c0,a0,a2,a1) D2(a0,a1,a3,a4) 
  // |-- [    1] --| W302(w,a2) += (    1.00000000) V2(c0,w,a4,a3) W301(c0,a3,a4,a2) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
  // Pref: 2
  orz::DTensor W301aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sc0));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no114_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO114_X0_TYPE1_ERI_C)
      (sa1, ia1, sc0, ic0, T2b.cptr(), W301aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  FC_FUNC(g_if_sigma_cooo_cooo_no114_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO114_X1_TYPE1_ERI_C)
    (sc0, ic0, V2_sym.cptr(), W301aaa_sigma_cooo_cooo.cptr(), W302ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  115] -- 
  // |-- [    0] --| W304(w,a4) += (    1.00000000) V2(w,c0,a4,a3) W303(c0,a3) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.50000000) D2(i,a4,j,k) W304(w,a4) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W304a_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_cooo_no115_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO115_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W303ca_sigma_cooo_cooo.cptr(), W304a_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no115_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO115_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W304a_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  116] -- 
  // |-- [    0] --| W306(w,a4,k,a2) += (    1.00000000) V2(w,c0,a4,a3) W305(c0,a3,k,a2) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.50000000) D2(i,a2,j,a4) W306(w,a4,k,a2) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W306aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_cooo_no116_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO116_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W305caaa_sigma_cooo_cooo.cptr(), W306aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no116_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO116_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W306aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  117] -- 
  // |-- [    0] --| W307(c0,a3,a4,a1) += (    1.00000000) T2(c0,a0,a2,a1) D2(a0,a2,a3,a4) 
  // |-- [    1] --| W308(w,a1) += (    1.00000000) V2(c0,w,a4,a3) W307(c0,a3,a4,a1) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W307aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sa1));
    FC_FUNC(g_if_sigma_cooo_cooo_no117_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO117_X0_TYPE1_ERI_C)
      (sa1, ia1, sc0, ic0, T2b.cptr(), W307aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no117_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO117_X1_TYPE1_ERI_C)
      (sa1, ia1, sc0, ic0, V2_sym.cptr(), W307aa_sigma_cooo_cooo.cptr(), W308ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  118] -- 
  // |-- [    0] --| W309(w,c0,a0,k) += (    1.00000000) V2(c0,w,a4,a3) D2(a0,k,a3,a4) 
  // |-- [    1] --| W310(c0,i,j,a0) += (    1.00000000) T2(c0,a0,a2,a1) D2(i,a2,j,a1) 
  // |-- [    2] --| S2(w,k,i,j) += (    0.50000000) W309(w,c0,a0,k) W310(c0,i,j,a0) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// ==> makeContractions: [1] Special care is taken for the contraction 1 <3, 2> 
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 1
  orz::DTensor W309caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sc0));
  FC_FUNC(g_if_sigma_cooo_cooo_no118_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO118_X0_TYPE1_ERI_C)
    (sc0, ic0, V2_sym.cptr(), W309caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W310aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sj));
    for(int sa1 = 0;sa1 < nir;++sa1){ 
    for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
      T2b = T2.get_amp2(ia1);
      FC_FUNC(g_if_sigma_cooo_cooo_no118_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO118_X1_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, T2b.cptr(), W310aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
    } // End ia1
    } // End sa1
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no118_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO118_X2_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, W309caa_sigma_cooo_cooo.cptr(), W310aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  119] -- 
  // |-- [    0] --| W312(w,a4,k,a1) += (    1.00000000) V2(w,c0,a4,a3) W311(c0,a3,k,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (   -0.33333333) D2(i,a1,j,a4) W312(w,a4,k,a1) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W312aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_cooo_no119_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO119_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W311caaa_sigma_cooo_cooo.cptr(), W312aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no119_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO119_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W312aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  120] -- 
  // |-- [    0] --| W314(w,a4,k,a1) += (    1.00000000) V2(w,c0,a4,a3) W313(c0,a3,k,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (   -0.16666667) D2(i,a4,j,a1) W314(w,a4,k,a1) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W314aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_cooo_no120_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO120_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W313caaa_sigma_cooo_cooo.cptr(), W314aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no120_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO120_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W314aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  121] -- 
  // |-- [    0] --| W316(w,a4,k,a1) += (    1.00000000) V2(w,c0,a4,a3) W315(c0,a3,k,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (   -0.16666667) D2(i,a1,j,a4) W316(w,a4,k,a1) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W316aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_cooo_no121_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO121_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W315caaa_sigma_cooo_cooo.cptr(), W316aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no121_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO121_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W316aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  122] -- 
  // |-- [    0] --| W318(w,a4,k,a1) += (    1.00000000) V2(w,c0,a4,a3) W317(c0,a3,k,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (   -0.33333333) D2(i,a4,j,a1) W318(w,a4,k,a1) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W318aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_cooo_no122_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO122_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W317caaa_sigma_cooo_cooo.cptr(), W318aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no122_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO122_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W318aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  123] -- 
  // |-- [    0] --| W319(w,c0,i,a2) += (    1.00000000) V2(c0,w,a4,a3) D2(a3,a4,i,a2) 
  // |-- [    1] --| W320(c0,j,k,a2) += (    1.00000000) T2(c0,a0,a2,a1) D2(a0,a1,j,k) 
  // |-- [    2] --| S2(w,k,i,j) += (   -1.00000000) W319(w,c0,i,a2) W320(c0,j,k,a2) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// ==> makeContractions: [1] Special care is taken for the contraction 1 <3, 2> 
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 1
  orz::DTensor W319caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sc0));
  FC_FUNC(g_if_sigma_cooo_cooo_no123_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO123_X0_TYPE1_ERI_C)
    (sc0, ic0, V2_sym.cptr(), W319caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W320aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sj));
    for(int sa1 = 0;sa1 < nir;++sa1){ 
    for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
      T2b = T2.get_amp2(ia1);
      FC_FUNC(g_if_sigma_cooo_cooo_no123_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO123_X1_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, T2b.cptr(), W320aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
    } // End ia1
    } // End sa1
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no123_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO123_X2_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, W319caa_sigma_cooo_cooo.cptr(), W320aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  124] -- 
  // |-- [    0] --| W322(w,a3,j,a2) += (    1.00000000) V2(w,c0,a4,a3) W321(c0,j,a4,a2) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.50000000) D2(a3,k,i,a2) W322(w,a3,j,a2) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W322aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no124_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO124_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, V2_sym.cptr(), W321caaa_sigma_cooo_cooo.cptr(), W322aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no124_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO124_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W322aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  125] -- 
  // |-- [    0] --| W324(w,a4,a3,j) += (    1.00000000) V2(w,c0,a4,a3) W323(c0,j) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.50000000) D2(a3,a4,i,k) W324(w,a4,a3,j) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W324aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no125_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO125_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, V2_sym.cptr(), W323ca_sigma_cooo_cooo.cptr(), W324aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no125_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO125_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W324aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  126] -- 
  // |-- [    0] --| W326(w,a4,i,a0) += (    1.00000000) V2(w,c0,a4,a3) W325(c0,a3,i,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.50000000) D2(a0,a4,j,k) W326(w,a4,i,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W326aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_cooo_no126_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO126_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W325caaa_sigma_cooo_cooo.cptr(), W326aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no126_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO126_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W326aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  127] -- 
  // |-- [    0] --| W327(w,c0,i,a1) += (    1.00000000) V2(c0,w,a4,a3) D2(a3,a4,i,a1) 
  // |-- [    1] --| W328(c0,j,k,a1) += (    1.00000000) T2(c0,a0,a2,a1) D2(a0,a2,j,k) 
  // |-- [    2] --| S2(w,k,i,j) += (    0.50000000) W327(w,c0,i,a1) W328(c0,j,k,a1) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W327ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sc0^sa1));
    FC_FUNC(g_if_sigma_cooo_cooo_no127_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO127_X0_TYPE1_ERI_C)
      (sa1, ia1, sc0, ic0, V2_sym.cptr(), W327ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      orz::DTensor W328a_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xa(symblockinfo, sc0^sj^sa1));
      FC_FUNC(g_if_sigma_cooo_cooo_no127_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO127_X1_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, T2b.cptr(), W328a_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
      FC_FUNC(g_if_sigma_cooo_cooo_no127_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO127_X2_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, W327ca_sigma_cooo_cooo.cptr(), W328a_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  128] -- 
  // |-- [    0] --| W330(w,a3,j,a1) += (    1.00000000) V2(w,c0,a4,a3) W329(c0,j,a4,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (   -0.33333333) D2(a3,a1,i,k) W330(w,a3,j,a1) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W330aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no128_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO128_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, V2_sym.cptr(), W329caaa_sigma_cooo_cooo.cptr(), W330aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no128_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO128_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W330aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  129] -- 
  // |-- [    0] --| W332(w,a3,j,a1) += (    1.00000000) V2(w,c0,a4,a3) W331(c0,j,a4,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (   -0.16666667) D2(a3,k,i,a1) W332(w,a3,j,a1) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W332aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no129_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO129_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, V2_sym.cptr(), W331caaa_sigma_cooo_cooo.cptr(), W332aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no129_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO129_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W332aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  130] -- 
  // |-- [    0] --| W334(w,a3,j,a1) += (    1.00000000) V2(w,c0,a4,a3) W333(c0,j,a4,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (   -0.16666667) D2(a3,a1,i,k) W334(w,a3,j,a1) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W334aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no130_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO130_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, V2_sym.cptr(), W333caaa_sigma_cooo_cooo.cptr(), W334aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no130_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO130_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W334aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  131] -- 
  // |-- [    0] --| W336(w,a3,j,a1) += (    1.00000000) V2(w,c0,a4,a3) W335(c0,j,a4,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (   -0.33333333) D2(a3,k,i,a1) W336(w,a3,j,a1) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W336aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no131_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO131_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, V2_sym.cptr(), W335caaa_sigma_cooo_cooo.cptr(), W336aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no131_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO131_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W336aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  132] -- 
  // |-- [    0] --| W338(w,a4,a3,i) += (    1.00000000) V2(w,c0,a4,a3) W337(c0,i) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) C2(a3,a4,j,k) W338(w,a4,a3,i) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W338aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_cooo_no132_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO132_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W337ca_sigma_cooo_cooo.cptr(), W338aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no132_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO132_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W338aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  133] -- 
  // |-- [    0] --| W340(w,a3,i,a2) += (    1.00000000) V2(w,c0,a4,a3) W339(c0,i,a4,a2) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) C2(a2,a3,k,j) W340(w,a3,i,a2) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W340aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_cooo_no133_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO133_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W339caaa_sigma_cooo_cooo.cptr(), W340aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no133_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO133_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W340aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  134] -- 
  // |-- [    0] --| W342(c0,i,k,a2) += (    1.00000000) T2(c0,a0,a2,a1) C2(a0,a1,i,k) 
  // |-- [    1] --| W341(w,c0,j,a2) += (    1.00000000) V2(c0,w,a4,a3) C2(a2,j,a4,a3) 
  // |-- [    2] --| S2(w,k,i,j) += (   -1.00000000) W341(w,c0,j,a2) W342(c0,i,k,a2) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W342aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sc0));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no134_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO134_X0_TYPE1_ERI_C)
      (sa1, ia1, sc0, ic0, T2b.cptr(), W342aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W341ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sc0^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no134_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO134_X1_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, V2_sym.cptr(), W341ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no134_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO134_X2_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, W341ca_sigma_cooo_cooo.cptr(), W342aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  135] -- 
  // |-- [    0] --| W344(w,a3,i,a1) += (    1.00000000) V2(w,c0,a4,a3) W343(c0,i,a4,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) C2(a1,a3,k,j) W344(w,a3,i,a1) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W344aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_cooo_no135_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO135_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W343caaa_sigma_cooo_cooo.cptr(), W344aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no135_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO135_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W344aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  136] -- 
  // |-- [    0] --| W346(c0,i,k,a1) += (    1.00000000) T2(c0,a0,a2,a1) C2(a0,k,i,a2) 
  // |-- [    1] --| W345(w,c0,j,a1) += (    1.00000000) V2(c0,w,a4,a3) C2(a1,j,a4,a3) 
  // |-- [    2] --| S2(w,k,i,j) += (   -1.00000000) W345(w,c0,j,a1) W346(c0,i,k,a1) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 0
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W346aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sa1));
    FC_FUNC(g_if_sigma_cooo_cooo_no136_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO136_X0_TYPE1_ERI_C)
      (sa1, ia1, sc0, ic0, T2b.cptr(), W346aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      orz::DTensor W345c_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, sc0^sj^sa1));
      FC_FUNC(g_if_sigma_cooo_cooo_no136_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO136_X1_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, V2_sym.cptr(), W345c_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_cooo_no136_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO136_X2_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, W345c_sigma_cooo_cooo.cptr(), W346aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  137] -- 
  // |-- [    0] --| W348(w,a4,j,a0) += (    1.00000000) V2(w,c0,a4,a3) W347(c0,a3,j,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.66666667) C2(a0,a4,i,k) W348(w,a4,j,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W348aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no137_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO137_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, V2_sym.cptr(), W347caaa_sigma_cooo_cooo.cptr(), W348aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no137_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO137_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W348aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  138] -- 
  // |-- [    0] --| W350(w,a4,j,a0) += (    1.00000000) V2(w,c0,a4,a3) W349(c0,a3,j,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.33333333) C2(a0,a4,i,k) W350(w,a4,j,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W350aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no138_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO138_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, V2_sym.cptr(), W349caaa_sigma_cooo_cooo.cptr(), W350aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no138_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO138_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W350aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  139] -- 
  // |-- [    0] --| W352(w,a4,j,a0) += (    1.00000000) V2(w,c0,a4,a3) W351(c0,a3,j,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.33333333) C2(a0,k,i,a4) W352(w,a4,j,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W352aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no139_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO139_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, V2_sym.cptr(), W351caaa_sigma_cooo_cooo.cptr(), W352aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no139_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO139_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W352aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  140] -- 
  // |-- [    0] --| W354(w,a4,j,a0) += (    1.00000000) V2(w,c0,a4,a3) W353(c0,a3,j,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.66666667) C2(a0,k,i,a4) W354(w,a4,j,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W354aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no140_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO140_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, V2_sym.cptr(), W353caaa_sigma_cooo_cooo.cptr(), W354aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no140_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO140_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W354aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  141] -- 
  // |-- [    0] --| W355(c0,a3,a4,a2) += (    1.00000000) T2(c0,a0,a2,a1) C2(a0,a1,a3,a4) 
  // |-- [    1] --| W356(w,a2) += (    1.00000000) V2(c0,w,a4,a3) W355(c0,a3,a4,a2) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
  // Pref: 2
  orz::DTensor W355aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sc0));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no141_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO141_X0_TYPE1_ERI_C)
      (sa1, ia1, sc0, ic0, T2b.cptr(), W355aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  FC_FUNC(g_if_sigma_cooo_cooo_no141_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO141_X1_TYPE1_ERI_C)
    (sc0, ic0, V2_sym.cptr(), W355aaa_sigma_cooo_cooo.cptr(), W356ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  142] -- 
  // |-- [    0] --| W358(w,a4) += (    1.00000000) V2(w,c0,a4,a3) W357(c0,a3) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) C2(a4,i,k,j) W358(w,a4) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W358a_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_cooo_no142_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO142_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W357ca_sigma_cooo_cooo.cptr(), W358a_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no142_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO142_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W358a_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  143] -- 
  // |-- [    0] --| W360(w,a4,k,a2) += (    1.00000000) V2(w,c0,a4,a3) W359(c0,a3,k,a2) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) C2(a2,i,a4,j) W360(w,a4,k,a2) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W360aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_cooo_no143_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO143_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W359caaa_sigma_cooo_cooo.cptr(), W360aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no143_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO143_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W360aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  144] -- 
  // |-- [    0] --| W361(c0,a3,a4,a1) += (    1.00000000) T2(c0,a0,a2,a1) C2(a0,a2,a3,a4) 
  // |-- [    1] --| W362(w,a1) += (    1.00000000) V2(c0,w,a4,a3) W361(c0,a3,a4,a1) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W361aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sa1));
    FC_FUNC(g_if_sigma_cooo_cooo_no144_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO144_X0_TYPE1_ERI_C)
      (sa1, ia1, sc0, ic0, T2b.cptr(), W361aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no144_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO144_X1_TYPE1_ERI_C)
      (sa1, ia1, sc0, ic0, V2_sym.cptr(), W361aa_sigma_cooo_cooo.cptr(), W362ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  145] -- 
  // |-- [    0] --| W363(w,c0,a0,k) += (    1.00000000) V2(c0,w,a4,a3) C2(a0,k,a3,a4) 
  // |-- [    1] --| W364(c0,i,j,a0) += (    1.00000000) T2(c0,a0,a2,a1) C2(a1,j,a2,i) 
  // |-- [    2] --| S2(w,k,i,j) += (   -1.00000000) W363(w,c0,a0,k) W364(c0,i,j,a0) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// ==> makeContractions: [1] Special care is taken for the contraction 1 <3, 2> 
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 1
  orz::DTensor W363caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sc0));
  FC_FUNC(g_if_sigma_cooo_cooo_no145_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO145_X0_TYPE1_ERI_C)
    (sc0, ic0, V2_sym.cptr(), W363caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W364aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sj));
    for(int sa1 = 0;sa1 < nir;++sa1){ 
    for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
      T2b = T2.get_amp2(ia1);
      FC_FUNC(g_if_sigma_cooo_cooo_no145_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO145_X1_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, T2b.cptr(), W364aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
    } // End ia1
    } // End sa1
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no145_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO145_X2_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, W363caa_sigma_cooo_cooo.cptr(), W364aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  146] -- 
  // |-- [    0] --| W366(w,a4,k,a1) += (    1.00000000) V2(w,c0,a4,a3) W365(c0,a3,k,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.66666667) C2(a1,i,a4,j) W366(w,a4,k,a1) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W366aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_cooo_no146_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO146_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W365caaa_sigma_cooo_cooo.cptr(), W366aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no146_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO146_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W366aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  147] -- 
  // |-- [    0] --| W368(w,a4,k,a1) += (    1.00000000) V2(w,c0,a4,a3) W367(c0,a3,k,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.33333333) C2(a1,j,a4,i) W368(w,a4,k,a1) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W368aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_cooo_no147_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO147_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W367caaa_sigma_cooo_cooo.cptr(), W368aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no147_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO147_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W368aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  148] -- 
  // |-- [    0] --| W370(w,a4,k,a1) += (    1.00000000) V2(w,c0,a4,a3) W369(c0,a3,k,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.33333333) C2(a1,i,a4,j) W370(w,a4,k,a1) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W370aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_cooo_no148_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO148_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W369caaa_sigma_cooo_cooo.cptr(), W370aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no148_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO148_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W370aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  149] -- 
  // |-- [    0] --| W372(w,a4,k,a1) += (    1.00000000) V2(w,c0,a4,a3) W371(c0,a3,k,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.66666667) C2(a1,j,a4,i) W372(w,a4,k,a1) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W372aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_cooo_no149_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO149_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W371caaa_sigma_cooo_cooo.cptr(), W372aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no149_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO149_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W372aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  150] -- 
  // |-- [    0] --| W373(w,c0,i,a2) += (    1.00000000) V2(c0,w,a4,a3) C2(a2,i,a4,a3) 
  // |-- [    1] --| W374(c0,j,k,a2) += (    1.00000000) T2(c0,a0,a2,a1) C2(a0,a1,j,k) 
  // |-- [    2] --| S2(w,k,i,j) += (    2.00000000) W373(w,c0,i,a2) W374(c0,j,k,a2) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// ==> makeContractions: [1] Special care is taken for the contraction 1 <3, 2> 
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 1
  orz::DTensor W373caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sc0));
  FC_FUNC(g_if_sigma_cooo_cooo_no150_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO150_X0_TYPE1_ERI_C)
    (sc0, ic0, V2_sym.cptr(), W373caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W374aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sj));
    for(int sa1 = 0;sa1 < nir;++sa1){ 
    for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
      T2b = T2.get_amp2(ia1);
      FC_FUNC(g_if_sigma_cooo_cooo_no150_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO150_X1_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, T2b.cptr(), W374aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
    } // End ia1
    } // End sa1
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no150_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO150_X2_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, W373caa_sigma_cooo_cooo.cptr(), W374aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  151] -- 
  // |-- [    0] --| W376(w,a3,j,a2) += (    1.00000000) V2(w,c0,a4,a3) W375(c0,j,a4,a2) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) C2(a2,i,k,a3) W376(w,a3,j,a2) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W376aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no151_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO151_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, V2_sym.cptr(), W375caaa_sigma_cooo_cooo.cptr(), W376aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no151_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO151_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W376aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  152] -- 
  // |-- [    0] --| W378(w,a4,a3,j) += (    1.00000000) V2(w,c0,a4,a3) W377(c0,j) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) C2(a3,a4,i,k) W378(w,a4,a3,j) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W378aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no152_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO152_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, V2_sym.cptr(), W377ca_sigma_cooo_cooo.cptr(), W378aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no152_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO152_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W378aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  153] -- 
  // |-- [    0] --| W380(w,a4,i,a0) += (    1.00000000) V2(w,c0,a4,a3) W379(c0,a3,i,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) C2(a0,a4,j,k) W380(w,a4,i,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W380aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_cooo_no153_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO153_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W379caaa_sigma_cooo_cooo.cptr(), W380aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no153_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO153_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W380aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  154] -- 
  // |-- [    0] --| W381(w,c0,i,a1) += (    1.00000000) V2(c0,w,a4,a3) C2(a1,i,a4,a3) 
  // |-- [    1] --| W382(c0,j,k,a1) += (    1.00000000) T2(c0,a0,a2,a1) C2(a0,a2,j,k) 
  // |-- [    2] --| S2(w,k,i,j) += (   -1.00000000) W381(w,c0,i,a1) W382(c0,j,k,a1) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W381ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sc0^sa1));
    FC_FUNC(g_if_sigma_cooo_cooo_no154_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO154_X0_TYPE1_ERI_C)
      (sa1, ia1, sc0, ic0, V2_sym.cptr(), W381ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      orz::DTensor W382a_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xa(symblockinfo, sc0^sj^sa1));
      FC_FUNC(g_if_sigma_cooo_cooo_no154_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO154_X1_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, T2b.cptr(), W382a_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
      FC_FUNC(g_if_sigma_cooo_cooo_no154_x2_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO154_X2_TYPE1_ERI_C)
        (sa1, ia1, sc0, ic0, sj, ij, W381ca_sigma_cooo_cooo.cptr(), W382a_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  155] -- 
  // |-- [    0] --| W384(w,a3,j,a1) += (    1.00000000) V2(w,c0,a4,a3) W383(c0,j,a4,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.66666667) C2(a1,a3,k,i) W384(w,a3,j,a1) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W384aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no155_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO155_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, V2_sym.cptr(), W383caaa_sigma_cooo_cooo.cptr(), W384aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no155_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO155_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W384aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  156] -- 
  // |-- [    0] --| W386(w,a3,j,a1) += (    1.00000000) V2(w,c0,a4,a3) W385(c0,j,a4,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.33333333) C2(a1,i,k,a3) W386(w,a3,j,a1) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W386aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no156_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO156_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, V2_sym.cptr(), W385caaa_sigma_cooo_cooo.cptr(), W386aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no156_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO156_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W386aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  157] -- 
  // |-- [    0] --| W388(w,a3,j,a1) += (    1.00000000) V2(w,c0,a4,a3) W387(c0,j,a4,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.33333333) C2(a1,a3,k,i) W388(w,a3,j,a1) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W388aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no157_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO157_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, V2_sym.cptr(), W387caaa_sigma_cooo_cooo.cptr(), W388aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no157_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO157_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W388aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [  158] -- 
  // |-- [    0] --| W390(w,a3,j,a1) += (    1.00000000) V2(w,c0,a4,a3) W389(c0,j,a4,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.66666667) C2(a1,i,k,a3) W390(w,a3,j,a1) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W390aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no158_x0_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO158_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, V2_sym.cptr(), W389caaa_sigma_cooo_cooo.cptr(), W390aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no158_x1_type1_eri_c,G_IF_SIGMA_COOO_COOO_NO158_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W390aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- Title : sigma_cooo_cooo
  //*-- Entering to take the type 2 contractions --*//

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    1.00000000) D1(j,k) W252(w,i) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no0_x0_type2_eri_c,G_IF_SIGMA_COOO_COOO_NO0_X0_TYPE2_ERI_C)
      (sj, ij, W252ca_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -0.50000000) D1(i,k) W274(w,j) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no1_x0_type2_eri_c,G_IF_SIGMA_COOO_COOO_NO1_X0_TYPE2_ERI_C)
      (sj, ij, W274ca_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -1.00000000) D2(i,a2,j,k) W302(w,a2) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no2_x0_type2_eri_c,G_IF_SIGMA_COOO_COOO_NO2_X0_TYPE2_ERI_C)
      (sj, ij, W302ca_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    0.50000000) D2(i,a1,j,k) W308(w,a1) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no3_x0_type2_eri_c,G_IF_SIGMA_COOO_COOO_NO3_X0_TYPE2_ERI_C)
      (sj, ij, W308ca_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    2.00000000) C2(a2,i,k,j) W356(w,a2) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no4_x0_type2_eri_c,G_IF_SIGMA_COOO_COOO_NO4_X0_TYPE2_ERI_C)
      (sj, ij, W356ca_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -1.00000000) C2(a1,i,k,j) W362(w,a1) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no5_x0_type2_eri_c,G_IF_SIGMA_COOO_COOO_NO5_X0_TYPE2_ERI_C)
      (sj, ij, W362ca_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
//-@type(2).contraction(end)
  orz::world().barrier();

  } // End Femto
  //*-- FEMTO ends --//*

//-@loadERI(c,end)

//-@loadERI(a,begin)
  //*-- FEMTO begins --//*
  // Label : eri_o
  {

//-@type(2).declaration(begin)
  // --  Title : sigma_cooo_cooo
  //  >> Intermediates for the external contractions are defined here << 
  orz::DTensor W95caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W96caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W99caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W100caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W109caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W110caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W420caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W28caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W30caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W31caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W33caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W52caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W53caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W54caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W56caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W57caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W59caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W78caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W79caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W80caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W81caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W86caaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W112ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W126ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W146ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W150ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W200ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W204ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
//-@type(2).declaration(end)

  //*-- Entering to take the type 0 contractions --*//
//-@type(0).contraction(begin)

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W95(c0,k,a3,i) += (    1.00000000) T2(c0,a0,a2,a1) D3(k,a3,a2,i,a1,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no0_x0_type0_eri_o,G_IF_SIGMA_COOO_COOO_NO0_X0_TYPE0_ERI_O)
      (sa1, ia1, T2b.cptr(), W95caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W96(c0,k,i,a3) += (    1.00000000) T2(c0,a0,a2,a1) D3(k,i,a2,a3,a1,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no1_x0_type0_eri_o,G_IF_SIGMA_COOO_COOO_NO1_X0_TYPE0_ERI_O)
      (sa1, ia1, T2b.cptr(), W96caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W99(c0,k,a2,i) += (    1.00000000) T2(c0,a0,a1,i) D2(k,a2,a1,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    FC_FUNC(g_if_sigma_cooo_cooo_no2_x0_type0_eri_o,G_IF_SIGMA_COOO_COOO_NO2_X0_TYPE0_ERI_O)
      (si, ii, T2b.cptr(), W99caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W100(c0,k,a2,i) += (    1.00000000) T2(c0,a0,a1,i) D2(k,a0,a1,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    FC_FUNC(g_if_sigma_cooo_cooo_no3_x0_type0_eri_o,G_IF_SIGMA_COOO_COOO_NO3_X0_TYPE0_ERI_O)
      (si, ii, T2b.cptr(), W100caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W109(c0,k,a2,i) += (    1.00000000) T2(c0,a0,i,a1) D2(k,a2,a1,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no4_x0_type0_eri_o,G_IF_SIGMA_COOO_COOO_NO4_X0_TYPE0_ERI_O)
      (sa1, ia1, T2b.cptr(), W109caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W110(c0,k,a2,i) += (    1.00000000) T2(c0,a0,i,a1) D2(k,a2,a1,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no5_x0_type0_eri_o,G_IF_SIGMA_COOO_COOO_NO5_X0_TYPE0_ERI_O)
      (sa1, ia1, T2b.cptr(), W110caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W420(w,k,a4,a3) += (    1.00000000) T2(w,a0,a2,a1) D3(k,a4,a2,a3,a1,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no6_x0_type0_eri_o,G_IF_SIGMA_COOO_COOO_NO6_X0_TYPE0_ERI_O)
      (sa1, ia1, T2b.cptr(), W420caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

//-@type(0).contraction(end)

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
  // -- Title : sigma_cooo_cooo

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W28(w,a0,a1,i) += (    1.00000000) V2(a2,c0,w,i) T2(a0,c0,a1,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(ia2);
  FC_FUNC(g_if_sigma_cooo_cooo_no0_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO0_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), V2_sym.cptr(), W28caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W29(k,a1,j,a0,i,a2) += (    1.00000000) V2(a2,a4,i,a3) D3(k,j,a4,a3,a1,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) T2(a0,w,a1,a2) W29(k,a1,j,a0,i,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(ia2);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W29aaaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sj^sa2));
    FC_FUNC(g_if_sigma_cooo_cooo_no1_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO1_X0_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, V2_sym.cptr(), W29aaaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no1_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO1_X1_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, T2b.cptr(), W29aaaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W30(w,a0,i,a2) += (    1.00000000) V2(a1,c0,w,a2) T2(a0,c0,i,a1) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(ia1);
  FC_FUNC(g_if_sigma_cooo_cooo_no2_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO2_X0_TYPE1_ERI_O)
    (sa1, ia1, T2b.cptr(), V2_sym.cptr(), W30caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W31(w,a0,i,a2) += (    1.00000000) V2(a1,a2,w,c0) T2(a0,c0,i,a1) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(ia1);
  FC_FUNC(g_if_sigma_cooo_cooo_no3_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO3_X0_TYPE1_ERI_O)
    (sa1, ia1, T2b.cptr(), V2_sym.cptr(), W31caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W32(k,j,a0,a1) += (    1.00000000) V2(a1,a3,a4,a2) D3(k,j,a4,a2,a3,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -0.50000000) T2(a0,w,i,a1) W32(k,j,a0,a1) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(ia1);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W32aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sa1));
    FC_FUNC(g_if_sigma_cooo_cooo_no4_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO4_X0_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, V2_sym.cptr(), W32aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no4_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO4_X1_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, T2b.cptr(), W32aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W33(w,a0,a1,i) += (    1.00000000) V2(a2,i,w,c0) T2(a0,c0,a1,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(ia2);
  FC_FUNC(g_if_sigma_cooo_cooo_no5_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO5_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), V2_sym.cptr(), W33caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W34(k,a1,j,a0,i,a3) += (    1.00000000) V2(a3,i,a4,a2) D3(k,j,a4,a2,a1,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) T2(a0,w,a1,a3) W34(k,a1,j,a0,i,a3) 
  int sa3(s_eri);
  int ia3(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(ia3);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W34aaaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sj^sa3));
    FC_FUNC(g_if_sigma_cooo_cooo_no6_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO6_X0_TYPE1_ERI_O)
      (sa3, ia3, sj, ij, V2_sym.cptr(), W34aaaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no6_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO6_X1_TYPE1_ERI_O)
      (sa3, ia3, sj, ij, T2b.cptr(), W34aaaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W35(w,a0,a1,j) += (    1.00000000) V2(j,w,c0,a2) T2(c0,a0,a2,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (   -2.00000000) D2(k,i,a1,a0) W35(w,a0,a1,j) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W35ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sa1^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no7_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO7_X0_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, T2b.cptr(), V2_sym.cptr(), W35ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no7_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO7_X1_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, W35ca_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| W36(k,a1,i,a0,j,a2) += (    1.00000000) V2(j,a3,a4,a2) D3(k,a3,a4,i,a1,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) T2(w,a0,a2,a1) W36(k,a1,i,a0,j,a2) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W36aaaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sa1^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no8_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO8_X0_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, V2_sym.cptr(), W36aaaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no8_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO8_X1_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, T2b.cptr(), W36aaaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W39(k,a0,i,a2) += (    1.00000000) V2(a2,a4,a3,a1) D3(k,a0,a4,i,a3,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (   -0.50000000) T2(w,a0,a2,j) W39(k,a0,i,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W39aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa2));
  FC_FUNC(g_if_sigma_cooo_cooo_no9_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO9_X0_TYPE1_ERI_O)
    (sa2, ia2, V2_sym.cptr(), W39aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no9_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO9_X1_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, T2b.cptr(), W39aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| W40(w,a0,a1,j) += (    1.00000000) V2(j,a2,w,c0) T2(c0,a0,a2,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) D2(k,i,a1,a0) W40(w,a0,a1,j) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W40ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sa1^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no10_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO10_X0_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, T2b.cptr(), V2_sym.cptr(), W40ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no10_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO10_X1_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, W40ca_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| W41(k,a1,i,a0,j,a3) += (    1.00000000) V2(j,a3,a4,a2) D3(k,i,a4,a2,a1,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) T2(w,a0,a3,a1) W41(k,a1,i,a0,j,a3) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W41aaaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sa1^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no11_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO11_X0_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, V2_sym.cptr(), W41aaaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no11_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO11_X1_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, T2b.cptr(), W41aaaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| W43(k,a0,i,a1) += (    1.00000000) V2(i,a2,a3,a1) D2(k,a0,a3,a2) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) T2(w,a0,a1,j) W43(k,a0,i,a1) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W43aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, si));
  FC_FUNC(g_if_sigma_cooo_cooo_no12_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO12_X0_TYPE1_ERI_O)
    (si, ii, V2_sym.cptr(), W43aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no12_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO12_X1_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), W43aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| W44(w,k,a2,a3) += (    1.00000000) T2(a0,w,a1,a3) D2(k,a2,a1,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) V2(a3,j,i,a2) W44(w,k,a2,a3) 
  int sa3(s_eri);
  int ia3(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia3);
  orz::DTensor W44caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa3));
  FC_FUNC(g_if_sigma_cooo_cooo_no13_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO13_X0_TYPE1_ERI_O)
    (sa3, ia3, T2b.cptr(), W44caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no13_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO13_X1_TYPE1_ERI_O)
      (sa3, ia3, sj, ij, V2_sym.cptr(), W44caa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   14] -- 
  // |-- [    0] --| W45(w,a0,i,j) += (    1.00000000) V2(j,w,c0,a1) T2(c0,a0,a1,i) 
  // |-- [    1] --| S2(w,k,i,j) += (   -2.00000000) D1(k,a0) W45(w,a0,i,j) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    orz::DTensor W45ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, si^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no14_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO14_X0_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), V2_sym.cptr(), W45ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no14_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO14_X1_TYPE1_ERI_O)
      (si, ii, sj, ij, W45ca_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   15] -- 
  // |-- [    0] --| W46(k,a0,j,a1) += (    1.00000000) V2(j,a2,a3,a1) D2(k,a2,a3,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) T2(w,a0,a1,i) W46(k,a0,j,a1) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W46aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sj));
  FC_FUNC(g_if_sigma_cooo_cooo_no15_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO15_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W46aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    FC_FUNC(g_if_sigma_cooo_cooo_no15_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO15_X1_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), W46aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   16] -- 
  // |-- [    0] --| W47(w,a0,i,j) += (    1.00000000) V2(j,a1,w,c0) T2(c0,a0,a1,i) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) D1(k,a0) W47(w,a0,i,j) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    orz::DTensor W47ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, si^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no16_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO16_X0_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), V2_sym.cptr(), W47ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no16_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO16_X1_TYPE1_ERI_O)
      (si, ii, sj, ij, W47ca_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   17] -- 
  // |-- [    0] --| W48(k,a0,j,a2) += (    1.00000000) V2(j,a2,a3,a1) D2(k,a0,a3,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) T2(w,a0,a2,i) W48(k,a0,j,a2) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W48aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sj));
  FC_FUNC(g_if_sigma_cooo_cooo_no17_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO17_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W48aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    FC_FUNC(g_if_sigma_cooo_cooo_no17_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO17_X1_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), W48aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   18] -- 
  // |-- [    0] --| W49(w,k,a3,a2) += (    1.00000000) T2(a0,w,a1,a2) D2(k,a3,a1,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) V2(a2,i,j,a3) W49(w,k,a3,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia2);
  orz::DTensor W49caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa2));
  FC_FUNC(g_if_sigma_cooo_cooo_no18_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO18_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), W49caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no18_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO18_X1_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, V2_sym.cptr(), W49caa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   19] -- 
  // |-- [    0] --| W51(k,a0,i,a2) += (    1.00000000) V2(i,a2,a3,a1) D2(k,a0,a3,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) T2(w,a0,a2,j) W51(k,a0,i,a2) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W51aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, si));
  FC_FUNC(g_if_sigma_cooo_cooo_no19_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO19_X0_TYPE1_ERI_O)
    (si, ii, V2_sym.cptr(), W51aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no19_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO19_X1_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), W51aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   20] -- 
  // |-- [    0] --| W52(w,a0,a1,a3) += (    1.00000000) V2(a2,c0,w,a3) T2(a0,c0,a1,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(ia2);
  FC_FUNC(g_if_sigma_cooo_cooo_no20_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO20_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), V2_sym.cptr(), W52caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   21] -- 
  // |-- [    0] --| W53(w,a0,a1,a3) += (    1.00000000) V2(a2,a3,w,c0) T2(a0,c0,a1,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(ia2);
  FC_FUNC(g_if_sigma_cooo_cooo_no21_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO21_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), V2_sym.cptr(), W53caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   22] -- 
  // |-- [    0] --| W54(w,a0,a1,i) += (    1.00000000) V2(a2,c0,w,i) T2(c0,a0,a1,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(ia2);
  FC_FUNC(g_if_sigma_cooo_cooo_no22_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO22_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), V2_sym.cptr(), W54caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   23] -- 
  // |-- [    0] --| W55(k,a1,j,a0,i,a2) += (    1.00000000) V2(a2,a4,i,a3) D3(k,j,a4,a0,a1,a3) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) T2(w,a0,a1,a2) W55(k,a1,j,a0,i,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(ia2);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W55aaaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sj^sa2));
    FC_FUNC(g_if_sigma_cooo_cooo_no23_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO23_X0_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, V2_sym.cptr(), W55aaaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no23_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO23_X1_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, T2b.cptr(), W55aaaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   24] -- 
  // |-- [    0] --| W56(w,a0,i,a2) += (    1.00000000) V2(a1,c0,w,a2) T2(c0,a0,i,a1) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(ia1);
  FC_FUNC(g_if_sigma_cooo_cooo_no24_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO24_X0_TYPE1_ERI_O)
    (sa1, ia1, T2b.cptr(), V2_sym.cptr(), W56caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   25] -- 
  // |-- [    0] --| W57(w,a0,i,a2) += (    1.00000000) V2(a1,a2,w,c0) T2(c0,a0,i,a1) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(ia1);
  FC_FUNC(g_if_sigma_cooo_cooo_no25_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO25_X0_TYPE1_ERI_O)
    (sa1, ia1, T2b.cptr(), V2_sym.cptr(), W57caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   26] -- 
  // |-- [    0] --| W58(k,j,a0,a1) += (    1.00000000) V2(a1,a3,a4,a2) D3(k,j,a4,a2,a3,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) T2(w,a0,i,a1) W58(k,j,a0,a1) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(ia1);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W58aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sa1));
    FC_FUNC(g_if_sigma_cooo_cooo_no26_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO26_X0_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, V2_sym.cptr(), W58aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no26_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO26_X1_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, T2b.cptr(), W58aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   27] -- 
  // |-- [    0] --| W59(w,a0,a1,i) += (    1.00000000) V2(a2,i,w,c0) T2(c0,a0,a1,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(ia2);
  FC_FUNC(g_if_sigma_cooo_cooo_no27_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO27_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), V2_sym.cptr(), W59caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   28] -- 
  // |-- [    0] --| W60(k,a1,j,a0,i,a3) += (    1.00000000) V2(a3,i,a4,a2) D3(k,j,a4,a2,a1,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) T2(w,a0,a1,a3) W60(k,a1,j,a0,i,a3) 
  int sa3(s_eri);
  int ia3(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(ia3);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W60aaaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sj^sa3));
    FC_FUNC(g_if_sigma_cooo_cooo_no28_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO28_X0_TYPE1_ERI_O)
      (sa3, ia3, sj, ij, V2_sym.cptr(), W60aaaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no28_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO28_X1_TYPE1_ERI_O)
      (sa3, ia3, sj, ij, T2b.cptr(), W60aaaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   29] -- 
  // |-- [    0] --| W61(w,a0,a1,j) += (    1.00000000) V2(j,w,c0,a2) T2(c0,a0,a1,a2) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) D2(k,i,a1,a0) W61(w,a0,a1,j) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W61caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
  for(int sa2 = 0;sa2 < nir;++sa2){ 
  for(int ia2 = symblockinfo.psym()(sa2,I_O,I_BEGIN);ia2 <= symblockinfo.psym()(sa2,I_O,I_END);++ia2){ 
    T2b = T2.get_amp2(ia2);
    FC_FUNC(g_if_sigma_cooo_cooo_no29_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO29_X0_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, T2b.cptr(), V2_sym.cptr(), W61caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a2, "active"] [notNeeded]
  } // End ia2
  } // End sa2
  FC_FUNC(g_if_sigma_cooo_cooo_no29_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO29_X1_TYPE1_ERI_O)
    (sj, ij, W61caa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   30] -- 
  // |-- [    0] --| W62(k,a1,a0,i,j,a2) += (    1.00000000) V2(j,a3,a4,a2) D3(k,a3,a4,a0,a1,i) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) T2(w,a0,a1,a2) W62(k,a1,a0,i,j,a2) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  for(int sa2 = 0;sa2 < nir;++sa2){ 
  for(int ia2 = symblockinfo.psym()(sa2,I_O,I_BEGIN);ia2 <= symblockinfo.psym()(sa2,I_O,I_END);++ia2){ 
    T2b = T2.get_amp2(ia2);
    orz::DTensor W62aaaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sj^sa2));
    FC_FUNC(g_if_sigma_cooo_cooo_no30_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO30_X0_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, V2_sym.cptr(), W62aaaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no30_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO30_X1_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, T2b.cptr(), W62aaaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a2, "active"] [notNeeded]
  } // End ia2
  } // End sa2
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   31] -- 
  // |-- [    0] --| W65(k,i,a0,a1) += (    1.00000000) V2(a1,a3,a4,a2) D3(k,i,a4,a2,a3,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -0.50000000) T2(w,a0,j,a1) W65(k,i,a0,a1) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia1);
  orz::DTensor W65aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa1));
  FC_FUNC(g_if_sigma_cooo_cooo_no31_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO31_X0_TYPE1_ERI_O)
    (sa1, ia1, V2_sym.cptr(), W65aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no31_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO31_X1_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, T2b.cptr(), W65aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   32] -- 
  // |-- [    0] --| W66(w,a0,a1,j) += (    1.00000000) V2(j,a2,w,c0) T2(c0,a0,a1,a2) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) D2(k,a0,a1,i) W66(w,a0,a1,j) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W66caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
  for(int sa2 = 0;sa2 < nir;++sa2){ 
  for(int ia2 = symblockinfo.psym()(sa2,I_O,I_BEGIN);ia2 <= symblockinfo.psym()(sa2,I_O,I_END);++ia2){ 
    T2b = T2.get_amp2(ia2);
    FC_FUNC(g_if_sigma_cooo_cooo_no32_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO32_X0_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, T2b.cptr(), V2_sym.cptr(), W66caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a2, "active"] [notNeeded]
  } // End ia2
  } // End sa2
  FC_FUNC(g_if_sigma_cooo_cooo_no32_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO32_X1_TYPE1_ERI_O)
    (sj, ij, W66caa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   33] -- 
  // |-- [    0] --| W67(k,a1,a0,i,j,a3) += (    1.00000000) V2(j,a3,a4,a2) D3(k,a0,a4,a2,a1,i) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) T2(w,a0,a1,a3) W67(k,a1,a0,i,j,a3) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  for(int sa3 = 0;sa3 < nir;++sa3){ 
  for(int ia3 = symblockinfo.psym()(sa3,I_O,I_BEGIN);ia3 <= symblockinfo.psym()(sa3,I_O,I_END);++ia3){ 
    T2b = T2.get_amp2(ia3);
    orz::DTensor W67aaaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sj^sa3));
    FC_FUNC(g_if_sigma_cooo_cooo_no33_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO33_X0_TYPE1_ERI_O)
      (sa3, ia3, sj, ij, V2_sym.cptr(), W67aaaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no33_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO33_X1_TYPE1_ERI_O)
      (sa3, ia3, sj, ij, T2b.cptr(), W67aaaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a3, "active"] [notNeeded]
  } // End ia3
  } // End sa3
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   34] -- 
  // |-- [    0] --| W69(k,a0,i,a1) += (    1.00000000) V2(i,a2,a3,a1) D2(k,a2,a3,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) T2(a0,w,a1,j) W69(k,a0,i,a1) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W69aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, si));
  FC_FUNC(g_if_sigma_cooo_cooo_no34_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO34_X0_TYPE1_ERI_O)
    (si, ii, V2_sym.cptr(), W69aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no34_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO34_X1_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), W69aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   35] -- 
  // |-- [    0] --| W70(w,k,a2,a3) += (    1.00000000) T2(w,a0,a1,a3) D2(k,a0,a1,a2) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) V2(a3,j,i,a2) W70(w,k,a2,a3) 
  int sa3(s_eri);
  int ia3(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia3);
  orz::DTensor W70caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa3));
  FC_FUNC(g_if_sigma_cooo_cooo_no35_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO35_X0_TYPE1_ERI_O)
    (sa3, ia3, T2b.cptr(), W70caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no35_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO35_X1_TYPE1_ERI_O)
      (sa3, ia3, sj, ij, V2_sym.cptr(), W70caa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   36] -- 
  // |-- [    0] --| W71(w,a0,i,j) += (    1.00000000) V2(j,w,c0,a1) T2(c0,a0,i,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) D1(k,a0) W71(w,a0,i,j) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W71caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no36_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO36_X0_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, T2b.cptr(), V2_sym.cptr(), W71caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  FC_FUNC(g_if_sigma_cooo_cooo_no36_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO36_X1_TYPE1_ERI_O)
    (sj, ij, W71caa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   37] -- 
  // |-- [    0] --| W72(k,a0,j,a1) += (    1.00000000) V2(j,a2,a3,a1) D2(k,a2,a3,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) T2(w,a0,i,a1) W72(k,a0,j,a1) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W72aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sa1));
    FC_FUNC(g_if_sigma_cooo_cooo_no37_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO37_X0_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, V2_sym.cptr(), W72aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no37_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO37_X1_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, T2b.cptr(), W72aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   38] -- 
  // |-- [    0] --| W73(w,a0,i,j) += (    1.00000000) V2(j,a1,w,c0) T2(c0,a0,i,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (   -2.00000000) D1(k,a0) W73(w,a0,i,j) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W73caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no38_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO38_X0_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, T2b.cptr(), V2_sym.cptr(), W73caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  FC_FUNC(g_if_sigma_cooo_cooo_no38_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO38_X1_TYPE1_ERI_O)
    (sj, ij, W73caa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   39] -- 
  // |-- [    0] --| W74(k,a0,j,a2) += (    1.00000000) V2(j,a2,a3,a1) D2(k,a0,a3,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) T2(w,a0,i,a2) W74(k,a0,j,a2) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  for(int sa2 = 0;sa2 < nir;++sa2){ 
  for(int ia2 = symblockinfo.psym()(sa2,I_O,I_BEGIN);ia2 <= symblockinfo.psym()(sa2,I_O,I_END);++ia2){ 
    T2b = T2.get_amp2(ia2);
    orz::DTensor W74aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sa2));
    FC_FUNC(g_if_sigma_cooo_cooo_no39_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO39_X0_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, V2_sym.cptr(), W74aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no39_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO39_X1_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, T2b.cptr(), W74aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a2, "active"] [notNeeded]
  } // End ia2
  } // End sa2
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   40] -- 
  // |-- [    0] --| W75(w,k,a3,a2) += (    1.00000000) T2(w,a0,a1,a2) D2(k,a3,a1,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) V2(a2,i,j,a3) W75(w,k,a3,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia2);
  orz::DTensor W75caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa2));
  FC_FUNC(g_if_sigma_cooo_cooo_no40_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO40_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), W75caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no40_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO40_X1_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, V2_sym.cptr(), W75caa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   41] -- 
  // |-- [    0] --| W77(k,a0,i,a2) += (    1.00000000) V2(i,a2,a3,a1) D2(k,a0,a3,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) T2(a0,w,a2,j) W77(k,a0,i,a2) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W77aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, si));
  FC_FUNC(g_if_sigma_cooo_cooo_no41_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO41_X0_TYPE1_ERI_O)
    (si, ii, V2_sym.cptr(), W77aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no41_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO41_X1_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), W77aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   42] -- 
  // |-- [    0] --| W78(w,a0,a1,a3) += (    1.00000000) V2(a2,c0,w,a3) T2(c0,a0,a1,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(ia2);
  FC_FUNC(g_if_sigma_cooo_cooo_no42_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO42_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), V2_sym.cptr(), W78caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   43] -- 
  // |-- [    0] --| W79(w,a0,a1,a3) += (    1.00000000) V2(a2,a3,w,c0) T2(c0,a0,a1,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(ia2);
  FC_FUNC(g_if_sigma_cooo_cooo_no43_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO43_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), V2_sym.cptr(), W79caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   44] -- 
  // |-- [    0] --| W80(w,a0,i,a3) += (    1.00000000) V2(a2,i,a3,a1) T2(a0,w,a1,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(ia2);
  FC_FUNC(g_if_sigma_cooo_cooo_no44_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO44_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), V2_sym.cptr(), W80caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   45] -- 
  // |-- [    0] --| W81(w,a0,i,a3) += (    1.00000000) V2(a2,i,a3,a1) T2(w,a0,a1,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(ia2);
  FC_FUNC(g_if_sigma_cooo_cooo_no45_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO45_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), V2_sym.cptr(), W81caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   46] -- 
  // |-- [    0] --| W82(w,a0,j,a3) += (    1.00000000) V2(j,a2,a3,a1) T2(w,a0,a2,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) D2(k,i,a3,a0) W82(w,a0,j,a3) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W82caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no46_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO46_X0_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, T2b.cptr(), V2_sym.cptr(), W82caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  FC_FUNC(g_if_sigma_cooo_cooo_no46_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO46_X1_TYPE1_ERI_O)
    (sj, ij, W82caa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   47] -- 
  // |-- [    0] --| W83(w,a0,j,a3) += (    1.00000000) V2(j,a2,a3,a1) T2(w,a0,a1,a2) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) D2(k,a0,a3,i) W83(w,a0,j,a3) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W83caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
  for(int sa2 = 0;sa2 < nir;++sa2){ 
  for(int ia2 = symblockinfo.psym()(sa2,I_O,I_BEGIN);ia2 <= symblockinfo.psym()(sa2,I_O,I_END);++ia2){ 
    T2b = T2.get_amp2(ia2);
    FC_FUNC(g_if_sigma_cooo_cooo_no47_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO47_X0_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, T2b.cptr(), V2_sym.cptr(), W83caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a2, "active"] [notNeeded]
  } // End ia2
  } // End sa2
  FC_FUNC(g_if_sigma_cooo_cooo_no47_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO47_X1_TYPE1_ERI_O)
    (sj, ij, W83caa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   48] -- 
  // |-- [    0] --| W84(w,a0,j,i) += (    1.00000000) V2(j,a2,i,a1) T2(w,a0,a1,a2) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) D1(k,a0) W84(w,a0,j,i) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W84caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
  for(int sa2 = 0;sa2 < nir;++sa2){ 
  for(int ia2 = symblockinfo.psym()(sa2,I_O,I_BEGIN);ia2 <= symblockinfo.psym()(sa2,I_O,I_END);++ia2){ 
    T2b = T2.get_amp2(ia2);
    FC_FUNC(g_if_sigma_cooo_cooo_no48_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO48_X0_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, T2b.cptr(), V2_sym.cptr(), W84caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a2, "active"] [notNeeded]
  } // End ia2
  } // End sa2
  FC_FUNC(g_if_sigma_cooo_cooo_no48_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO48_X1_TYPE1_ERI_O)
    (sj, ij, W84caa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   49] -- 
  // |-- [    0] --| W85(w,a0,j,i) += (    1.00000000) V2(j,a2,i,a1) T2(w,a0,a2,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) D1(k,a0) W85(w,a0,j,i) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W85caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no49_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO49_X0_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, T2b.cptr(), V2_sym.cptr(), W85caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  FC_FUNC(g_if_sigma_cooo_cooo_no49_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO49_X1_TYPE1_ERI_O)
    (sj, ij, W85caa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   50] -- 
  // |-- [    0] --| W86(w,a0,a4,a3) += (    1.00000000) V2(a2,a4,a3,a1) T2(a0,w,a1,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(ia2);
  FC_FUNC(g_if_sigma_cooo_cooo_no50_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO50_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), V2_sym.cptr(), W86caaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   51] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    1.00000000) V2(j,a3,w,c0) W95(c0,k,a3,i) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  FC_FUNC(g_if_sigma_cooo_cooo_no51_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO51_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W95caaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   52] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    1.00000000) V2(j,w,c0,a3) W96(c0,k,i,a3) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  FC_FUNC(g_if_sigma_cooo_cooo_no52_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO52_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W96caaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   53] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    1.00000000) V2(j,a2,w,c0) W99(c0,k,a2,i) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  FC_FUNC(g_if_sigma_cooo_cooo_no53_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO53_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W99caaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   54] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    1.00000000) V2(j,w,c0,a2) W100(c0,k,a2,i) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  FC_FUNC(g_if_sigma_cooo_cooo_no54_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO54_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W100caaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   55] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -2.00000000) V2(j,a2,w,c0) W109(c0,k,a2,i) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  FC_FUNC(g_if_sigma_cooo_cooo_no55_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO55_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W109caaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   56] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    1.00000000) V2(j,w,c0,a2) W110(c0,k,a2,i) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  FC_FUNC(g_if_sigma_cooo_cooo_no56_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO56_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W110caaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   57] -- 
  // |-- [    0] --| W111(c0,a3,i,a4) += (    1.00000000) T2(c0,a0,a2,a1) D3(a0,a1,a3,a2,i,a4) 
  // |-- [    1] --| W112(w,i) += (    1.00000000) V2(a4,w,c0,a3) W111(c0,a3,i,a4) 
  int sa4(s_eri);
  int ia4(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
  // Pref: 2
  orz::DTensor W111caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa4));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no57_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO57_X0_TYPE1_ERI_O)
      (sa1, ia1, sa4, ia4, T2b.cptr(), W111caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  FC_FUNC(g_if_sigma_cooo_cooo_no57_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO57_X1_TYPE1_ERI_O)
    (sa4, ia4, V2_sym.cptr(), W111caa_sigma_cooo_cooo.cptr(), W112ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   58] -- 
  // |-- [    0] --| W125(c0,a3,j,a4) += (    1.00000000) T2(c0,a0,a2,a1) D3(a0,a1,a3,a2,j,a4) 
  // |-- [    1] --| W126(w,j) += (    1.00000000) V2(a4,w,c0,a3) W125(c0,a3,j,a4) 
  int sa4(s_eri);
  int ia4(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
  // Pref: 2
  orz::DTensor W125caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa4));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no58_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO58_X0_TYPE1_ERI_O)
      (sa1, ia1, sa4, ia4, T2b.cptr(), W125caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  FC_FUNC(g_if_sigma_cooo_cooo_no58_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO58_X1_TYPE1_ERI_O)
    (sa4, ia4, V2_sym.cptr(), W125caa_sigma_cooo_cooo.cptr(), W126ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   59] -- 
  // |-- [    0] --| W133(w,c0,k,a4) += (    1.00000000) V2(a4,w,c0,a3) D1(a3,k) 
  // |-- [    1] --| W134(c0,i,j,a4) += (    1.00000000) T2(c0,a0,a2,a1) D3(a0,a1,i,a4,j,a2) 
  // |-- [    2] --| S2(w,k,i,j) += (   -0.50000000) W133(w,c0,k,a4) W134(c0,i,j,a4) 
  int sa4(s_eri);
  int ia4(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// ==> makeContractions: [1] Special care is taken for the contraction 1 <3, 2> 
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 1
  orz::DTensor W133cca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcca(symblockinfo, sa4));
  FC_FUNC(g_if_sigma_cooo_cooo_no59_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO59_X0_TYPE1_ERI_O)
    (sa4, ia4, V2_sym.cptr(), W133cca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W134ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sa4));
    for(int sa1 = 0;sa1 < nir;++sa1){ 
    for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
      T2b = T2.get_amp2(ia1);
      FC_FUNC(g_if_sigma_cooo_cooo_no59_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO59_X1_TYPE1_ERI_O)
        (sa1, ia1, sa4, ia4, sj, ij, T2b.cptr(), W134ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
    } // End ia1
    } // End sa1
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no59_x2_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO59_X2_TYPE1_ERI_O)
      (sa4, ia4, sj, ij, W133cca_sigma_cooo_cooo.cptr(), W134ca_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   60] -- 
  // |-- [    0] --| W145(c0,a3,a4,a2) += (    1.00000000) T2(c0,a0,a2,a1) D2(a0,a1,a3,a4) 
  // |-- [    1] --| W146(w,a2) += (    1.00000000) V2(a4,w,c0,a3) W145(c0,a3,a4,a2) 
  int sa4(s_eri);
  int ia4(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
  // Pref: 2
  orz::DTensor W145caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa4));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no60_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO60_X0_TYPE1_ERI_O)
      (sa1, ia1, sa4, ia4, T2b.cptr(), W145caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  FC_FUNC(g_if_sigma_cooo_cooo_no60_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO60_X1_TYPE1_ERI_O)
    (sa4, ia4, V2_sym.cptr(), W145caa_sigma_cooo_cooo.cptr(), W146ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   61] -- 
  // |-- [    0] --| W149(c0,a3,a4,a1) += (    1.00000000) T2(c0,a0,a2,a1) D2(a0,a4,a3,a2) 
  // |-- [    1] --| W150(w,a1) += (    1.00000000) V2(a4,w,c0,a3) W149(c0,a3,a4,a1) 
  int sa4(s_eri);
  int ia4(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W149ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sa4^sa1));
    FC_FUNC(g_if_sigma_cooo_cooo_no61_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO61_X0_TYPE1_ERI_O)
      (sa1, ia1, sa4, ia4, T2b.cptr(), W149ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no61_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO61_X1_TYPE1_ERI_O)
      (sa1, ia1, sa4, ia4, V2_sym.cptr(), W149ca_sigma_cooo_cooo.cptr(), W150ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   62] -- 
  // |-- [    0] --| W199(c0,a3,a4,a2) += (    1.00000000) T2(c0,a0,a2,a1) C2(a0,a1,a3,a4) 
  // |-- [    1] --| W200(w,a2) += (    1.00000000) V2(a4,w,c0,a3) W199(c0,a3,a4,a2) 
  int sa4(s_eri);
  int ia4(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
  // Pref: 2
  orz::DTensor W199caa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa4));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no62_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO62_X0_TYPE1_ERI_O)
      (sa1, ia1, sa4, ia4, T2b.cptr(), W199caa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  FC_FUNC(g_if_sigma_cooo_cooo_no62_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO62_X1_TYPE1_ERI_O)
    (sa4, ia4, V2_sym.cptr(), W199caa_sigma_cooo_cooo.cptr(), W200ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   63] -- 
  // |-- [    0] --| W203(c0,a3,a4,a1) += (    1.00000000) T2(c0,a0,a2,a1) C2(a0,a4,a3,a2) 
  // |-- [    1] --| W204(w,a1) += (    1.00000000) V2(a4,w,c0,a3) W203(c0,a3,a4,a1) 
  int sa4(s_eri);
  int ia4(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W203ca_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sa4^sa1));
    FC_FUNC(g_if_sigma_cooo_cooo_no63_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO63_X0_TYPE1_ERI_O)
      (sa1, ia1, sa4, ia4, T2b.cptr(), W203ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no63_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO63_X1_TYPE1_ERI_O)
      (sa1, ia1, sa4, ia4, V2_sym.cptr(), W203ca_sigma_cooo_cooo.cptr(), W204ca_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   64] -- 
  // |-- [    0] --| W391(k,a0,j,a1) += (    1.00000000) V2(a1,a3,a4,a2) D3(k,j,a4,a2,a0,a3) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) T2(w,a1,i,a0) W391(k,a0,j,a1) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      orz::DTensor W391a_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xa(symblockinfo, sa0^sj^sa1));
      FC_FUNC(g_if_sigma_cooo_cooo_no64_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO64_X0_TYPE1_ERI_O)
        (sa0, ia0, sa1, ia1, sj, ij, V2_sym.cptr(), W391a_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_cooo_no64_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO64_X1_TYPE1_ERI_O)
        (sa0, ia0, sa1, ia1, sj, ij, T2b.cptr(), W391a_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   65] -- 
  // |-- [    0] --| W392(k,a1,a0,j,i,a2) += (    1.00000000) V2(i,a3,a4,a2) D3(k,j,a1,a4,a0,a3) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.50000000) T2(w,a2,a0,a1) W392(k,a1,a0,j,i,a2) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      orz::DTensor W392aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa1^sj^si));
      FC_FUNC(g_if_sigma_cooo_cooo_no65_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO65_X0_TYPE1_ERI_O)
        (sa1, ia1, si, ii, sj, ij, V2_sym.cptr(), W392aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_cooo_no65_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO65_X1_TYPE1_ERI_O)
        (sa1, ia1, si, ii, sj, ij, T2b.cptr(), W392aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   66] -- 
  // |-- [    0] --| W393(k,a0,j,a1) += (    1.00000000) V2(a1,a3,a4,a2) D3(k,j,a4,a2,a0,a3) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.50000000) T2(w,a1,a0,i) W393(k,a0,j,a1) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W393aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sa1));
    FC_FUNC(g_if_sigma_cooo_cooo_no66_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO66_X0_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, V2_sym.cptr(), W393aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    for(int si = 0;si < nir;++si){ 
    for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
      T2b = T2.get_amp2(ii);
      FC_FUNC(g_if_sigma_cooo_cooo_no66_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO66_X1_TYPE1_ERI_O)
        (sa1, ia1, si, ii, sj, ij, T2b.cptr(), W393aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    } // End ii
    } // End si
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   67] -- 
  // |-- [    0] --| W394(k,a0,i,a1) += (    1.00000000) V2(a1,a3,a4,a2) D3(k,i,a4,a2,a0,a3) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.50000000) T2(a1,w,a0,j) W394(k,a0,i,a1) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W394aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa1));
  FC_FUNC(g_if_sigma_cooo_cooo_no67_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO67_X0_TYPE1_ERI_O)
    (sa1, ia1, V2_sym.cptr(), W394aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no67_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO67_X1_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, T2b.cptr(), W394aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   68] -- 
  // |-- [    0] --| W395(k,a1,a0,i,j,a2) += (    1.00000000) V2(j,a3,a4,a2) D3(k,a3,a1,i,a0,a4) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.50000000) T2(w,a2,a1,a0) W395(k,a1,a0,i,j,a2) 
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
    orz::DTensor W395aaaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sa0^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no68_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO68_X0_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, V2_sym.cptr(), W395aaaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no68_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO68_X1_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, T2b.cptr(), W395aaaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   69] -- 
  // |-- [    0] --| W396(k,a0,i,a1) += (    1.00000000) V2(a1,a3,a4,a2) D3(k,a3,a4,a2,a0,i) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.50000000) T2(w,a1,a0,j) W396(k,a0,i,a1) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W396aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa1));
  FC_FUNC(g_if_sigma_cooo_cooo_no69_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO69_X0_TYPE1_ERI_O)
    (sa1, ia1, V2_sym.cptr(), W396aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no69_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO69_X1_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, T2b.cptr(), W396aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   70] -- 
  // |-- [    0] --| W397(k,a0,j,a1) += (    1.00000000) V2(j,a2,a3,a1) D2(k,a2,a0,a3) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) T2(w,a1,i,a0) W397(k,a0,j,a1) 
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
    orz::DTensor W397aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa0^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no70_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO70_X0_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, V2_sym.cptr(), W397aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no70_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO70_X1_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, T2b.cptr(), W397aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   71] -- 
  // |-- [    0] --| W398(k,a0) += (    1.00000000) V2(a0,a2,a3,a1) D2(k,a2,a3,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) T2(w,a0,i,j) W398(k,a0) 
  int sa0(s_eri);
  int ia0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W398a_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xa(symblockinfo, sa0));
  FC_FUNC(g_if_sigma_cooo_cooo_no71_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO71_X0_TYPE1_ERI_O)
    (sa0, ia0, V2_sym.cptr(), W398a_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no71_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO71_X1_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, T2b.cptr(), W398a_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   72] -- 
  // |-- [    0] --| W399(k,a0,i,a1) += (    1.00000000) V2(i,a2,a3,a1) D2(k,a2,a0,a3) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.50000000) T2(a1,w,a0,j) W399(k,a0,i,a1) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W399aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, si));
  FC_FUNC(g_if_sigma_cooo_cooo_no72_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO72_X0_TYPE1_ERI_O)
    (si, ii, V2_sym.cptr(), W399aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no72_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO72_X1_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), W399aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   73] -- 
  // |-- [    0] --| W400(k,a0,i,a1) += (    1.00000000) V2(i,a2,a3,a1) D2(k,a3,a0,a2) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.50000000) T2(w,a1,a0,j) W400(k,a0,i,a1) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W400aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, si));
  FC_FUNC(g_if_sigma_cooo_cooo_no73_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO73_X0_TYPE1_ERI_O)
    (si, ii, V2_sym.cptr(), W400aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no73_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO73_X1_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), W400aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   74] -- 
  // |-- [    0] --| W401(k,a0) += (    1.00000000) V2(a0,a2,a3,a1) D2(k,a2,a3,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.50000000) T2(a0,w,i,j) W401(k,a0) 
  int sa0(s_eri);
  int ia0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W401a_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xa(symblockinfo, sa0));
  FC_FUNC(g_if_sigma_cooo_cooo_no74_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO74_X0_TYPE1_ERI_O)
    (sa0, ia0, V2_sym.cptr(), W401a_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no74_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO74_X1_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, T2b.cptr(), W401a_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   75] -- 
  // |-- [    0] --| W402(k,a0,j,a1) += (    1.00000000) V2(j,a2,a3,a1) D2(k,a2,a0,a3) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.50000000) T2(w,a1,a0,i) W402(k,a0,j,a1) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W402aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sj));
  FC_FUNC(g_if_sigma_cooo_cooo_no75_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO75_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W402aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    FC_FUNC(g_if_sigma_cooo_cooo_no75_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO75_X1_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), W402aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   76] -- 
  // |-- [    0] --| W416(k,a1,a0,i) += (    1.00000000) V2(i,a3,a4,a2) D3(k,a3,a4,a2,a1,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -0.50000000) T2(a0,w,a1,j) W416(k,a1,a0,i) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W416aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, si));
  FC_FUNC(g_if_sigma_cooo_cooo_no76_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO76_X0_TYPE1_ERI_O)
    (si, ii, V2_sym.cptr(), W416aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no76_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO76_X1_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), W416aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   77] -- 
  // |-- [    0] --| W417(k,a1,a0,i) += (    1.00000000) V2(i,a3,a4,a2) D3(k,a0,a4,a2,a1,a3) 
  // |-- [    1] --| S2(w,k,i,j) += (   -0.50000000) T2(w,a0,a1,j) W417(k,a1,a0,i) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W417aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, si));
  FC_FUNC(g_if_sigma_cooo_cooo_no77_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO77_X0_TYPE1_ERI_O)
    (si, ii, V2_sym.cptr(), W417aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no77_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO77_X1_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), W417aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   78] -- 
  // |-- [    0] --| W418(k,a1,a0,j) += (    1.00000000) V2(j,a3,a4,a2) D3(k,a3,a4,a2,a1,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) T2(w,a0,i,a1) W418(k,a1,a0,j) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W418aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa1^sj));
    FC_FUNC(g_if_sigma_cooo_cooo_no78_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO78_X0_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, V2_sym.cptr(), W418aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no78_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO78_X1_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, T2b.cptr(), W418aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   79] -- 
  // |-- [    0] --| W419(k,a1,a0,j) += (    1.00000000) V2(j,a3,a4,a2) D3(k,a3,a4,a2,a1,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -0.50000000) T2(w,a0,a1,i) W419(k,a1,a0,j) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W419aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sj));
  FC_FUNC(g_if_sigma_cooo_cooo_no79_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO79_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W419aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    FC_FUNC(g_if_sigma_cooo_cooo_no79_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO79_X1_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), W419aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   80] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -0.50000000) V2(j,a4,i,a3) W420(w,k,a4,a3) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  FC_FUNC(g_if_sigma_cooo_cooo_no80_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO80_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W420caaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   81] -- 
  // |-- [    0] --| W421(j,a1,a0,k) += (    1.00000000) V2(k,a3,a4,a2) D3(j,a3,a4,a2,a1,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) T2(w,a1,i,a0) W421(j,a1,a0,k) 
  int sk(s_eri);
  int ik(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      orz::DTensor W421a_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xa(symblockinfo, sj^sa0^sk));
      FC_FUNC(g_if_sigma_cooo_cooo_no81_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO81_X0_TYPE1_ERI_O)
        (sa0, ia0, sj, ij, sk, ik, V2_sym.cptr(), W421a_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_cooo_no81_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO81_X1_TYPE1_ERI_O)
        (sa0, ia0, sj, ij, sk, ik, T2b.cptr(), W421a_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   82] -- 
  // |-- [    0] --| W422(j,i,a1,a0,k,a2) += (    1.00000000) V2(a2,a4,k,a3) D3(j,a3,i,a4,a1,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.50000000) T2(a1,w,a0,a2) W422(j,i,a1,a0,k,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(ia2);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W422aaaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sj^sa2));
    FC_FUNC(g_if_sigma_cooo_cooo_no82_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO82_X0_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, V2_sym.cptr(), W422aaaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no82_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO82_X1_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, T2b.cptr(), W422aaaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   83] -- 
  // |-- [    0] --| W423(i,a1,a0,k) += (    1.00000000) V2(k,a3,a4,a2) D3(i,a3,a4,a2,a1,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.50000000) T2(a1,w,a0,j) W423(i,a1,a0,k) 
  int sk(s_eri);
  int ik(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W423aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sk));
  FC_FUNC(g_if_sigma_cooo_cooo_no83_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO83_X0_TYPE1_ERI_O)
    (sk, ik, V2_sym.cptr(), W423aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no83_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO83_X1_TYPE1_ERI_O)
      (sj, ij, sk, ik, T2b.cptr(), W423aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   84] -- 
  // |-- [    0] --| W424(j,a1,a0,k) += (    1.00000000) V2(k,a3,a4,a2) D3(j,a3,a4,a2,a1,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.50000000) T2(w,a1,a0,i) W424(j,a1,a0,k) 
  int sk(s_eri);
  int ik(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W424aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sk));
    FC_FUNC(g_if_sigma_cooo_cooo_no84_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO84_X0_TYPE1_ERI_O)
      (sj, ij, sk, ik, V2_sym.cptr(), W424aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    for(int si = 0;si < nir;++si){ 
    for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
      T2b = T2.get_amp2(ii);
      FC_FUNC(g_if_sigma_cooo_cooo_no84_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO84_X1_TYPE1_ERI_O)
        (si, ii, sj, ij, sk, ik, T2b.cptr(), W424aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    } // End ii
    } // End si
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   85] -- 
  // |-- [    0] --| W425(j,i,a1,a0,k,a2) += (    1.00000000) V2(a2,a4,k,a3) D3(j,a3,i,a0,a1,a4) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.50000000) T2(w,a1,a0,a2) W425(j,i,a1,a0,k,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(ia2);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W425aaaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sj^sa2));
    FC_FUNC(g_if_sigma_cooo_cooo_no85_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO85_X0_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, V2_sym.cptr(), W425aaaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no85_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO85_X1_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, T2b.cptr(), W425aaaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   86] -- 
  // |-- [    0] --| W426(i,a1,a0,k) += (    1.00000000) V2(k,a3,a4,a2) D3(i,a0,a4,a2,a1,a3) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.50000000) T2(w,a1,a0,j) W426(i,a1,a0,k) 
  int sk(s_eri);
  int ik(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W426aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sk));
  FC_FUNC(g_if_sigma_cooo_cooo_no86_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO86_X0_TYPE1_ERI_O)
    (sk, ik, V2_sym.cptr(), W426aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no86_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO86_X1_TYPE1_ERI_O)
      (sj, ij, sk, ik, T2b.cptr(), W426aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   87] -- 
  // |-- [    0] --| W427(j,a0,k,a1) += (    1.00000000) V2(a1,a3,k,a2) D2(j,a2,a0,a3) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) T2(w,a0,i,a1) W427(j,a0,k,a1) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(ia1);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W427aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sa1));
    FC_FUNC(g_if_sigma_cooo_cooo_no87_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO87_X0_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, V2_sym.cptr(), W427aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no87_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO87_X1_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, T2b.cptr(), W427aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   88] -- 
  // |-- [    0] --| W428(a0,k) += (    1.00000000) V2(k,a2,a3,a1) D2(a3,a1,a2,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) T2(w,a0,i,j) W428(a0,k) 
  int sk(s_eri);
  int ik(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W428a_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xa(symblockinfo, sk));
  FC_FUNC(g_if_sigma_cooo_cooo_no88_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO88_X0_TYPE1_ERI_O)
    (sk, ik, V2_sym.cptr(), W428a_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no88_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO88_X1_TYPE1_ERI_O)
      (sj, ij, sk, ik, T2b.cptr(), W428a_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   89] -- 
  // |-- [    0] --| W429(j,a0,k,a1) += (    1.00000000) V2(a1,a3,k,a2) D2(j,a2,a0,a3) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.50000000) T2(a0,w,i,a1) W429(j,a0,k,a1) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(ia1);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W429aa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sa1));
    FC_FUNC(g_if_sigma_cooo_cooo_no89_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO89_X0_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, V2_sym.cptr(), W429aa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_cooo_no89_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO89_X1_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, T2b.cptr(), W429aa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   90] -- 
  // |-- [    0] --| W430(i,a0,k,a1) += (    1.00000000) V2(k,a2,a3,a1) D2(i,a3,a0,a2) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.50000000) T2(w,a0,a1,j) W430(i,a0,k,a1) 
  int sk(s_eri);
  int ik(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W430aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sk));
  FC_FUNC(g_if_sigma_cooo_cooo_no90_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO90_X0_TYPE1_ERI_O)
    (sk, ik, V2_sym.cptr(), W430aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no90_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO90_X1_TYPE1_ERI_O)
      (sj, ij, sk, ik, T2b.cptr(), W430aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   91] -- 
  // |-- [    0] --| W431(a0,k) += (    1.00000000) V2(k,a2,a3,a1) D2(a3,a1,a2,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.50000000) T2(a0,w,i,j) W431(a0,k) 
  int sk(s_eri);
  int ik(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W431a_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xa(symblockinfo, sk));
  FC_FUNC(g_if_sigma_cooo_cooo_no91_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO91_X0_TYPE1_ERI_O)
    (sk, ik, V2_sym.cptr(), W431a_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no91_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO91_X1_TYPE1_ERI_O)
      (sj, ij, sk, ik, T2b.cptr(), W431a_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   92] -- 
  // |-- [    0] --| W432(i,a0,k,a1) += (    1.00000000) V2(k,a2,a3,a1) D2(i,a2,a0,a3) 
  // |-- [    1] --| S2(w,k,i,j) += (    0.50000000) T2(a0,w,a1,j) W432(i,a0,k,a1) 
  int sk(s_eri);
  int ik(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W432aaa_sigma_cooo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sk));
  FC_FUNC(g_if_sigma_cooo_cooo_no92_x0_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO92_X0_TYPE1_ERI_O)
    (sk, ik, V2_sym.cptr(), W432aaa_sigma_cooo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no92_x1_type1_eri_o,G_IF_SIGMA_COOO_COOO_NO92_X1_TYPE1_ERI_O)
      (sj, ij, sk, ik, T2b.cptr(), W432aaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- Title : sigma_cooo_cooo
  //*-- Entering to take the type 2 contractions --*//

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    4.00000000) D2(k,j,a1,a0) W28(w,a0,a1,i) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no0_x0_type2_eri_o,G_IF_SIGMA_COOO_COOO_NO0_X0_TYPE2_ERI_O)
      (sj, ij, W28caaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -2.00000000) D2(k,j,a2,a0) W30(w,a0,i,a2) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no1_x0_type2_eri_o,G_IF_SIGMA_COOO_COOO_NO1_X0_TYPE2_ERI_O)
      (sj, ij, W30caaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    1.00000000) D2(k,j,a2,a0) W31(w,a0,i,a2) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no2_x0_type2_eri_o,G_IF_SIGMA_COOO_COOO_NO2_X0_TYPE2_ERI_O)
      (sj, ij, W31caaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -2.00000000) D2(k,j,a1,a0) W33(w,a0,a1,i) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no3_x0_type2_eri_o,G_IF_SIGMA_COOO_COOO_NO3_X0_TYPE2_ERI_O)
      (sj, ij, W33caaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -2.00000000) D3(k,j,a3,i,a1,a0) W52(w,a0,a1,a3) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no4_x0_type2_eri_o,G_IF_SIGMA_COOO_COOO_NO4_X0_TYPE2_ERI_O)
      (sj, ij, W52caaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    1.00000000) D3(k,j,a3,i,a1,a0) W53(w,a0,a1,a3) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no5_x0_type2_eri_o,G_IF_SIGMA_COOO_COOO_NO5_X0_TYPE2_ERI_O)
      (sj, ij, W53caaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -2.00000000) D2(k,j,a1,a0) W54(w,a0,a1,i) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no6_x0_type2_eri_o,G_IF_SIGMA_COOO_COOO_NO6_X0_TYPE2_ERI_O)
      (sj, ij, W54caaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    1.00000000) D2(k,j,a2,a0) W56(w,a0,i,a2) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no7_x0_type2_eri_o,G_IF_SIGMA_COOO_COOO_NO7_X0_TYPE2_ERI_O)
      (sj, ij, W56caaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -2.00000000) D2(k,j,a2,a0) W57(w,a0,i,a2) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no8_x0_type2_eri_o,G_IF_SIGMA_COOO_COOO_NO8_X0_TYPE2_ERI_O)
      (sj, ij, W57caaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    1.00000000) D2(k,j,a1,a0) W59(w,a0,a1,i) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no9_x0_type2_eri_o,G_IF_SIGMA_COOO_COOO_NO9_X0_TYPE2_ERI_O)
      (sj, ij, W59caaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    1.00000000) D3(k,j,a3,i,a1,a0) W78(w,a0,a1,a3) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no10_x0_type2_eri_o,G_IF_SIGMA_COOO_COOO_NO10_X0_TYPE2_ERI_O)
      (sj, ij, W78caaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    1.00000000) D3(k,j,a3,a0,a1,i) W79(w,a0,a1,a3) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no11_x0_type2_eri_o,G_IF_SIGMA_COOO_COOO_NO11_X0_TYPE2_ERI_O)
      (sj, ij, W79caaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    2.00000000) D2(k,j,a3,a0) W80(w,a0,i,a3) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no12_x0_type2_eri_o,G_IF_SIGMA_COOO_COOO_NO12_X0_TYPE2_ERI_O)
      (sj, ij, W80caaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -1.00000000) D2(k,j,a3,a0) W81(w,a0,i,a3) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no13_x0_type2_eri_o,G_IF_SIGMA_COOO_COOO_NO13_X0_TYPE2_ERI_O)
      (sj, ij, W81caaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   14] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -0.50000000) D3(k,j,a4,i,a3,a0) W86(w,a0,a4,a3) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no14_x0_type2_eri_o,G_IF_SIGMA_COOO_COOO_NO14_X0_TYPE2_ERI_O)
      (sj, ij, W86caaa_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   15] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    1.00000000) D1(j,k) W112(w,i) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no15_x0_type2_eri_o,G_IF_SIGMA_COOO_COOO_NO15_X0_TYPE2_ERI_O)
      (sj, ij, W112ca_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   16] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -0.50000000) D1(i,k) W126(w,j) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no16_x0_type2_eri_o,G_IF_SIGMA_COOO_COOO_NO16_X0_TYPE2_ERI_O)
      (sj, ij, W126ca_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   17] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    0.50000000) D2(i,a2,j,k) W146(w,a2) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no17_x0_type2_eri_o,G_IF_SIGMA_COOO_COOO_NO17_X0_TYPE2_ERI_O)
      (sj, ij, W146ca_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   18] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    0.50000000) D2(i,a1,j,k) W150(w,a1) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no18_x0_type2_eri_o,G_IF_SIGMA_COOO_COOO_NO18_X0_TYPE2_ERI_O)
      (sj, ij, W150ca_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   19] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -1.00000000) C2(a2,i,k,j) W200(w,a2) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no19_x0_type2_eri_o,G_IF_SIGMA_COOO_COOO_NO19_X0_TYPE2_ERI_O)
      (sj, ij, W200ca_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   20] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -1.00000000) C2(a1,i,k,j) W204(w,a1) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no20_x0_type2_eri_o,G_IF_SIGMA_COOO_COOO_NO20_X0_TYPE2_ERI_O)
      (sj, ij, W204ca_sigma_cooo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
//-@type(2).contraction(end)
  orz::world().barrier();

  } // End Femto
  //*-- FEMTO ends --//*

//-@loadERI(a,end)

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
  // -- Title : sigma_cooo_cooo

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -0.50000000) T2(w,a0,a3,a1) C5(a0,a1,j,k,i,a3) 
  int sa3(s_eri);
  int ia3(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_cooo_no0_x0_type1_d4c_o,G_IF_SIGMA_COOO_COOO_NO0_X0_TYPE1_D4C_O)
        (sa1, ia1, sa3, ia3, sj, ij, C5.cptr(), T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -0.50000000) T2(w,a0,a1,a2) C5(i,a1,j,k,a0,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(ia2);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_cooo_no1_x0_type1_d4c_o,G_IF_SIGMA_COOO_COOO_NO1_X0_TYPE1_D4C_O)
      (sa2, ia2, sj, ij, C5.cptr(), T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    0.50000000) T2(w,a2,a1,a0) C5(k,j,a1,i,a0,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_cooo_no2_x0_type1_d4c_o,G_IF_SIGMA_COOO_COOO_NO2_X0_TYPE1_D4C_O)
        (sa0, ia0, sa2, ia2, sj, ij, C5.cptr(), T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -0.50000000) T2(w,a0,a2,a1) C5(a1,a0,k,j,a2,i) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_cooo_no3_x0_type1_d4c_o,G_IF_SIGMA_COOO_COOO_NO3_X0_TYPE1_D4C_O)
        (sa1, ia1, si, ii, sj, ij, C5.cptr(), T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -0.50000000) T2(w,a0,a2,a1) C5(a1,a0,a2,i,k,j) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_cooo_no4_x0_type1_d4c_o,G_IF_SIGMA_COOO_COOO_NO4_X0_TYPE1_D4C_O)
      (sa1, ia1, sj, ij, C5.cptr(), T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    0.50000000) T2(w,a2,a1,a0) C5(a2,a0,i,a1,j,k) 
  int sk(s_eri);
  int ik(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_cooo_no5_x0_type1_d4c_o,G_IF_SIGMA_COOO_COOO_NO5_X0_TYPE1_D4C_O)
        (sa0, ia0, sj, ij, sk, ik, C5.cptr(), T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
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
  orz::mr::Femto::my_timer.push_back(boost::make_tuple("sigma_cooo_cooo", time_sigma.elapsed_cputime(), time_sigma.elapsed_wallclocktime()));
  orz::mr::Femto::file_timing << "* " << boost::format("%20s : %10.7f %10.7f ") % orz::mr::Femto::my_timer.back().get<0>() % orz::mr::Femto::my_timer.back().get<1>() % orz::mr::Femto::my_timer.back().get<2>() << endl;
  flush(orz::mr::Femto::file_timing);

  return retval; 
} 
