                                                                                
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
#include <sci/icmr/Femto/elems/c_sigma_cooo_ccoo.h>                                  
                                                                                
using std::cout;                                                                
using std::endl;                                                                
                                                                                
#define FLOPCOUNT                                                               
                                                                                
using orz::mr::Femto::Fc0;                                                      
using orz::mr::Femto::h1_int;                                                   
using orz::mr::Femto::h6_int;                                                   
                                                                                
//  ___________                __               
//  \_   _____/____    _____ _/  |_  ____      
//   |    __)_/ __ \  /     \\   __\/  _ \ 
//   |     \ \  ___/ |  Y Y  \|  | (  <_> )  
//   \___  /  \___  >|__|_|  /|__|  \____/   
//       \/       \/       \/                

//                                   Generated date : Sun Apr 20 10:26:17 2014

                                                                                
// ***************************************************************************  
// orz::mr::mrci                                                                
// ***************************************************************************  
									    /*!        
   @brief CT input                                                              
									     */        
                                                                                
orz::mr::BareAmpPack orz::mr::Femto::sigma_cooo_ccoo(const orz::mr::Input &ctinp,                                    
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
  // |-- [    0] --| W0(w,a2,a1,a0) += (    1.00000000) T2(w,c0,a1,a0) Fc1(c0,a2) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) D3(k,j,a1,i,a0,a2) W0(w,a2,a1,a0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W0caa_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa0));
    FC_FUNC(g_if_sigma_cooo_ccoo_no0_x0_type0_noeri,G_IF_SIGMA_COOO_CCOO_NO0_X0_TYPE0_NOERI)
      (sa0, ia0, T2b.cptr(), W0caa_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_ccoo_no0_x1_type0_noeri,G_IF_SIGMA_COOO_CCOO_NO0_X1_TYPE0_NOERI)
        (sa0, ia0, sj, ij, W0caa_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- No. [    1] -- 
  // |-- [    0] --| W1(w,a1,i,a0) += (    1.00000000) T2(w,c0,i,a0) Fc1(c0,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (   -4.00000000) D2(k,j,a0,a1) W1(w,a1,i,a0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W1caa_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa0));
    FC_FUNC(g_if_sigma_cooo_ccoo_no1_x0_type0_noeri,G_IF_SIGMA_COOO_CCOO_NO1_X0_TYPE0_NOERI)
      (sa0, ia0, T2b.cptr(), W1caa_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_ccoo_no1_x1_type0_noeri,G_IF_SIGMA_COOO_CCOO_NO1_X1_TYPE0_NOERI)
        (sa0, ia0, sj, ij, W1caa_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // |-- [    0] --| W2(w,a1,a0,i) += (    1.00000000) T2(w,c0,a0,i) Fc1(c0,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) D2(k,j,a0,a1) W2(w,a1,a0,i) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    orz::DTensor W2caa_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, si));
    FC_FUNC(g_if_sigma_cooo_ccoo_no2_x0_type0_noeri,G_IF_SIGMA_COOO_CCOO_NO2_X0_TYPE0_NOERI)
      (si, ii, T2b.cptr(), W2caa_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_ccoo_no2_x1_type0_noeri,G_IF_SIGMA_COOO_CCOO_NO2_X1_TYPE0_NOERI)
        (si, ii, sj, ij, W2caa_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- No. [    3] -- 
  // |-- [    0] --| W3(w,a1,j,a0) += (    1.00000000) T2(c0,w,a0,j) Fc1(c0,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) D2(k,i,a0,a1) W3(w,a1,j,a0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W3caa_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccoo_no3_x0_type0_noeri,G_IF_SIGMA_COOO_CCOO_NO3_X0_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W3caa_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_ccoo_no3_x1_type0_noeri,G_IF_SIGMA_COOO_CCOO_NO3_X1_TYPE0_NOERI)
      (sj, ij, W3caa_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W4(w,a1,a0,j) += (    1.00000000) T2(w,c0,a0,j) Fc1(c0,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) D2(k,a1,a0,i) W4(w,a1,a0,j) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W4caa_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccoo_no4_x0_type0_noeri,G_IF_SIGMA_COOO_CCOO_NO4_X0_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W4caa_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_ccoo_no4_x1_type0_noeri,G_IF_SIGMA_COOO_CCOO_NO4_X1_TYPE0_NOERI)
      (sj, ij, W4caa_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W5(c0,k) += (    1.00000000) D1(k,a0) Fc1(c0,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) T2(c0,w,i,j) W5(c0,k) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W5ca_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  FC_FUNC(g_if_sigma_cooo_ccoo_no5_x0_type0_noeri,G_IF_SIGMA_COOO_CCOO_NO5_X0_TYPE0_NOERI)
    (W5ca_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccoo_no5_x1_type0_noeri,G_IF_SIGMA_COOO_CCOO_NO5_X1_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W5ca_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W6(c0,k) += (    1.00000000) D1(k,a0) Fc1(c0,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -4.00000000) T2(w,c0,i,j) W6(c0,k) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W6ca_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  FC_FUNC(g_if_sigma_cooo_ccoo_no6_x0_type0_noeri,G_IF_SIGMA_COOO_CCOO_NO6_X0_TYPE0_NOERI)
    (W6ca_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccoo_no6_x1_type0_noeri,G_IF_SIGMA_COOO_CCOO_NO6_X1_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W6ca_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W7(w,a0) += (    1.00000000) T2(w,c0,a1,a0) Fc1(c0,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) D2(k,j,a0,i) W7(w,a0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W7c_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xc(symblockinfo, sa0));
    FC_FUNC(g_if_sigma_cooo_ccoo_no7_x0_type0_noeri,G_IF_SIGMA_COOO_CCOO_NO7_X0_TYPE0_NOERI)
      (sa0, ia0, T2b.cptr(), W7c_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_ccoo_no7_x1_type0_noeri,G_IF_SIGMA_COOO_CCOO_NO7_X1_TYPE0_NOERI)
        (sa0, ia0, sj, ij, W7c_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- No. [    8] -- 
  // |-- [    0] --| W8(w,a0) += (    1.00000000) T2(w,c0,a0,a1) Fc1(c0,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (   -4.00000000) D2(k,j,a0,i) W8(w,a0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W8ca_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_ccoo_no8_x0_type0_noeri,G_IF_SIGMA_COOO_CCOO_NO8_X0_TYPE0_NOERI)
      (sa1, ia1, T2b.cptr(), W8ca_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccoo_no8_x1_type0_noeri,G_IF_SIGMA_COOO_CCOO_NO8_X1_TYPE0_NOERI)
      (sj, ij, W8ca_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W9(w,i) += (    1.00000000) T2(w,c0,i,a0) Fc1(c0,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    8.00000000) D1(k,j) W9(w,i) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W9ca_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_cooo_ccoo_no9_x0_type0_noeri,G_IF_SIGMA_COOO_CCOO_NO9_X0_TYPE0_NOERI)
      (sa0, ia0, T2b.cptr(), W9ca_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccoo_no9_x1_type0_noeri,G_IF_SIGMA_COOO_CCOO_NO9_X1_TYPE0_NOERI)
      (sj, ij, W9ca_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| W10(w,i) += (    1.00000000) T2(w,c0,a0,i) Fc1(c0,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -4.00000000) D1(k,j) W10(w,i) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    orz::DTensor W10c_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xc(symblockinfo, si));
    FC_FUNC(g_if_sigma_cooo_ccoo_no10_x0_type0_noeri,G_IF_SIGMA_COOO_CCOO_NO10_X0_TYPE0_NOERI)
      (si, ii, T2b.cptr(), W10c_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_ccoo_no10_x1_type0_noeri,G_IF_SIGMA_COOO_CCOO_NO10_X1_TYPE0_NOERI)
        (si, ii, sj, ij, W10c_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // |-- [    0] --| W11(w,j) += (    1.00000000) T2(c0,w,a0,j) Fc1(c0,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -4.00000000) D1(k,i) W11(w,j) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W11c_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xc(symblockinfo, sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccoo_no11_x0_type0_noeri,G_IF_SIGMA_COOO_CCOO_NO11_X0_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W11c_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_ccoo_no11_x1_type0_noeri,G_IF_SIGMA_COOO_CCOO_NO11_X1_TYPE0_NOERI)
      (sj, ij, W11c_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| W12(w,j) += (    1.00000000) T2(w,c0,a0,j) Fc1(c0,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) D1(k,i) W12(w,j) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W12c_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xc(symblockinfo, sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccoo_no12_x0_type0_noeri,G_IF_SIGMA_COOO_CCOO_NO12_X0_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W12c_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_ccoo_no12_x1_type0_noeri,G_IF_SIGMA_COOO_CCOO_NO12_X1_TYPE0_NOERI)
      (sj, ij, W12c_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- Title : sigma_cooo_ccoo

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W13(w,a0,a1,a2) += (    1.00000000) V2(w,c0,c1,a2) T2(c1,c0,a0,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (   -2.00000000) D3(k,j,a1,i,a0,a2) W13(w,a0,a1,a2) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W13aa_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sa1));
    FC_FUNC(g_if_sigma_cooo_ccoo_no0_x0_type1_eri_c,G_IF_SIGMA_COOO_CCOO_NO0_X0_TYPE1_ERI_C)
      (sa1, ia1, sw, iw, T2b.cptr(), V2_sym.cptr(), W13aa_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_ccoo_no0_x1_type1_eri_c,G_IF_SIGMA_COOO_CCOO_NO0_X1_TYPE1_ERI_C)
        (sa1, ia1, sj, ij, sw, iw, W13aa_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // |-- [    0] --| W15(w,a0,i,a1) += (    1.00000000) V2(w,c0,c1,a1) T2(c1,c0,a0,i) 
  // |-- [    1] --| S2(w,k,i,j) += (    4.00000000) D2(k,j,a0,a1) W15(w,a0,i,a1) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    orz::DTensor W15aa_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^si));
    FC_FUNC(g_if_sigma_cooo_ccoo_no1_x0_type1_eri_c,G_IF_SIGMA_COOO_CCOO_NO1_X0_TYPE1_ERI_C)
      (si, ii, sw, iw, T2b.cptr(), V2_sym.cptr(), W15aa_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_ccoo_no1_x1_type1_eri_c,G_IF_SIGMA_COOO_CCOO_NO1_X1_TYPE1_ERI_C)
        (si, ii, sj, ij, sw, iw, W15aa_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W16(w,i,a0,a1) += (    1.00000000) V2(w,c0,c1,a1) T2(c1,c0,i,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -2.00000000) D2(k,j,a0,a1) W16(w,i,a0,a1) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W16aa_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sa0));
    FC_FUNC(g_if_sigma_cooo_ccoo_no2_x0_type1_eri_c,G_IF_SIGMA_COOO_CCOO_NO2_X0_TYPE1_ERI_C)
      (sa0, ia0, sw, iw, T2b.cptr(), V2_sym.cptr(), W16aa_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_ccoo_no2_x1_type1_eri_c,G_IF_SIGMA_COOO_CCOO_NO2_X1_TYPE1_ERI_C)
        (sa0, ia0, sj, ij, sw, iw, W16aa_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // |-- [    0] --| W17(c0,k,a0,j) += (    1.00000000) V2(c0,a2,a3,a1) D3(k,j,a3,a1,a0,a2) 
  // |-- [    1] --| S2(w,k,i,j) += (   -4.00000000) T2(w,c0,i,a0) W17(c0,k,a0,j) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      orz::DTensor W17a_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xa(symblockinfo, sc0^sa0^sj));
      FC_FUNC(g_if_sigma_cooo_ccoo_no3_x0_type1_eri_c,G_IF_SIGMA_COOO_CCOO_NO3_X0_TYPE1_ERI_C)
        (sa0, ia0, sc0, ic0, sj, ij, V2_sym.cptr(), W17a_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_ccoo_no3_x1_type1_eri_c,G_IF_SIGMA_COOO_CCOO_NO3_X1_TYPE1_ERI_C)
        (sa0, ia0, sc0, ic0, sj, ij, T2b.cptr(), W17a_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W18(c0,k,a0,j) += (    1.00000000) V2(c0,a2,a3,a1) D3(k,j,a3,a1,a0,a2) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) T2(w,c0,a0,i) W18(c0,k,a0,j) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W18aa_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sj));
    FC_FUNC(g_if_sigma_cooo_ccoo_no4_x0_type1_eri_c,G_IF_SIGMA_COOO_CCOO_NO4_X0_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, V2_sym.cptr(), W18aa_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    for(int si = 0;si < nir;++si){ 
    for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
      T2b = T2.get_amp2(ii);
      FC_FUNC(g_if_sigma_cooo_ccoo_no4_x1_type1_eri_c,G_IF_SIGMA_COOO_CCOO_NO4_X1_TYPE1_ERI_C)
        (sc0, ic0, si, ii, sj, ij, T2b.cptr(), W18aa_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    } // End ii
    } // End si
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W20(w,a0,j,a1) += (    1.00000000) V2(w,c0,c1,a1) T2(c1,c0,a0,j) 
  // |-- [    1] --| S2(w,k,i,j) += (   -2.00000000) D2(k,i,a0,a1) W20(w,a0,j,a1) 
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
    orz::DTensor W20aa_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccoo_no5_x0_type1_eri_c,G_IF_SIGMA_COOO_CCOO_NO5_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), V2_sym.cptr(), W20aa_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_ccoo_no5_x1_type1_eri_c,G_IF_SIGMA_COOO_CCOO_NO5_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W20aa_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W21(w,j,a0,a1) += (    1.00000000) V2(w,c0,c1,a1) T2(c0,c1,a0,j) 
  // |-- [    1] --| S2(w,k,i,j) += (   -2.00000000) D2(k,a1,a0,i) W21(w,j,a0,a1) 
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
    orz::DTensor W21aa_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccoo_no6_x0_type1_eri_c,G_IF_SIGMA_COOO_CCOO_NO6_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), V2_sym.cptr(), W21aa_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_ccoo_no6_x1_type1_eri_c,G_IF_SIGMA_COOO_CCOO_NO6_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W21aa_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W22(c0,k,a0,i) += (    1.00000000) V2(c0,a2,a3,a1) D3(k,i,a3,a1,a0,a2) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) T2(c0,w,a0,j) W22(c0,k,a0,i) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W22aaa_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sc0));
  FC_FUNC(g_if_sigma_cooo_ccoo_no7_x0_type1_eri_c,G_IF_SIGMA_COOO_CCOO_NO7_X0_TYPE1_ERI_C)
    (sc0, ic0, V2_sym.cptr(), W22aaa_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccoo_no7_x1_type1_eri_c,G_IF_SIGMA_COOO_CCOO_NO7_X1_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, T2b.cptr(), W22aaa_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| W23(c0,k,a0,i) += (    1.00000000) V2(c0,a2,a3,a1) D3(k,a2,a3,a1,a0,i) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) T2(w,c0,a0,j) W23(c0,k,a0,i) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W23aaa_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sc0));
  FC_FUNC(g_if_sigma_cooo_ccoo_no8_x0_type1_eri_c,G_IF_SIGMA_COOO_CCOO_NO8_X0_TYPE1_ERI_C)
    (sc0, ic0, V2_sym.cptr(), W23aaa_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccoo_no8_x1_type1_eri_c,G_IF_SIGMA_COOO_CCOO_NO8_X1_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, T2b.cptr(), W23aaa_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W28(w,c1,c0,k) += (    1.00000000) V2(w,c0,c1,a0) D1(k,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -2.00000000) T2(c1,c0,i,j) W28(w,c1,c0,k) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W28cca_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xcca(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_ccoo_no9_x0_type1_eri_c,G_IF_SIGMA_COOO_CCOO_NO9_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W28cca_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccoo_no9_x1_type1_eri_c,G_IF_SIGMA_COOO_CCOO_NO9_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), W28cca_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| W29(w,c1,c0,k) += (    1.00000000) V2(w,c0,c1,a0) D1(k,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    4.00000000) T2(c0,c1,i,j) W29(w,c1,c0,k) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W29cca_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xcca(symblockinfo, sw));
  FC_FUNC(g_if_sigma_cooo_ccoo_no10_x0_type1_eri_c,G_IF_SIGMA_COOO_CCOO_NO10_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W29cca_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccoo_no10_x1_type1_eri_c,G_IF_SIGMA_COOO_CCOO_NO10_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), W29cca_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| W30(c0,k) += (    1.00000000) V2(c0,a1,a2,a0) D2(k,a1,a2,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -4.00000000) T2(w,c0,i,j) W30(c0,k) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W30a_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xa(symblockinfo, sc0));
  FC_FUNC(g_if_sigma_cooo_ccoo_no11_x0_type1_eri_c,G_IF_SIGMA_COOO_CCOO_NO11_X0_TYPE1_ERI_C)
    (sc0, ic0, V2_sym.cptr(), W30a_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccoo_no11_x1_type1_eri_c,G_IF_SIGMA_COOO_CCOO_NO11_X1_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, T2b.cptr(), W30a_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| W31(c0,k) += (    1.00000000) V2(c0,a1,a2,a0) D2(k,a1,a2,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) T2(c0,w,i,j) W31(c0,k) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W31a_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xa(symblockinfo, sc0));
  FC_FUNC(g_if_sigma_cooo_ccoo_no12_x0_type1_eri_c,G_IF_SIGMA_COOO_CCOO_NO12_X0_TYPE1_ERI_C)
    (sc0, ic0, V2_sym.cptr(), W31a_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccoo_no12_x1_type1_eri_c,G_IF_SIGMA_COOO_CCOO_NO12_X1_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, T2b.cptr(), W31a_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| W50(w,j) += (    1.00000000) V2(w,c0,c1,a0) T2(c1,c0,a0,j) 
  // |-- [    1] --| S2(w,k,i,j) += (    4.00000000) D1(k,i) W50(w,j) 
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
    double W50_sigma_cooo_ccoo(0);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccoo_no13_x0_type1_eri_c,G_IF_SIGMA_COOO_CCOO_NO13_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), V2_sym.cptr(), &W50_sigma_cooo_ccoo, nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_ccoo_no13_x1_type1_eri_c,G_IF_SIGMA_COOO_CCOO_NO13_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, &W50_sigma_cooo_ccoo, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   14] -- 
  // |-- [    0] --| W51(w,j) += (    1.00000000) V2(w,c0,c1,a0) T2(c0,c1,a0,j) 
  // |-- [    1] --| S2(w,k,i,j) += (   -2.00000000) D1(k,i) W51(w,j) 
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
    double W51_sigma_cooo_ccoo(0);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccoo_no14_x0_type1_eri_c,G_IF_SIGMA_COOO_CCOO_NO14_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), V2_sym.cptr(), &W51_sigma_cooo_ccoo, nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_ccoo_no14_x1_type1_eri_c,G_IF_SIGMA_COOO_CCOO_NO14_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, &W51_sigma_cooo_ccoo, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   15] -- 
  // |-- [    0] --| W52(c0,k,i,a0) += (    1.00000000) V2(c0,a1,a2,a0) D2(k,i,a2,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) T2(c0,w,a0,j) W52(c0,k,i,a0) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W52aaa_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sc0));
  FC_FUNC(g_if_sigma_cooo_ccoo_no15_x0_type1_eri_c,G_IF_SIGMA_COOO_CCOO_NO15_X0_TYPE1_ERI_C)
    (sc0, ic0, V2_sym.cptr(), W52aaa_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccoo_no15_x1_type1_eri_c,G_IF_SIGMA_COOO_CCOO_NO15_X1_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, T2b.cptr(), W52aaa_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   16] -- 
  // |-- [    0] --| W53(c0,k,i,a0) += (    1.00000000) V2(c0,a1,a2,a0) D2(k,a1,a2,i) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) T2(w,c0,a0,j) W53(c0,k,i,a0) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W53aaa_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sc0));
  FC_FUNC(g_if_sigma_cooo_ccoo_no16_x0_type1_eri_c,G_IF_SIGMA_COOO_CCOO_NO16_X0_TYPE1_ERI_C)
    (sc0, ic0, V2_sym.cptr(), W53aaa_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccoo_no16_x1_type1_eri_c,G_IF_SIGMA_COOO_CCOO_NO16_X1_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, T2b.cptr(), W53aaa_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   17] -- 
  // |-- [    0] --| W54(c0,k,i,a1) += (    1.00000000) V2(c0,a1,a2,a0) D2(k,i,a2,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -4.00000000) T2(c0,w,a1,j) W54(c0,k,i,a1) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W54aaa_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sc0));
  FC_FUNC(g_if_sigma_cooo_ccoo_no17_x0_type1_eri_c,G_IF_SIGMA_COOO_CCOO_NO17_X0_TYPE1_ERI_C)
    (sc0, ic0, V2_sym.cptr(), W54aaa_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccoo_no17_x1_type1_eri_c,G_IF_SIGMA_COOO_CCOO_NO17_X1_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, T2b.cptr(), W54aaa_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   18] -- 
  // |-- [    0] --| W55(c0,k,i,a1) += (    1.00000000) V2(c0,a1,a2,a0) D2(k,i,a2,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) T2(w,c0,a1,j) W55(c0,k,i,a1) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W55aaa_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sc0));
  FC_FUNC(g_if_sigma_cooo_ccoo_no18_x0_type1_eri_c,G_IF_SIGMA_COOO_CCOO_NO18_X0_TYPE1_ERI_C)
    (sc0, ic0, V2_sym.cptr(), W55aaa_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccoo_no18_x1_type1_eri_c,G_IF_SIGMA_COOO_CCOO_NO18_X1_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, T2b.cptr(), W55aaa_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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

//-@type(2).declaration(begin)
  // --  Title : sigma_cooo_ccoo
  //  >> Intermediates for the external contractions are defined here << 
  orz::DTensor W32ca_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W33ca_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W34caaa_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W35caaa_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W36caaa_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W37caaa_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W38caaa_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W39caaa_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W40ca_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W41ca_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W46caaa_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W47caaa_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W66ca_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W67ca_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W68ca_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W69ca_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
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
  // -- Title : sigma_cooo_ccoo

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W14(c0,k,a1,a0,j,i) += (    1.00000000) V2(i,a3,c0,a2) D3(k,j,a1,a3,a0,a2) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) T2(w,c0,a1,a0) W14(c0,k,a1,a0,j,i) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      orz::DTensor W14caa_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa0^sj^si));
      FC_FUNC(g_if_sigma_cooo_ccoo_no0_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO0_X0_TYPE1_ERI_O)
        (sa0, ia0, si, ii, sj, ij, V2_sym.cptr(), W14caa_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_ccoo_no0_x1_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO0_X1_TYPE1_ERI_O)
        (sa0, ia0, si, ii, sj, ij, T2b.cptr(), W14caa_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // |-- [    0] --| W19(c0,k,a1,a0,i,j) += (    1.00000000) V2(j,a3,c0,a2) D3(k,a3,a1,i,a0,a2) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) T2(w,c0,a1,a0) W19(c0,k,a1,a0,i,j) 
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
    orz::DTensor W19caaa_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, sa0^sj));
    FC_FUNC(g_if_sigma_cooo_ccoo_no1_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO1_X0_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, V2_sym.cptr(), W19caaa_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_ccoo_no1_x1_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO1_X1_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, T2b.cptr(), W19caaa_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W24(c0,k,a0,i) += (    1.00000000) V2(i,a2,c0,a1) D2(k,a2,a0,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) T2(c0,w,a0,j) W24(c0,k,a0,i) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W24caa_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, si));
  FC_FUNC(g_if_sigma_cooo_ccoo_no2_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO2_X0_TYPE1_ERI_O)
    (si, ii, V2_sym.cptr(), W24caa_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccoo_no2_x1_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO2_X1_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), W24caa_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W25(c0,k,a0,i) += (    1.00000000) V2(i,a2,c0,a1) D2(k,a1,a0,a2) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) T2(w,c0,a0,j) W25(c0,k,a0,i) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W25caa_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, si));
  FC_FUNC(g_if_sigma_cooo_ccoo_no3_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO3_X0_TYPE1_ERI_O)
    (si, ii, V2_sym.cptr(), W25caa_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccoo_no3_x1_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO3_X1_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), W25caa_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W26(c0,k,a0,j) += (    1.00000000) V2(j,a2,c0,a1) D2(k,a2,a0,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (   -4.00000000) T2(w,c0,i,a0) W26(c0,k,a0,j) 
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
    orz::DTensor W26ca_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xca(symblockinfo, sa0^sj));
    FC_FUNC(g_if_sigma_cooo_ccoo_no4_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO4_X0_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, V2_sym.cptr(), W26ca_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_ccoo_no4_x1_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO4_X1_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, T2b.cptr(), W26ca_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W27(c0,k,a0,j) += (    1.00000000) V2(j,a2,c0,a1) D2(k,a2,a0,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) T2(w,c0,a0,i) W27(c0,k,a0,j) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W27caa_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
  FC_FUNC(g_if_sigma_cooo_ccoo_no5_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO5_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W27caa_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    FC_FUNC(g_if_sigma_cooo_ccoo_no5_x1_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO5_X1_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), W27caa_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W32(w,a0) += (    1.00000000) V2(a1,c1,w,c0) T2(c1,c0,a0,a1) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(ia1);
  FC_FUNC(g_if_sigma_cooo_ccoo_no6_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO6_X0_TYPE1_ERI_O)
    (sa1, ia1, T2b.cptr(), V2_sym.cptr(), W32ca_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W33(w,a0) += (    1.00000000) V2(a1,c1,w,c0) T2(c0,c1,a0,a1) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(ia1);
  FC_FUNC(g_if_sigma_cooo_ccoo_no7_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO7_X0_TYPE1_ERI_O)
    (sa1, ia1, T2b.cptr(), V2_sym.cptr(), W33ca_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| W34(w,a0,a3,a2) += (    1.00000000) V2(a1,a3,c0,a2) T2(c0,w,a0,a1) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(ia1);
  FC_FUNC(g_if_sigma_cooo_ccoo_no8_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO8_X0_TYPE1_ERI_O)
    (sa1, ia1, T2b.cptr(), V2_sym.cptr(), W34caaa_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W35(w,a0,a3,a2) += (    1.00000000) V2(a1,a3,c0,a2) T2(w,c0,a0,a1) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(ia1);
  FC_FUNC(g_if_sigma_cooo_ccoo_no9_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO9_X0_TYPE1_ERI_O)
    (sa1, ia1, T2b.cptr(), V2_sym.cptr(), W35caaa_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| W36(w,a0,a3,a1) += (    1.00000000) V2(a2,c0,a3,a1) T2(c0,w,a0,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(ia2);
  FC_FUNC(g_if_sigma_cooo_ccoo_no10_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO10_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), V2_sym.cptr(), W36caaa_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| W37(w,a0,a3,a1) += (    1.00000000) V2(a2,c0,a3,a1) T2(w,c0,a0,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(ia2);
  FC_FUNC(g_if_sigma_cooo_ccoo_no11_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO11_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), V2_sym.cptr(), W37caaa_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| W38(w,a0,i,a2) += (    1.00000000) V2(a1,c0,i,a2) T2(c0,w,a0,a1) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(ia1);
  FC_FUNC(g_if_sigma_cooo_ccoo_no12_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO12_X0_TYPE1_ERI_O)
    (sa1, ia1, T2b.cptr(), V2_sym.cptr(), W38caaa_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| W39(w,a0,i,a2) += (    1.00000000) V2(a1,c0,i,a2) T2(w,c0,a0,a1) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(ia1);
  FC_FUNC(g_if_sigma_cooo_ccoo_no13_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO13_X0_TYPE1_ERI_O)
    (sa1, ia1, T2b.cptr(), V2_sym.cptr(), W39caaa_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   14] -- 
  // |-- [    0] --| W40(w,i) += (    1.00000000) V2(a0,c1,w,c0) T2(c0,c1,i,a0) 
  int sa0(s_eri);
  int ia0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(ia0);
  FC_FUNC(g_if_sigma_cooo_ccoo_no14_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO14_X0_TYPE1_ERI_O)
    (sa0, ia0, T2b.cptr(), V2_sym.cptr(), W40ca_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   15] -- 
  // |-- [    0] --| W41(w,i) += (    1.00000000) V2(a0,c1,w,c0) T2(c1,c0,i,a0) 
  int sa0(s_eri);
  int ia0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(ia0);
  FC_FUNC(g_if_sigma_cooo_ccoo_no15_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO15_X0_TYPE1_ERI_O)
    (sa0, ia0, T2b.cptr(), V2_sym.cptr(), W41ca_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   16] -- 
  // |-- [    0] --| W42(c0,k,j,a0) += (    1.00000000) V2(a0,a2,c0,a1) D2(k,j,a2,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (   -4.00000000) T2(w,c0,i,a0) W42(c0,k,j,a0) 
  int sa0(s_eri);
  int ia0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(ia0);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W42ca_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sa0));
    FC_FUNC(g_if_sigma_cooo_ccoo_no16_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO16_X0_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, V2_sym.cptr(), W42ca_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_ccoo_no16_x1_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO16_X1_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, T2b.cptr(), W42ca_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   17] -- 
  // |-- [    0] --| W43(c0,k,j,a0) += (    1.00000000) V2(a0,a2,c0,a1) D2(k,j,a2,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) T2(c0,w,i,a0) W43(c0,k,j,a0) 
  int sa0(s_eri);
  int ia0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(ia0);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W43ca_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sa0));
    FC_FUNC(g_if_sigma_cooo_ccoo_no17_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO17_X0_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, V2_sym.cptr(), W43ca_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_ccoo_no17_x1_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO17_X1_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, T2b.cptr(), W43ca_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   18] -- 
  // |-- [    0] --| W44(c0,k,j,a1) += (    1.00000000) V2(a1,c0,a2,a0) D2(k,j,a2,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    8.00000000) T2(w,c0,i,a1) W44(c0,k,j,a1) 
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
    orz::DTensor W44ca_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sa1));
    FC_FUNC(g_if_sigma_cooo_ccoo_no18_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO18_X0_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, V2_sym.cptr(), W44ca_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_ccoo_no18_x1_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO18_X1_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, T2b.cptr(), W44ca_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   19] -- 
  // |-- [    0] --| W45(c0,k,j,a1) += (    1.00000000) V2(a1,c0,a2,a0) D2(k,j,a2,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -4.00000000) T2(c0,w,i,a1) W45(c0,k,j,a1) 
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
    orz::DTensor W45ca_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sa1));
    FC_FUNC(g_if_sigma_cooo_ccoo_no19_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO19_X0_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, V2_sym.cptr(), W45ca_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_ccoo_no19_x1_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO19_X1_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, T2b.cptr(), W45ca_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   20] -- 
  // |-- [    0] --| W46(w,a0,i,a1) += (    1.00000000) V2(a2,i,c0,a1) T2(c0,w,a0,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(ia2);
  FC_FUNC(g_if_sigma_cooo_ccoo_no20_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO20_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), V2_sym.cptr(), W46caaa_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   21] -- 
  // |-- [    0] --| W47(w,a0,i,a1) += (    1.00000000) V2(a2,i,c0,a1) T2(w,c0,a0,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(ia2);
  FC_FUNC(g_if_sigma_cooo_ccoo_no21_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO21_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), V2_sym.cptr(), W47caaa_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   22] -- 
  // |-- [    0] --| W48(w,a0,j,a2) += (    1.00000000) V2(j,a2,c0,a1) T2(w,c0,a1,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) D2(k,a2,a0,i) W48(w,a0,j,a2) 
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
    orz::DTensor W48ca_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xca(symblockinfo, sa0^sj));
    FC_FUNC(g_if_sigma_cooo_ccoo_no22_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO22_X0_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, T2b.cptr(), V2_sym.cptr(), W48ca_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_ccoo_no22_x1_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO22_X1_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, W48ca_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   23] -- 
  // |-- [    0] --| W49(w,a0,j,a2) += (    1.00000000) V2(j,a2,c0,a1) T2(w,c0,a0,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (   -4.00000000) D2(k,a2,a0,i) W49(w,a0,j,a2) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W49caa_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_ccoo_no23_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO23_X0_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, T2b.cptr(), V2_sym.cptr(), W49caa_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  FC_FUNC(g_if_sigma_cooo_ccoo_no23_x1_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO23_X1_TYPE1_ERI_O)
    (sj, ij, W49caa_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   24] -- 
  // |-- [    0] --| W56(w,a0,j,a1) += (    1.00000000) V2(j,a2,c0,a1) T2(w,c0,a2,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) D2(k,i,a0,a1) W56(w,a0,j,a1) 
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
    orz::DTensor W56ca_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xca(symblockinfo, sa0^sj));
    FC_FUNC(g_if_sigma_cooo_ccoo_no24_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO24_X0_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, T2b.cptr(), V2_sym.cptr(), W56ca_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_ccoo_no24_x1_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO24_X1_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, W56ca_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   25] -- 
  // |-- [    0] --| W57(w,a0,j,a1) += (    1.00000000) V2(j,a2,c0,a1) T2(w,c0,a0,a2) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) D2(k,a1,a0,i) W57(w,a0,j,a1) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W57caa_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
  for(int sa2 = 0;sa2 < nir;++sa2){ 
  for(int ia2 = symblockinfo.psym()(sa2,I_O,I_BEGIN);ia2 <= symblockinfo.psym()(sa2,I_O,I_END);++ia2){ 
    T2b = T2.get_amp2(ia2);
    FC_FUNC(g_if_sigma_cooo_ccoo_no25_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO25_X0_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, T2b.cptr(), V2_sym.cptr(), W57caa_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
  // --> @[a2, "active"] [notNeeded]
  } // End ia2
  } // End sa2
  FC_FUNC(g_if_sigma_cooo_ccoo_no25_x1_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO25_X1_TYPE1_ERI_O)
    (sj, ij, W57caa_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   26] -- 
  // |-- [    0] --| W58(c0,k,i,a0) += (    1.00000000) V2(i,a1,c0,a0) D1(k,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (   -4.00000000) T2(c0,w,a0,j) W58(c0,k,i,a0) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W58caa_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, si));
  FC_FUNC(g_if_sigma_cooo_ccoo_no26_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO26_X0_TYPE1_ERI_O)
    (si, ii, V2_sym.cptr(), W58caa_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccoo_no26_x1_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO26_X1_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), W58caa_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   27] -- 
  // |-- [    0] --| W59(c0,k,i,a0) += (    1.00000000) V2(i,a1,c0,a0) D1(k,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) T2(w,c0,a0,j) W59(c0,k,i,a0) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W59caa_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, si));
  FC_FUNC(g_if_sigma_cooo_ccoo_no27_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO27_X0_TYPE1_ERI_O)
    (si, ii, V2_sym.cptr(), W59caa_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccoo_no27_x1_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO27_X1_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), W59caa_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   28] -- 
  // |-- [    0] --| W60(c0,k,j,a0) += (    1.00000000) V2(j,a1,c0,a0) D1(k,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (    8.00000000) T2(w,c0,i,a0) W60(c0,k,j,a0) 
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
    orz::DTensor W60ca_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sa0));
    FC_FUNC(g_if_sigma_cooo_ccoo_no28_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO28_X0_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, V2_sym.cptr(), W60ca_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_ccoo_no28_x1_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO28_X1_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, T2b.cptr(), W60ca_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   29] -- 
  // |-- [    0] --| W61(c0,k,j,a0) += (    1.00000000) V2(j,a1,c0,a0) D1(k,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (   -4.00000000) T2(w,c0,a0,i) W61(c0,k,j,a0) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W61caa_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
  FC_FUNC(g_if_sigma_cooo_ccoo_no29_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO29_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W61caa_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    FC_FUNC(g_if_sigma_cooo_ccoo_no29_x1_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO29_X1_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), W61caa_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   30] -- 
  // |-- [    0] --| W62(c0,k,j,a1) += (    1.00000000) V2(j,a1,c0,a0) D1(k,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -4.00000000) T2(w,c0,i,a1) W62(c0,k,j,a1) 
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
    orz::DTensor W62ca_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sa1));
    FC_FUNC(g_if_sigma_cooo_ccoo_no30_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO30_X0_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, V2_sym.cptr(), W62ca_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_ccoo_no30_x1_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO30_X1_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, T2b.cptr(), W62ca_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   31] -- 
  // |-- [    0] --| W63(c0,k,j,a1) += (    1.00000000) V2(j,a1,c0,a0) D1(k,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) T2(w,c0,a1,i) W63(c0,k,j,a1) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W63caa_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
  FC_FUNC(g_if_sigma_cooo_ccoo_no31_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO31_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W63caa_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    FC_FUNC(g_if_sigma_cooo_ccoo_no31_x1_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO31_X1_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), W63caa_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   32] -- 
  // |-- [    0] --| W64(c0,k,i,a1) += (    1.00000000) V2(i,a1,c0,a0) D1(k,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -4.00000000) T2(w,c0,a1,j) W64(c0,k,i,a1) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W64caa_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, si));
  FC_FUNC(g_if_sigma_cooo_ccoo_no32_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO32_X0_TYPE1_ERI_O)
    (si, ii, V2_sym.cptr(), W64caa_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccoo_no32_x1_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO32_X1_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), W64caa_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   33] -- 
  // |-- [    0] --| W65(c0,k,i,a1) += (    1.00000000) V2(i,a1,c0,a0) D1(k,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) T2(c0,w,a1,j) W65(c0,k,i,a1) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W65caa_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, si));
  FC_FUNC(g_if_sigma_cooo_ccoo_no33_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO33_X0_TYPE1_ERI_O)
    (si, ii, V2_sym.cptr(), W65caa_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccoo_no33_x1_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO33_X1_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), W65caa_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   34] -- 
  // |-- [    0] --| W66(w,a2) += (    1.00000000) V2(a1,c0,a2,a0) T2(c0,w,a0,a1) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(ia1);
  FC_FUNC(g_if_sigma_cooo_ccoo_no34_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO34_X0_TYPE1_ERI_O)
    (sa1, ia1, T2b.cptr(), V2_sym.cptr(), W66ca_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   35] -- 
  // |-- [    0] --| W67(w,a2) += (    1.00000000) V2(a1,c0,a2,a0) T2(w,c0,a0,a1) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(ia1);
  FC_FUNC(g_if_sigma_cooo_ccoo_no35_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO35_X0_TYPE1_ERI_O)
    (sa1, ia1, T2b.cptr(), V2_sym.cptr(), W67ca_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   36] -- 
  // |-- [    0] --| W68(w,i) += (    1.00000000) V2(a1,i,c0,a0) T2(c0,w,a0,a1) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(ia1);
  FC_FUNC(g_if_sigma_cooo_ccoo_no36_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO36_X0_TYPE1_ERI_O)
    (sa1, ia1, T2b.cptr(), V2_sym.cptr(), W68ca_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   37] -- 
  // |-- [    0] --| W69(w,i) += (    1.00000000) V2(a1,i,c0,a0) T2(w,c0,a0,a1) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(ia1);
  FC_FUNC(g_if_sigma_cooo_ccoo_no37_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO37_X0_TYPE1_ERI_O)
    (sa1, ia1, T2b.cptr(), V2_sym.cptr(), W69ca_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   38] -- 
  // |-- [    0] --| W70(w,j) += (    1.00000000) V2(j,a1,c0,a0) T2(w,c0,a1,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -4.00000000) D1(k,i) W70(w,j) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W70c_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xc(symblockinfo, sj));
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_cooo_ccoo_no38_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO38_X0_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, T2b.cptr(), V2_sym.cptr(), W70c_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  FC_FUNC(g_if_sigma_cooo_ccoo_no38_x1_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO38_X1_TYPE1_ERI_O)
    (sj, ij, W70c_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   39] -- 
  // |-- [    0] --| W71(w,j) += (    1.00000000) V2(j,a1,c0,a0) T2(w,c0,a0,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) D1(k,i) W71(w,j) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W71c_sigma_cooo_ccoo(orz::mr::sizeof_sympack_Xc(symblockinfo, sj));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_cooo_ccoo_no39_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO39_X0_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, T2b.cptr(), V2_sym.cptr(), W71c_sigma_cooo_ccoo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  FC_FUNC(g_if_sigma_cooo_ccoo_no39_x1_type1_eri_o,G_IF_SIGMA_COOO_CCOO_NO39_X1_TYPE1_ERI_O)
    (sj, ij, W71c_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

//-@type(1).contraction(end)

  } // End myrank
  }
  }
  orz::world().barrier();

//-@type(2).contraction(begin)
  // -- Title : sigma_cooo_ccoo
  //*-- Entering to take the type 2 contractions --*//

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -2.00000000) D2(k,j,a0,i) W32(w,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccoo_no0_x0_type2_eri_o,G_IF_SIGMA_COOO_CCOO_NO0_X0_TYPE2_ERI_O)
      (sj, ij, W32ca_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    4.00000000) D2(k,j,a0,i) W33(w,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccoo_no1_x0_type2_eri_o,G_IF_SIGMA_COOO_CCOO_NO1_X0_TYPE2_ERI_O)
      (sj, ij, W33ca_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    2.00000000) D3(k,j,a3,i,a0,a2) W34(w,a0,a3,a2) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccoo_no2_x0_type2_eri_o,G_IF_SIGMA_COOO_CCOO_NO2_X0_TYPE2_ERI_O)
      (sj, ij, W34caaa_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    2.00000000) D3(k,j,a3,a2,a0,i) W35(w,a0,a3,a2) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccoo_no3_x0_type2_eri_o,G_IF_SIGMA_COOO_CCOO_NO3_X0_TYPE2_ERI_O)
      (sj, ij, W35caaa_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    2.00000000) D3(k,j,a3,a1,a0,i) W36(w,a0,a3,a1) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccoo_no4_x0_type2_eri_o,G_IF_SIGMA_COOO_CCOO_NO4_X0_TYPE2_ERI_O)
      (sj, ij, W36caaa_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -4.00000000) D3(k,j,a3,a1,a0,i) W37(w,a0,a3,a1) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccoo_no5_x0_type2_eri_o,G_IF_SIGMA_COOO_CCOO_NO5_X0_TYPE2_ERI_O)
      (sj, ij, W37caaa_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    2.00000000) D2(k,j,a0,a2) W38(w,a0,i,a2) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccoo_no6_x0_type2_eri_o,G_IF_SIGMA_COOO_CCOO_NO6_X0_TYPE2_ERI_O)
      (sj, ij, W38caaa_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -4.00000000) D2(k,j,a0,a2) W39(w,a0,i,a2) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccoo_no7_x0_type2_eri_o,G_IF_SIGMA_COOO_CCOO_NO7_X0_TYPE2_ERI_O)
      (sj, ij, W39caaa_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -8.00000000) D1(k,j) W40(w,i) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccoo_no8_x0_type2_eri_o,G_IF_SIGMA_COOO_CCOO_NO8_X0_TYPE2_ERI_O)
      (sj, ij, W40ca_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    4.00000000) D1(k,j) W41(w,i) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccoo_no9_x0_type2_eri_o,G_IF_SIGMA_COOO_CCOO_NO9_X0_TYPE2_ERI_O)
      (sj, ij, W41ca_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -4.00000000) D2(k,j,a0,a1) W46(w,a0,i,a1) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccoo_no10_x0_type2_eri_o,G_IF_SIGMA_COOO_CCOO_NO10_X0_TYPE2_ERI_O)
      (sj, ij, W46caaa_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    2.00000000) D2(k,j,a0,a1) W47(w,a0,i,a1) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccoo_no11_x0_type2_eri_o,G_IF_SIGMA_COOO_CCOO_NO11_X0_TYPE2_ERI_O)
      (sj, ij, W47caaa_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    2.00000000) D2(k,j,a2,i) W66(w,a2) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccoo_no12_x0_type2_eri_o,G_IF_SIGMA_COOO_CCOO_NO12_X0_TYPE2_ERI_O)
      (sj, ij, W66ca_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -4.00000000) D2(k,j,a2,i) W67(w,a2) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccoo_no13_x0_type2_eri_o,G_IF_SIGMA_COOO_CCOO_NO13_X0_TYPE2_ERI_O)
      (sj, ij, W67ca_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   14] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    8.00000000) D1(k,j) W68(w,i) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccoo_no14_x0_type2_eri_o,G_IF_SIGMA_COOO_CCOO_NO14_X0_TYPE2_ERI_O)
      (sj, ij, W68ca_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   15] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -4.00000000) D1(k,j) W69(w,i) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccoo_no15_x0_type2_eri_o,G_IF_SIGMA_COOO_CCOO_NO15_X0_TYPE2_ERI_O)
      (sj, ij, W69ca_sigma_cooo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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

//-@loadD4C(c,begin)
  //*-- FEMTO begins --//*
  // Label : d4c_c
  {


  for(int s_eri = 0;s_eri < nir;++s_eri){ 
  for(int i_eri = symblockinfo.psym()(s_eri,I_C,I_BEGIN);i_eri <= symblockinfo.psym()(s_eri,I_C,I_END);++i_eri){ 
  if(hintmo.iproc_havingimo()[i_eri] == myrank) {           
  orz::DTensor C5;
  orz::LoadBin(ctinp.dir()/(format("D4C_g[%d]")%i_eri).str()) >> C5;

  //*-- Entering to take the type 1 contractions --*//
//-@type(1).contraction(begin)
  // -- Title : sigma_cooo_ccoo

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    2.00000000) T2(w,c0,a1,a0) C5(k,j,a1,i,a0,c0) 
  int sc0(s_eri);
  int ic0(i_eri);
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
      FC_FUNC(g_if_sigma_cooo_ccoo_no0_x0_type1_d4c_c,G_IF_SIGMA_COOO_CCOO_NO0_X0_TYPE1_D4C_C)
        (sa0, ia0, sc0, ic0, sj, ij, C5.cptr(), T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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

//-@loadD4C(c,end)

//-@D4C.contractions(end)

  // Do timing!
  orz::mr::Femto::my_timer.push_back(boost::make_tuple("sigma_cooo_ccoo", time_sigma.elapsed_cputime(), time_sigma.elapsed_wallclocktime()));
  orz::mr::Femto::file_timing << "* " << boost::format("%20s : %10.7f %10.7f ") % orz::mr::Femto::my_timer.back().get<0>() % orz::mr::Femto::my_timer.back().get<1>() % orz::mr::Femto::my_timer.back().get<2>() << endl;
  flush(orz::mr::Femto::file_timing);

  return retval; 
} 
