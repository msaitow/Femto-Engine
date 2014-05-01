                                                                                
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
#include <sci/icmr/Femto/elems/c_sigma_ccoo_cooo.h>                                  
                                                                                
using std::cout;                                                                
using std::endl;                                                                
                                                                                
#define FLOPCOUNT                                                               
                                                                                
using orz::mr::Femto::Fc0;                                                      
using orz::mr::Femto::h1_int;                                                   
using orz::mr::Femto::h6_int;                                                   
                                                                                
//                                                              
//   _______________                                  ______    
//  |          |                 .'. .`. `````|`````.~      ~.  
//  |______    |______         .'   `   `.    |    |          | 
//  |          |             .'           `.  |    |          | 
//  |          |___________.'               `.|     `.______.'  
//                                                              

//                                   Generated date : Sun Apr 20 10:26:21 2014

                                                                                
// ***************************************************************************  
// orz::mr::mrci                                                                
// ***************************************************************************  
									    /*!        
   @brief CT input                                                              
									     */        
                                                                                
orz::mr::BareAmpPack orz::mr::Femto::sigma_ccoo_cooo(const orz::mr::Input &ctinp,                                    
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
  // |-- [    0] --| W0(x,j,i,a3) += (    1.00000000) T2(x,a2,a1,a0) D3(j,a1,i,a3,a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) Fc1(w,a3) W0(x,j,i,a3) 
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W0caa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
    for(int sa0 = 0;sa0 < nir;++sa0){ 
    for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
      T2b = T2.get_amp2(ia0);
      FC_FUNC(g_if_sigma_ccoo_cooo_no0_x0_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO0_X0_TYPE0_NOERI)
        (sa0, ia0, sj, ij, T2b.cptr(), W0caa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
    } // End ia0
    } // End sa0
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no0_x1_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO0_X1_TYPE0_NOERI)
      (sj, ij, W0caa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W1(w,j,i,a3) += (    1.00000000) T2(w,a2,a1,a0) D3(j,a3,i,a1,a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) Fc1(x,a3) W1(w,j,i,a3) 
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W1caa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
    for(int sa0 = 0;sa0 < nir;++sa0){ 
    for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
      T2b = T2.get_amp2(ia0);
      FC_FUNC(g_if_sigma_ccoo_cooo_no1_x0_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO1_X0_TYPE0_NOERI)
        (sa0, ia0, sj, ij, T2b.cptr(), W1caa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
    } // End ia0
    } // End sa0
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no1_x1_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO1_X1_TYPE0_NOERI)
      (sj, ij, W1caa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W2(x,j) += (    1.00000000) T2(x,a2,a1,a0) D2(j,a1,a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) Fc1(w,i) W2(x,j) 
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W2c_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, sj));
    for(int sa0 = 0;sa0 < nir;++sa0){ 
    for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
      T2b = T2.get_amp2(ia0);
      FC_FUNC(g_if_sigma_ccoo_cooo_no2_x0_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO2_X0_TYPE0_NOERI)
        (sa0, ia0, sj, ij, T2b.cptr(), W2c_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
    } // End ia0
    } // End sa0
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no2_x1_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO2_X1_TYPE0_NOERI)
      (sj, ij, W2c_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W3(w,j) += (    1.00000000) T2(w,a2,a1,a0) D2(j,a1,a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) Fc1(x,i) W3(w,j) 
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W3c_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, sj));
    for(int sa0 = 0;sa0 < nir;++sa0){ 
    for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
      T2b = T2.get_amp2(ia0);
      FC_FUNC(g_if_sigma_ccoo_cooo_no3_x0_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO3_X0_TYPE0_NOERI)
        (sa0, ia0, sj, ij, T2b.cptr(), W3c_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
    } // End ia0
    } // End sa0
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no3_x1_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO3_X1_TYPE0_NOERI)
      (sj, ij, W3c_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W4(x,j,a2,i) += (    1.00000000) T2(x,a1,i,a0) D2(j,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) Fc1(w,a2) W4(x,j,a2,i) 
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W4caa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
    for(int sa0 = 0;sa0 < nir;++sa0){ 
    for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
      T2b = T2.get_amp2(ia0);
      FC_FUNC(g_if_sigma_ccoo_cooo_no4_x0_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO4_X0_TYPE0_NOERI)
        (sa0, ia0, sj, ij, T2b.cptr(), W4caa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
    } // End ia0
    } // End sa0
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no4_x1_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO4_X1_TYPE0_NOERI)
      (sj, ij, W4caa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W5(w,j,a2,i) += (    1.00000000) T2(w,a1,i,a0) D2(j,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) Fc1(x,a2) W5(w,j,a2,i) 
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W5caa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
    for(int sa0 = 0;sa0 < nir;++sa0){ 
    for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
      T2b = T2.get_amp2(ia0);
      FC_FUNC(g_if_sigma_ccoo_cooo_no5_x0_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO5_X0_TYPE0_NOERI)
        (sa0, ia0, sj, ij, T2b.cptr(), W5caa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
    } // End ia0
    } // End sa0
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no5_x1_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO5_X1_TYPE0_NOERI)
      (sj, ij, W5caa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W6(x,j,a2,i) += (    1.00000000) T2(x,a1,a0,i) D2(j,a0,a1,a2) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) Fc1(w,a2) W6(x,j,a2,i) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      orz::DTensor W6ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^si));
      FC_FUNC(g_if_sigma_ccoo_cooo_no6_x0_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO6_X0_TYPE0_NOERI)
        (si, ii, sj, ij, T2b.cptr(), W6ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_ccoo_cooo_no6_x1_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO6_X1_TYPE0_NOERI)
        (si, ii, sj, ij, W6ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- No. [    7] -- 
  // |-- [    0] --| W7(w,j,a2,i) += (    1.00000000) T2(w,a1,a0,i) D2(j,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) Fc1(x,a2) W7(w,j,a2,i) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      orz::DTensor W7ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^si));
      FC_FUNC(g_if_sigma_ccoo_cooo_no7_x0_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO7_X0_TYPE0_NOERI)
        (si, ii, sj, ij, T2b.cptr(), W7ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_ccoo_cooo_no7_x1_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO7_X1_TYPE0_NOERI)
        (si, ii, sj, ij, W7ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- No. [    8] -- 
  // |-- [    0] --| W8(x,i) += (    1.00000000) T2(x,a2,a1,a0) D2(i,a1,a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) Fc1(w,j) W8(x,i) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W8ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccoo_cooo_no8_x0_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO8_X0_TYPE0_NOERI)
      (sa0, ia0, T2b.cptr(), W8ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no8_x1_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO8_X1_TYPE0_NOERI)
      (sj, ij, W8ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W9(w,i) += (    1.00000000) T2(w,a2,a1,a0) D2(i,a1,a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) Fc1(x,j) W9(w,i) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W9ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccoo_cooo_no9_x0_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO9_X0_TYPE0_NOERI)
      (sa0, ia0, T2b.cptr(), W9ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no9_x1_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO9_X1_TYPE0_NOERI)
      (sj, ij, W9ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| W10(x,i,a2,j) += (    1.00000000) T2(a1,x,a0,j) D2(i,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) Fc1(w,a2) W10(x,i,a2,j) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W10caa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no10_x0_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO10_X0_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W10caa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_cooo_no10_x1_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO10_X1_TYPE0_NOERI)
      (sj, ij, W10caa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| W11(w,i,a2,j) += (    1.00000000) T2(a1,w,a0,j) D2(i,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) Fc1(x,a2) W11(w,i,a2,j) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W11caa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no11_x0_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO11_X0_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W11caa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_cooo_no11_x1_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO11_X1_TYPE0_NOERI)
      (sj, ij, W11caa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| W12(x,i,a2,j) += (    1.00000000) T2(x,a1,a0,j) D2(i,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) Fc1(w,a2) W12(x,i,a2,j) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W12caa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no12_x0_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO12_X0_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W12caa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_cooo_no12_x1_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO12_X1_TYPE0_NOERI)
      (sj, ij, W12caa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| W13(w,i,a2,j) += (    1.00000000) T2(w,a1,a0,j) D2(i,a0,a1,a2) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) Fc1(x,a2) W13(w,i,a2,j) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W13caa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no13_x0_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO13_X0_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W13caa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_cooo_no13_x1_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO13_X1_TYPE0_NOERI)
      (sj, ij, W13caa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   14] -- 
  // |-- [    0] --| W14(x,j) += (    1.00000000) T2(a1,x,a0,j) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    4.00000000) Fc1(w,i) W14(x,j) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W14c_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no14_x0_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO14_X0_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W14c_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_cooo_no14_x1_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO14_X1_TYPE0_NOERI)
      (sj, ij, W14c_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   15] -- 
  // |-- [    0] --| W15(w,j) += (    1.00000000) T2(a1,w,a0,j) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) Fc1(x,i) W15(w,j) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W15c_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no15_x0_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO15_X0_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W15c_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_cooo_no15_x1_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO15_X1_TYPE0_NOERI)
      (sj, ij, W15c_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   16] -- 
  // |-- [    0] --| W16(x,j) += (    1.00000000) T2(x,a1,a0,j) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) Fc1(w,i) W16(x,j) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W16c_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no16_x0_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO16_X0_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W16c_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_cooo_no16_x1_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO16_X1_TYPE0_NOERI)
      (sj, ij, W16c_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   17] -- 
  // |-- [    0] --| W17(w,j) += (    1.00000000) T2(w,a1,a0,j) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) Fc1(x,i) W17(w,j) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W17c_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no17_x0_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO17_X0_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W17c_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_cooo_no17_x1_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO17_X1_TYPE0_NOERI)
      (sj, ij, W17c_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   18] -- 
  // |-- [    0] --| W18(w,i) += (    1.00000000) T2(w,a1,i,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    4.00000000) Fc1(x,j) W18(w,i) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W18ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccoo_cooo_no18_x0_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO18_X0_TYPE0_NOERI)
      (sa0, ia0, T2b.cptr(), W18ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no18_x1_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO18_X1_TYPE0_NOERI)
      (sj, ij, W18ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   19] -- 
  // |-- [    0] --| W19(x,i) += (    1.00000000) T2(x,a1,i,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) Fc1(w,j) W19(x,i) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W19ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccoo_cooo_no19_x0_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO19_X0_TYPE0_NOERI)
      (sa0, ia0, T2b.cptr(), W19ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no19_x1_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO19_X1_TYPE0_NOERI)
      (sj, ij, W19ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   20] -- 
  // |-- [    0] --| W20(w,a0) += (    1.00000000) D1(a1,a0) Fc1(w,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(x,a0,i,j) W20(w,a0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W20ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  FC_FUNC(g_if_sigma_ccoo_cooo_no20_x0_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO20_X0_TYPE0_NOERI)
    (W20ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no20_x1_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO20_X1_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W20ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   21] -- 
  // |-- [    0] --| W21(x,a0) += (    1.00000000) D1(a1,a0) Fc1(x,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(w,a0,i,j) W21(x,a0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W21ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  FC_FUNC(g_if_sigma_ccoo_cooo_no21_x0_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO21_X0_TYPE0_NOERI)
    (W21ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no21_x1_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO21_X1_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W21ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   22] -- 
  // |-- [    0] --| W22(x,i) += (    1.00000000) T2(x,a1,a0,i) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) Fc1(w,j) W22(x,i) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    orz::DTensor W22c_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, si));
    FC_FUNC(g_if_sigma_ccoo_cooo_no22_x0_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO22_X0_TYPE0_NOERI)
      (si, ii, T2b.cptr(), W22c_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_ccoo_cooo_no22_x1_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO22_X1_TYPE0_NOERI)
        (si, ii, sj, ij, W22c_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- No. [   23] -- 
  // |-- [    0] --| W23(w,i) += (    1.00000000) T2(w,a1,a0,i) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) Fc1(x,j) W23(w,i) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    orz::DTensor W23c_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, si));
    FC_FUNC(g_if_sigma_ccoo_cooo_no23_x0_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO23_X0_TYPE0_NOERI)
      (si, ii, T2b.cptr(), W23c_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_ccoo_cooo_no23_x1_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO23_X1_TYPE0_NOERI)
        (si, ii, sj, ij, W23c_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- No. [   24] -- 
  // |-- [    0] --| W24(w,a0) += (    1.00000000) D1(a1,a0) Fc1(w,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(a0,x,i,j) W24(w,a0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W24ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  FC_FUNC(g_if_sigma_ccoo_cooo_no24_x0_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO24_X0_TYPE0_NOERI)
    (W24ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no24_x1_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO24_X1_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W24ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   25] -- 
  // |-- [    0] --| W25(x,a0) += (    1.00000000) D1(a1,a0) Fc1(x,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(a0,w,i,j) W25(x,a0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W25ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  FC_FUNC(g_if_sigma_ccoo_cooo_no25_x0_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO25_X0_TYPE0_NOERI)
    (W25ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no25_x1_type0_noeri,G_IF_SIGMA_CCOO_COOO_NO25_X1_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W25ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- Title : sigma_ccoo_cooo

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W27(c0,j,i,a3) += (    1.00000000) T2(c0,a2,a1,a0) D3(j,a1,i,a3,a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -1.00000000) V2(c0,x,w,a3) W27(c0,j,i,a3) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <3, 2> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W27aa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sj));
    for(int sa0 = 0;sa0 < nir;++sa0){ 
    for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
      T2b = T2.get_amp2(ia0);
      FC_FUNC(g_if_sigma_ccoo_cooo_no0_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO0_X0_TYPE1_ERI_C)
        (sa0, ia0, sc0, ic0, sj, ij, T2b.cptr(), W27aa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
    } // End ia0
    } // End sa0
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no0_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO0_X1_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, V2_sym.cptr(), W27aa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W28(c0,j) += (    1.00000000) T2(c0,a2,a1,a0) D2(j,a1,a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -1.00000000) V2(c0,w,x,i) W28(c0,j) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <3, 2> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    double W28_sigma_ccoo_cooo(0);
    for(int sa0 = 0;sa0 < nir;++sa0){ 
    for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
      T2b = T2.get_amp2(ia0);
      FC_FUNC(g_if_sigma_ccoo_cooo_no1_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO1_X0_TYPE1_ERI_C)
        (sa0, ia0, sc0, ic0, sj, ij, T2b.cptr(), &W28_sigma_ccoo_cooo, nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
    } // End ia0
    } // End sa0
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no1_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO1_X1_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, V2_sym.cptr(), &W28_sigma_ccoo_cooo, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W29(c0,j) += (    1.00000000) T2(c0,a2,a1,a0) D2(j,a1,a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) V2(c0,x,w,i) W29(c0,j) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <3, 2> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    double W29_sigma_ccoo_cooo(0);
    for(int sa0 = 0;sa0 < nir;++sa0){ 
    for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
      T2b = T2.get_amp2(ia0);
      FC_FUNC(g_if_sigma_ccoo_cooo_no2_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO2_X0_TYPE1_ERI_C)
        (sa0, ia0, sc0, ic0, sj, ij, T2b.cptr(), &W29_sigma_ccoo_cooo, nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
    } // End ia0
    } // End sa0
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no2_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO2_X1_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, V2_sym.cptr(), &W29_sigma_ccoo_cooo, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W34(c0,j,a2,i) += (    1.00000000) T2(c0,a1,i,a0) D2(j,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -1.00000000) V2(c0,x,w,a2) W34(c0,j,a2,i) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <3, 2> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W34aa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sj));
    for(int sa0 = 0;sa0 < nir;++sa0){ 
    for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
      T2b = T2.get_amp2(ia0);
      FC_FUNC(g_if_sigma_ccoo_cooo_no3_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO3_X0_TYPE1_ERI_C)
        (sa0, ia0, sc0, ic0, sj, ij, T2b.cptr(), W34aa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
    } // End ia0
    } // End sa0
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no3_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO3_X1_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, V2_sym.cptr(), W34aa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W36(w,j,a1,a0) += (    1.00000000) V2(w,a3,a4,a2) D3(j,a3,a4,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(x,a1,i,a0) W36(w,j,a1,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      orz::DTensor W36a_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xa(symblockinfo, sw^sj^sa0));
      FC_FUNC(g_if_sigma_ccoo_cooo_no4_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO4_X0_TYPE1_ERI_C)
        (sa0, ia0, sj, ij, sw, iw, V2_sym.cptr(), W36a_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_ccoo_cooo_no4_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO4_X1_TYPE1_ERI_C)
        (sa0, ia0, sj, ij, sw, iw, T2b.cptr(), W36a_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W37(x,j,a1,a0) += (    1.00000000) V2(x,a3,a4,a2) D3(j,a3,a4,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(w,a1,i,a0) W37(x,j,a1,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      orz::DTensor W37a_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xa(symblockinfo, sx^sj^sa0));
      FC_FUNC(g_if_sigma_ccoo_cooo_no5_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO5_X0_TYPE1_ERI_C)
        (sa0, ia0, sj, ij, sx, ix, V2_sym.cptr(), W37a_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_ccoo_cooo_no5_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO5_X1_TYPE1_ERI_C)
        (sa0, ia0, sj, ij, sx, ix, T2b.cptr(), W37a_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W39(c0,j,a2,i) += (    1.00000000) T2(c0,a1,a0,i) D2(j,a0,a1,a2) 
  // |-- [    1] --| S2(w,x,i,j) += (   -1.00000000) V2(c0,x,w,a2) W39(c0,j,a2,i) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      orz::DTensor W39a_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xa(symblockinfo, sc0^sj^si));
      FC_FUNC(g_if_sigma_ccoo_cooo_no6_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO6_X0_TYPE1_ERI_C)
        (sc0, ic0, si, ii, sj, ij, T2b.cptr(), W39a_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_ccoo_cooo_no6_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO6_X1_TYPE1_ERI_C)
        (sc0, ic0, si, ii, sj, ij, V2_sym.cptr(), W39a_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W40(x,j,a1,a0) += (    1.00000000) V2(x,a3,a4,a2) D3(j,a3,a4,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(w,a1,a0,i) W40(x,j,a1,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W40aa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sx^sj));
    FC_FUNC(g_if_sigma_ccoo_cooo_no7_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO7_X0_TYPE1_ERI_C)
      (sj, ij, sx, ix, V2_sym.cptr(), W40aa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    for(int si = 0;si < nir;++si){ 
    for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
      T2b = T2.get_amp2(ii);
      FC_FUNC(g_if_sigma_ccoo_cooo_no7_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO7_X1_TYPE1_ERI_C)
        (si, ii, sj, ij, sx, ix, T2b.cptr(), W40aa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    } // End ii
    } // End si
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| W41(w,j,a1,a0) += (    1.00000000) V2(w,a3,a4,a2) D3(j,a0,a4,a2,a1,a3) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(x,a1,a0,i) W41(w,j,a1,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W41aa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    FC_FUNC(g_if_sigma_ccoo_cooo_no8_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO8_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, V2_sym.cptr(), W41aa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    for(int si = 0;si < nir;++si){ 
    for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
      T2b = T2.get_amp2(ii);
      FC_FUNC(g_if_sigma_ccoo_cooo_no8_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO8_X1_TYPE1_ERI_C)
        (si, ii, sj, ij, sw, iw, T2b.cptr(), W41aa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    } // End ii
    } // End si
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W48(c0,i,a2,j) += (    1.00000000) T2(a1,c0,a0,j) D2(i,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) V2(c0,x,w,a2) W48(c0,i,a2,j) 
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
    orz::DTensor W48aa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no9_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO9_X0_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, T2b.cptr(), W48aa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_cooo_no9_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO9_X1_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, V2_sym.cptr(), W48aa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| W50(w,i,a1,a0) += (    1.00000000) V2(w,a3,a4,a2) D3(i,a3,a4,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(a1,x,a0,j) W50(w,i,a1,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W50aaa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_cooo_no10_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO10_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W50aaa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no10_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO10_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), W50aaa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| W51(x,i,a1,a0) += (    1.00000000) V2(x,a3,a4,a2) D3(i,a3,a4,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(a1,w,a0,j) W51(x,i,a1,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W51aaa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_cooo_no11_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO11_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W51aaa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no11_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO11_X1_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), W51aaa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| W53(c0,i,a2,j) += (    1.00000000) T2(c0,a1,a0,j) D2(i,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -1.00000000) V2(c0,x,w,a2) W53(c0,i,a2,j) 
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
    orz::DTensor W53aa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no12_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO12_X0_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, T2b.cptr(), W53aa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_cooo_no12_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO12_X1_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, V2_sym.cptr(), W53aa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| W54(x,i,a1,a0) += (    1.00000000) V2(x,a3,a4,a2) D3(i,a0,a4,a2,a1,a3) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(w,a1,a0,j) W54(x,i,a1,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W54aaa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_cooo_no13_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO13_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W54aaa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no13_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO13_X1_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), W54aaa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   14] -- 
  // |-- [    0] --| W55(w,i,a1,a0) += (    1.00000000) V2(w,a3,a4,a2) D3(i,a3,a4,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(x,a1,a0,j) W55(w,i,a1,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W55aaa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_cooo_no14_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO14_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W55aaa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no14_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO14_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), W55aaa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   15] -- 
  // |-- [    0] --| W60(c0,j) += (    1.00000000) T2(a1,c0,a0,j) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -4.00000000) V2(c0,x,w,i) W60(c0,j) 
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
    double W60_sigma_ccoo_cooo(0);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no15_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO15_X0_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, T2b.cptr(), &W60_sigma_ccoo_cooo, nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_cooo_no15_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO15_X1_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, V2_sym.cptr(), &W60_sigma_ccoo_cooo, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   16] -- 
  // |-- [    0] --| W61(c0,j) += (    1.00000000) T2(a1,c0,a0,j) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) V2(c0,w,x,i) W61(c0,j) 
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
    double W61_sigma_ccoo_cooo(0);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no16_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO16_X0_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, T2b.cptr(), &W61_sigma_ccoo_cooo, nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_cooo_no16_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO16_X1_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, V2_sym.cptr(), &W61_sigma_ccoo_cooo, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   17] -- 
  // |-- [    0] --| W62(w,a1,a0,i) += (    1.00000000) V2(w,i,a3,a2) D2(a3,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    4.00000000) T2(a1,x,a0,j) W62(w,a1,a0,i) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W62aaa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_cooo_no17_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO17_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W62aaa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no17_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO17_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), W62aaa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   18] -- 
  // |-- [    0] --| W63(x,a1,a0,i) += (    1.00000000) V2(x,i,a3,a2) D2(a3,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(a1,w,a0,j) W63(x,a1,a0,i) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W63aaa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_cooo_no18_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO18_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W63aaa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no18_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO18_X1_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), W63aaa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   19] -- 
  // |-- [    0] --| W64(w,a1,a0,i) += (    1.00000000) V2(w,a3,i,a2) D2(a3,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(a0,x,a1,j) W64(w,a1,a0,i) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W64aaa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_cooo_no19_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO19_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W64aaa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no19_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO19_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), W64aaa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   20] -- 
  // |-- [    0] --| W65(x,a1,a0,i) += (    1.00000000) V2(x,a3,i,a2) D2(a3,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(a0,w,a1,j) W65(x,a1,a0,i) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W65aaa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_cooo_no20_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO20_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W65aaa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no20_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO20_X1_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), W65aaa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   21] -- 
  // |-- [    0] --| W66(c0,j) += (    1.00000000) T2(c0,a1,a0,j) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -1.00000000) V2(c0,w,x,i) W66(c0,j) 
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
    double W66_sigma_ccoo_cooo(0);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no21_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO21_X0_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, T2b.cptr(), &W66_sigma_ccoo_cooo, nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_cooo_no21_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO21_X1_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, V2_sym.cptr(), &W66_sigma_ccoo_cooo, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   22] -- 
  // |-- [    0] --| W67(c0,j) += (    1.00000000) T2(c0,a1,a0,j) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) V2(c0,x,w,i) W67(c0,j) 
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
    double W67_sigma_ccoo_cooo(0);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no22_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO22_X0_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, T2b.cptr(), &W67_sigma_ccoo_cooo, nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_cooo_no22_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO22_X1_TYPE1_ERI_C)
      (sc0, ic0, sj, ij, V2_sym.cptr(), &W67_sigma_ccoo_cooo, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   23] -- 
  // |-- [    0] --| W68(x,a1,a0,i) += (    1.00000000) V2(x,i,a3,a2) D2(a3,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(w,a1,a0,j) W68(x,a1,a0,i) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W68aaa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_cooo_no23_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO23_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W68aaa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no23_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO23_X1_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), W68aaa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   24] -- 
  // |-- [    0] --| W69(w,a1,a0,i) += (    1.00000000) V2(w,i,a3,a2) D2(a3,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(x,a1,a0,j) W69(w,a1,a0,i) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W69aaa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_cooo_no24_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO24_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W69aaa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no24_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO24_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), W69aaa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   25] -- 
  // |-- [    0] --| W70(x,a1,a0,i) += (    1.00000000) V2(x,a3,i,a2) D2(a3,a0,a1,a2) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(w,a0,a1,j) W70(x,a1,a0,i) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W70aaa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_cooo_no25_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO25_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W70aaa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no25_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO25_X1_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), W70aaa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   26] -- 
  // |-- [    0] --| W71(w,a1,a0,i) += (    1.00000000) V2(w,a3,i,a2) D2(a3,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(x,a0,a1,j) W71(w,a1,a0,i) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W71aaa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_cooo_no26_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO26_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W71aaa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no26_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO26_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), W71aaa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   27] -- 
  // |-- [    0] --| W78(x,w,c0,a0) += (    1.00000000) V2(x,c0,w,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -1.00000000) T2(c0,a0,i,j) W78(x,w,c0,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W78cca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xcca(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_cooo_no27_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO27_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W78cca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no27_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO27_X1_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), W78cca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   28] -- 
  // |-- [    0] --| W79(x,w,c0,a0) += (    1.00000000) V2(x,a1,w,c0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) T2(c0,a0,i,j) W79(x,w,c0,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W79cca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xcca(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_cooo_no28_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO28_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W79cca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no28_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO28_X1_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), W79cca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   29] -- 
  // |-- [    0] --| W80(w,a0) += (    1.00000000) V2(w,a2,a3,a1) D2(a3,a1,a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(x,a0,i,j) W80(w,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W80a_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_cooo_no29_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO29_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W80a_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no29_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO29_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), W80a_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   30] -- 
  // |-- [    0] --| W81(x,a0) += (    1.00000000) V2(x,a2,a3,a1) D2(a3,a1,a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(w,a0,i,j) W81(x,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W81a_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_cooo_no30_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO30_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W81a_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no30_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO30_X1_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), W81a_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   31] -- 
  // |-- [    0] --| W88(x,w,c0,a0) += (    1.00000000) V2(x,c0,w,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) T2(a0,c0,i,j) W88(x,w,c0,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W88cca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xcca(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_cooo_no31_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO31_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W88cca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no31_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO31_X1_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), W88cca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   32] -- 
  // |-- [    0] --| W89(x,w,c0,a0) += (    1.00000000) V2(x,a1,w,c0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -1.00000000) T2(a0,c0,i,j) W89(x,w,c0,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W89cca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xcca(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_cooo_no32_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO32_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W89cca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no32_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO32_X1_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), W89cca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   33] -- 
  // |-- [    0] --| W90(w,a0) += (    1.00000000) V2(w,a2,a3,a1) D2(a3,a1,a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(a0,x,i,j) W90(w,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W90a_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_cooo_no33_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO33_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W90a_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no33_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO33_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), W90a_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   34] -- 
  // |-- [    0] --| W91(x,a0) += (    1.00000000) V2(x,a2,a3,a1) D2(a3,a1,a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(a0,w,i,j) W91(x,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W91a_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_cooo_no34_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO34_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W91a_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no34_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO34_X1_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), W91a_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   35] -- 
  // |-- [    0] --| W102(x,i,a0,a1) += (    1.00000000) V2(x,a2,a3,a1) D2(i,a3,a0,a2) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(w,a0,a1,j) W102(x,i,a0,a1) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W102aaa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_cooo_no35_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO35_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W102aaa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no35_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO35_X1_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), W102aaa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   36] -- 
  // |-- [    0] --| W103(w,i,a0,a1) += (    1.00000000) V2(w,a2,a3,a1) D2(i,a2,a0,a3) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(x,a0,a1,j) W103(w,i,a0,a1) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W103aaa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_cooo_no36_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO36_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W103aaa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no36_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO36_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), W103aaa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   37] -- 
  // |-- [    0] --| W106(x,a0,i,a1) += (    1.00000000) V2(x,i,a2,a1) D1(a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(w,a0,a1,j) W106(x,a0,i,a1) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W106aaa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_cooo_no37_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO37_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W106aaa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no37_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO37_X1_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), W106aaa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   38] -- 
  // |-- [    0] --| W107(w,a0,i,a1) += (    1.00000000) V2(w,i,a2,a1) D1(a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(x,a0,a1,j) W107(w,a0,i,a1) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W107aaa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_cooo_no38_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO38_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W107aaa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no38_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO38_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), W107aaa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   39] -- 
  // |-- [    0] --| W116(w,a0,i,a1) += (    1.00000000) V2(w,a2,i,a1) D1(a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(x,a0,a1,j) W116(w,a0,i,a1) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W116aaa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_cooo_no39_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO39_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W116aaa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no39_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO39_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), W116aaa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   40] -- 
  // |-- [    0] --| W117(x,a0,i,a1) += (    1.00000000) V2(x,a2,i,a1) D1(a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(w,a0,a1,j) W117(x,a0,i,a1) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W117aaa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_cooo_no40_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO40_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W117aaa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no40_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO40_X1_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), W117aaa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   41] -- 
  // |-- [    0] --| W128(w,i,a0,a1) += (    1.00000000) V2(w,a2,a3,a1) D2(i,a2,a0,a3) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(a0,x,a1,j) W128(w,i,a0,a1) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W128aaa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_cooo_no41_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO41_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W128aaa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no41_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO41_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), W128aaa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   42] -- 
  // |-- [    0] --| W129(x,i,a0,a1) += (    1.00000000) V2(x,a2,a3,a1) D2(i,a2,a0,a3) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(a0,w,a1,j) W129(x,i,a0,a1) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W129aaa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_cooo_no42_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO42_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W129aaa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no42_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO42_X1_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), W129aaa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   43] -- 
  // |-- [    0] --| W132(w,a0,i,a1) += (    1.00000000) V2(w,i,a2,a1) D1(a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    4.00000000) T2(a0,x,a1,j) W132(w,a0,i,a1) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W132aaa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_cooo_no43_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO43_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W132aaa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no43_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO43_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), W132aaa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   44] -- 
  // |-- [    0] --| W133(x,a0,i,a1) += (    1.00000000) V2(x,i,a2,a1) D1(a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(a0,w,a1,j) W133(x,a0,i,a1) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W133aaa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_cooo_no44_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO44_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W133aaa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no44_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO44_X1_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), W133aaa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   45] -- 
  // |-- [    0] --| W142(w,a0,i,a1) += (    1.00000000) V2(w,a2,i,a1) D1(a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(a0,x,a1,j) W142(w,a0,i,a1) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W142aaa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_cooo_no45_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO45_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W142aaa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no45_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO45_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), W142aaa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   46] -- 
  // |-- [    0] --| W143(x,a0,i,a1) += (    1.00000000) V2(x,a2,i,a1) D1(a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(a0,w,a1,j) W143(x,a0,i,a1) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W143aaa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_cooo_no46_x0_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO46_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W143aaa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no46_x1_type1_eri_c,G_IF_SIGMA_CCOO_COOO_NO46_X1_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), W143aaa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // --  Title : sigma_ccoo_cooo
  //  >> Intermediates for the external contractions are defined here << 
  orz::DTensor W42ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W43ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W44caaa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W45caaa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W46caaa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W47caaa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W56ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W57ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W58ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W59ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W72ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W73ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W82ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W83ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
//-@type(2).declaration(end)

  //*-- Entering to take the type 0 contractions --*//
//-@type(0).contraction(begin)

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W42(c0,i) += (    1.00000000) T2(c0,a2,a1,a0) D2(i,a1,a2,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccoo_cooo_no0_x0_type0_eri_o,G_IF_SIGMA_CCOO_COOO_NO0_X0_TYPE0_ERI_O)
      (sa0, ia0, T2b.cptr(), W42ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W43(c0,i) += (    1.00000000) T2(c0,a2,a1,a0) D2(i,a1,a2,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccoo_cooo_no1_x0_type0_eri_o,G_IF_SIGMA_CCOO_COOO_NO1_X0_TYPE0_ERI_O)
      (sa0, ia0, T2b.cptr(), W43ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W44(w,i,a4,a3) += (    1.00000000) T2(w,a2,a1,a0) D3(i,a1,a4,a3,a2,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccoo_cooo_no2_x0_type0_eri_o,G_IF_SIGMA_CCOO_COOO_NO2_X0_TYPE0_ERI_O)
      (sa0, ia0, T2b.cptr(), W44caaa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W45(x,i,a4,a3) += (    1.00000000) T2(x,a2,a1,a0) D3(i,a1,a4,a3,a2,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccoo_cooo_no3_x0_type0_eri_o,G_IF_SIGMA_CCOO_COOO_NO3_X0_TYPE0_ERI_O)
      (sa0, ia0, T2b.cptr(), W45caaa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W46(w,i,a3,a4) += (    1.00000000) T2(w,a2,a1,a0) D3(i,a1,a3,a4,a2,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccoo_cooo_no4_x0_type0_eri_o,G_IF_SIGMA_CCOO_COOO_NO4_X0_TYPE0_ERI_O)
      (sa0, ia0, T2b.cptr(), W46caaa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W47(x,i,a3,a4) += (    1.00000000) T2(x,a2,a1,a0) D3(i,a4,a3,a1,a2,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccoo_cooo_no5_x0_type0_eri_o,G_IF_SIGMA_CCOO_COOO_NO5_X0_TYPE0_ERI_O)
      (sa0, ia0, T2b.cptr(), W47caaa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W56(w,a3) += (    1.00000000) T2(w,a2,a1,a0) D2(a3,a1,a2,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccoo_cooo_no6_x0_type0_eri_o,G_IF_SIGMA_CCOO_COOO_NO6_X0_TYPE0_ERI_O)
      (sa0, ia0, T2b.cptr(), W56ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W57(x,a3) += (    1.00000000) T2(x,a2,a1,a0) D2(a3,a1,a2,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccoo_cooo_no7_x0_type0_eri_o,G_IF_SIGMA_CCOO_COOO_NO7_X0_TYPE0_ERI_O)
      (sa0, ia0, T2b.cptr(), W57ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| W58(w,a3) += (    1.00000000) T2(w,a2,a1,a0) D2(a3,a1,a2,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccoo_cooo_no8_x0_type0_eri_o,G_IF_SIGMA_CCOO_COOO_NO8_X0_TYPE0_ERI_O)
      (sa0, ia0, T2b.cptr(), W58ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W59(x,a3) += (    1.00000000) T2(x,a2,a1,a0) D2(a3,a1,a2,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccoo_cooo_no9_x0_type0_eri_o,G_IF_SIGMA_CCOO_COOO_NO9_X0_TYPE0_ERI_O)
      (sa0, ia0, T2b.cptr(), W59ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| W72(c0,i) += (    1.00000000) T2(c0,a1,i,a0) D1(a1,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccoo_cooo_no10_x0_type0_eri_o,G_IF_SIGMA_CCOO_COOO_NO10_X0_TYPE0_ERI_O)
      (sa0, ia0, T2b.cptr(), W72ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| W73(c0,i) += (    1.00000000) T2(c0,a1,i,a0) D1(a1,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccoo_cooo_no11_x0_type0_eri_o,G_IF_SIGMA_CCOO_COOO_NO11_X0_TYPE0_ERI_O)
      (sa0, ia0, T2b.cptr(), W73ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| W82(c0,i) += (    1.00000000) T2(c0,a1,a0,i) D1(a1,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    FC_FUNC(g_if_sigma_ccoo_cooo_no12_x0_type0_eri_o,G_IF_SIGMA_CCOO_COOO_NO12_X0_TYPE0_ERI_O)
      (si, ii, T2b.cptr(), W82ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| W83(c0,i) += (    1.00000000) T2(c0,a1,a0,i) D1(a1,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    FC_FUNC(g_if_sigma_ccoo_cooo_no13_x0_type0_eri_o,G_IF_SIGMA_CCOO_COOO_NO13_X0_TYPE0_ERI_O)
      (si, ii, T2b.cptr(), W83ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
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
  // -- Title : sigma_ccoo_cooo

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W26(c0,j,i,a3) += (    1.00000000) T2(c0,a2,a1,a0) D3(j,a3,i,a1,a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -1.00000000) V2(a3,x,w,c0) W26(c0,j,i,a3) 
  int sa3(s_eri);
  int ia3(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <3, 2> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W26ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sa3));
    for(int sa0 = 0;sa0 < nir;++sa0){ 
    for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
      T2b = T2.get_amp2(ia0);
      FC_FUNC(g_if_sigma_ccoo_cooo_no0_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO0_X0_TYPE1_ERI_O)
        (sa0, ia0, sa3, ia3, sj, ij, T2b.cptr(), W26ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
    } // End ia0
    } // End sa0
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no0_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO0_X1_TYPE1_ERI_O)
      (sa3, ia3, sj, ij, V2_sym.cptr(), W26ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W30(w,j,a4,a3) += (    1.00000000) T2(w,a2,a1,a0) D3(j,a1,a4,a3,a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) V2(a4,a3,x,i) W30(w,j,a4,a3) 
  int sa4(s_eri);
  int ia4(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <3, 2> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W30ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sa4));
    for(int sa0 = 0;sa0 < nir;++sa0){ 
    for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
      T2b = T2.get_amp2(ia0);
      FC_FUNC(g_if_sigma_ccoo_cooo_no1_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO1_X0_TYPE1_ERI_O)
        (sa0, ia0, sa4, ia4, sj, ij, T2b.cptr(), W30ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
    } // End ia0
    } // End sa0
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no1_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO1_X1_TYPE1_ERI_O)
      (sa4, ia4, sj, ij, V2_sym.cptr(), W30ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W31(x,j,a4,a3) += (    1.00000000) T2(x,a2,a1,a0) D3(j,a1,a4,a3,a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) V2(a4,a3,w,i) W31(x,j,a4,a3) 
  int sa4(s_eri);
  int ia4(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <3, 2> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W31ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sa4));
    for(int sa0 = 0;sa0 < nir;++sa0){ 
    for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
      T2b = T2.get_amp2(ia0);
      FC_FUNC(g_if_sigma_ccoo_cooo_no2_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO2_X0_TYPE1_ERI_O)
        (sa0, ia0, sa4, ia4, sj, ij, T2b.cptr(), W31ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
    } // End ia0
    } // End sa0
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no2_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO2_X1_TYPE1_ERI_O)
      (sa4, ia4, sj, ij, V2_sym.cptr(), W31ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W32(w,j,a3,a4) += (    1.00000000) T2(w,a2,a1,a0) D3(j,a4,a3,a1,a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) V2(a4,x,i,a3) W32(w,j,a3,a4) 
  int sa4(s_eri);
  int ia4(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <3, 2> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W32ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sa4));
    for(int sa0 = 0;sa0 < nir;++sa0){ 
    for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
      T2b = T2.get_amp2(ia0);
      FC_FUNC(g_if_sigma_ccoo_cooo_no3_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO3_X0_TYPE1_ERI_O)
        (sa0, ia0, sa4, ia4, sj, ij, T2b.cptr(), W32ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
    } // End ia0
    } // End sa0
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no3_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO3_X1_TYPE1_ERI_O)
      (sa4, ia4, sj, ij, V2_sym.cptr(), W32ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W33(x,j,a3,a4) += (    1.00000000) T2(x,a2,a1,a0) D3(j,a1,a3,a4,a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) V2(a4,w,i,a3) W33(x,j,a3,a4) 
  int sa4(s_eri);
  int ia4(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <3, 2> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W33ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sa4));
    for(int sa0 = 0;sa0 < nir;++sa0){ 
    for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
      T2b = T2.get_amp2(ia0);
      FC_FUNC(g_if_sigma_ccoo_cooo_no4_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO4_X0_TYPE1_ERI_O)
        (sa0, ia0, sa4, ia4, sj, ij, T2b.cptr(), W33ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
    } // End ia0
    } // End sa0
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no4_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO4_X1_TYPE1_ERI_O)
      (sa4, ia4, sj, ij, V2_sym.cptr(), W33ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W35(c0,j,a2,i) += (    1.00000000) T2(c0,a1,i,a0) D2(j,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) V2(a2,x,w,c0) W35(c0,j,a2,i) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <3, 2> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W35ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sa2));
    for(int sa0 = 0;sa0 < nir;++sa0){ 
    for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
      T2b = T2.get_amp2(ia0);
      FC_FUNC(g_if_sigma_ccoo_cooo_no5_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO5_X0_TYPE1_ERI_O)
        (sa0, ia0, sa2, ia2, sj, ij, T2b.cptr(), W35ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
    } // End ia0
    } // End sa0
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no5_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO5_X1_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, V2_sym.cptr(), W35ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W38(c0,j,a2,i) += (    1.00000000) T2(c0,a1,a0,i) D2(j,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -1.00000000) V2(a2,x,w,c0) W38(c0,j,a2,i) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      orz::DTensor W38c_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, sj^sa2^si));
      FC_FUNC(g_if_sigma_ccoo_cooo_no6_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO6_X0_TYPE1_ERI_O)
        (sa2, ia2, si, ii, sj, ij, T2b.cptr(), W38c_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_ccoo_cooo_no6_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO6_X1_TYPE1_ERI_O)
        (sa2, ia2, si, ii, sj, ij, V2_sym.cptr(), W38c_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    2.00000000) V2(j,x,w,c0) W42(c0,i) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  FC_FUNC(g_if_sigma_ccoo_cooo_no7_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO7_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W42ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -1.00000000) V2(j,w,x,c0) W43(c0,i) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  FC_FUNC(g_if_sigma_ccoo_cooo_no8_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO8_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W43ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -2.00000000) V2(j,x,a4,a3) W44(w,i,a4,a3) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  FC_FUNC(g_if_sigma_ccoo_cooo_no9_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO9_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W44caaa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    1.00000000) V2(j,w,a4,a3) W45(x,i,a4,a3) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  FC_FUNC(g_if_sigma_ccoo_cooo_no10_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO10_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W45caaa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    1.00000000) V2(j,a3,x,a4) W46(w,i,a3,a4) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  FC_FUNC(g_if_sigma_ccoo_cooo_no11_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO11_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W46caaa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    1.00000000) V2(j,a3,w,a4) W47(x,i,a3,a4) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  FC_FUNC(g_if_sigma_ccoo_cooo_no12_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO12_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W47caaa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| W49(c0,i,a2,j) += (    1.00000000) T2(a1,c0,a0,j) D2(i,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -1.00000000) V2(a2,x,w,c0) W49(c0,i,a2,j) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W49ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sa2^sj));
    FC_FUNC(g_if_sigma_ccoo_cooo_no13_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO13_X0_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, T2b.cptr(), W49ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no13_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO13_X1_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, V2_sym.cptr(), W49ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   14] -- 
  // |-- [    0] --| W52(c0,i,a2,j) += (    1.00000000) T2(c0,a1,a0,j) D2(i,a0,a1,a2) 
  // |-- [    1] --| S2(w,x,i,j) += (   -1.00000000) V2(a2,x,w,c0) W52(c0,i,a2,j) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W52ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sa2^sj));
    FC_FUNC(g_if_sigma_ccoo_cooo_no14_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO14_X0_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, T2b.cptr(), W52ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no14_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO14_X1_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, V2_sym.cptr(), W52ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   15] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -2.00000000) V2(j,x,i,a3) W56(w,a3) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  FC_FUNC(g_if_sigma_ccoo_cooo_no15_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO15_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W56ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   16] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    1.00000000) V2(j,w,i,a3) W57(x,a3) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  FC_FUNC(g_if_sigma_ccoo_cooo_no16_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO16_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W57ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   17] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    1.00000000) V2(j,a3,x,i) W58(w,a3) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  FC_FUNC(g_if_sigma_ccoo_cooo_no17_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO17_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W58ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   18] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -2.00000000) V2(j,a3,w,i) W59(x,a3) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  FC_FUNC(g_if_sigma_ccoo_cooo_no18_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO18_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W59ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   19] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    2.00000000) V2(j,w,x,c0) W72(c0,i) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  FC_FUNC(g_if_sigma_ccoo_cooo_no19_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO19_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W72ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   20] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -4.00000000) V2(j,x,w,c0) W73(c0,i) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  FC_FUNC(g_if_sigma_ccoo_cooo_no20_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO20_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W73ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   21] -- 
  // |-- [    0] --| W74(x,a1,a0,j) += (    1.00000000) V2(j,x,a3,a2) D2(a3,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    4.00000000) T2(w,a1,i,a0) W74(x,a1,a0,j) 
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
    orz::DTensor W74ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sa0^sj));
    FC_FUNC(g_if_sigma_ccoo_cooo_no21_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO21_X0_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, V2_sym.cptr(), W74ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_cooo_no21_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO21_X1_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, T2b.cptr(), W74ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   22] -- 
  // |-- [    0] --| W75(w,a1,a0,j) += (    1.00000000) V2(j,w,a3,a2) D2(a3,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(x,a1,i,a0) W75(w,a1,a0,j) 
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
    orz::DTensor W75ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sa0^sj));
    FC_FUNC(g_if_sigma_ccoo_cooo_no22_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO22_X0_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, V2_sym.cptr(), W75ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_cooo_no22_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO22_X1_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, T2b.cptr(), W75ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   23] -- 
  // |-- [    0] --| W76(w,a1,a0,j) += (    1.00000000) V2(j,a2,w,a3) D2(a3,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(x,a0,i,a1) W76(w,a1,a0,j) 
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
    orz::DTensor W76ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sa1^sj));
    FC_FUNC(g_if_sigma_ccoo_cooo_no23_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO23_X0_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, V2_sym.cptr(), W76ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_cooo_no23_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO23_X1_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, T2b.cptr(), W76ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   24] -- 
  // |-- [    0] --| W77(x,a1,a0,j) += (    1.00000000) V2(j,a2,x,a3) D2(a3,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(w,a0,i,a1) W77(x,a1,a0,j) 
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
    orz::DTensor W77ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sa1^sj));
    FC_FUNC(g_if_sigma_ccoo_cooo_no24_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO24_X0_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, V2_sym.cptr(), W77ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_cooo_no24_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO24_X1_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, T2b.cptr(), W77ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   25] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    2.00000000) V2(j,x,w,c0) W82(c0,i) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  FC_FUNC(g_if_sigma_ccoo_cooo_no25_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO25_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W82ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   26] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -1.00000000) V2(j,w,x,c0) W83(c0,i) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  FC_FUNC(g_if_sigma_ccoo_cooo_no26_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO26_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W83ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   27] -- 
  // |-- [    0] --| W84(x,a1,a0,j) += (    1.00000000) V2(j,x,a3,a2) D2(a3,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(w,a1,a0,i) W84(x,a1,a0,j) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W84caa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
  FC_FUNC(g_if_sigma_ccoo_cooo_no27_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO27_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W84caa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    FC_FUNC(g_if_sigma_ccoo_cooo_no27_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO27_X1_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), W84caa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   28] -- 
  // |-- [    0] --| W85(w,a1,a0,j) += (    1.00000000) V2(j,w,a3,a2) D2(a3,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(x,a1,a0,i) W85(w,a1,a0,j) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W85caa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
  FC_FUNC(g_if_sigma_ccoo_cooo_no28_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO28_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W85caa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    FC_FUNC(g_if_sigma_ccoo_cooo_no28_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO28_X1_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), W85caa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   29] -- 
  // |-- [    0] --| W86(x,a1,a0,j) += (    1.00000000) V2(j,a2,x,a3) D2(a3,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(w,a0,a1,i) W86(x,a1,a0,j) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W86caa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
  FC_FUNC(g_if_sigma_ccoo_cooo_no29_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO29_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W86caa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    FC_FUNC(g_if_sigma_ccoo_cooo_no29_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO29_X1_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), W86caa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   30] -- 
  // |-- [    0] --| W87(w,a1,a0,j) += (    1.00000000) V2(j,a2,w,a3) D2(a3,a0,a1,a2) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(x,a0,a1,i) W87(w,a1,a0,j) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W87caa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
  FC_FUNC(g_if_sigma_ccoo_cooo_no30_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO30_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W87caa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    FC_FUNC(g_if_sigma_ccoo_cooo_no30_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO30_X1_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), W87caa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   31] -- 
  // |-- [    0] --| W92(x,j,i,a1,a0,a2) += (    1.00000000) V2(a2,a4,x,a3) D3(j,a3,i,a4,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(a1,w,a0,a2) W92(x,j,i,a1,a0,a2) 
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
    orz::DTensor W92caaa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, sj^sa2));
    FC_FUNC(g_if_sigma_ccoo_cooo_no31_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO31_X0_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, V2_sym.cptr(), W92caaa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_cooo_no31_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO31_X1_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, T2b.cptr(), W92caaa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   32] -- 
  // |-- [    0] --| W93(w,j,i,a1,a0,a2) += (    1.00000000) V2(a2,a4,w,a3) D3(j,a4,i,a3,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(a1,x,a0,a2) W93(w,j,i,a1,a0,a2) 
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
    orz::DTensor W93caaa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, sj^sa2));
    FC_FUNC(g_if_sigma_ccoo_cooo_no32_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO32_X0_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, V2_sym.cptr(), W93caaa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_cooo_no32_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO32_X1_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, T2b.cptr(), W93caaa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   33] -- 
  // |-- [    0] --| W94(w,j,a3,a2) += (    1.00000000) T2(a1,w,a0,a2) D2(j,a3,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) V2(a2,a3,x,i) W94(w,j,a3,a2) 
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
    orz::DTensor W94ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sa2));
    FC_FUNC(g_if_sigma_ccoo_cooo_no33_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO33_X0_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, T2b.cptr(), W94ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_cooo_no33_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO33_X1_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, V2_sym.cptr(), W94ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   34] -- 
  // |-- [    0] --| W95(x,j,a3,a2) += (    1.00000000) T2(a1,x,a0,a2) D2(j,a3,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) V2(a2,a3,w,i) W95(x,j,a3,a2) 
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
    orz::DTensor W95ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sa2));
    FC_FUNC(g_if_sigma_ccoo_cooo_no34_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO34_X0_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, T2b.cptr(), W95ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_cooo_no34_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO34_X1_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, V2_sym.cptr(), W95ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   35] -- 
  // |-- [    0] --| W96(x,j,a0,a1) += (    1.00000000) V2(a1,a3,x,a2) D2(j,a2,a0,a3) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(a0,w,i,a1) W96(x,j,a0,a1) 
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
    orz::DTensor W96ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sa1));
    FC_FUNC(g_if_sigma_ccoo_cooo_no35_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO35_X0_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, V2_sym.cptr(), W96ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_cooo_no35_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO35_X1_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, T2b.cptr(), W96ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   36] -- 
  // |-- [    0] --| W97(w,j,a0,a1) += (    1.00000000) V2(a1,a3,w,a2) D2(j,a3,a0,a2) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(a0,x,i,a1) W97(w,j,a0,a1) 
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
    orz::DTensor W97ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sa1));
    FC_FUNC(g_if_sigma_ccoo_cooo_no36_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO36_X0_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, V2_sym.cptr(), W97ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_cooo_no36_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO36_X1_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, T2b.cptr(), W97ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   37] -- 
  // |-- [    0] --| W98(x,j,a3,a2) += (    1.00000000) T2(a1,x,a0,a2) D2(j,a3,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) V2(a2,i,w,a3) W98(x,j,a3,a2) 
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
    orz::DTensor W98ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sa2));
    FC_FUNC(g_if_sigma_ccoo_cooo_no37_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO37_X0_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, T2b.cptr(), W98ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_cooo_no37_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO37_X1_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, V2_sym.cptr(), W98ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   38] -- 
  // |-- [    0] --| W99(w,j,a3,a2) += (    1.00000000) T2(a1,w,a0,a2) D2(j,a3,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) V2(a2,i,x,a3) W99(w,j,a3,a2) 
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
    orz::DTensor W99ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sa2));
    FC_FUNC(g_if_sigma_ccoo_cooo_no38_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO38_X0_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, T2b.cptr(), W99ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_cooo_no38_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO38_X1_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, V2_sym.cptr(), W99ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   39] -- 
  // |-- [    0] --| W100(w,i,a3,a2) += (    1.00000000) T2(a1,w,a0,a2) D2(i,a3,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) V2(a2,a3,x,j) W100(w,i,a3,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia2);
  orz::DTensor W100caa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa2));
  FC_FUNC(g_if_sigma_ccoo_cooo_no39_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO39_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), W100caa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no39_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO39_X1_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, V2_sym.cptr(), W100caa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   40] -- 
  // |-- [    0] --| W101(x,i,a3,a2) += (    1.00000000) T2(a1,x,a0,a2) D2(i,a3,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) V2(a2,a3,w,j) W101(x,i,a3,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia2);
  orz::DTensor W101caa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa2));
  FC_FUNC(g_if_sigma_ccoo_cooo_no40_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO40_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), W101caa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no40_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO40_X1_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, V2_sym.cptr(), W101caa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   41] -- 
  // |-- [    0] --| W104(x,i,a3,a2) += (    1.00000000) T2(a1,x,a0,a2) D2(i,a3,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) V2(a2,j,w,a3) W104(x,i,a3,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia2);
  orz::DTensor W104caa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa2));
  FC_FUNC(g_if_sigma_ccoo_cooo_no41_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO41_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), W104caa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no41_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO41_X1_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, V2_sym.cptr(), W104caa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   42] -- 
  // |-- [    0] --| W105(w,i,a3,a2) += (    1.00000000) T2(a1,w,a0,a2) D2(i,a3,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) V2(a2,j,x,a3) W105(w,i,a3,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia2);
  orz::DTensor W105caa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa2));
  FC_FUNC(g_if_sigma_ccoo_cooo_no42_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO42_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), W105caa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no42_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO42_X1_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, V2_sym.cptr(), W105caa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   43] -- 
  // |-- [    0] --| W108(x,a2) += (    1.00000000) T2(a1,x,a0,a2) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    4.00000000) V2(a2,j,w,i) W108(x,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia2);
  orz::DTensor W108c_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, sa2));
  FC_FUNC(g_if_sigma_ccoo_cooo_no43_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO43_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), W108c_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no43_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO43_X1_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, V2_sym.cptr(), W108c_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   44] -- 
  // |-- [    0] --| W109(w,a2) += (    1.00000000) T2(a1,w,a0,a2) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) V2(a2,j,x,i) W109(w,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia2);
  orz::DTensor W109c_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, sa2));
  FC_FUNC(g_if_sigma_ccoo_cooo_no44_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO44_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), W109c_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no44_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO44_X1_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, V2_sym.cptr(), W109c_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   45] -- 
  // |-- [    0] --| W110(x,a0,j,a1) += (    1.00000000) V2(j,x,a2,a1) D1(a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(w,a0,a1,i) W110(x,a0,j,a1) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W110caa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
  FC_FUNC(g_if_sigma_ccoo_cooo_no45_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO45_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W110caa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    FC_FUNC(g_if_sigma_ccoo_cooo_no45_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO45_X1_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), W110caa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   46] -- 
  // |-- [    0] --| W111(w,a0,j,a1) += (    1.00000000) V2(j,w,a2,a1) D1(a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(x,a0,a1,i) W111(w,a0,j,a1) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W111caa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
  FC_FUNC(g_if_sigma_ccoo_cooo_no46_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO46_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W111caa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    FC_FUNC(g_if_sigma_ccoo_cooo_no46_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO46_X1_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), W111caa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   47] -- 
  // |-- [    0] --| W112(w,a0,j,a1) += (    1.00000000) V2(j,a1,w,a2) D1(a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(x,a0,a1,i) W112(w,a0,j,a1) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W112caa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
  FC_FUNC(g_if_sigma_ccoo_cooo_no47_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO47_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W112caa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    FC_FUNC(g_if_sigma_ccoo_cooo_no47_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO47_X1_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), W112caa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   48] -- 
  // |-- [    0] --| W113(x,a0,j,a1) += (    1.00000000) V2(j,a1,x,a2) D1(a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(w,a0,a1,i) W113(x,a0,j,a1) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W113caa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
  FC_FUNC(g_if_sigma_ccoo_cooo_no48_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO48_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W113caa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    FC_FUNC(g_if_sigma_ccoo_cooo_no48_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO48_X1_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), W113caa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   49] -- 
  // |-- [    0] --| W114(w,a2) += (    1.00000000) T2(a1,w,a0,a2) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    4.00000000) V2(a2,i,x,j) W114(w,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia2);
  orz::DTensor W114c_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, sa2));
  FC_FUNC(g_if_sigma_ccoo_cooo_no49_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO49_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), W114c_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no49_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO49_X1_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, V2_sym.cptr(), W114c_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   50] -- 
  // |-- [    0] --| W115(x,a2) += (    1.00000000) T2(a1,x,a0,a2) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) V2(a2,i,w,j) W115(x,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia2);
  orz::DTensor W115c_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, sa2));
  FC_FUNC(g_if_sigma_ccoo_cooo_no50_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO50_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), W115c_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no50_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO50_X1_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, V2_sym.cptr(), W115c_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   51] -- 
  // |-- [    0] --| W118(x,j,i,a1,a0,a2) += (    1.00000000) V2(a2,a4,x,a3) D3(j,a3,i,a0,a1,a4) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(w,a1,a0,a2) W118(x,j,i,a1,a0,a2) 
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
    orz::DTensor W118caaa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, sj^sa2));
    FC_FUNC(g_if_sigma_ccoo_cooo_no51_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO51_X0_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, V2_sym.cptr(), W118caaa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_cooo_no51_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO51_X1_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, T2b.cptr(), W118caaa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   52] -- 
  // |-- [    0] --| W119(w,j,i,a1,a0,a2) += (    1.00000000) V2(a2,a4,w,a3) D3(j,a0,i,a3,a1,a4) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(x,a1,a0,a2) W119(w,j,i,a1,a0,a2) 
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
    orz::DTensor W119caaa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, sj^sa2));
    FC_FUNC(g_if_sigma_ccoo_cooo_no52_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO52_X0_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, V2_sym.cptr(), W119caaa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_cooo_no52_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO52_X1_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, T2b.cptr(), W119caaa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   53] -- 
  // |-- [    0] --| W120(w,j,a3,a2) += (    1.00000000) T2(w,a1,a0,a2) D2(j,a0,a1,a3) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) V2(a2,a3,x,i) W120(w,j,a3,a2) 
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
    orz::DTensor W120ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sa2));
    FC_FUNC(g_if_sigma_ccoo_cooo_no53_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO53_X0_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, T2b.cptr(), W120ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_cooo_no53_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO53_X1_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, V2_sym.cptr(), W120ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   54] -- 
  // |-- [    0] --| W121(x,j,a3,a2) += (    1.00000000) T2(x,a1,a0,a2) D2(j,a0,a1,a3) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) V2(a2,a3,w,i) W121(x,j,a3,a2) 
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
    orz::DTensor W121ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sa2));
    FC_FUNC(g_if_sigma_ccoo_cooo_no54_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO54_X0_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, T2b.cptr(), W121ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_cooo_no54_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO54_X1_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, V2_sym.cptr(), W121ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   55] -- 
  // |-- [    0] --| W122(w,j,a0,a1) += (    1.00000000) V2(a1,a3,w,a2) D2(j,a2,a0,a3) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(x,a0,i,a1) W122(w,j,a0,a1) 
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
    orz::DTensor W122ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sa1));
    FC_FUNC(g_if_sigma_ccoo_cooo_no55_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO55_X0_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, V2_sym.cptr(), W122ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_cooo_no55_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO55_X1_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, T2b.cptr(), W122ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   56] -- 
  // |-- [    0] --| W123(x,j,a0,a1) += (    1.00000000) V2(a1,a3,x,a2) D2(j,a2,a0,a3) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(w,a0,i,a1) W123(x,j,a0,a1) 
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
    orz::DTensor W123ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sa1));
    FC_FUNC(g_if_sigma_ccoo_cooo_no56_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO56_X0_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, V2_sym.cptr(), W123ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_cooo_no56_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO56_X1_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, T2b.cptr(), W123ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   57] -- 
  // |-- [    0] --| W124(w,j,a3,a2) += (    1.00000000) T2(w,a1,a0,a2) D2(j,a3,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) V2(a2,i,x,a3) W124(w,j,a3,a2) 
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
    orz::DTensor W124ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sa2));
    FC_FUNC(g_if_sigma_ccoo_cooo_no57_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO57_X0_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, T2b.cptr(), W124ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_cooo_no57_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO57_X1_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, V2_sym.cptr(), W124ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   58] -- 
  // |-- [    0] --| W125(x,j,a3,a2) += (    1.00000000) T2(x,a1,a0,a2) D2(j,a0,a1,a3) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) V2(a2,i,w,a3) W125(x,j,a3,a2) 
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
    orz::DTensor W125ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sa2));
    FC_FUNC(g_if_sigma_ccoo_cooo_no58_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO58_X0_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, T2b.cptr(), W125ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_cooo_no58_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO58_X1_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, V2_sym.cptr(), W125ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   59] -- 
  // |-- [    0] --| W126(w,i,a3,a2) += (    1.00000000) T2(w,a1,a0,a2) D2(i,a0,a1,a3) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) V2(a2,a3,x,j) W126(w,i,a3,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia2);
  orz::DTensor W126caa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa2));
  FC_FUNC(g_if_sigma_ccoo_cooo_no59_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO59_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), W126caa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no59_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO59_X1_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, V2_sym.cptr(), W126caa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   60] -- 
  // |-- [    0] --| W127(x,i,a3,a2) += (    1.00000000) T2(x,a1,a0,a2) D2(i,a0,a1,a3) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) V2(a2,a3,w,j) W127(x,i,a3,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia2);
  orz::DTensor W127caa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa2));
  FC_FUNC(g_if_sigma_ccoo_cooo_no60_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO60_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), W127caa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no60_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO60_X1_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, V2_sym.cptr(), W127caa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   61] -- 
  // |-- [    0] --| W130(w,i,a3,a2) += (    1.00000000) T2(w,a1,a0,a2) D2(i,a0,a1,a3) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) V2(a2,j,x,a3) W130(w,i,a3,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia2);
  orz::DTensor W130caa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa2));
  FC_FUNC(g_if_sigma_ccoo_cooo_no61_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO61_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), W130caa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no61_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO61_X1_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, V2_sym.cptr(), W130caa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   62] -- 
  // |-- [    0] --| W131(x,i,a3,a2) += (    1.00000000) T2(x,a1,a0,a2) D2(i,a3,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) V2(a2,j,w,a3) W131(x,i,a3,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia2);
  orz::DTensor W131caa_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa2));
  FC_FUNC(g_if_sigma_ccoo_cooo_no62_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO62_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), W131caa_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no62_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO62_X1_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, V2_sym.cptr(), W131caa_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   63] -- 
  // |-- [    0] --| W134(w,a2) += (    1.00000000) T2(w,a1,a0,a2) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) V2(a2,j,x,i) W134(w,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia2);
  orz::DTensor W134c_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, sa2));
  FC_FUNC(g_if_sigma_ccoo_cooo_no63_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO63_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), W134c_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no63_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO63_X1_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, V2_sym.cptr(), W134c_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   64] -- 
  // |-- [    0] --| W135(x,a2) += (    1.00000000) T2(x,a1,a0,a2) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) V2(a2,j,w,i) W135(x,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia2);
  orz::DTensor W135c_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, sa2));
  FC_FUNC(g_if_sigma_ccoo_cooo_no64_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO64_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), W135c_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no64_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO64_X1_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, V2_sym.cptr(), W135c_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   65] -- 
  // |-- [    0] --| W136(x,a0,j,a1) += (    1.00000000) V2(j,x,a2,a1) D1(a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    4.00000000) T2(w,a0,i,a1) W136(x,a0,j,a1) 
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
    orz::DTensor W136ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sa1));
    FC_FUNC(g_if_sigma_ccoo_cooo_no65_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO65_X0_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, V2_sym.cptr(), W136ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_cooo_no65_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO65_X1_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, T2b.cptr(), W136ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   66] -- 
  // |-- [    0] --| W137(w,a0,j,a1) += (    1.00000000) V2(j,w,a2,a1) D1(a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(x,a0,i,a1) W137(w,a0,j,a1) 
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
    orz::DTensor W137ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sa1));
    FC_FUNC(g_if_sigma_ccoo_cooo_no66_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO66_X0_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, V2_sym.cptr(), W137ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_cooo_no66_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO66_X1_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, T2b.cptr(), W137ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   67] -- 
  // |-- [    0] --| W138(w,a0,j,a1) += (    1.00000000) V2(j,a1,w,a2) D1(a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) T2(x,a0,i,a1) W138(w,a0,j,a1) 
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
    orz::DTensor W138ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sa1));
    FC_FUNC(g_if_sigma_ccoo_cooo_no67_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO67_X0_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, V2_sym.cptr(), W138ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_cooo_no67_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO67_X1_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, T2b.cptr(), W138ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   68] -- 
  // |-- [    0] --| W139(x,a0,j,a1) += (    1.00000000) V2(j,a1,x,a2) D1(a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(w,a0,i,a1) W139(x,a0,j,a1) 
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
    orz::DTensor W139ca_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sa1));
    FC_FUNC(g_if_sigma_ccoo_cooo_no68_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO68_X0_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, V2_sym.cptr(), W139ca_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_cooo_no68_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO68_X1_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, T2b.cptr(), W139ca_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   69] -- 
  // |-- [    0] --| W140(w,a2) += (    1.00000000) T2(w,a1,a0,a2) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) V2(a2,i,x,j) W140(w,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia2);
  orz::DTensor W140c_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, sa2));
  FC_FUNC(g_if_sigma_ccoo_cooo_no69_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO69_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), W140c_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no69_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO69_X1_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, V2_sym.cptr(), W140c_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   70] -- 
  // |-- [    0] --| W141(x,a2) += (    1.00000000) T2(x,a1,a0,a2) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    1.00000000) V2(a2,i,w,j) W141(x,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia2);
  orz::DTensor W141c_sigma_ccoo_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, sa2));
  FC_FUNC(g_if_sigma_ccoo_cooo_no70_x0_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO70_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), W141c_sigma_ccoo_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_cooo_no70_x1_type1_eri_o,G_IF_SIGMA_CCOO_COOO_NO70_X1_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, V2_sym.cptr(), W141c_sigma_ccoo_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- Title : sigma_ccoo_cooo

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    1.00000000) T2(w,a2,a1,a0) C5(a2,a0,i,a1,j,x) 
  int sx(s_eri);
  int ix(i_eri);
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
      FC_FUNC(g_if_sigma_ccoo_cooo_no0_x0_type1_d4c_c,G_IF_SIGMA_CCOO_COOO_NO0_X0_TYPE1_D4C_C)
        (sa0, ia0, sj, ij, sx, ix, C5.cptr(), T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // |-- [    0] --| S2(w,x,i,j) += (    1.00000000) T2(x,a2,a1,a0) C5(a2,a0,j,a1,i,w) 
  int sw(s_eri);
  int iw(i_eri);
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
      FC_FUNC(g_if_sigma_ccoo_cooo_no1_x0_type1_d4c_c,G_IF_SIGMA_CCOO_COOO_NO1_X0_TYPE1_D4C_C)
        (sa0, ia0, sj, ij, sw, iw, C5.cptr(), T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  orz::mr::Femto::my_timer.push_back(boost::make_tuple("sigma_ccoo_cooo", time_sigma.elapsed_cputime(), time_sigma.elapsed_wallclocktime()));
  orz::mr::Femto::file_timing << "* " << boost::format("%20s : %10.7f %10.7f ") % orz::mr::Femto::my_timer.back().get<0>() % orz::mr::Femto::my_timer.back().get<1>() % orz::mr::Femto::my_timer.back().get<2>() << endl;
  flush(orz::mr::Femto::file_timing);

  return retval; 
} 
