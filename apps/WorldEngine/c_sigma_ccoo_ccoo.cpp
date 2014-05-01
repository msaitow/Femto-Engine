                                                                                
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
#include <sci/icmr/Femto/elems/c_sigma_ccoo_ccoo.h>                                  
                                                                                
using std::cout;                                                                
using std::endl;                                                                
                                                                                
#define FLOPCOUNT                                                               
                                                                                
using orz::mr::Femto::Fc0;                                                      
using orz::mr::Femto::h1_int;                                                   
using orz::mr::Femto::h6_int;                                                   
                                                                                
//    o__ __o__/_                            o                         
//   <|    v                                <|>                        
//   < >                                    < >                        
//    |         o__  __o   \o__ __o__ __o    |        o__ __o         
//    o__/_    /v      |>   |     |     |>   o__/_   /v     v\        
//    |       />      //   / \   / \   / \   |      />       <\    
//   <o>      \o    o/     \o/   \o/   \o/   |      \         /   
//    |        v\  /v __o   |     |     |    o       o       o        
//   / \        <\/> __/>  / \   / \   / \   <\__    <\__ __/>  

//                                   Generated date : Sun Apr 20 10:26:22 2014

                                                                                
// ***************************************************************************  
// orz::mr::mrci                                                                
// ***************************************************************************  
									    /*!        
   @brief CT input                                                              
									     */        
                                                                                
orz::mr::BareAmpPack orz::mr::Femto::sigma_ccoo_ccoo(const orz::mr::Input &ctinp,                                    
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
  // |-- [    0] --| S2(w,x,i,j) += (    2.00000000) Fc0 T2(x,w,a1,a0) D2(j,a1,i,a0) 
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
      FC_FUNC(g_if_sigma_ccoo_ccoo_no0_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO0_X0_TYPE0_NOERI)
        (sa0, ia0, sj, ij, &Fc0, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // |-- [    0] --| S2(x,w,j,i) += (   -4.00000000) Fc0 T2(x,w,a0,i) D1(j,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    S2b = orz::DTensor(retval.namps_iamp()[ii]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no1_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO1_X0_TYPE0_NOERI)
      (si, ii, &Fc0, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    retval.acc_amp2(ii, S2b);
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| S2(x,w,j,i) += (    2.00000000) Fc0 T2(w,x,a0,i) D1(j,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    S2b = orz::DTensor(retval.namps_iamp()[ii]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no2_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO2_X0_TYPE0_NOERI)
      (si, ii, &Fc0, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    retval.acc_amp2(ii, S2b);
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    2.00000000) Fc0 T2(x,w,a0,j) D1(i,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no3_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO3_X0_TYPE0_NOERI)
      (sj, ij, &Fc0, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -4.00000000) Fc0 T2(w,x,a0,j) D1(i,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no4_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO4_X0_TYPE0_NOERI)
      (sj, ij, &Fc0, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    8.00000000) Fc0 T2(w,x,i,j) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no5_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO5_X0_TYPE0_NOERI)
      (sj, ij, &Fc0, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -4.00000000) Fc0 T2(x,w,i,j) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no6_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO6_X0_TYPE0_NOERI)
      (sj, ij, &Fc0, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W0(x,c0,j,i) += (    1.00000000) T2(x,c0,a1,a0) D2(j,a1,i,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) Fc1(w,c0) W0(x,c0,j,i) 
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W0cca_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcca(symblockinfo, sj));
    for(int sa0 = 0;sa0 < nir;++sa0){ 
    for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
      T2b = T2.get_amp2(ia0);
      FC_FUNC(g_if_sigma_ccoo_ccoo_no7_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO7_X0_TYPE0_NOERI)
        (sa0, ia0, sj, ij, T2b.cptr(), W0cca_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
    } // End ia0
    } // End sa0
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no7_x1_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO7_X1_TYPE0_NOERI)
      (sj, ij, W0cca_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| W1(w,c0,j,i) += (    1.00000000) T2(w,c0,a0,a1) D2(j,a1,i,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) Fc1(x,c0) W1(w,c0,j,i) 
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W1cca_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcca(symblockinfo, sj));
    for(int sa1 = 0;sa1 < nir;++sa1){ 
    for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
      T2b = T2.get_amp2(ia1);
      FC_FUNC(g_if_sigma_ccoo_ccoo_no8_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO8_X0_TYPE0_NOERI)
        (sa1, ia1, sj, ij, T2b.cptr(), W1cca_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
    } // End ia1
    } // End sa1
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no8_x1_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO8_X1_TYPE0_NOERI)
      (sj, ij, W1cca_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W2(j,i,a1,a0) += (    1.00000000) D3(j,a1,i,a0,a3,a2) Fc1(a3,a2) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) T2(x,w,a1,a0) W2(j,i,a1,a0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      orz::DTensor W2aa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sa0));
      FC_FUNC(g_if_sigma_ccoo_ccoo_no9_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO9_X0_TYPE0_NOERI)
        (sa0, ia0, sj, ij, W2aa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_ccoo_ccoo_no9_x1_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO9_X1_TYPE0_NOERI)
        (sa0, ia0, sj, ij, T2b.cptr(), W2aa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- No. [   10] -- 
  // |-- [    0] --| W3(j,a1,a0,i) += (    1.00000000) D2(j,a1,a2,a0) Fc1(i,a2) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) T2(x,w,a1,a0) W3(j,a1,a0,i) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      orz::DTensor W3aa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sa0));
      FC_FUNC(g_if_sigma_ccoo_ccoo_no10_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO10_X0_TYPE0_NOERI)
        (sa0, ia0, sj, ij, W3aa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_ccoo_ccoo_no10_x1_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO10_X1_TYPE0_NOERI)
        (sa0, ia0, sj, ij, T2b.cptr(), W3aa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // |-- [    0] --| W4(w,c0,j,i) += (    1.00000000) T2(c0,w,a0,i) D1(j,a0) 
  // |-- [    1] --| S2(x,w,j,i) += (    4.00000000) Fc1(x,c0) W4(w,c0,j,i) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    orz::DTensor W4cca_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcca(symblockinfo, si));
    S2b = orz::DTensor(retval.namps_iamp()[ii]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no11_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO11_X0_TYPE0_NOERI)
      (si, ii, T2b.cptr(), W4cca_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no11_x1_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO11_X1_TYPE0_NOERI)
      (si, ii, W4cca_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    retval.acc_amp2(ii, S2b);
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| W5(x,c0,j,i) += (    1.00000000) T2(c0,x,a0,i) D1(j,a0) 
  // |-- [    1] --| S2(x,w,j,i) += (   -2.00000000) Fc1(w,c0) W5(x,c0,j,i) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    orz::DTensor W5cca_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcca(symblockinfo, si));
    S2b = orz::DTensor(retval.namps_iamp()[ii]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no12_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO12_X0_TYPE0_NOERI)
      (si, ii, T2b.cptr(), W5cca_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no12_x1_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO12_X1_TYPE0_NOERI)
      (si, ii, W5cca_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    retval.acc_amp2(ii, S2b);
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| W6(x,c0,j,i) += (    1.00000000) T2(x,c0,a0,i) D1(j,a0) 
  // |-- [    1] --| S2(x,w,j,i) += (    4.00000000) Fc1(w,c0) W6(x,c0,j,i) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    orz::DTensor W6cca_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcca(symblockinfo, si));
    S2b = orz::DTensor(retval.namps_iamp()[ii]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no13_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO13_X0_TYPE0_NOERI)
      (si, ii, T2b.cptr(), W6cca_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no13_x1_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO13_X1_TYPE0_NOERI)
      (si, ii, W6cca_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    retval.acc_amp2(ii, S2b);
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   14] -- 
  // |-- [    0] --| W7(w,c0,j,i) += (    1.00000000) T2(w,c0,a0,i) D1(j,a0) 
  // |-- [    1] --| S2(x,w,j,i) += (   -2.00000000) Fc1(x,c0) W7(w,c0,j,i) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    orz::DTensor W7cca_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcca(symblockinfo, si));
    S2b = orz::DTensor(retval.namps_iamp()[ii]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no14_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO14_X0_TYPE0_NOERI)
      (si, ii, T2b.cptr(), W7cca_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no14_x1_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO14_X1_TYPE0_NOERI)
      (si, ii, W7cca_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    retval.acc_amp2(ii, S2b);
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   15] -- 
  // |-- [    0] --| W8(j,a0) += (    1.00000000) D2(j,a0,a2,a1) Fc1(a2,a1) 
  // |-- [    1] --| S2(x,w,j,i) += (   -4.00000000) T2(x,w,a0,i) W8(j,a0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W8aa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, 0));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no15_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO15_X0_TYPE0_NOERI)
    (W8aa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    S2b = orz::DTensor(retval.namps_iamp()[ii]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no15_x1_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO15_X1_TYPE0_NOERI)
      (si, ii, T2b.cptr(), W8aa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    retval.acc_amp2(ii, S2b);
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   16] -- 
  // |-- [    0] --| W9(j,a0) += (    1.00000000) D2(j,a0,a2,a1) Fc1(a2,a1) 
  // |-- [    1] --| S2(x,w,j,i) += (    2.00000000) T2(w,x,a0,i) W9(j,a0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W9aa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, 0));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no16_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO16_X0_TYPE0_NOERI)
    (W9aa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    S2b = orz::DTensor(retval.namps_iamp()[ii]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no16_x1_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO16_X1_TYPE0_NOERI)
      (si, ii, T2b.cptr(), W9aa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    retval.acc_amp2(ii, S2b);
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   17] -- 
  // |-- [    0] --| W10(i,a1,a0,j) += (    1.00000000) D2(i,a1,a2,a0) Fc1(j,a2) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) T2(x,w,a0,a1) W10(i,a1,a0,j) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      orz::DTensor W10aa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa1^sj));
      FC_FUNC(g_if_sigma_ccoo_ccoo_no17_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO17_X0_TYPE0_NOERI)
        (sa1, ia1, sj, ij, W10aa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
      FC_FUNC(g_if_sigma_ccoo_ccoo_no17_x1_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO17_X1_TYPE0_NOERI)
        (sa1, ia1, sj, ij, T2b.cptr(), W10aa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- No. [   18] -- 
  // |-- [    0] --| W11(w,c0,i,j) += (    1.00000000) T2(c0,w,a0,j) D1(i,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) Fc1(x,c0) W11(w,c0,i,j) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W11cca_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcca(symblockinfo, sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no18_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO18_X0_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W11cca_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no18_x1_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO18_X1_TYPE0_NOERI)
      (sj, ij, W11cca_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   19] -- 
  // |-- [    0] --| W12(x,c0,i,j) += (    1.00000000) T2(c0,x,a0,j) D1(i,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    4.00000000) Fc1(w,c0) W12(x,c0,i,j) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W12cca_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcca(symblockinfo, sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no19_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO19_X0_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W12cca_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no19_x1_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO19_X1_TYPE0_NOERI)
      (sj, ij, W12cca_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   20] -- 
  // |-- [    0] --| W13(x,c0,i,j) += (    1.00000000) T2(x,c0,a0,j) D1(i,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) Fc1(w,c0) W13(x,c0,i,j) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W13cca_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcca(symblockinfo, sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no20_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO20_X0_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W13cca_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no20_x1_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO20_X1_TYPE0_NOERI)
      (sj, ij, W13cca_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   21] -- 
  // |-- [    0] --| W14(w,c0,i,j) += (    1.00000000) T2(w,c0,a0,j) D1(i,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    4.00000000) Fc1(x,c0) W14(w,c0,i,j) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W14cca_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcca(symblockinfo, sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no21_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO21_X0_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W14cca_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no21_x1_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO21_X1_TYPE0_NOERI)
      (sj, ij, W14cca_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   22] -- 
  // |-- [    0] --| W15(i,a0) += (    1.00000000) D2(i,a0,a2,a1) Fc1(a2,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) T2(x,w,a0,j) W15(i,a0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W15aa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, 0));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no22_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO22_X0_TYPE0_NOERI)
    (W15aa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no22_x1_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO22_X1_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W15aa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   23] -- 
  // |-- [    0] --| W16(i,a0) += (    1.00000000) D2(i,a0,a2,a1) Fc1(a2,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (   -4.00000000) T2(w,x,a0,j) W16(i,a0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W16aa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, 0));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no23_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO23_X0_TYPE0_NOERI)
    (W16aa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no23_x1_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO23_X1_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W16aa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   24] -- 
  // |-- [    0] --| W17(a0,i) += (    1.00000000) D1(a1,a0) Fc1(i,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) T2(x,w,a0,j) W17(a0,i) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W17aa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, 0));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no24_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO24_X0_TYPE0_NOERI)
    (W17aa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no24_x1_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO24_X1_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W17aa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   25] -- 
  // |-- [    0] --| W18(a0,i) += (    1.00000000) D1(a1,a0) Fc1(i,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (   -4.00000000) T2(w,x,a0,j) W18(a0,i) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W18aa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, 0));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no25_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO25_X0_TYPE0_NOERI)
    (W18aa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no25_x1_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO25_X1_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W18aa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   26] -- 
  // |-- [    0] --| W19(a0,j) += (    1.00000000) D1(a1,a0) Fc1(j,a1) 
  // |-- [    1] --| S2(x,w,j,i) += (   -4.00000000) T2(x,w,a0,i) W19(a0,j) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W19aa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, 0));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no26_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO26_X0_TYPE0_NOERI)
    (W19aa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    S2b = orz::DTensor(retval.namps_iamp()[ii]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no26_x1_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO26_X1_TYPE0_NOERI)
      (si, ii, T2b.cptr(), W19aa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    retval.acc_amp2(ii, S2b);
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   27] -- 
  // |-- [    0] --| W20(a0,j) += (    1.00000000) D1(a1,a0) Fc1(j,a1) 
  // |-- [    1] --| S2(x,w,j,i) += (    2.00000000) T2(w,x,a0,i) W20(a0,j) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W20aa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, 0));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no27_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO27_X0_TYPE0_NOERI)
    (W20aa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    S2b = orz::DTensor(retval.namps_iamp()[ii]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no27_x1_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO27_X1_TYPE0_NOERI)
      (si, ii, T2b.cptr(), W20aa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    retval.acc_amp2(ii, S2b);
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   28] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -8.00000000) T2(w,c0,i,j) Fc1(x,c0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no28_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO28_X0_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   29] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    4.00000000) T2(x,c0,i,j) Fc1(w,c0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no29_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO29_X0_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   30] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    4.00000000) T2(c0,w,i,j) Fc1(x,c0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no30_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO30_X0_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   31] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -8.00000000) T2(c0,x,i,j) Fc1(w,c0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no31_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO31_X0_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   32] -- 
  // |-- [    0] --| W21() += (    1.00000000) D1(a1,a0) Fc1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    8.00000000) T2(w,x,i,j) W21() 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  double W21_sigma_ccoo_ccoo(0);
  FC_FUNC(g_if_sigma_ccoo_ccoo_no32_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO32_X0_TYPE0_NOERI)
    (&W21_sigma_ccoo_ccoo, nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no32_x1_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO32_X1_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), &W21_sigma_ccoo_ccoo, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   33] -- 
  // |-- [    0] --| W22() += (    1.00000000) D1(a1,a0) Fc1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -4.00000000) T2(x,w,i,j) W22() 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  double W22_sigma_ccoo_ccoo(0);
  FC_FUNC(g_if_sigma_ccoo_ccoo_no33_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO33_X0_TYPE0_NOERI)
    (&W22_sigma_ccoo_ccoo, nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no33_x1_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO33_X1_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), &W22_sigma_ccoo_ccoo, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   34] -- 
  // |-- [    0] --| W23(j,i,a0,a1) += (    1.00000000) D2(j,a0,i,a2) Fc1(a2,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) T2(x,w,a0,a1) W23(j,i,a0,a1) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      orz::DTensor W23aa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sa1));
      FC_FUNC(g_if_sigma_ccoo_ccoo_no34_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO34_X0_TYPE0_NOERI)
        (sa1, ia1, sj, ij, W23aa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_ccoo_ccoo_no34_x1_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO34_X1_TYPE0_NOERI)
        (sa1, ia1, sj, ij, T2b.cptr(), W23aa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- No. [   35] -- 
  // |-- [    0] --| W24(j,i,a0,a1) += (    1.00000000) D2(j,a2,i,a0) Fc1(a2,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) T2(x,w,a1,a0) W24(j,i,a0,a1) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      orz::DTensor W24aa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sa0));
      FC_FUNC(g_if_sigma_ccoo_ccoo_no35_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO35_X0_TYPE0_NOERI)
        (sa0, ia0, sj, ij, W24aa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_ccoo_ccoo_no35_x1_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO35_X1_TYPE0_NOERI)
        (sa0, ia0, sj, ij, T2b.cptr(), W24aa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- No. [   36] -- 
  // |-- [    0] --| W25(j,a0) += (    1.00000000) D1(j,a1) Fc1(a1,a0) 
  // |-- [    1] --| S2(x,w,j,i) += (   -4.00000000) T2(x,w,a0,i) W25(j,a0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W25aa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, 0));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no36_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO36_X0_TYPE0_NOERI)
    (W25aa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    S2b = orz::DTensor(retval.namps_iamp()[ii]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no36_x1_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO36_X1_TYPE0_NOERI)
      (si, ii, T2b.cptr(), W25aa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    retval.acc_amp2(ii, S2b);
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   37] -- 
  // |-- [    0] --| W26(j,a0) += (    1.00000000) D1(j,a1) Fc1(a1,a0) 
  // |-- [    1] --| S2(x,w,j,i) += (    2.00000000) T2(w,x,a0,i) W26(j,a0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W26aa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, 0));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no37_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO37_X0_TYPE0_NOERI)
    (W26aa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    S2b = orz::DTensor(retval.namps_iamp()[ii]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no37_x1_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO37_X1_TYPE0_NOERI)
      (si, ii, T2b.cptr(), W26aa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    retval.acc_amp2(ii, S2b);
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   38] -- 
  // |-- [    0] --| W27(x,w,j,a1) += (    1.00000000) T2(x,w,a0,a1) D1(j,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -4.00000000) Fc1(i,a1) W27(x,w,j,a1) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      orz::DTensor W27cc_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcc(symblockinfo, sj^sa1));
      FC_FUNC(g_if_sigma_ccoo_ccoo_no38_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO38_X0_TYPE0_NOERI)
        (sa1, ia1, sj, ij, T2b.cptr(), W27cc_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
      FC_FUNC(g_if_sigma_ccoo_ccoo_no38_x1_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO38_X1_TYPE0_NOERI)
        (sa1, ia1, sj, ij, W27cc_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- No. [   39] -- 
  // |-- [    0] --| W28(x,w,j,a1) += (    1.00000000) T2(x,w,a1,a0) D1(j,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) Fc1(i,a1) W28(x,w,j,a1) 
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W28cca_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcca(symblockinfo, sj));
    for(int sa0 = 0;sa0 < nir;++sa0){ 
    for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
      T2b = T2.get_amp2(ia0);
      FC_FUNC(g_if_sigma_ccoo_ccoo_no39_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO39_X0_TYPE0_NOERI)
        (sa0, ia0, sj, ij, T2b.cptr(), W28cca_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
    } // End ia0
    } // End sa0
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no39_x1_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO39_X1_TYPE0_NOERI)
      (sj, ij, W28cca_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   40] -- 
  // |-- [    0] --| W29(i,a0) += (    1.00000000) D1(i,a1) Fc1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) T2(x,w,a0,j) W29(i,a0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W29aa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, 0));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no40_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO40_X0_TYPE0_NOERI)
    (W29aa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no40_x1_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO40_X1_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W29aa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   41] -- 
  // |-- [    0] --| W30(i,a0) += (    1.00000000) D1(i,a1) Fc1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -4.00000000) T2(w,x,a0,j) W30(i,a0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W30aa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, 0));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no41_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO41_X0_TYPE0_NOERI)
    (W30aa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no41_x1_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO41_X1_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W30aa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   42] -- 
  // |-- [    0] --| W31(x,w,i,a1) += (    1.00000000) T2(x,w,a0,a1) D1(i,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) Fc1(j,a1) W31(x,w,i,a1) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W31cca_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcca(symblockinfo, sa1));
    FC_FUNC(g_if_sigma_ccoo_ccoo_no42_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO42_X0_TYPE0_NOERI)
      (sa1, ia1, T2b.cptr(), W31cca_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_ccoo_ccoo_no42_x1_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO42_X1_TYPE0_NOERI)
        (sa1, ia1, sj, ij, W31cca_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- No. [   43] -- 
  // |-- [    0] --| W32(x,w,i,a1) += (    1.00000000) T2(x,w,a1,a0) D1(i,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -4.00000000) Fc1(j,a1) W32(x,w,i,a1) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W32ccaa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xccaa(symblockinfo, 0));
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no43_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO43_X0_TYPE0_NOERI)
      (sa0, ia0, T2b.cptr(), W32ccaa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no43_x1_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO43_X1_TYPE0_NOERI)
      (sj, ij, W32ccaa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   44] -- 
  // |-- [    0] --| S2(x,w,j,i) += (    8.00000000) T2(x,w,a0,i) Fc1(j,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    S2b = orz::DTensor(retval.namps_iamp()[ii]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no44_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO44_X0_TYPE0_NOERI)
      (si, ii, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    retval.acc_amp2(ii, S2b);
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   45] -- 
  // |-- [    0] --| S2(x,w,j,i) += (   -4.00000000) T2(w,x,a0,i) Fc1(j,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    S2b = orz::DTensor(retval.namps_iamp()[ii]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no45_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO45_X0_TYPE0_NOERI)
      (si, ii, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    retval.acc_amp2(ii, S2b);
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   46] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    8.00000000) T2(w,x,a0,j) Fc1(i,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no46_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO46_X0_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   47] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -4.00000000) T2(x,w,a0,j) Fc1(i,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no47_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO47_X0_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   48] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    1.00000000) T2(x,w,a1,a0) C4(a1,j,a0,i) 
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
      FC_FUNC(g_if_sigma_ccoo_ccoo_no48_x0_type0_noeri,G_IF_SIGMA_CCOO_CCOO_NO48_X0_TYPE0_NOERI)
        (sa0, ia0, sj, ij, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
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
  // -- Title : sigma_ccoo_ccoo

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W33(x,w,a1,a0) += (    1.00000000) V2(x,c1,w,c0) T2(c1,c0,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) D2(j,a1,i,a0) W33(x,w,a1,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W33ca_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xca(symblockinfo, sx^sa0));
    FC_FUNC(g_if_sigma_ccoo_ccoo_no0_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO0_X0_TYPE1_ERI_C)
      (sa0, ia0, sx, ix, T2b.cptr(), V2_sym.cptr(), W33ca_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_ccoo_ccoo_no0_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO0_X1_TYPE1_ERI_C)
        (sa0, ia0, sj, ij, sx, ix, W33ca_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // |-- [    0] --| W34(x,c0,j,i,a1,a0) += (    1.00000000) V2(x,a3,c0,a2) D3(j,a3,i,a1,a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(w,c0,a1,a0) W34(x,c0,j,i,a1,a0) 
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
      orz::DTensor W34caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sx^sj^sa0));
      FC_FUNC(g_if_sigma_ccoo_ccoo_no1_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO1_X0_TYPE1_ERI_C)
        (sa0, ia0, sj, ij, sx, ix, V2_sym.cptr(), W34caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_ccoo_ccoo_no1_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO1_X1_TYPE1_ERI_C)
        (sa0, ia0, sj, ij, sx, ix, T2b.cptr(), W34caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // |-- [    0] --| W35(w,c0,j,i,a1,a0) += (    1.00000000) V2(w,a3,c0,a2) D3(j,a1,i,a3,a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(x,c0,a1,a0) W35(w,c0,j,i,a1,a0) 
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
      orz::DTensor W35caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw^sj^sa0));
      FC_FUNC(g_if_sigma_ccoo_ccoo_no2_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO2_X0_TYPE1_ERI_C)
        (sa0, ia0, sj, ij, sw, iw, V2_sym.cptr(), W35caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_ccoo_ccoo_no2_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO2_X1_TYPE1_ERI_C)
        (sa0, ia0, sj, ij, sw, iw, T2b.cptr(), W35caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // |-- [    0] --| W36(x,c0,j,i,a1,a0) += (    1.00000000) V2(x,c0,a3,a2) D3(j,a1,i,a0,a3,a2) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(w,c0,a0,a1) W36(x,c0,j,i,a1,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      orz::DTensor W36caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sx^sj^sa1));
      FC_FUNC(g_if_sigma_ccoo_ccoo_no3_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO3_X0_TYPE1_ERI_C)
        (sa1, ia1, sj, ij, sx, ix, V2_sym.cptr(), W36caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_ccoo_ccoo_no3_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO3_X1_TYPE1_ERI_C)
        (sa1, ia1, sj, ij, sx, ix, T2b.cptr(), W36caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // |-- [    0] --| W37(w,c0,j,i,a1,a0) += (    1.00000000) V2(w,c0,a3,a2) D3(j,a1,i,a0,a3,a2) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(x,c0,a1,a0) W37(w,c0,j,i,a1,a0) 
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
      orz::DTensor W37caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw^sj^sa0));
      FC_FUNC(g_if_sigma_ccoo_ccoo_no4_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO4_X0_TYPE1_ERI_C)
        (sa0, ia0, sj, ij, sw, iw, V2_sym.cptr(), W37caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_ccoo_ccoo_no4_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO4_X1_TYPE1_ERI_C)
        (sa0, ia0, sj, ij, sw, iw, T2b.cptr(), W37caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // |-- [    0] --| W43(x,w,a0,i) += (    1.00000000) V2(x,c1,w,c0) T2(c1,c0,a0,i) 
  // |-- [    1] --| S2(x,w,j,i) += (   -4.00000000) D1(j,a0) W43(x,w,a0,i) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    orz::DTensor W43ca_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xca(symblockinfo, sx^si));
    S2b = orz::DTensor(retval.namps_iamp()[ii]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no5_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO5_X0_TYPE1_ERI_C)
      (si, ii, sx, ix, T2b.cptr(), V2_sym.cptr(), W43ca_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no5_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO5_X1_TYPE1_ERI_C)
      (si, ii, sx, ix, W43ca_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    retval.acc_amp2(ii, S2b);
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W44(x,w,i,a0) += (    1.00000000) V2(x,c1,w,c0) T2(c0,c1,a0,i) 
  // |-- [    1] --| S2(x,w,j,i) += (    2.00000000) D1(j,a0) W44(x,w,i,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    orz::DTensor W44ca_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xca(symblockinfo, sx^si));
    S2b = orz::DTensor(retval.namps_iamp()[ii]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no6_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO6_X0_TYPE1_ERI_C)
      (si, ii, sx, ix, T2b.cptr(), V2_sym.cptr(), W44ca_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no6_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO6_X1_TYPE1_ERI_C)
      (si, ii, sx, ix, W44ca_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    retval.acc_amp2(ii, S2b);
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W45(w,c0,j,a0) += (    1.00000000) V2(w,a2,c0,a1) D2(j,a2,a1,a0) 
  // |-- [    1] --| S2(x,w,j,i) += (   -2.00000000) T2(c0,x,a0,i) W45(w,c0,j,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W45caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no7_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO7_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W45caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    S2b = orz::DTensor(retval.namps_iamp()[ii]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no7_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO7_X1_TYPE1_ERI_C)
      (si, ii, sw, iw, T2b.cptr(), W45caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    retval.acc_amp2(ii, S2b);
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| W46(x,c0,j,a0) += (    1.00000000) V2(x,a2,c0,a1) D2(j,a2,a1,a0) 
  // |-- [    1] --| S2(x,w,j,i) += (    4.00000000) T2(c0,w,a0,i) W46(x,c0,j,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W46caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no8_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO8_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W46caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    S2b = orz::DTensor(retval.namps_iamp()[ii]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no8_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO8_X1_TYPE1_ERI_C)
      (si, ii, sx, ix, T2b.cptr(), W46caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    retval.acc_amp2(ii, S2b);
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W47(w,c0,j,a0) += (    1.00000000) V2(w,a2,c0,a1) D2(j,a0,a1,a2) 
  // |-- [    1] --| S2(x,w,j,i) += (   -2.00000000) T2(x,c0,a0,i) W47(w,c0,j,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W47caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no9_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO9_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W47caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    S2b = orz::DTensor(retval.namps_iamp()[ii]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no9_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO9_X1_TYPE1_ERI_C)
      (si, ii, sw, iw, T2b.cptr(), W47caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    retval.acc_amp2(ii, S2b);
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| W48(x,c0,j,a0) += (    1.00000000) V2(x,a2,c0,a1) D2(j,a2,a1,a0) 
  // |-- [    1] --| S2(x,w,j,i) += (   -2.00000000) T2(w,c0,a0,i) W48(x,c0,j,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W48caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no10_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO10_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W48caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    S2b = orz::DTensor(retval.namps_iamp()[ii]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no10_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO10_X1_TYPE1_ERI_C)
      (si, ii, sx, ix, T2b.cptr(), W48caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    retval.acc_amp2(ii, S2b);
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| W49(w,c0,j,a0) += (    1.00000000) V2(w,c0,a2,a1) D2(j,a0,a2,a1) 
  // |-- [    1] --| S2(x,w,j,i) += (   -2.00000000) T2(c0,x,a0,i) W49(w,c0,j,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W49caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no11_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO11_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W49caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    S2b = orz::DTensor(retval.namps_iamp()[ii]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no11_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO11_X1_TYPE1_ERI_C)
      (si, ii, sw, iw, T2b.cptr(), W49caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    retval.acc_amp2(ii, S2b);
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| W50(x,c0,j,a0) += (    1.00000000) V2(x,c0,a2,a1) D2(j,a0,a2,a1) 
  // |-- [    1] --| S2(x,w,j,i) += (    4.00000000) T2(c0,w,a0,i) W50(x,c0,j,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W50caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no12_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO12_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W50caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    S2b = orz::DTensor(retval.namps_iamp()[ii]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no12_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO12_X1_TYPE1_ERI_C)
      (si, ii, sx, ix, T2b.cptr(), W50caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    retval.acc_amp2(ii, S2b);
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| W51(w,c0,j,a0) += (    1.00000000) V2(w,c0,a2,a1) D2(j,a0,a2,a1) 
  // |-- [    1] --| S2(x,w,j,i) += (    4.00000000) T2(x,c0,a0,i) W51(w,c0,j,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W51caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no13_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO13_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W51caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    S2b = orz::DTensor(retval.namps_iamp()[ii]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no13_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO13_X1_TYPE1_ERI_C)
      (si, ii, sw, iw, T2b.cptr(), W51caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    retval.acc_amp2(ii, S2b);
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   14] -- 
  // |-- [    0] --| W52(x,c0,j,a0) += (    1.00000000) V2(x,c0,a2,a1) D2(j,a0,a2,a1) 
  // |-- [    1] --| S2(x,w,j,i) += (   -2.00000000) T2(w,c0,a0,i) W52(x,c0,j,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W52caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no14_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO14_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W52caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    S2b = orz::DTensor(retval.namps_iamp()[ii]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no14_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO14_X1_TYPE1_ERI_C)
      (si, ii, sx, ix, T2b.cptr(), W52caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    retval.acc_amp2(ii, S2b);
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   15] -- 
  // |-- [    0] --| W60(x,w,a0,j) += (    1.00000000) V2(x,c1,w,c0) T2(c1,c0,a0,j) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) D1(i,a0) W60(x,w,a0,j) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W60ca_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xca(symblockinfo, sx^sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no15_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO15_X0_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), V2_sym.cptr(), W60ca_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no15_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO15_X1_TYPE1_ERI_C)
      (sj, ij, sx, ix, W60ca_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   16] -- 
  // |-- [    0] --| W61(x,w,j,a0) += (    1.00000000) V2(x,c1,w,c0) T2(c0,c1,a0,j) 
  // |-- [    1] --| S2(w,x,i,j) += (   -4.00000000) D1(i,a0) W61(x,w,j,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W61ca_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xca(symblockinfo, sx^sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no16_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO16_X0_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), V2_sym.cptr(), W61ca_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no16_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO16_X1_TYPE1_ERI_C)
      (sj, ij, sx, ix, W61ca_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   17] -- 
  // |-- [    0] --| W62(w,c0,i,a0) += (    1.00000000) V2(w,a2,c0,a1) D2(i,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    4.00000000) T2(c0,x,a0,j) W62(w,c0,i,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W62caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no17_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO17_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W62caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no17_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO17_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), W62caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   18] -- 
  // |-- [    0] --| W63(x,c0,i,a0) += (    1.00000000) V2(x,a2,c0,a1) D2(i,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(c0,w,a0,j) W63(x,c0,i,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W63caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no18_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO18_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W63caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no18_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO18_X1_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), W63caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   19] -- 
  // |-- [    0] --| W64(w,c0,i,a0) += (    1.00000000) V2(w,a2,c0,a1) D2(i,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(x,c0,a0,j) W64(w,c0,i,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W64caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no19_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO19_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W64caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no19_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO19_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), W64caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   20] -- 
  // |-- [    0] --| W65(x,c0,i,a0) += (    1.00000000) V2(x,a2,c0,a1) D2(i,a0,a1,a2) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(w,c0,a0,j) W65(x,c0,i,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W65caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no20_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO20_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W65caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no20_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO20_X1_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), W65caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   21] -- 
  // |-- [    0] --| W66(w,c0,i,a0) += (    1.00000000) V2(w,c0,a2,a1) D2(i,a0,a2,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (    4.00000000) T2(c0,x,a0,j) W66(w,c0,i,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W66caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no21_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO21_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W66caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no21_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO21_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), W66caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   22] -- 
  // |-- [    0] --| W67(x,c0,i,a0) += (    1.00000000) V2(x,c0,a2,a1) D2(i,a0,a2,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(c0,w,a0,j) W67(x,c0,i,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W67caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no22_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO22_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W67caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no22_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO22_X1_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), W67caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   23] -- 
  // |-- [    0] --| W68(w,c0,i,a0) += (    1.00000000) V2(w,c0,a2,a1) D2(i,a0,a2,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(x,c0,a0,j) W68(w,c0,i,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W68caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no23_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO23_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W68caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no23_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO23_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), W68caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   24] -- 
  // |-- [    0] --| W69(x,c0,i,a0) += (    1.00000000) V2(x,c0,a2,a1) D2(i,a0,a2,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (    4.00000000) T2(w,c0,a0,j) W69(x,c0,i,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W69caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no24_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO24_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W69caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no24_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO24_X1_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), W69caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   25] -- 
  // |-- [    0] --| W73(w,c0,a0,i) += (    1.00000000) V2(w,i,c0,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -8.00000000) T2(c0,x,a0,j) W73(w,c0,a0,i) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W73caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no25_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO25_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W73caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no25_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO25_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), W73caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   26] -- 
  // |-- [    0] --| W74(x,c0,a0,i) += (    1.00000000) V2(x,i,c0,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    4.00000000) T2(c0,w,a0,j) W74(x,c0,a0,i) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W74caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no26_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO26_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W74caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no26_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO26_X1_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), W74caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   27] -- 
  // |-- [    0] --| W75(w,c0,a0,i) += (    1.00000000) V2(w,i,c0,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    4.00000000) T2(x,c0,a0,j) W75(w,c0,a0,i) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W75caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no27_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO27_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W75caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no27_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO27_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), W75caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   28] -- 
  // |-- [    0] --| W76(x,c0,a0,i) += (    1.00000000) V2(x,i,c0,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(w,c0,a0,j) W76(x,c0,a0,i) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W76caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no28_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO28_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W76caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no28_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO28_X1_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), W76caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   29] -- 
  // |-- [    0] --| W77(w,c0,a0,i) += (    1.00000000) V2(w,c0,i,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    4.00000000) T2(c0,x,a0,j) W77(w,c0,a0,i) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W77caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no29_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO29_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W77caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no29_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO29_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), W77caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   30] -- 
  // |-- [    0] --| W78(x,c0,a0,i) += (    1.00000000) V2(x,c0,i,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(c0,w,a0,j) W78(x,c0,a0,i) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W78caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no30_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO30_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W78caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no30_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO30_X1_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), W78caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   31] -- 
  // |-- [    0] --| W79(w,c0,a0,i) += (    1.00000000) V2(w,c0,i,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(x,c0,a0,j) W79(w,c0,a0,i) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W79caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no31_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO31_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W79caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no31_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO31_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), W79caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   32] -- 
  // |-- [    0] --| W80(x,c0,a0,i) += (    1.00000000) V2(x,c0,i,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    4.00000000) T2(w,c0,a0,j) W80(x,c0,a0,i) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W80caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no32_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO32_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W80caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no32_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO32_X1_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), W80caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   33] -- 
  // |-- [    0] --| W83(x,c0,a0,j) += (    1.00000000) V2(x,j,c0,a1) D1(a1,a0) 
  // |-- [    1] --| S2(x,w,j,i) += (   -8.00000000) T2(c0,w,a0,i) W83(x,c0,a0,j) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W83caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no33_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO33_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W83caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    S2b = orz::DTensor(retval.namps_iamp()[ii]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no33_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO33_X1_TYPE1_ERI_C)
      (si, ii, sx, ix, T2b.cptr(), W83caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    retval.acc_amp2(ii, S2b);
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   34] -- 
  // |-- [    0] --| W84(w,c0,a0,j) += (    1.00000000) V2(w,j,c0,a1) D1(a1,a0) 
  // |-- [    1] --| S2(x,w,j,i) += (    4.00000000) T2(c0,x,a0,i) W84(w,c0,a0,j) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W84caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no34_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO34_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W84caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    S2b = orz::DTensor(retval.namps_iamp()[ii]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no34_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO34_X1_TYPE1_ERI_C)
      (si, ii, sw, iw, T2b.cptr(), W84caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    retval.acc_amp2(ii, S2b);
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   35] -- 
  // |-- [    0] --| W85(w,c0,a0,j) += (    1.00000000) V2(w,j,c0,a1) D1(a1,a0) 
  // |-- [    1] --| S2(x,w,j,i) += (   -2.00000000) T2(x,c0,a0,i) W85(w,c0,a0,j) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W85caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no35_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO35_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W85caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    S2b = orz::DTensor(retval.namps_iamp()[ii]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no35_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO35_X1_TYPE1_ERI_C)
      (si, ii, sw, iw, T2b.cptr(), W85caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    retval.acc_amp2(ii, S2b);
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   36] -- 
  // |-- [    0] --| W86(x,c0,a0,j) += (    1.00000000) V2(x,j,c0,a1) D1(a1,a0) 
  // |-- [    1] --| S2(x,w,j,i) += (    4.00000000) T2(w,c0,a0,i) W86(x,c0,a0,j) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W86caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no36_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO36_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W86caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    S2b = orz::DTensor(retval.namps_iamp()[ii]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no36_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO36_X1_TYPE1_ERI_C)
      (si, ii, sx, ix, T2b.cptr(), W86caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    retval.acc_amp2(ii, S2b);
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   37] -- 
  // |-- [    0] --| W87(w,c0,a0,j) += (    1.00000000) V2(w,c0,j,a1) D1(a1,a0) 
  // |-- [    1] --| S2(x,w,j,i) += (   -2.00000000) T2(c0,x,a0,i) W87(w,c0,a0,j) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W87caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no37_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO37_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W87caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    S2b = orz::DTensor(retval.namps_iamp()[ii]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no37_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO37_X1_TYPE1_ERI_C)
      (si, ii, sw, iw, T2b.cptr(), W87caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    retval.acc_amp2(ii, S2b);
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   38] -- 
  // |-- [    0] --| W88(x,c0,a0,j) += (    1.00000000) V2(x,c0,j,a1) D1(a1,a0) 
  // |-- [    1] --| S2(x,w,j,i) += (    4.00000000) T2(c0,w,a0,i) W88(x,c0,a0,j) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W88caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no38_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO38_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W88caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    S2b = orz::DTensor(retval.namps_iamp()[ii]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no38_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO38_X1_TYPE1_ERI_C)
      (si, ii, sx, ix, T2b.cptr(), W88caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    retval.acc_amp2(ii, S2b);
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   39] -- 
  // |-- [    0] --| W89(w,c0,a0,j) += (    1.00000000) V2(w,c0,j,a1) D1(a1,a0) 
  // |-- [    1] --| S2(x,w,j,i) += (    4.00000000) T2(x,c0,a0,i) W89(w,c0,a0,j) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W89caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no39_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO39_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W89caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    S2b = orz::DTensor(retval.namps_iamp()[ii]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no39_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO39_X1_TYPE1_ERI_C)
      (si, ii, sw, iw, T2b.cptr(), W89caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    retval.acc_amp2(ii, S2b);
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   40] -- 
  // |-- [    0] --| W90(x,c0,a0,j) += (    1.00000000) V2(x,c0,j,a1) D1(a1,a0) 
  // |-- [    1] --| S2(x,w,j,i) += (   -2.00000000) T2(w,c0,a0,i) W90(x,c0,a0,j) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W90caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no40_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO40_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W90caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    S2b = orz::DTensor(retval.namps_iamp()[ii]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no40_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO40_X1_TYPE1_ERI_C)
      (si, ii, sx, ix, T2b.cptr(), W90caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    retval.acc_amp2(ii, S2b);
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   41] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    8.00000000) V2(x,c1,w,c0) T2(c0,c1,i,j) 
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
    FC_FUNC(g_if_sigma_ccoo_ccoo_no41_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO41_X0_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   42] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -4.00000000) V2(x,c0,w,c1) T2(c0,c1,i,j) 
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
    FC_FUNC(g_if_sigma_ccoo_ccoo_no42_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO42_X0_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   43] -- 
  // |-- [    0] --| W93(x,c0) += (    1.00000000) V2(x,a1,c0,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    4.00000000) T2(w,c0,i,j) W93(x,c0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W93c_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xc(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no43_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO43_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W93c_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no43_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO43_X1_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), W93c_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   44] -- 
  // |-- [    0] --| W94(w,c0) += (    1.00000000) V2(w,a1,c0,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(x,c0,i,j) W94(w,c0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W94c_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xc(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no44_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO44_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W94c_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no44_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO44_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), W94c_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   45] -- 
  // |-- [    0] --| W95(x,c0) += (    1.00000000) V2(x,a1,c0,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(c0,w,i,j) W95(x,c0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W95c_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xc(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no45_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO45_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W95c_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no45_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO45_X1_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), W95c_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   46] -- 
  // |-- [    0] --| W96(w,c0) += (    1.00000000) V2(w,a1,c0,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    4.00000000) T2(c0,x,i,j) W96(w,c0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W96c_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xc(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no46_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO46_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W96c_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no46_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO46_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), W96c_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   47] -- 
  // |-- [    0] --| W97(x,c0) += (    1.00000000) V2(x,c0,a1,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -8.00000000) T2(w,c0,i,j) W97(x,c0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W97c_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xc(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no47_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO47_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W97c_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no47_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO47_X1_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), W97c_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   48] -- 
  // |-- [    0] --| W98(w,c0) += (    1.00000000) V2(w,c0,a1,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    4.00000000) T2(x,c0,i,j) W98(w,c0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W98c_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xc(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no48_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO48_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W98c_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no48_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO48_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), W98c_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   49] -- 
  // |-- [    0] --| W99(x,c0) += (    1.00000000) V2(x,c0,a1,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    4.00000000) T2(c0,w,i,j) W99(x,c0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W99c_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xc(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no49_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO49_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W99c_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no49_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO49_X1_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), W99c_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   50] -- 
  // |-- [    0] --| W100(w,c0) += (    1.00000000) V2(w,c0,a1,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -8.00000000) T2(c0,x,i,j) W100(w,c0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W100c_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xc(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no50_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO50_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W100c_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no50_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO50_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), W100c_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   51] -- 
  // |-- [    0] --| W119(w,c0,j,a0) += (    1.00000000) V2(w,a1,c0,a0) D1(j,a1) 
  // |-- [    1] --| S2(x,w,j,i) += (    4.00000000) T2(c0,x,a0,i) W119(w,c0,j,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W119caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no51_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO51_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W119caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    S2b = orz::DTensor(retval.namps_iamp()[ii]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no51_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO51_X1_TYPE1_ERI_C)
      (si, ii, sw, iw, T2b.cptr(), W119caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    retval.acc_amp2(ii, S2b);
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   52] -- 
  // |-- [    0] --| W120(x,c0,j,a0) += (    1.00000000) V2(x,a1,c0,a0) D1(j,a1) 
  // |-- [    1] --| S2(x,w,j,i) += (   -8.00000000) T2(c0,w,a0,i) W120(x,c0,j,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W120caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no52_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO52_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W120caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    S2b = orz::DTensor(retval.namps_iamp()[ii]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no52_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO52_X1_TYPE1_ERI_C)
      (si, ii, sx, ix, T2b.cptr(), W120caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    retval.acc_amp2(ii, S2b);
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   53] -- 
  // |-- [    0] --| W121(w,c0,j,a0) += (    1.00000000) V2(w,a1,c0,a0) D1(j,a1) 
  // |-- [    1] --| S2(x,w,j,i) += (   -2.00000000) T2(x,c0,a0,i) W121(w,c0,j,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W121caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no53_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO53_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W121caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    S2b = orz::DTensor(retval.namps_iamp()[ii]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no53_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO53_X1_TYPE1_ERI_C)
      (si, ii, sw, iw, T2b.cptr(), W121caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    retval.acc_amp2(ii, S2b);
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   54] -- 
  // |-- [    0] --| W122(x,c0,j,a0) += (    1.00000000) V2(x,a1,c0,a0) D1(j,a1) 
  // |-- [    1] --| S2(x,w,j,i) += (    4.00000000) T2(w,c0,a0,i) W122(x,c0,j,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W122caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no54_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO54_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W122caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    S2b = orz::DTensor(retval.namps_iamp()[ii]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no54_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO54_X1_TYPE1_ERI_C)
      (si, ii, sx, ix, T2b.cptr(), W122caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    retval.acc_amp2(ii, S2b);
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   55] -- 
  // |-- [    0] --| W123(w,c0,j,a0) += (    1.00000000) V2(w,c0,a1,a0) D1(j,a1) 
  // |-- [    1] --| S2(x,w,j,i) += (   -2.00000000) T2(c0,x,a0,i) W123(w,c0,j,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W123caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no55_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO55_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W123caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    S2b = orz::DTensor(retval.namps_iamp()[ii]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no55_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO55_X1_TYPE1_ERI_C)
      (si, ii, sw, iw, T2b.cptr(), W123caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    retval.acc_amp2(ii, S2b);
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   56] -- 
  // |-- [    0] --| W124(x,c0,j,a0) += (    1.00000000) V2(x,c0,a1,a0) D1(j,a1) 
  // |-- [    1] --| S2(x,w,j,i) += (    4.00000000) T2(c0,w,a0,i) W124(x,c0,j,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W124caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no56_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO56_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W124caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    S2b = orz::DTensor(retval.namps_iamp()[ii]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no56_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO56_X1_TYPE1_ERI_C)
      (si, ii, sx, ix, T2b.cptr(), W124caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    retval.acc_amp2(ii, S2b);
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   57] -- 
  // |-- [    0] --| W125(w,c0,j,a0) += (    1.00000000) V2(w,c0,a1,a0) D1(j,a1) 
  // |-- [    1] --| S2(x,w,j,i) += (    4.00000000) T2(x,c0,a0,i) W125(w,c0,j,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W125caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no57_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO57_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W125caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    S2b = orz::DTensor(retval.namps_iamp()[ii]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no57_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO57_X1_TYPE1_ERI_C)
      (si, ii, sw, iw, T2b.cptr(), W125caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    retval.acc_amp2(ii, S2b);
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   58] -- 
  // |-- [    0] --| W126(x,c0,j,a0) += (    1.00000000) V2(x,c0,a1,a0) D1(j,a1) 
  // |-- [    1] --| S2(x,w,j,i) += (   -2.00000000) T2(w,c0,a0,i) W126(x,c0,j,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W126caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no58_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO58_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W126caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    S2b = orz::DTensor(retval.namps_iamp()[ii]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no58_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO58_X1_TYPE1_ERI_C)
      (si, ii, sx, ix, T2b.cptr(), W126caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    retval.acc_amp2(ii, S2b);
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   59] -- 
  // |-- [    0] --| W141(w,c0,i,a0) += (    1.00000000) V2(w,a1,c0,a0) D1(i,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (   -8.00000000) T2(c0,x,a0,j) W141(w,c0,i,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W141caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no59_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO59_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W141caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no59_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO59_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), W141caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   60] -- 
  // |-- [    0] --| W142(x,c0,i,a0) += (    1.00000000) V2(x,a1,c0,a0) D1(i,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (    4.00000000) T2(c0,w,a0,j) W142(x,c0,i,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W142caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no60_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO60_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W142caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no60_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO60_X1_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), W142caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   61] -- 
  // |-- [    0] --| W143(w,c0,i,a0) += (    1.00000000) V2(w,a1,c0,a0) D1(i,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (    4.00000000) T2(x,c0,a0,j) W143(w,c0,i,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W143caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no61_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO61_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W143caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no61_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO61_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), W143caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   62] -- 
  // |-- [    0] --| W144(x,c0,i,a0) += (    1.00000000) V2(x,a1,c0,a0) D1(i,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(w,c0,a0,j) W144(x,c0,i,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W144caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no62_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO62_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W144caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no62_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO62_X1_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), W144caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   63] -- 
  // |-- [    0] --| W145(w,c0,i,a0) += (    1.00000000) V2(w,c0,a1,a0) D1(i,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (    4.00000000) T2(c0,x,a0,j) W145(w,c0,i,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W145caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no63_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO63_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W145caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no63_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO63_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), W145caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   64] -- 
  // |-- [    0] --| W146(x,c0,i,a0) += (    1.00000000) V2(x,c0,a1,a0) D1(i,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(c0,w,a0,j) W146(x,c0,i,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W146caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no64_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO64_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W146caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no64_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO64_X1_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), W146caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   65] -- 
  // |-- [    0] --| W147(w,c0,i,a0) += (    1.00000000) V2(w,c0,a1,a0) D1(i,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) T2(x,c0,a0,j) W147(w,c0,i,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W147caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no65_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO65_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W147caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no65_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO65_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), W147caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   66] -- 
  // |-- [    0] --| W148(x,c0,i,a0) += (    1.00000000) V2(x,c0,a1,a0) D1(i,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (    4.00000000) T2(w,c0,a0,j) W148(x,c0,i,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W148caa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no66_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO66_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W148caa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no66_x1_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO66_X1_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), W148caa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   67] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   16.00000000) V2(w,i,c0,a0) T2(c0,x,a0,j) 
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
    FC_FUNC(g_if_sigma_ccoo_ccoo_no67_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO67_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   68] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -8.00000000) V2(x,i,c0,a0) T2(c0,w,a0,j) 
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
    FC_FUNC(g_if_sigma_ccoo_ccoo_no68_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO68_X0_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   69] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -8.00000000) V2(w,i,c0,a0) T2(x,c0,a0,j) 
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
    FC_FUNC(g_if_sigma_ccoo_ccoo_no69_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO69_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   70] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    4.00000000) V2(x,i,c0,a0) T2(w,c0,a0,j) 
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
    FC_FUNC(g_if_sigma_ccoo_ccoo_no70_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO70_X0_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   71] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -8.00000000) V2(x,c0,i,a0) T2(w,c0,a0,j) 
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
    FC_FUNC(g_if_sigma_ccoo_ccoo_no71_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO71_X0_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   72] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    4.00000000) V2(w,c0,i,a0) T2(x,c0,a0,j) 
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
    FC_FUNC(g_if_sigma_ccoo_ccoo_no72_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO72_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   73] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    4.00000000) V2(x,c0,i,a0) T2(c0,w,a0,j) 
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
    FC_FUNC(g_if_sigma_ccoo_ccoo_no73_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO73_X0_TYPE1_ERI_C)
      (sj, ij, sx, ix, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   74] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -8.00000000) V2(w,c0,i,a0) T2(c0,x,a0,j) 
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
    FC_FUNC(g_if_sigma_ccoo_ccoo_no74_x0_type1_eri_c,G_IF_SIGMA_CCOO_CCOO_NO74_X0_TYPE1_ERI_C)
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

//-@type(2).declaration(begin)
  // --  Title : sigma_ccoo_ccoo
  //  >> Intermediates for the external contractions are defined here << 
  orz::DTensor W38ccaa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xccaa(symblockinfo, 0));
  orz::DTensor W39ccaa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xccaa(symblockinfo, 0));
  orz::DTensor W40ccaa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xccaa(symblockinfo, 0));
  orz::DTensor W41ccaa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xccaa(symblockinfo, 0));
  orz::DTensor W55ccaa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xccaa(symblockinfo, 0));
  orz::DTensor W56ccaa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xccaa(symblockinfo, 0));
  orz::DTensor W57ccaa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xccaa(symblockinfo, 0));
  orz::DTensor W58ccaa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xccaa(symblockinfo, 0));
  orz::DTensor W53aa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, 0));
  orz::DTensor W54aa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, 0));
  orz::DTensor W70aa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, 0));
  orz::DTensor W71aa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, 0));
  double W101_sigma_ccoo_ccoo(0);
  double W102_sigma_ccoo_ccoo(0);
  orz::DTensor W103ccaa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xccaa(symblockinfo, 0));
  orz::DTensor W104ccaa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xccaa(symblockinfo, 0));
  orz::DTensor W105ccaa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xccaa(symblockinfo, 0));
  orz::DTensor W106ccaa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xccaa(symblockinfo, 0));
  orz::DTensor W107ccaa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xccaa(symblockinfo, 0));
  orz::DTensor W108ccaa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xccaa(symblockinfo, 0));
  orz::DTensor W109ccaa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xccaa(symblockinfo, 0));
  orz::DTensor W110ccaa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xccaa(symblockinfo, 0));
//-@type(2).declaration(end)

  //*-- Entering to take the type 0 contractions --*//
//-@type(0).contraction(begin)

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W38(w,c0,j,a2) += (    1.00000000) T2(w,c0,a1,a0) D2(j,a1,a2,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no0_x0_type0_eri_o,G_IF_SIGMA_CCOO_CCOO_NO0_X0_TYPE0_ERI_O)
      (sa0, ia0, T2b.cptr(), W38ccaa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W39(x,c0,j,a2) += (    1.00000000) T2(x,c0,a1,a0) D2(j,a1,a2,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no1_x0_type0_eri_o,G_IF_SIGMA_CCOO_CCOO_NO1_X0_TYPE0_ERI_O)
      (sa0, ia0, T2b.cptr(), W39ccaa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W40(w,c0,j,a2) += (    1.00000000) T2(w,c0,a0,a1) D2(j,a1,a2,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no2_x0_type0_eri_o,G_IF_SIGMA_CCOO_CCOO_NO2_X0_TYPE0_ERI_O)
      (sa1, ia1, T2b.cptr(), W40ccaa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W41(x,c0,j,a2) += (    1.00000000) T2(x,c0,a1,a0) D2(j,a1,a2,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no3_x0_type0_eri_o,G_IF_SIGMA_CCOO_CCOO_NO3_X0_TYPE0_ERI_O)
      (sa0, ia0, T2b.cptr(), W41ccaa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W55(w,c0,i,a2) += (    1.00000000) T2(w,c0,a1,a0) D2(i,a1,a2,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no4_x0_type0_eri_o,G_IF_SIGMA_CCOO_CCOO_NO4_X0_TYPE0_ERI_O)
      (sa0, ia0, T2b.cptr(), W55ccaa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W56(x,c0,i,a2) += (    1.00000000) T2(x,c0,a1,a0) D2(i,a1,a2,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no5_x0_type0_eri_o,G_IF_SIGMA_CCOO_CCOO_NO5_X0_TYPE0_ERI_O)
      (sa0, ia0, T2b.cptr(), W56ccaa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W57(w,c0,i,a2) += (    1.00000000) T2(w,c0,a1,a0) D2(i,a1,a2,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no6_x0_type0_eri_o,G_IF_SIGMA_CCOO_CCOO_NO6_X0_TYPE0_ERI_O)
      (sa0, ia0, T2b.cptr(), W57ccaa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W58(x,c0,i,a2) += (    1.00000000) T2(x,c0,a0,a1) D2(i,a1,a2,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no7_x0_type0_eri_o,G_IF_SIGMA_CCOO_CCOO_NO7_X0_TYPE0_ERI_O)
      (sa1, ia1, T2b.cptr(), W58ccaa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
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
  // -- Title : sigma_ccoo_ccoo

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| S2(x,w,j,i) += (   -2.00000000) V2(i,x,c0,a2) W38(w,c0,j,a2) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ii]);
  FC_FUNC(g_if_sigma_ccoo_ccoo_no0_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO0_X0_TYPE1_ERI_O)
    (si, ii, V2_sym.cptr(), W38ccaa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ii, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| S2(x,w,j,i) += (    4.00000000) V2(i,w,c0,a2) W39(x,c0,j,a2) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ii]);
  FC_FUNC(g_if_sigma_ccoo_ccoo_no1_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO1_X0_TYPE1_ERI_O)
    (si, ii, V2_sym.cptr(), W39ccaa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ii, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| S2(x,w,j,i) += (   -2.00000000) V2(i,a2,x,c0) W40(w,c0,j,a2) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ii]);
  FC_FUNC(g_if_sigma_ccoo_ccoo_no2_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO2_X0_TYPE1_ERI_O)
    (si, ii, V2_sym.cptr(), W40ccaa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ii, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| S2(x,w,j,i) += (   -2.00000000) V2(i,a2,w,c0) W41(x,c0,j,a2) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ii]);
  FC_FUNC(g_if_sigma_ccoo_ccoo_no3_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO3_X0_TYPE1_ERI_O)
    (si, ii, V2_sym.cptr(), W41ccaa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ii, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W42(j,a1,a0,i) += (    1.00000000) V2(i,a3,a4,a2) D3(j,a1,a4,a2,a3,a0) 
  // |-- [    1] --| S2(x,w,j,i) += (    2.00000000) T2(x,w,a1,a0) W42(j,a1,a0,i) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ii]);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W42aa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa0^si));
    FC_FUNC(g_if_sigma_ccoo_ccoo_no4_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO4_X0_TYPE1_ERI_O)
      (sa0, ia0, si, ii, V2_sym.cptr(), W42aa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no4_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO4_X1_TYPE1_ERI_O)
      (sa0, ia0, si, ii, T2b.cptr(), W42aa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ii, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W53(j,a0) += (    1.00000000) V2(a4,a2,a3,a1) D3(j,a0,a4,a2,a3,a1) 
  int sa4(s_eri);
  int ia4(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  FC_FUNC(g_if_sigma_ccoo_ccoo_no5_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO5_X0_TYPE1_ERI_O)
    (sa4, ia4, V2_sym.cptr(), W53aa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W54(j,a0) += (    1.00000000) V2(a4,a2,a3,a1) D3(j,a0,a4,a2,a3,a1) 
  int sa4(s_eri);
  int ia4(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  FC_FUNC(g_if_sigma_ccoo_ccoo_no6_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO6_X0_TYPE1_ERI_O)
    (sa4, ia4, V2_sym.cptr(), W54aa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    4.00000000) V2(j,x,c0,a2) W55(w,c0,i,a2) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  FC_FUNC(g_if_sigma_ccoo_ccoo_no7_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO7_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W55ccaa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -2.00000000) V2(j,w,c0,a2) W56(x,c0,i,a2) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  FC_FUNC(g_if_sigma_ccoo_ccoo_no8_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO8_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W56ccaa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -2.00000000) V2(j,a2,x,c0) W57(w,c0,i,a2) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  FC_FUNC(g_if_sigma_ccoo_ccoo_no9_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO9_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W57ccaa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -2.00000000) V2(j,a2,w,c0) W58(x,c0,i,a2) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  FC_FUNC(g_if_sigma_ccoo_ccoo_no10_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO10_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W58ccaa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| W59(i,a1,a0,j) += (    1.00000000) V2(j,a3,a4,a2) D3(i,a1,a4,a2,a3,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) T2(x,w,a0,a1) W59(i,a1,a0,j) 
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
    orz::DTensor W59aa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa1^sj));
    FC_FUNC(g_if_sigma_ccoo_ccoo_no11_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO11_X0_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, V2_sym.cptr(), W59aa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no11_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO11_X1_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, T2b.cptr(), W59aa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| W70(i,a0) += (    1.00000000) V2(a4,a2,a3,a1) D3(i,a0,a4,a2,a3,a1) 
  int sa4(s_eri);
  int ia4(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  FC_FUNC(g_if_sigma_ccoo_ccoo_no12_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO12_X0_TYPE1_ERI_O)
    (sa4, ia4, V2_sym.cptr(), W70aa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| W71(i,a0) += (    1.00000000) V2(a4,a2,a3,a1) D3(i,a0,a4,a2,a3,a1) 
  int sa4(s_eri);
  int ia4(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  FC_FUNC(g_if_sigma_ccoo_ccoo_no13_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO13_X0_TYPE1_ERI_O)
    (sa4, ia4, V2_sym.cptr(), W71aa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   14] -- 
  // |-- [    0] --| W72(a1,a0,j,i) += (    1.00000000) V2(j,a3,i,a2) D2(a3,a1,a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) T2(x,w,a1,a0) W72(a1,a0,j,i) 
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
    orz::DTensor W72aa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa0^sj));
    FC_FUNC(g_if_sigma_ccoo_ccoo_no14_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO14_X0_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, V2_sym.cptr(), W72aa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no14_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO14_X1_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, T2b.cptr(), W72aa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   15] -- 
  // |-- [    0] --| W81(a0,i) += (    1.00000000) V2(i,a2,a3,a1) D2(a3,a1,a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -4.00000000) T2(w,x,a0,j) W81(a0,i) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W81a_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xa(symblockinfo, si));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no15_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO15_X0_TYPE1_ERI_O)
    (si, ii, V2_sym.cptr(), W81a_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no15_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO15_X1_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), W81a_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   16] -- 
  // |-- [    0] --| W82(a0,i) += (    1.00000000) V2(i,a2,a3,a1) D2(a3,a1,a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) T2(x,w,a0,j) W82(a0,i) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W82a_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xa(symblockinfo, si));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no16_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO16_X0_TYPE1_ERI_O)
    (si, ii, V2_sym.cptr(), W82a_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no16_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO16_X1_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), W82a_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   17] -- 
  // |-- [    0] --| W91(a0,j) += (    1.00000000) V2(j,a2,a3,a1) D2(a3,a1,a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) T2(x,w,i,a0) W91(a0,j) 
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
    double W91_sigma_ccoo_ccoo(0);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no17_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO17_X0_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, V2_sym.cptr(), &W91_sigma_ccoo_ccoo, nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no17_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO17_X1_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, T2b.cptr(), &W91_sigma_ccoo_ccoo, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   18] -- 
  // |-- [    0] --| W92(a0,j) += (    1.00000000) V2(j,a2,a3,a1) D2(a3,a1,a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -4.00000000) T2(x,w,a0,i) W92(a0,j) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W92a_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xa(symblockinfo, sj));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no18_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO18_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W92a_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no18_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO18_X1_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), W92a_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   19] -- 
  // |-- [    0] --| W101() += (    1.00000000) V2(a3,a1,a2,a0) D2(a3,a1,a2,a0) 
  int sa3(s_eri);
  int ia3(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  FC_FUNC(g_if_sigma_ccoo_ccoo_no19_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO19_X0_TYPE1_ERI_O)
    (sa3, ia3, V2_sym.cptr(), &W101_sigma_ccoo_ccoo, nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   20] -- 
  // |-- [    0] --| W102() += (    1.00000000) V2(a3,a1,a2,a0) D2(a3,a1,a2,a0) 
  int sa3(s_eri);
  int ia3(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  FC_FUNC(g_if_sigma_ccoo_ccoo_no20_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO20_X0_TYPE1_ERI_O)
    (sa3, ia3, V2_sym.cptr(), &W102_sigma_ccoo_ccoo, nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   21] -- 
  // |-- [    0] --| W103(w,x,a0,a2) += (    1.00000000) V2(a1,c0,x,a2) T2(c0,w,a0,a1) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(ia1);
  FC_FUNC(g_if_sigma_ccoo_ccoo_no21_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO21_X0_TYPE1_ERI_O)
    (sa1, ia1, T2b.cptr(), V2_sym.cptr(), W103ccaa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   22] -- 
  // |-- [    0] --| W104(x,w,a0,a2) += (    1.00000000) V2(a1,c0,w,a2) T2(c0,x,a0,a1) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(ia1);
  FC_FUNC(g_if_sigma_ccoo_ccoo_no22_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO22_X0_TYPE1_ERI_O)
    (sa1, ia1, T2b.cptr(), V2_sym.cptr(), W104ccaa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   23] -- 
  // |-- [    0] --| W105(x,w,a0,a2) += (    1.00000000) V2(a1,c0,w,a2) T2(x,c0,a0,a1) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(ia1);
  FC_FUNC(g_if_sigma_ccoo_ccoo_no23_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO23_X0_TYPE1_ERI_O)
    (sa1, ia1, T2b.cptr(), V2_sym.cptr(), W105ccaa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   24] -- 
  // |-- [    0] --| W106(w,x,a0,a2) += (    1.00000000) V2(a1,c0,x,a2) T2(w,c0,a0,a1) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(ia1);
  FC_FUNC(g_if_sigma_ccoo_ccoo_no24_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO24_X0_TYPE1_ERI_O)
    (sa1, ia1, T2b.cptr(), V2_sym.cptr(), W106ccaa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   25] -- 
  // |-- [    0] --| W107(w,x,a0,a2) += (    1.00000000) V2(a1,a2,x,c0) T2(c0,w,a0,a1) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(ia1);
  FC_FUNC(g_if_sigma_ccoo_ccoo_no25_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO25_X0_TYPE1_ERI_O)
    (sa1, ia1, T2b.cptr(), V2_sym.cptr(), W107ccaa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   26] -- 
  // |-- [    0] --| W108(x,w,a0,a2) += (    1.00000000) V2(a1,a2,w,c0) T2(c0,x,a0,a1) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(ia1);
  FC_FUNC(g_if_sigma_ccoo_ccoo_no26_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO26_X0_TYPE1_ERI_O)
    (sa1, ia1, T2b.cptr(), V2_sym.cptr(), W108ccaa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   27] -- 
  // |-- [    0] --| W109(x,w,a0,a2) += (    1.00000000) V2(a1,a2,w,c0) T2(x,c0,a0,a1) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(ia1);
  FC_FUNC(g_if_sigma_ccoo_ccoo_no27_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO27_X0_TYPE1_ERI_O)
    (sa1, ia1, T2b.cptr(), V2_sym.cptr(), W109ccaa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   28] -- 
  // |-- [    0] --| W110(w,x,a0,a2) += (    1.00000000) V2(a1,a2,x,c0) T2(w,c0,a0,a1) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(ia1);
  FC_FUNC(g_if_sigma_ccoo_ccoo_no28_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO28_X0_TYPE1_ERI_O)
    (sa1, ia1, T2b.cptr(), V2_sym.cptr(), W110ccaa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   29] -- 
  // |-- [    0] --| W111(j,i,a0,a1) += (    1.00000000) V2(a1,a3,a4,a2) D3(j,a3,i,a0,a4,a2) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) T2(w,x,a0,a1) W111(j,i,a0,a1) 
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
    orz::DTensor W111aa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sa1));
    FC_FUNC(g_if_sigma_ccoo_ccoo_no29_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO29_X0_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, V2_sym.cptr(), W111aa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no29_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO29_X1_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, T2b.cptr(), W111aa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   30] -- 
  // |-- [    0] --| W112(j,i,a0,a1) += (    1.00000000) V2(a1,a3,a4,a2) D3(j,a0,i,a3,a4,a2) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) T2(x,w,a0,a1) W112(j,i,a0,a1) 
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
    orz::DTensor W112aa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sa1));
    FC_FUNC(g_if_sigma_ccoo_ccoo_no30_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO30_X0_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, V2_sym.cptr(), W112aa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no30_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO30_X1_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, T2b.cptr(), W112aa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   31] -- 
  // |-- [    0] --| W113(w,x,a0,i) += (    1.00000000) V2(i,x,c0,a1) T2(w,c0,a1,a0) 
  // |-- [    1] --| S2(x,w,j,i) += (   -2.00000000) D1(j,a0) W113(w,x,a0,i) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ii]);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W113cc_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcc(symblockinfo, sa0^si));
    FC_FUNC(g_if_sigma_ccoo_ccoo_no31_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO31_X0_TYPE1_ERI_O)
      (sa0, ia0, si, ii, T2b.cptr(), V2_sym.cptr(), W113cc_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no31_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO31_X1_TYPE1_ERI_O)
      (sa0, ia0, si, ii, W113cc_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ii, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   32] -- 
  // |-- [    0] --| W114(x,w,a0,i) += (    1.00000000) V2(i,w,c0,a1) T2(x,c0,a1,a0) 
  // |-- [    1] --| S2(x,w,j,i) += (    4.00000000) D1(j,a0) W114(x,w,a0,i) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ii]);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W114cc_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcc(symblockinfo, sa0^si));
    FC_FUNC(g_if_sigma_ccoo_ccoo_no32_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO32_X0_TYPE1_ERI_O)
      (sa0, ia0, si, ii, T2b.cptr(), V2_sym.cptr(), W114cc_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no32_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO32_X1_TYPE1_ERI_O)
      (sa0, ia0, si, ii, W114cc_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ii, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   33] -- 
  // |-- [    0] --| W115(x,w,a0,i) += (    1.00000000) V2(i,w,c0,a1) T2(x,c0,a0,a1) 
  // |-- [    1] --| S2(x,w,j,i) += (   -8.00000000) D1(j,a0) W115(x,w,a0,i) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ii]);
  orz::DTensor W115cca_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcca(symblockinfo, si));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no33_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO33_X0_TYPE1_ERI_O)
      (sa1, ia1, si, ii, T2b.cptr(), V2_sym.cptr(), W115cca_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  FC_FUNC(g_if_sigma_ccoo_ccoo_no33_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO33_X1_TYPE1_ERI_O)
    (si, ii, W115cca_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ii, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   34] -- 
  // |-- [    0] --| W116(w,x,a0,i) += (    1.00000000) V2(i,x,c0,a1) T2(w,c0,a0,a1) 
  // |-- [    1] --| S2(x,w,j,i) += (    4.00000000) D1(j,a0) W116(w,x,a0,i) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ii]);
  orz::DTensor W116cca_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcca(symblockinfo, si));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no34_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO34_X0_TYPE1_ERI_O)
      (sa1, ia1, si, ii, T2b.cptr(), V2_sym.cptr(), W116cca_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  FC_FUNC(g_if_sigma_ccoo_ccoo_no34_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO34_X1_TYPE1_ERI_O)
    (si, ii, W116cca_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ii, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   35] -- 
  // |-- [    0] --| W117(j,a0,i,a1) += (    1.00000000) V2(i,a2,a3,a1) D2(j,a3,a2,a0) 
  // |-- [    1] --| S2(x,w,j,i) += (    2.00000000) T2(x,w,a1,a0) W117(j,a0,i,a1) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ii]);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W117aa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa0^si));
    FC_FUNC(g_if_sigma_ccoo_ccoo_no35_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO35_X0_TYPE1_ERI_O)
      (sa0, ia0, si, ii, V2_sym.cptr(), W117aa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no35_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO35_X1_TYPE1_ERI_O)
      (sa0, ia0, si, ii, T2b.cptr(), W117aa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ii, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   36] -- 
  // |-- [    0] --| W118(j,a0,i,a1) += (    1.00000000) V2(i,a2,a3,a1) D2(j,a0,a2,a3) 
  // |-- [    1] --| S2(x,w,j,i) += (    2.00000000) T2(x,w,a0,a1) W118(j,a0,i,a1) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ii]);
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W118aa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, si^sa1));
    FC_FUNC(g_if_sigma_ccoo_ccoo_no36_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO36_X0_TYPE1_ERI_O)
      (sa1, ia1, si, ii, V2_sym.cptr(), W118aa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no36_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO36_X1_TYPE1_ERI_O)
      (sa1, ia1, si, ii, T2b.cptr(), W118aa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  retval.acc_amp2(ii, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   37] -- 
  // |-- [    0] --| W127(j,a0) += (    1.00000000) V2(a0,a2,a3,a1) D2(j,a2,a3,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) T2(x,w,i,a0) W127(j,a0) 
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
    double W127_sigma_ccoo_ccoo(0);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no37_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO37_X0_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, V2_sym.cptr(), &W127_sigma_ccoo_ccoo, nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no37_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO37_X1_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, T2b.cptr(), &W127_sigma_ccoo_ccoo, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   38] -- 
  // |-- [    0] --| W128(j,a0) += (    1.00000000) V2(a0,a2,a3,a1) D2(j,a2,a3,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (   -4.00000000) T2(w,x,i,a0) W128(j,a0) 
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
    double W128_sigma_ccoo_ccoo(0);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no38_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO38_X0_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, V2_sym.cptr(), &W128_sigma_ccoo_ccoo, nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no38_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO38_X1_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, T2b.cptr(), &W128_sigma_ccoo_ccoo, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   39] -- 
  // |-- [    0] --| W129(x,w,a0,i) += (    1.00000000) V2(i,a1,w,c0) T2(x,c0,a1,a0) 
  // |-- [    1] --| S2(x,w,j,i) += (   -2.00000000) D1(j,a0) W129(x,w,a0,i) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ii]);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W129cc_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcc(symblockinfo, sa0^si));
    FC_FUNC(g_if_sigma_ccoo_ccoo_no39_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO39_X0_TYPE1_ERI_O)
      (sa0, ia0, si, ii, T2b.cptr(), V2_sym.cptr(), W129cc_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no39_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO39_X1_TYPE1_ERI_O)
      (sa0, ia0, si, ii, W129cc_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ii, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   40] -- 
  // |-- [    0] --| W130(w,x,a0,i) += (    1.00000000) V2(i,a1,x,c0) T2(w,c0,a1,a0) 
  // |-- [    1] --| S2(x,w,j,i) += (    4.00000000) D1(j,a0) W130(w,x,a0,i) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ii]);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W130cc_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcc(symblockinfo, sa0^si));
    FC_FUNC(g_if_sigma_ccoo_ccoo_no40_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO40_X0_TYPE1_ERI_O)
      (sa0, ia0, si, ii, T2b.cptr(), V2_sym.cptr(), W130cc_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no40_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO40_X1_TYPE1_ERI_O)
      (sa0, ia0, si, ii, W130cc_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ii, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   41] -- 
  // |-- [    0] --| W131(x,w,a0,i) += (    1.00000000) V2(i,a1,w,c0) T2(x,c0,a0,a1) 
  // |-- [    1] --| S2(x,w,j,i) += (    4.00000000) D1(j,a0) W131(x,w,a0,i) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ii]);
  orz::DTensor W131cca_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcca(symblockinfo, si));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no41_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO41_X0_TYPE1_ERI_O)
      (sa1, ia1, si, ii, T2b.cptr(), V2_sym.cptr(), W131cca_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  FC_FUNC(g_if_sigma_ccoo_ccoo_no41_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO41_X1_TYPE1_ERI_O)
    (si, ii, W131cca_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ii, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   42] -- 
  // |-- [    0] --| W132(w,x,a0,i) += (    1.00000000) V2(i,a1,x,c0) T2(w,c0,a0,a1) 
  // |-- [    1] --| S2(x,w,j,i) += (   -2.00000000) D1(j,a0) W132(w,x,a0,i) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ii]);
  orz::DTensor W132cca_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcca(symblockinfo, si));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no42_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO42_X0_TYPE1_ERI_O)
      (sa1, ia1, si, ii, T2b.cptr(), V2_sym.cptr(), W132cca_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  FC_FUNC(g_if_sigma_ccoo_ccoo_no42_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO42_X1_TYPE1_ERI_O)
    (si, ii, W132cca_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ii, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   43] -- 
  // |-- [    0] --| W133(j,a0,i,a2) += (    1.00000000) V2(i,a2,a3,a1) D2(j,a0,a3,a1) 
  // |-- [    1] --| S2(x,w,j,i) += (    2.00000000) T2(x,w,a2,a0) W133(j,a0,i,a2) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ii]);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W133aa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa0^si));
    FC_FUNC(g_if_sigma_ccoo_ccoo_no43_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO43_X0_TYPE1_ERI_O)
      (sa0, ia0, si, ii, V2_sym.cptr(), W133aa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no43_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO43_X1_TYPE1_ERI_O)
      (sa0, ia0, si, ii, T2b.cptr(), W133aa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ii, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   44] -- 
  // |-- [    0] --| W134(j,a0,i,a2) += (    1.00000000) V2(i,a2,a3,a1) D2(j,a0,a3,a1) 
  // |-- [    1] --| S2(x,w,j,i) += (   -4.00000000) T2(x,w,a0,a2) W134(j,a0,i,a2) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ii]);
  for(int sa2 = 0;sa2 < nir;++sa2){ 
  for(int ia2 = symblockinfo.psym()(sa2,I_O,I_BEGIN);ia2 <= symblockinfo.psym()(sa2,I_O,I_END);++ia2){ 
    T2b = T2.get_amp2(ia2);
    orz::DTensor W134aa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, si^sa2));
    FC_FUNC(g_if_sigma_ccoo_ccoo_no44_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO44_X0_TYPE1_ERI_O)
      (sa2, ia2, si, ii, V2_sym.cptr(), W134aa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no44_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO44_X1_TYPE1_ERI_O)
      (sa2, ia2, si, ii, T2b.cptr(), W134aa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a2, "active"] [notNeeded]
  } // End ia2
  } // End sa2
  retval.acc_amp2(ii, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   45] -- 
  // |-- [    0] --| W135(w,x,a0,j) += (    1.00000000) V2(j,x,c0,a1) T2(w,c0,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    4.00000000) D1(i,a0) W135(w,x,a0,j) 
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
    orz::DTensor W135cc_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcc(symblockinfo, sa0^sj));
    FC_FUNC(g_if_sigma_ccoo_ccoo_no45_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO45_X0_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, T2b.cptr(), V2_sym.cptr(), W135cc_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no45_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO45_X1_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, W135cc_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   46] -- 
  // |-- [    0] --| W136(x,w,a0,j) += (    1.00000000) V2(j,w,c0,a1) T2(x,c0,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) D1(i,a0) W136(x,w,a0,j) 
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
    orz::DTensor W136cc_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcc(symblockinfo, sa0^sj));
    FC_FUNC(g_if_sigma_ccoo_ccoo_no46_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO46_X0_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, T2b.cptr(), V2_sym.cptr(), W136cc_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no46_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO46_X1_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, W136cc_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   47] -- 
  // |-- [    0] --| W137(x,w,a0,j) += (    1.00000000) V2(j,w,c0,a1) T2(x,c0,a0,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (    4.00000000) D1(i,a0) W137(x,w,a0,j) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W137cca_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcca(symblockinfo, sj));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no47_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO47_X0_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, T2b.cptr(), V2_sym.cptr(), W137cca_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  FC_FUNC(g_if_sigma_ccoo_ccoo_no47_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO47_X1_TYPE1_ERI_O)
    (sj, ij, W137cca_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   48] -- 
  // |-- [    0] --| W138(w,x,a0,j) += (    1.00000000) V2(j,x,c0,a1) T2(w,c0,a0,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (   -8.00000000) D1(i,a0) W138(w,x,a0,j) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W138cca_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcca(symblockinfo, sj));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no48_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO48_X0_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, T2b.cptr(), V2_sym.cptr(), W138cca_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  FC_FUNC(g_if_sigma_ccoo_ccoo_no48_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO48_X1_TYPE1_ERI_O)
    (sj, ij, W138cca_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   49] -- 
  // |-- [    0] --| W139(i,a0,j,a1) += (    1.00000000) V2(j,a2,a3,a1) D2(i,a0,a2,a3) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) T2(x,w,a1,a0) W139(i,a0,j,a1) 
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
    orz::DTensor W139aa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa0^sj));
    FC_FUNC(g_if_sigma_ccoo_ccoo_no49_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO49_X0_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, V2_sym.cptr(), W139aa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no49_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO49_X1_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, T2b.cptr(), W139aa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   50] -- 
  // |-- [    0] --| W140(i,a0,j,a1) += (    1.00000000) V2(j,a2,a3,a1) D2(i,a3,a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) T2(x,w,a0,a1) W140(i,a0,j,a1) 
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
    orz::DTensor W140aa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sa1));
    FC_FUNC(g_if_sigma_ccoo_ccoo_no50_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO50_X0_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, V2_sym.cptr(), W140aa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no50_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO50_X1_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, T2b.cptr(), W140aa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   51] -- 
  // |-- [    0] --| W149(i,a0) += (    1.00000000) V2(a0,a2,a3,a1) D2(i,a2,a3,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (   -4.00000000) T2(x,w,j,a0) W149(i,a0) 
  int sa0(s_eri);
  int ia0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia0);
  orz::DTensor W149a_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xa(symblockinfo, sa0));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no51_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO51_X0_TYPE1_ERI_O)
    (sa0, ia0, V2_sym.cptr(), W149a_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no51_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO51_X1_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, T2b.cptr(), W149a_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   52] -- 
  // |-- [    0] --| W150(i,a0) += (    1.00000000) V2(a0,a2,a3,a1) D2(i,a2,a3,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) T2(x,w,a0,j) W150(i,a0) 
  int sa0(s_eri);
  int ia0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W150a_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xa(symblockinfo, sa0));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no52_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO52_X0_TYPE1_ERI_O)
    (sa0, ia0, V2_sym.cptr(), W150a_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no52_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO52_X1_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, T2b.cptr(), W150a_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   53] -- 
  // |-- [    0] --| W151(x,w,a0,j) += (    1.00000000) V2(j,a1,w,c0) T2(x,c0,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    4.00000000) D1(i,a0) W151(x,w,a0,j) 
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
    orz::DTensor W151cc_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcc(symblockinfo, sa0^sj));
    FC_FUNC(g_if_sigma_ccoo_ccoo_no53_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO53_X0_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, T2b.cptr(), V2_sym.cptr(), W151cc_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no53_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO53_X1_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, W151cc_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   54] -- 
  // |-- [    0] --| W152(w,x,a0,j) += (    1.00000000) V2(j,a1,x,c0) T2(w,c0,a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) D1(i,a0) W152(w,x,a0,j) 
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
    orz::DTensor W152cc_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcc(symblockinfo, sa0^sj));
    FC_FUNC(g_if_sigma_ccoo_ccoo_no54_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO54_X0_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, T2b.cptr(), V2_sym.cptr(), W152cc_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no54_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO54_X1_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, W152cc_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   55] -- 
  // |-- [    0] --| W153(x,w,a0,j) += (    1.00000000) V2(j,a1,w,c0) T2(x,c0,a0,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (   -2.00000000) D1(i,a0) W153(x,w,a0,j) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W153cca_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcca(symblockinfo, sj));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no55_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO55_X0_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, T2b.cptr(), V2_sym.cptr(), W153cca_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  FC_FUNC(g_if_sigma_ccoo_ccoo_no55_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO55_X1_TYPE1_ERI_O)
    (sj, ij, W153cca_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   56] -- 
  // |-- [    0] --| W154(w,x,a0,j) += (    1.00000000) V2(j,a1,x,c0) T2(w,c0,a0,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (    4.00000000) D1(i,a0) W154(w,x,a0,j) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W154cca_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xcca(symblockinfo, sj));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no56_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO56_X0_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, T2b.cptr(), V2_sym.cptr(), W154cca_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  FC_FUNC(g_if_sigma_ccoo_ccoo_no56_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO56_X1_TYPE1_ERI_O)
    (sj, ij, W154cca_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   57] -- 
  // |-- [    0] --| W155(i,a0,j,a2) += (    1.00000000) V2(j,a2,a3,a1) D2(i,a0,a3,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (   -4.00000000) T2(x,w,a2,a0) W155(i,a0,j,a2) 
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
    orz::DTensor W155aa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa0^sj));
    FC_FUNC(g_if_sigma_ccoo_ccoo_no57_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO57_X0_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, V2_sym.cptr(), W155aa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no57_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO57_X1_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, T2b.cptr(), W155aa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   58] -- 
  // |-- [    0] --| W156(i,a0,j,a2) += (    1.00000000) V2(j,a2,a3,a1) D2(i,a0,a3,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) T2(x,w,a0,a2) W156(i,a0,j,a2) 
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
    orz::DTensor W156aa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sa2));
    FC_FUNC(g_if_sigma_ccoo_ccoo_no58_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO58_X0_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, V2_sym.cptr(), W156aa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no58_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO58_X1_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, T2b.cptr(), W156aa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a2, "active"] [notNeeded]
  } // End ia2
  } // End sa2
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   59] -- 
  // |-- [    0] --| W157(i,a0) += (    1.00000000) V2(i,a1,a2,a0) D1(a2,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (   -4.00000000) T2(w,x,a0,j) W157(i,a0) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W157a_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xa(symblockinfo, si));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no59_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO59_X0_TYPE1_ERI_O)
    (si, ii, V2_sym.cptr(), W157a_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no59_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO59_X1_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), W157a_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   60] -- 
  // |-- [    0] --| W158(i,a0) += (    1.00000000) V2(i,a1,a2,a0) D1(a2,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) T2(x,w,a0,j) W158(i,a0) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W158a_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xa(symblockinfo, si));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no60_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO60_X0_TYPE1_ERI_O)
    (si, ii, V2_sym.cptr(), W158a_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no60_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO60_X1_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), W158a_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   61] -- 
  // |-- [    0] --| W159(a0,j,i,a2) += (    1.00000000) V2(j,a2,i,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -4.00000000) T2(x,w,a2,a0) W159(a0,j,i,a2) 
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
    orz::DTensor W159aa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa0^sj));
    FC_FUNC(g_if_sigma_ccoo_ccoo_no61_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO61_X0_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, V2_sym.cptr(), W159aa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no61_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO61_X1_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, T2b.cptr(), W159aa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   62] -- 
  // |-- [    0] --| W160(a0,j,i,a2) += (    1.00000000) V2(j,a2,i,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) T2(x,w,a0,a2) W160(a0,j,i,a2) 
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
    orz::DTensor W160aa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sa2));
    FC_FUNC(g_if_sigma_ccoo_ccoo_no62_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO62_X0_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, V2_sym.cptr(), W160aa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no62_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO62_X1_TYPE1_ERI_O)
      (sa2, ia2, sj, ij, T2b.cptr(), W160aa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a2, "active"] [notNeeded]
  } // End ia2
  } // End sa2
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   63] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   16.00000000) V2(j,x,c0,a0) T2(w,c0,i,a0) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no63_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO63_X0_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   64] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -8.00000000) V2(j,w,c0,a0) T2(x,c0,i,a0) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no64_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO64_X0_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   65] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    4.00000000) V2(j,w,c0,a0) T2(x,c0,a0,i) 
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
    FC_FUNC(g_if_sigma_ccoo_ccoo_no65_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO65_X0_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   66] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -8.00000000) V2(j,x,c0,a0) T2(w,c0,a0,i) 
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
    FC_FUNC(g_if_sigma_ccoo_ccoo_no66_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO66_X0_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   67] -- 
  // |-- [    0] --| W161(j,a0) += (    1.00000000) V2(j,a1,a2,a0) D1(a2,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) T2(x,w,i,a0) W161(j,a0) 
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
    double W161_sigma_ccoo_ccoo(0);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no67_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO67_X0_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, V2_sym.cptr(), &W161_sigma_ccoo_ccoo, nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no67_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO67_X1_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, T2b.cptr(), &W161_sigma_ccoo_ccoo, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   68] -- 
  // |-- [    0] --| W162(j,a0) += (    1.00000000) V2(j,a1,a2,a0) D1(a2,a1) 
  // |-- [    1] --| S2(w,x,i,j) += (   -4.00000000) T2(x,w,a0,i) W162(j,a0) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W162a_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xa(symblockinfo, sj));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no68_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO68_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W162a_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no68_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO68_X1_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), W162a_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   69] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -8.00000000) V2(j,a0,x,c0) T2(w,c0,i,a0) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no69_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO69_X0_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   70] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    4.00000000) V2(j,a0,w,c0) T2(x,c0,i,a0) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no70_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO70_X0_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   71] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    4.00000000) V2(j,a0,x,c0) T2(w,c0,a0,i) 
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
    FC_FUNC(g_if_sigma_ccoo_ccoo_no71_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO71_X0_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   72] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -8.00000000) V2(j,a0,w,c0) T2(x,c0,a0,i) 
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
    FC_FUNC(g_if_sigma_ccoo_ccoo_no72_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO72_X0_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   73] -- 
  // |-- [    0] --| W163(j,a1) += (    1.00000000) V2(j,a1,a2,a0) D1(a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    8.00000000) T2(x,w,a1,i) W163(j,a1) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W163a_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xa(symblockinfo, sj));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no73_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO73_X0_TYPE1_ERI_O)
    (sj, ij, V2_sym.cptr(), W163a_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no73_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO73_X1_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), W163a_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   74] -- 
  // |-- [    0] --| W164(j,a1) += (    1.00000000) V2(j,a1,a2,a0) D1(a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -4.00000000) T2(x,w,i,a1) W164(j,a1) 
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
    double W164_sigma_ccoo_ccoo(0);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no74_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO74_X0_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, V2_sym.cptr(), &W164_sigma_ccoo_ccoo, nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no74_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO74_X1_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, T2b.cptr(), &W164_sigma_ccoo_ccoo, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   75] -- 
  // |-- [    0] --| W165(a0,j,i,a1) += (    1.00000000) V2(j,a2,i,a1) D1(a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) T2(x,w,a1,a0) W165(a0,j,i,a1) 
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
    orz::DTensor W165aa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa0^sj));
    FC_FUNC(g_if_sigma_ccoo_ccoo_no75_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO75_X0_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, V2_sym.cptr(), W165aa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no75_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO75_X1_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, T2b.cptr(), W165aa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   76] -- 
  // |-- [    0] --| W166(a0,j,i,a1) += (    1.00000000) V2(j,a2,i,a1) D1(a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -4.00000000) T2(x,w,a0,a1) W166(a0,j,i,a1) 
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
    orz::DTensor W166aa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sa1));
    FC_FUNC(g_if_sigma_ccoo_ccoo_no76_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO76_X0_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, V2_sym.cptr(), W166aa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no76_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO76_X1_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, T2b.cptr(), W166aa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   77] -- 
  // |-- [    0] --| W167(i,a1) += (    1.00000000) V2(i,a1,a2,a0) D1(a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    8.00000000) T2(w,x,a1,j) W167(i,a1) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W167a_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xa(symblockinfo, si));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no77_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO77_X0_TYPE1_ERI_O)
    (si, ii, V2_sym.cptr(), W167a_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no77_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO77_X1_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), W167a_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   78] -- 
  // |-- [    0] --| W168(i,a1) += (    1.00000000) V2(i,a1,a2,a0) D1(a2,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -4.00000000) T2(x,w,a1,j) W168(i,a1) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W168a_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xa(symblockinfo, si));
  FC_FUNC(g_if_sigma_ccoo_ccoo_no78_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO78_X0_TYPE1_ERI_O)
    (si, ii, V2_sym.cptr(), W168a_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no78_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO78_X1_TYPE1_ERI_O)
      (si, ii, sj, ij, T2b.cptr(), W168a_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   79] -- 
  // |-- [    0] --| W169(j,i,a1,a0) += (    1.00000000) V2(a1,a3,a2,a0) D2(j,a3,i,a2) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) T2(w,x,a0,a1) W169(j,i,a1,a0) 
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
    orz::DTensor W169aa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sa1));
    FC_FUNC(g_if_sigma_ccoo_ccoo_no79_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO79_X0_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, V2_sym.cptr(), W169aa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no79_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO79_X1_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, T2b.cptr(), W169aa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   80] -- 
  // |-- [    0] --| W170(j,i,a1,a0) += (    1.00000000) V2(i,a1,a2,a0) D1(j,a2) 
  // |-- [    1] --| S2(x,w,j,i) += (    2.00000000) T2(x,w,a1,a0) W170(j,i,a1,a0) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ii]);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W170aa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, si^sa0));
    FC_FUNC(g_if_sigma_ccoo_ccoo_no80_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO80_X0_TYPE1_ERI_O)
      (sa0, ia0, si, ii, V2_sym.cptr(), W170aa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no80_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO80_X1_TYPE1_ERI_O)
      (sa0, ia0, si, ii, T2b.cptr(), W170aa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ii, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   81] -- 
  // |-- [    0] --| W171(j,i,a1,a0) += (    1.00000000) V2(i,a1,a2,a0) D1(j,a2) 
  // |-- [    1] --| S2(x,w,j,i) += (   -4.00000000) T2(x,w,a0,a1) W171(j,i,a1,a0) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ii]);
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W171aa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, si^sa1));
    FC_FUNC(g_if_sigma_ccoo_ccoo_no81_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO81_X0_TYPE1_ERI_O)
      (sa1, ia1, si, ii, V2_sym.cptr(), W171aa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no81_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO81_X1_TYPE1_ERI_O)
      (sa1, ia1, si, ii, T2b.cptr(), W171aa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  retval.acc_amp2(ii, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   82] -- 
  // |-- [    0] --| W172(i,j,a1,a0) += (    1.00000000) V2(j,a1,a2,a0) D1(i,a2) 
  // |-- [    1] --| S2(w,x,i,j) += (   -4.00000000) T2(x,w,a1,a0) W172(i,j,a1,a0) 
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
    orz::DTensor W172aa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sa0));
    FC_FUNC(g_if_sigma_ccoo_ccoo_no82_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO82_X0_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, V2_sym.cptr(), W172aa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no82_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO82_X1_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, T2b.cptr(), W172aa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   83] -- 
  // |-- [    0] --| W173(i,j,a1,a0) += (    1.00000000) V2(j,a1,a2,a0) D1(i,a2) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) T2(x,w,a0,a1) W173(i,j,a1,a0) 
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
    orz::DTensor W173aa_sigma_ccoo_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sa1));
    FC_FUNC(g_if_sigma_ccoo_ccoo_no83_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO83_X0_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, V2_sym.cptr(), W173aa_sigma_ccoo_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no83_x1_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO83_X1_TYPE1_ERI_O)
      (sa1, ia1, sj, ij, T2b.cptr(), W173aa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   84] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    8.00000000) V2(j,a1,i,a0) T2(x,w,a1,a0) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no84_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO84_X0_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   85] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -4.00000000) V2(j,a0,i,a1) T2(x,w,a1,a0) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no85_x0_type1_eri_o,G_IF_SIGMA_CCOO_CCOO_NO85_X0_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

//-@type(1).contraction(end)

  } // End myrank
  }
  }
  orz::world().barrier();

//-@type(2).contraction(begin)
  // -- Title : sigma_ccoo_ccoo
  //*-- Entering to take the type 2 contractions --*//

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| S2(x,w,j,i) += (    1.00000000) T2(w,x,a0,i) W53(j,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    S2b = orz::DTensor(retval.namps_iamp()[ii]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no0_x0_type2_eri_o,G_IF_SIGMA_CCOO_CCOO_NO0_X0_TYPE2_ERI_O)
      (si, ii, T2b.cptr(), W53aa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    retval.acc_amp2(ii, S2b);
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| S2(x,w,j,i) += (   -2.00000000) T2(x,w,a0,i) W54(j,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    S2b = orz::DTensor(retval.namps_iamp()[ii]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no1_x0_type2_eri_o,G_IF_SIGMA_CCOO_CCOO_NO1_X0_TYPE2_ERI_O)
      (si, ii, T2b.cptr(), W54aa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    retval.acc_amp2(ii, S2b);
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -2.00000000) T2(w,x,a0,j) W70(i,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no2_x0_type2_eri_o,G_IF_SIGMA_CCOO_CCOO_NO2_X0_TYPE2_ERI_O)
      (sj, ij, T2b.cptr(), W70aa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    1.00000000) T2(x,w,a0,j) W71(i,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no3_x0_type2_eri_o,G_IF_SIGMA_CCOO_CCOO_NO3_X0_TYPE2_ERI_O)
      (sj, ij, T2b.cptr(), W71aa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    4.00000000) T2(w,x,i,j) W101() 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no4_x0_type2_eri_o,G_IF_SIGMA_CCOO_CCOO_NO4_X0_TYPE2_ERI_O)
      (sj, ij, T2b.cptr(), &W101_sigma_ccoo_ccoo, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -2.00000000) T2(x,w,i,j) W102() 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no5_x0_type2_eri_o,G_IF_SIGMA_CCOO_CCOO_NO5_X0_TYPE2_ERI_O)
      (sj, ij, T2b.cptr(), &W102_sigma_ccoo_ccoo, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -2.00000000) D2(j,a2,i,a0) W103(w,x,a0,a2) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no6_x0_type2_eri_o,G_IF_SIGMA_CCOO_CCOO_NO6_X0_TYPE2_ERI_O)
      (sj, ij, W103ccaa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -2.00000000) D2(j,a0,i,a2) W104(x,w,a0,a2) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no7_x0_type2_eri_o,G_IF_SIGMA_CCOO_CCOO_NO7_X0_TYPE2_ERI_O)
      (sj, ij, W104ccaa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    4.00000000) D2(j,a0,i,a2) W105(x,w,a0,a2) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no8_x0_type2_eri_o,G_IF_SIGMA_CCOO_CCOO_NO8_X0_TYPE2_ERI_O)
      (sj, ij, W105ccaa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    4.00000000) D2(j,a2,i,a0) W106(w,x,a0,a2) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no9_x0_type2_eri_o,G_IF_SIGMA_CCOO_CCOO_NO9_X0_TYPE2_ERI_O)
      (sj, ij, W106ccaa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -2.00000000) D2(j,a0,i,a2) W107(w,x,a0,a2) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no10_x0_type2_eri_o,G_IF_SIGMA_CCOO_CCOO_NO10_X0_TYPE2_ERI_O)
      (sj, ij, W107ccaa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -2.00000000) D2(j,a2,i,a0) W108(x,w,a0,a2) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no11_x0_type2_eri_o,G_IF_SIGMA_CCOO_CCOO_NO11_X0_TYPE2_ERI_O)
      (sj, ij, W108ccaa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -2.00000000) D2(j,a0,i,a2) W109(x,w,a0,a2) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no12_x0_type2_eri_o,G_IF_SIGMA_CCOO_CCOO_NO12_X0_TYPE2_ERI_O)
      (sj, ij, W109ccaa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -2.00000000) D2(j,a2,i,a0) W110(w,x,a0,a2) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccoo_no13_x0_type2_eri_o,G_IF_SIGMA_CCOO_CCOO_NO13_X0_TYPE2_ERI_O)
      (sj, ij, W110ccaa_sigma_ccoo_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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

//-@D4C.contractions(end)

  // Do timing!
  orz::mr::Femto::my_timer.push_back(boost::make_tuple("sigma_ccoo_ccoo", time_sigma.elapsed_cputime(), time_sigma.elapsed_wallclocktime()));
  orz::mr::Femto::file_timing << "* " << boost::format("%20s : %10.7f %10.7f ") % orz::mr::Femto::my_timer.back().get<0>() % orz::mr::Femto::my_timer.back().get<1>() % orz::mr::Femto::my_timer.back().get<2>() << endl;
  flush(orz::mr::Femto::file_timing);

  return retval; 
} 
