                                                                                
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
#include <sci/icmr/Femto/elems/c_sigma_ocov_ooov.h>                                  
                                                                                
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

//                                   Generated date : Sun Apr 20 10:26:10 2014

                                                                                
// ***************************************************************************  
// orz::mr::mrci                                                                
// ***************************************************************************  
									    /*!        
   @brief CT input                                                              
									     */        
                                                                                
orz::mr::BareAmpPack orz::mr::Femto::sigma_ocov_ooov(const orz::mr::Input &ctinp,                                    
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
  // |-- [    0] --| W0(j,i,a3,a) += (    1.00000000) T2(a1,a2,a0,a) D3(j,i,a2,a3,a1,a0) 
  // |-- [    1] --| S2(i,w,j,a) += (   -1.00000000) Fc1(w,a3) W0(j,i,a3,a) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W0aaa_sigma_ocov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ocov_ooov_no0_x0_type0_noeri,G_IF_SIGMA_OCOV_OOOV_NO0_X0_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W0aaa_sigma_ocov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ocov_ooov_no0_x1_type0_noeri,G_IF_SIGMA_OCOV_OOOV_NO0_X1_TYPE0_NOERI)
      (sa, ia, W0aaa_sigma_ocov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W1(i,a) += (    1.00000000) T2(a0,a1,a2,a) D2(i,a1,a2,a0) 
  // |-- [    1] --| S2(i,w,j,a) += (   -1.00000000) Fc1(w,j) W1(i,a) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W1a_sigma_ocov_ooov(orz::mr::sizeof_sympack_Xa(symblockinfo, sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ocov_ooov_no1_x0_type0_noeri,G_IF_SIGMA_OCOV_OOOV_NO1_X0_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W1a_sigma_ocov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ocov_ooov_no1_x1_type0_noeri,G_IF_SIGMA_OCOV_OOOV_NO1_X1_TYPE0_NOERI)
      (sa, ia, W1a_sigma_ocov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W2(i,a2,j,a) += (    1.00000000) T2(a1,a0,j,a) D2(i,a1,a2,a0) 
  // |-- [    1] --| S2(i,w,j,a) += (   -1.00000000) Fc1(w,a2) W2(i,a2,j,a) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W2aaa_sigma_ocov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ocov_ooov_no2_x0_type0_noeri,G_IF_SIGMA_OCOV_OOOV_NO2_X0_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W2aaa_sigma_ocov_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ocov_ooov_no2_x1_type0_noeri,G_IF_SIGMA_OCOV_OOOV_NO2_X1_TYPE0_NOERI)
      (sa, ia, W2aaa_sigma_ocov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- Title : sigma_ocov_ooov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W5(w,i,a1,a0) += (    1.00000000) V2(w,a3,a4,a2) D3(i,a1,a4,a2,a3,a0) 
  // |-- [    1] --| S2(i,w,j,a) += (   -1.00000000) T2(a1,a0,j,a) W5(w,i,a1,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W5aaa_sigma_ocov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ocov_ooov_no0_x0_type1_eri_c,G_IF_SIGMA_OCOV_OOOV_NO0_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W5aaa_sigma_ocov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ocov_ooov_no0_x1_type1_eri_c,G_IF_SIGMA_OCOV_OOOV_NO0_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), W5aaa_sigma_ocov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W6(w,j,a1,a0,i,a2) += (    1.00000000) V2(w,a3,a4,a2) D3(j,i,a1,a4,a0,a3) 
  // |-- [    1] --| S2(i,w,j,a) += (   -1.00000000) T2(a1,a0,a2,a) W6(w,j,a1,a0,i,a2) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W6aaaaa_sigma_ocov_ooov(orz::mr::sizeof_sympack_Xaaaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ocov_ooov_no1_x0_type1_eri_c,G_IF_SIGMA_OCOV_OOOV_NO1_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W6aaaaa_sigma_ocov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ocov_ooov_no1_x1_type1_eri_c,G_IF_SIGMA_OCOV_OOOV_NO1_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), W6aaaaa_sigma_ocov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- Title : sigma_ocov_ooov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W3(i,a4,a3,a) += (    1.00000000) T2(a0,a1,a2,a) D3(i,a1,a4,a3,a2,a0) 
  // |-- [    1] --| S2(i,w,j,a) += (   -1.00000000) V2(a4,a3,w,j) W3(i,a4,a3,a) 
  int sa4(s_eri);
  int ia4(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W3aa_sigma_ocov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa4^sa));
    FC_FUNC(g_if_sigma_ocov_ooov_no0_x0_type1_eri_o,G_IF_SIGMA_OCOV_OOOV_NO0_X0_TYPE1_ERI_O)
      (sa, ia, sa4, ia4, T2b.cptr(), W3aa_sigma_ocov_ooov.cptr(), nir, nsym, psym, &flops);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ocov_ooov_no0_x1_type1_eri_o,G_IF_SIGMA_OCOV_OOOV_NO0_X1_TYPE1_ERI_O)
      (sa, ia, sa4, ia4, V2_sym.cptr(), W3aa_sigma_ocov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W4(i,a4,a3,a) += (    1.00000000) T2(a0,a1,a2,a) D3(i,a3,a4,a1,a2,a0) 
  // |-- [    1] --| S2(i,w,j,a) += (   -1.00000000) V2(a4,w,j,a3) W4(i,a4,a3,a) 
  int sa4(s_eri);
  int ia4(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W4aa_sigma_ocov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa4^sa));
    FC_FUNC(g_if_sigma_ocov_ooov_no1_x0_type1_eri_o,G_IF_SIGMA_OCOV_OOOV_NO1_X0_TYPE1_ERI_O)
      (sa, ia, sa4, ia4, T2b.cptr(), W4aa_sigma_ocov_ooov.cptr(), nir, nsym, psym, &flops);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ocov_ooov_no1_x1_type1_eri_o,G_IF_SIGMA_OCOV_OOOV_NO1_X1_TYPE1_ERI_O)
      (sa, ia, sa4, ia4, V2_sym.cptr(), W4aa_sigma_ocov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W7(i,a3,a2,a) += (    1.00000000) T2(a1,a0,a,a2) D2(i,a1,a3,a0) 
  // |-- [    1] --| S2(i,w,j,a) += (   -1.00000000) V2(a2,a3,w,j) W7(i,a3,a2,a) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  T2b = T2.get_amp2(ia2);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    orz::DTensor W7aa_sigma_ocov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa2^sa));
    FC_FUNC(g_if_sigma_ocov_ooov_no2_x0_type1_eri_o,G_IF_SIGMA_OCOV_OOOV_NO2_X0_TYPE1_ERI_O)
      (sa, ia, sa2, ia2, T2b.cptr(), W7aa_sigma_ocov_ooov.cptr(), nir, nsym, psym, &flops);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ocov_ooov_no2_x1_type1_eri_o,G_IF_SIGMA_OCOV_OOOV_NO2_X1_TYPE1_ERI_O)
      (sa, ia, sa2, ia2, V2_sym.cptr(), W7aa_sigma_ocov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W8(i,a3,a2,a) += (    1.00000000) T2(a1,a0,a2,a) D2(i,a1,a3,a0) 
  // |-- [    1] --| S2(i,w,j,a) += (   -1.00000000) V2(a3,w,j,a2) W8(i,a3,a2,a) 
  int sa3(s_eri);
  int ia3(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W8aa_sigma_ocov_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa3^sa));
    FC_FUNC(g_if_sigma_ocov_ooov_no3_x0_type1_eri_o,G_IF_SIGMA_OCOV_OOOV_NO3_X0_TYPE1_ERI_O)
      (sa, ia, sa3, ia3, T2b.cptr(), W8aa_sigma_ocov_ooov.cptr(), nir, nsym, psym, &flops);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ocov_ooov_no3_x1_type1_eri_o,G_IF_SIGMA_OCOV_OOOV_NO3_X1_TYPE1_ERI_O)
      (sa, ia, sa3, ia3, V2_sym.cptr(), W8aa_sigma_ocov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- Title : sigma_ocov_ooov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W9(j,i,a3,v0) += (    1.00000000) T2(a1,a2,a0,v0) D3(j,i,a2,a3,a1,a0) 
  // |-- [    1] --| S2(i,w,j,a) += (    2.00000000) V2(v0,a3,w,a) W9(j,i,a3,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W9aaa_sigma_ocov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ocov_ooov_no0_x0_type1_eri_v,G_IF_SIGMA_OCOV_OOOV_NO0_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W9aaa_sigma_ocov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ocov_ooov_no0_x1_type1_eri_v,G_IF_SIGMA_OCOV_OOOV_NO0_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W9aaa_sigma_ocov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W10(j,i,a3,v0) += (    1.00000000) T2(a1,a2,a0,v0) D3(j,i,a2,a3,a1,a0) 
  // |-- [    1] --| S2(i,w,j,a) += (   -1.00000000) V2(v0,a,w,a3) W10(j,i,a3,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W10aaa_sigma_ocov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ocov_ooov_no1_x0_type1_eri_v,G_IF_SIGMA_OCOV_OOOV_NO1_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W10aaa_sigma_ocov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ocov_ooov_no1_x1_type1_eri_v,G_IF_SIGMA_OCOV_OOOV_NO1_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W10aaa_sigma_ocov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W11(i,v0) += (    1.00000000) T2(a0,a1,a2,v0) D2(i,a1,a2,a0) 
  // |-- [    1] --| S2(i,w,j,a) += (   -1.00000000) V2(v0,a,w,j) W11(i,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W11a_sigma_ocov_ooov(orz::mr::sizeof_sympack_Xa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ocov_ooov_no2_x0_type1_eri_v,G_IF_SIGMA_OCOV_OOOV_NO2_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W11a_sigma_ocov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ocov_ooov_no2_x1_type1_eri_v,G_IF_SIGMA_OCOV_OOOV_NO2_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W11a_sigma_ocov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W12(i,a2,j,v0) += (    1.00000000) T2(a1,a0,j,v0) D2(i,a1,a2,a0) 
  // |-- [    1] --| S2(i,w,j,a) += (    2.00000000) V2(v0,a2,w,a) W12(i,a2,j,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W12aaa_sigma_ocov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ocov_ooov_no3_x0_type1_eri_v,G_IF_SIGMA_OCOV_OOOV_NO3_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W12aaa_sigma_ocov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ocov_ooov_no3_x1_type1_eri_v,G_IF_SIGMA_OCOV_OOOV_NO3_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W12aaa_sigma_ocov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W13(i,a2,j,v0) += (    1.00000000) T2(a1,a0,j,v0) D2(i,a1,a2,a0) 
  // |-- [    1] --| S2(i,w,j,a) += (   -1.00000000) V2(v0,a,w,a2) W13(i,a2,j,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W13aaa_sigma_ocov_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ocov_ooov_no4_x0_type1_eri_v,G_IF_SIGMA_OCOV_OOOV_NO4_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W13aaa_sigma_ocov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ocov_ooov_no4_x1_type1_eri_v,G_IF_SIGMA_OCOV_OOOV_NO4_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W13aaa_sigma_ocov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W14(i,v0) += (    1.00000000) T2(a0,a1,a2,v0) D2(i,a1,a2,a0) 
  // |-- [    1] --| S2(i,w,j,a) += (    2.00000000) V2(v0,j,w,a) W14(i,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W14a_sigma_ocov_ooov(orz::mr::sizeof_sympack_Xa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ocov_ooov_no5_x0_type1_eri_v,G_IF_SIGMA_OCOV_OOOV_NO5_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W14a_sigma_ocov_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ocov_ooov_no5_x1_type1_eri_v,G_IF_SIGMA_OCOV_OOOV_NO5_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, V2_sym.cptr(), W14a_sigma_ocov_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- Title : sigma_ocov_ooov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| S2(i,w,j,a) += (   -1.00000000) T2(a1,a2,a0,a) C5(a1,a0,j,i,a2,w) 
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
    FC_FUNC(g_if_sigma_ocov_ooov_no0_x0_type1_d4c_c,G_IF_SIGMA_OCOV_OOOV_NO0_X0_TYPE1_D4C_C)
      (sa, ia, sw, iw, C5.cptr(), T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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

//-@loadD4C(c,end)

//-@D4C.contractions(end)

  // Do timing!
  orz::mr::Femto::my_timer.push_back(boost::make_tuple("sigma_ocov_ooov", time_sigma.elapsed_cputime(), time_sigma.elapsed_wallclocktime()));
  orz::mr::Femto::file_timing << "* " << boost::format("%20s : %10.7f %10.7f ") % orz::mr::Femto::my_timer.back().get<0>() % orz::mr::Femto::my_timer.back().get<1>() % orz::mr::Femto::my_timer.back().get<2>() << endl;
  flush(orz::mr::Femto::file_timing);

  return retval; 
} 
