                                                                                
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
#include <sci/icmr/Femto/elems/c_sigma_ccov_cooo.h>                                  
                                                                                
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
                                                                                
orz::mr::BareAmpPack orz::mr::Femto::sigma_ccov_cooo(const orz::mr::Input &ctinp,                                    
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
  // |-- [    0] --| W0(x,i) += (    1.00000000) T2(x,a2,a1,a0) D2(i,a1,a2,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) Fc1(w,a) W0(x,i) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W0ca_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccov_cooo_no0_x0_type0_noeri,G_IF_SIGMA_CCOV_COOO_NO0_X0_TYPE0_NOERI)
      (sa0, ia0, T2b.cptr(), W0ca_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_cooo_no0_x1_type0_noeri,G_IF_SIGMA_CCOV_COOO_NO0_X1_TYPE0_NOERI)
      (sa, ia, W0ca_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W1(w,i) += (    1.00000000) T2(w,a2,a1,a0) D2(i,a1,a2,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) Fc1(x,a) W1(w,i) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W1ca_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccov_cooo_no1_x0_type0_noeri,G_IF_SIGMA_CCOV_COOO_NO1_X0_TYPE0_NOERI)
      (sa0, ia0, T2b.cptr(), W1ca_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_cooo_no1_x1_type0_noeri,G_IF_SIGMA_CCOV_COOO_NO1_X1_TYPE0_NOERI)
      (sa, ia, W1ca_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W2(w,i) += (    1.00000000) T2(w,a1,i,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    4.00000000) Fc1(x,a) W2(w,i) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W2ca_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccov_cooo_no2_x0_type0_noeri,G_IF_SIGMA_CCOV_COOO_NO2_X0_TYPE0_NOERI)
      (sa0, ia0, T2b.cptr(), W2ca_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_cooo_no2_x1_type0_noeri,G_IF_SIGMA_CCOV_COOO_NO2_X1_TYPE0_NOERI)
      (sa, ia, W2ca_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W3(x,i) += (    1.00000000) T2(x,a1,i,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) Fc1(w,a) W3(x,i) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W3ca_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccov_cooo_no3_x0_type0_noeri,G_IF_SIGMA_CCOV_COOO_NO3_X0_TYPE0_NOERI)
      (sa0, ia0, T2b.cptr(), W3ca_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_cooo_no3_x1_type0_noeri,G_IF_SIGMA_CCOV_COOO_NO3_X1_TYPE0_NOERI)
      (sa, ia, W3ca_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W4(x,i) += (    1.00000000) T2(x,a1,a0,i) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) Fc1(w,a) W4(x,i) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    orz::DTensor W4c_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, si));
    FC_FUNC(g_if_sigma_ccov_cooo_no4_x0_type0_noeri,G_IF_SIGMA_CCOV_COOO_NO4_X0_TYPE0_NOERI)
      (si, ii, T2b.cptr(), W4c_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      S2b = orz::DTensor(retval.namps_iamp()[ia]);
      FC_FUNC(g_if_sigma_ccov_cooo_no4_x1_type0_noeri,G_IF_SIGMA_CCOV_COOO_NO4_X1_TYPE0_NOERI)
        (sa, ia, si, ii, W4c_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- No. [    5] -- 
  // |-- [    0] --| W5(w,i) += (    1.00000000) T2(w,a1,a0,i) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) Fc1(x,a) W5(w,i) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    orz::DTensor W5c_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, si));
    FC_FUNC(g_if_sigma_ccov_cooo_no5_x0_type0_noeri,G_IF_SIGMA_CCOV_COOO_NO5_X0_TYPE0_NOERI)
      (si, ii, T2b.cptr(), W5c_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      S2b = orz::DTensor(retval.namps_iamp()[ia]);
      FC_FUNC(g_if_sigma_ccov_cooo_no5_x1_type0_noeri,G_IF_SIGMA_CCOV_COOO_NO5_X1_TYPE0_NOERI)
        (sa, ia, si, ii, W5c_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
      retval.acc_amp2(ia, S2b);
    } // End ia
    } // End sa
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
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
  // -- Title : sigma_ccov_cooo

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W28(w,i,a3,a2) += (    1.00000000) T2(a1,w,a0,a2) D2(i,a3,a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) V2(a2,a3,x,a) W28(w,i,a3,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia2);
  orz::DTensor W28caa_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa2));
  FC_FUNC(g_if_sigma_ccov_cooo_no0_x0_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO0_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), W28caa_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_cooo_no0_x1_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO0_X1_TYPE1_ERI_O)
      (sa, ia, sa2, ia2, V2_sym.cptr(), W28caa_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W29(x,i,a3,a2) += (    1.00000000) T2(a1,x,a0,a2) D2(i,a3,a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) V2(a2,a3,w,a) W29(x,i,a3,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia2);
  orz::DTensor W29caa_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa2));
  FC_FUNC(g_if_sigma_ccov_cooo_no1_x0_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO1_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), W29caa_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_cooo_no1_x1_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO1_X1_TYPE1_ERI_O)
      (sa, ia, sa2, ia2, V2_sym.cptr(), W29caa_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W30(x,i,a2,a3) += (    1.00000000) T2(a1,x,a0,a3) D2(i,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) V2(a3,a,w,a2) W30(x,i,a2,a3) 
  int sa3(s_eri);
  int ia3(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia3);
  orz::DTensor W30caa_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa3));
  FC_FUNC(g_if_sigma_ccov_cooo_no2_x0_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO2_X0_TYPE1_ERI_O)
    (sa3, ia3, T2b.cptr(), W30caa_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_cooo_no2_x1_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO2_X1_TYPE1_ERI_O)
      (sa, ia, sa3, ia3, V2_sym.cptr(), W30caa_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W31(w,i,a2,a3) += (    1.00000000) T2(a1,w,a0,a3) D2(i,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) V2(a3,a,x,a2) W31(w,i,a2,a3) 
  int sa3(s_eri);
  int ia3(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia3);
  orz::DTensor W31caa_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa3));
  FC_FUNC(g_if_sigma_ccov_cooo_no3_x0_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO3_X0_TYPE1_ERI_O)
    (sa3, ia3, T2b.cptr(), W31caa_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_cooo_no3_x1_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO3_X1_TYPE1_ERI_O)
      (sa, ia, sa3, ia3, V2_sym.cptr(), W31caa_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W32(x,a2) += (    1.00000000) T2(a1,x,a0,a2) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    4.00000000) V2(a2,a,w,i) W32(x,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia2);
  orz::DTensor W32c_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, sa2));
  FC_FUNC(g_if_sigma_ccov_cooo_no4_x0_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO4_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), W32c_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_cooo_no4_x1_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO4_X1_TYPE1_ERI_O)
      (sa, ia, sa2, ia2, V2_sym.cptr(), W32c_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W33(w,a2) += (    1.00000000) T2(a1,w,a0,a2) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) V2(a2,a,x,i) W33(w,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia2);
  orz::DTensor W33c_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, sa2));
  FC_FUNC(g_if_sigma_ccov_cooo_no5_x0_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO5_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), W33c_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_cooo_no5_x1_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO5_X1_TYPE1_ERI_O)
      (sa, ia, sa2, ia2, V2_sym.cptr(), W33c_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W36(w,a0,a2,a) += (    1.00000000) V2(a2,a,w,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) T2(a0,x,i,a2) W36(w,a0,a2,a) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(ia2);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W36ca_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sa2^sa));
    FC_FUNC(g_if_sigma_ccov_cooo_no6_x0_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO6_X0_TYPE1_ERI_O)
      (sa, ia, sa2, ia2, V2_sym.cptr(), W36ca_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_cooo_no6_x1_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO6_X1_TYPE1_ERI_O)
      (sa, ia, sa2, ia2, T2b.cptr(), W36ca_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W37(x,a0,a2,a) += (    1.00000000) V2(a2,a,x,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) T2(a0,w,i,a2) W37(x,a0,a2,a) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(ia2);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W37ca_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sa2^sa));
    FC_FUNC(g_if_sigma_ccov_cooo_no7_x0_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO7_X0_TYPE1_ERI_O)
      (sa, ia, sa2, ia2, V2_sym.cptr(), W37ca_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_cooo_no7_x1_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO7_X1_TYPE1_ERI_O)
      (sa, ia, sa2, ia2, T2b.cptr(), W37ca_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| W38(w,a2) += (    1.00000000) T2(a1,w,a0,a2) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    4.00000000) V2(a2,i,x,a) W38(w,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia2);
  orz::DTensor W38c_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, sa2));
  FC_FUNC(g_if_sigma_ccov_cooo_no8_x0_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO8_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), W38c_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_cooo_no8_x1_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO8_X1_TYPE1_ERI_O)
      (sa, ia, sa2, ia2, V2_sym.cptr(), W38c_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W39(x,a2) += (    1.00000000) T2(a1,x,a0,a2) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) V2(a2,i,w,a) W39(x,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia2);
  orz::DTensor W39c_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, sa2));
  FC_FUNC(g_if_sigma_ccov_cooo_no9_x0_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO9_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), W39c_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_cooo_no9_x1_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO9_X1_TYPE1_ERI_O)
      (sa, ia, sa2, ia2, V2_sym.cptr(), W39c_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| W40(w,i,a3,a2) += (    1.00000000) T2(w,a1,a0,a2) D2(i,a0,a1,a3) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) V2(a2,a3,x,a) W40(w,i,a3,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia2);
  orz::DTensor W40caa_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa2));
  FC_FUNC(g_if_sigma_ccov_cooo_no10_x0_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO10_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), W40caa_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_cooo_no10_x1_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO10_X1_TYPE1_ERI_O)
      (sa, ia, sa2, ia2, V2_sym.cptr(), W40caa_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| W41(x,i,a3,a2) += (    1.00000000) T2(x,a1,a0,a2) D2(i,a0,a1,a3) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) V2(a2,a3,w,a) W41(x,i,a3,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia2);
  orz::DTensor W41caa_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa2));
  FC_FUNC(g_if_sigma_ccov_cooo_no11_x0_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO11_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), W41caa_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_cooo_no11_x1_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO11_X1_TYPE1_ERI_O)
      (sa, ia, sa2, ia2, V2_sym.cptr(), W41caa_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| W42(w,i,a2,a3) += (    1.00000000) T2(w,a1,a0,a3) D2(i,a0,a1,a2) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) V2(a3,a,x,a2) W42(w,i,a2,a3) 
  int sa3(s_eri);
  int ia3(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia3);
  orz::DTensor W42caa_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa3));
  FC_FUNC(g_if_sigma_ccov_cooo_no12_x0_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO12_X0_TYPE1_ERI_O)
    (sa3, ia3, T2b.cptr(), W42caa_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_cooo_no12_x1_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO12_X1_TYPE1_ERI_O)
      (sa, ia, sa3, ia3, V2_sym.cptr(), W42caa_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| W43(x,i,a2,a3) += (    1.00000000) T2(x,a1,a0,a3) D2(i,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) V2(a3,a,w,a2) W43(x,i,a2,a3) 
  int sa3(s_eri);
  int ia3(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia3);
  orz::DTensor W43caa_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa3));
  FC_FUNC(g_if_sigma_ccov_cooo_no13_x0_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO13_X0_TYPE1_ERI_O)
    (sa3, ia3, T2b.cptr(), W43caa_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_cooo_no13_x1_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO13_X1_TYPE1_ERI_O)
      (sa, ia, sa3, ia3, V2_sym.cptr(), W43caa_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   14] -- 
  // |-- [    0] --| W44(w,a2) += (    1.00000000) T2(w,a1,a0,a2) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) V2(a2,a,x,i) W44(w,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia2);
  orz::DTensor W44c_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, sa2));
  FC_FUNC(g_if_sigma_ccov_cooo_no14_x0_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO14_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), W44c_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_cooo_no14_x1_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO14_X1_TYPE1_ERI_O)
      (sa, ia, sa2, ia2, V2_sym.cptr(), W44c_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   15] -- 
  // |-- [    0] --| W45(x,a2) += (    1.00000000) T2(x,a1,a0,a2) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) V2(a2,a,w,i) W45(x,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia2);
  orz::DTensor W45c_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, sa2));
  FC_FUNC(g_if_sigma_ccov_cooo_no15_x0_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO15_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), W45c_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_cooo_no15_x1_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO15_X1_TYPE1_ERI_O)
      (sa, ia, sa2, ia2, V2_sym.cptr(), W45c_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   16] -- 
  // |-- [    0] --| W48(w,a0,a2,a) += (    1.00000000) V2(a2,a,w,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) T2(x,a0,i,a2) W48(w,a0,a2,a) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(ia2);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W48ca_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sa2^sa));
    FC_FUNC(g_if_sigma_ccov_cooo_no16_x0_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO16_X0_TYPE1_ERI_O)
      (sa, ia, sa2, ia2, V2_sym.cptr(), W48ca_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_cooo_no16_x1_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO16_X1_TYPE1_ERI_O)
      (sa, ia, sa2, ia2, T2b.cptr(), W48ca_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   17] -- 
  // |-- [    0] --| W49(x,a0,a2,a) += (    1.00000000) V2(a2,a,x,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) T2(w,a0,i,a2) W49(x,a0,a2,a) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(ia2);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W49ca_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sa2^sa));
    FC_FUNC(g_if_sigma_ccov_cooo_no17_x0_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO17_X0_TYPE1_ERI_O)
      (sa, ia, sa2, ia2, V2_sym.cptr(), W49ca_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_cooo_no17_x1_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO17_X1_TYPE1_ERI_O)
      (sa, ia, sa2, ia2, T2b.cptr(), W49ca_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   18] -- 
  // |-- [    0] --| W50(w,a2) += (    1.00000000) T2(w,a1,a0,a2) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) V2(a2,i,x,a) W50(w,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia2);
  orz::DTensor W50c_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, sa2));
  FC_FUNC(g_if_sigma_ccov_cooo_no18_x0_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO18_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), W50c_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_cooo_no18_x1_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO18_X1_TYPE1_ERI_O)
      (sa, ia, sa2, ia2, V2_sym.cptr(), W50c_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   19] -- 
  // |-- [    0] --| W51(x,a2) += (    1.00000000) T2(x,a1,a0,a2) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) V2(a2,i,w,a) W51(x,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia2);
  orz::DTensor W51c_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, sa2));
  FC_FUNC(g_if_sigma_ccov_cooo_no19_x0_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO19_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), W51c_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_cooo_no19_x1_type1_eri_o,G_IF_SIGMA_CCOV_COOO_NO19_X1_TYPE1_ERI_O)
      (sa, ia, sa2, ia2, V2_sym.cptr(), W51c_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // --  Title : sigma_ccov_cooo
  //  >> Intermediates for the external contractions are defined here << 
  orz::DTensor W6ca_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W7ca_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W8caaa_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W9caaa_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W10caaa_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W11caaa_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W12ca_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W13ca_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W14ca_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W15ca_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W16ca_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W17ca_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W22ca_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W23ca_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
//-@type(2).declaration(end)

  //*-- Entering to take the type 0 contractions --*//
//-@type(0).contraction(begin)

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W6(c0,i) += (    1.00000000) T2(c0,a2,a1,a0) D2(i,a1,a2,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccov_cooo_no0_x0_type0_eri_v,G_IF_SIGMA_CCOV_COOO_NO0_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W6ca_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W7(c0,i) += (    1.00000000) T2(c0,a2,a1,a0) D2(i,a1,a2,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccov_cooo_no1_x0_type0_eri_v,G_IF_SIGMA_CCOV_COOO_NO1_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W7ca_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W8(w,i,a4,a3) += (    1.00000000) T2(w,a2,a1,a0) D3(i,a1,a4,a3,a2,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccov_cooo_no2_x0_type0_eri_v,G_IF_SIGMA_CCOV_COOO_NO2_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W8caaa_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W9(x,i,a4,a3) += (    1.00000000) T2(x,a2,a1,a0) D3(i,a1,a4,a3,a2,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccov_cooo_no3_x0_type0_eri_v,G_IF_SIGMA_CCOV_COOO_NO3_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W9caaa_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W10(w,i,a4,a3) += (    1.00000000) T2(w,a2,a1,a0) D3(i,a1,a4,a3,a2,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccov_cooo_no4_x0_type0_eri_v,G_IF_SIGMA_CCOV_COOO_NO4_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W10caaa_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W11(x,i,a4,a3) += (    1.00000000) T2(x,a2,a1,a0) D3(i,a3,a4,a1,a2,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccov_cooo_no5_x0_type0_eri_v,G_IF_SIGMA_CCOV_COOO_NO5_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W11caaa_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W12(w,a3) += (    1.00000000) T2(w,a2,a1,a0) D2(a3,a1,a2,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccov_cooo_no6_x0_type0_eri_v,G_IF_SIGMA_CCOV_COOO_NO6_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W12ca_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W13(x,a3) += (    1.00000000) T2(x,a2,a1,a0) D2(a3,a1,a2,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccov_cooo_no7_x0_type0_eri_v,G_IF_SIGMA_CCOV_COOO_NO7_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W13ca_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| W14(w,a3) += (    1.00000000) T2(w,a2,a1,a0) D2(a3,a1,a2,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccov_cooo_no8_x0_type0_eri_v,G_IF_SIGMA_CCOV_COOO_NO8_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W14ca_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W15(x,a3) += (    1.00000000) T2(x,a2,a1,a0) D2(a3,a1,a2,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccov_cooo_no9_x0_type0_eri_v,G_IF_SIGMA_CCOV_COOO_NO9_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W15ca_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| W16(c0,i) += (    1.00000000) T2(c0,a1,i,a0) D1(a1,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccov_cooo_no10_x0_type0_eri_v,G_IF_SIGMA_CCOV_COOO_NO10_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W16ca_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| W17(c0,i) += (    1.00000000) T2(c0,a1,i,a0) D1(a1,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccov_cooo_no11_x0_type0_eri_v,G_IF_SIGMA_CCOV_COOO_NO11_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W17ca_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| W22(c0,i) += (    1.00000000) T2(c0,a1,a0,i) D1(a1,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    FC_FUNC(g_if_sigma_ccov_cooo_no12_x0_type0_eri_v,G_IF_SIGMA_CCOV_COOO_NO12_X0_TYPE0_ERI_V)
      (si, ii, T2b.cptr(), W22ca_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| W23(c0,i) += (    1.00000000) T2(c0,a1,a0,i) D1(a1,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    FC_FUNC(g_if_sigma_ccov_cooo_no13_x0_type0_eri_v,G_IF_SIGMA_CCOV_COOO_NO13_X0_TYPE0_ERI_V)
      (si, ii, T2b.cptr(), W23ca_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
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
  // -- Title : sigma_ccov_cooo

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| S2(w,x,i,a) += (    2.00000000) V2(a,x,w,c0) W6(c0,i) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  FC_FUNC(g_if_sigma_ccov_cooo_no0_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO0_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W6ca_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| S2(w,x,i,a) += (   -1.00000000) V2(a,w,x,c0) W7(c0,i) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  FC_FUNC(g_if_sigma_ccov_cooo_no1_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO1_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W7ca_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| S2(w,x,i,a) += (   -2.00000000) V2(a,x,a4,a3) W8(w,i,a4,a3) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  FC_FUNC(g_if_sigma_ccov_cooo_no2_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO2_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W8caaa_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| S2(w,x,i,a) += (    1.00000000) V2(a,w,a4,a3) W9(x,i,a4,a3) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  FC_FUNC(g_if_sigma_ccov_cooo_no3_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO3_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W9caaa_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| S2(w,x,i,a) += (    1.00000000) V2(a,a4,x,a3) W10(w,i,a4,a3) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  FC_FUNC(g_if_sigma_ccov_cooo_no4_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO4_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W10caaa_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| S2(w,x,i,a) += (    1.00000000) V2(a,a4,w,a3) W11(x,i,a4,a3) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  FC_FUNC(g_if_sigma_ccov_cooo_no5_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO5_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W11caaa_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| S2(w,x,i,a) += (   -2.00000000) V2(a,x,i,a3) W12(w,a3) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  FC_FUNC(g_if_sigma_ccov_cooo_no6_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO6_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W12ca_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| S2(w,x,i,a) += (    1.00000000) V2(a,w,i,a3) W13(x,a3) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  FC_FUNC(g_if_sigma_ccov_cooo_no7_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO7_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W13ca_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| S2(w,x,i,a) += (    1.00000000) V2(a,a3,x,i) W14(w,a3) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  FC_FUNC(g_if_sigma_ccov_cooo_no8_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO8_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W14ca_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| S2(w,x,i,a) += (   -2.00000000) V2(a,a3,w,i) W15(x,a3) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  FC_FUNC(g_if_sigma_ccov_cooo_no9_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO9_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W15ca_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| S2(w,x,i,a) += (    2.00000000) V2(a,w,x,c0) W16(c0,i) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  FC_FUNC(g_if_sigma_ccov_cooo_no10_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO10_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W16ca_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| S2(w,x,i,a) += (   -4.00000000) V2(a,x,w,c0) W17(c0,i) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  FC_FUNC(g_if_sigma_ccov_cooo_no11_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO11_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W17ca_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| W18(x,a1,a0,a) += (    1.00000000) V2(a,x,a3,a2) D2(a3,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    4.00000000) T2(w,a1,i,a0) W18(x,a1,a0,a) 
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
    orz::DTensor W18ca_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_ccov_cooo_no12_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO12_X0_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, V2_sym.cptr(), W18ca_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_cooo_no12_x1_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO12_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), W18ca_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| W19(w,a1,a0,a) += (    1.00000000) V2(a,w,a3,a2) D2(a3,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) T2(x,a1,i,a0) W19(w,a1,a0,a) 
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
    orz::DTensor W19ca_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_ccov_cooo_no13_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO13_X0_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, V2_sym.cptr(), W19ca_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_cooo_no13_x1_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO13_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), W19ca_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   14] -- 
  // |-- [    0] --| W20(w,a1,a0,a) += (    1.00000000) V2(a,a3,w,a2) D2(a3,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) T2(x,a1,i,a0) W20(w,a1,a0,a) 
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
    orz::DTensor W20ca_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_ccov_cooo_no14_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO14_X0_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, V2_sym.cptr(), W20ca_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_cooo_no14_x1_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO14_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), W20ca_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   15] -- 
  // |-- [    0] --| W21(x,a1,a0,a) += (    1.00000000) V2(a,a3,x,a2) D2(a3,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) T2(w,a1,i,a0) W21(x,a1,a0,a) 
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
    orz::DTensor W21ca_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_ccov_cooo_no15_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO15_X0_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, V2_sym.cptr(), W21ca_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_cooo_no15_x1_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO15_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), W21ca_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   16] -- 
  // |-- [    0] --| S2(w,x,i,a) += (    2.00000000) V2(a,x,w,c0) W22(c0,i) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  FC_FUNC(g_if_sigma_ccov_cooo_no16_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO16_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W22ca_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   17] -- 
  // |-- [    0] --| S2(w,x,i,a) += (   -1.00000000) V2(a,w,x,c0) W23(c0,i) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  FC_FUNC(g_if_sigma_ccov_cooo_no17_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO17_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W23ca_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   18] -- 
  // |-- [    0] --| W24(x,a1,a0,a) += (    1.00000000) V2(a,x,a3,a2) D2(a3,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) T2(w,a1,a0,i) W24(x,a1,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W24caa_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ccov_cooo_no18_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO18_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W24caa_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    FC_FUNC(g_if_sigma_ccov_cooo_no18_x1_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO18_X1_TYPE1_ERI_V)
      (sa, ia, si, ii, T2b.cptr(), W24caa_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   19] -- 
  // |-- [    0] --| W25(w,a1,a0,a) += (    1.00000000) V2(a,w,a3,a2) D2(a3,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) T2(x,a1,a0,i) W25(w,a1,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W25caa_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ccov_cooo_no19_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO19_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W25caa_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    FC_FUNC(g_if_sigma_ccov_cooo_no19_x1_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO19_X1_TYPE1_ERI_V)
      (sa, ia, si, ii, T2b.cptr(), W25caa_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   20] -- 
  // |-- [    0] --| W26(x,a1,a0,a) += (    1.00000000) V2(a,a3,x,a2) D2(a3,a2,a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) T2(w,a1,a0,i) W26(x,a1,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W26caa_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ccov_cooo_no20_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO20_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W26caa_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    FC_FUNC(g_if_sigma_ccov_cooo_no20_x1_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO20_X1_TYPE1_ERI_V)
      (sa, ia, si, ii, T2b.cptr(), W26caa_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   21] -- 
  // |-- [    0] --| W27(w,a1,a0,a) += (    1.00000000) V2(a,a3,w,a2) D2(a3,a0,a1,a2) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) T2(x,a1,a0,i) W27(w,a1,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W27caa_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ccov_cooo_no21_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO21_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W27caa_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    FC_FUNC(g_if_sigma_ccov_cooo_no21_x1_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO21_X1_TYPE1_ERI_V)
      (sa, ia, si, ii, T2b.cptr(), W27caa_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   22] -- 
  // |-- [    0] --| W34(x,a0,a1,a) += (    1.00000000) V2(a,x,a2,a1) D1(a2,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) T2(w,a0,a1,i) W34(x,a0,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W34caa_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ccov_cooo_no22_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO22_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W34caa_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    FC_FUNC(g_if_sigma_ccov_cooo_no22_x1_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO22_X1_TYPE1_ERI_V)
      (sa, ia, si, ii, T2b.cptr(), W34caa_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   23] -- 
  // |-- [    0] --| W35(w,a0,a1,a) += (    1.00000000) V2(a,w,a2,a1) D1(a2,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) T2(x,a0,a1,i) W35(w,a0,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W35caa_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ccov_cooo_no23_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO23_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W35caa_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    FC_FUNC(g_if_sigma_ccov_cooo_no23_x1_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO23_X1_TYPE1_ERI_V)
      (sa, ia, si, ii, T2b.cptr(), W35caa_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   24] -- 
  // |-- [    0] --| W46(x,a0,a1,a) += (    1.00000000) V2(a,x,a2,a1) D1(a2,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    4.00000000) T2(w,a0,i,a1) W46(x,a0,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W46ca_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sa1^sa));
    FC_FUNC(g_if_sigma_ccov_cooo_no24_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO24_X0_TYPE1_ERI_V)
      (sa, ia, sa1, ia1, V2_sym.cptr(), W46ca_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_cooo_no24_x1_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO24_X1_TYPE1_ERI_V)
      (sa, ia, sa1, ia1, T2b.cptr(), W46ca_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   25] -- 
  // |-- [    0] --| W47(w,a0,a1,a) += (    1.00000000) V2(a,w,a2,a1) D1(a2,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) T2(x,a0,i,a1) W47(w,a0,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W47ca_sigma_ccov_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sa1^sa));
    FC_FUNC(g_if_sigma_ccov_cooo_no25_x0_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO25_X0_TYPE1_ERI_V)
      (sa, ia, sa1, ia1, V2_sym.cptr(), W47ca_sigma_ccov_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_cooo_no25_x1_type1_eri_v,G_IF_SIGMA_CCOV_COOO_NO25_X1_TYPE1_ERI_V)
      (sa, ia, sa1, ia1, T2b.cptr(), W47ca_sigma_ccov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
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
  orz::mr::Femto::my_timer.push_back(boost::make_tuple("sigma_ccov_cooo", time_sigma.elapsed_cputime(), time_sigma.elapsed_wallclocktime()));
  orz::mr::Femto::file_timing << "* " << boost::format("%20s : %10.7f %10.7f ") % orz::mr::Femto::my_timer.back().get<0>() % orz::mr::Femto::my_timer.back().get<1>() % orz::mr::Femto::my_timer.back().get<2>() << endl;
  flush(orz::mr::Femto::file_timing);

  return retval; 
} 
