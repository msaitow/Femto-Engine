                                                                                
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
#include <sci/icmr/Femto/elems/c_sigma_ccov_ccoo.h>                                  
                                                                                
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
                                                                                
orz::mr::BareAmpPack orz::mr::Femto::sigma_ccov_ccoo(const orz::mr::Input &ctinp,                                    
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
  // |-- [    0] --| W0(x,w,i,a2) += (    1.00000000) T2(x,w,a0,a1) D2(i,a1,a2,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    2.00000000) Fc1(a2,a) W0(x,w,i,a2) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W0ccaa_sigma_ccov_ccoo(orz::mr::sizeof_sympack_Xccaa(symblockinfo, 0));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_ccov_ccoo_no0_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOO_NO0_X0_TYPE0_NOERI)
      (sa1, ia1, T2b.cptr(), W0ccaa_sigma_ccov_ccoo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccoo_no0_x1_type0_noeri,G_IF_SIGMA_CCOV_CCOO_NO0_X1_TYPE0_NOERI)
      (sa, ia, W0ccaa_sigma_ccov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W1(a0,a) += (    1.00000000) D1(a1,a0) Fc1(a1,a) 
  // |-- [    1] --| S2(w,x,i,a) += (   -4.00000000) T2(x,w,a0,i) W1(a0,a) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    orz::DTensor W1a_sigma_ccov_ccoo(orz::mr::sizeof_sympack_Xa(symblockinfo, sa));
    FC_FUNC(g_if_sigma_ccov_ccoo_no1_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOO_NO1_X0_TYPE0_NOERI)
      (sa, ia, W1a_sigma_ccov_ccoo.cptr(), nir, nsym, psym, &flops);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    for(int si = 0;si < nir;++si){ 
    for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
      T2b = T2.get_amp2(ii);
      FC_FUNC(g_if_sigma_ccov_ccoo_no1_x1_type0_noeri,G_IF_SIGMA_CCOV_CCOO_NO1_X1_TYPE0_NOERI)
        (sa, ia, si, ii, T2b.cptr(), W1a_sigma_ccov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
    } // End ii
    } // End si
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W2(a0,a) += (    1.00000000) D1(a1,a0) Fc1(a1,a) 
  // |-- [    1] --| S2(w,x,i,a) += (    2.00000000) T2(x,w,i,a0) W2(a0,a) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      S2b = orz::DTensor(retval.namps_iamp()[ia]);
      double W2_sigma_ccov_ccoo(0);
      FC_FUNC(g_if_sigma_ccov_ccoo_no2_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOO_NO2_X0_TYPE0_NOERI)
        (sa, ia, sa0, ia0, &W2_sigma_ccov_ccoo, nir, nsym, psym, &flops);
      FC_FUNC(g_if_sigma_ccov_ccoo_no2_x1_type0_noeri,G_IF_SIGMA_CCOV_CCOO_NO2_X1_TYPE0_NOERI)
        (sa, ia, sa0, ia0, T2b.cptr(), &W2_sigma_ccov_ccoo, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
      retval.acc_amp2(ia, S2b);
    } // End ia
    } // End sa
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W3(x,w,i,a1) += (    1.00000000) T2(x,w,a0,a1) D1(i,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    2.00000000) Fc1(a1,a) W3(x,w,i,a1) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W3cca_sigma_ccov_ccoo(orz::mr::sizeof_sympack_Xcca(symblockinfo, sa1));
    FC_FUNC(g_if_sigma_ccov_ccoo_no3_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOO_NO3_X0_TYPE0_NOERI)
      (sa1, ia1, T2b.cptr(), W3cca_sigma_ccov_ccoo.cptr(), nir, nsym, psym, &flops);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      S2b = orz::DTensor(retval.namps_iamp()[ia]);
      FC_FUNC(g_if_sigma_ccov_ccoo_no3_x1_type0_noeri,G_IF_SIGMA_CCOV_CCOO_NO3_X1_TYPE0_NOERI)
        (sa, ia, sa1, ia1, W3cca_sigma_ccov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- No. [    4] -- 
  // |-- [    0] --| W4(x,w,i,a1) += (    1.00000000) T2(x,w,a1,a0) D1(i,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -4.00000000) Fc1(a1,a) W4(x,w,i,a1) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W4ccaa_sigma_ccov_ccoo(orz::mr::sizeof_sympack_Xccaa(symblockinfo, 0));
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccov_ccoo_no4_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOO_NO4_X0_TYPE0_NOERI)
      (sa0, ia0, T2b.cptr(), W4ccaa_sigma_ccov_ccoo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccoo_no4_x1_type0_noeri,G_IF_SIGMA_CCOV_CCOO_NO4_X1_TYPE0_NOERI)
      (sa, ia, W4ccaa_sigma_ccov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| S2(w,x,i,a) += (    8.00000000) T2(x,w,a0,i) Fc1(a0,a) 
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
      FC_FUNC(g_if_sigma_ccov_ccoo_no5_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOO_NO5_X0_TYPE0_NOERI)
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
  // -- No. [    6] -- 
  // |-- [    0] --| S2(w,x,i,a) += (   -4.00000000) T2(x,w,i,a0) Fc1(a0,a) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      S2b = orz::DTensor(retval.namps_iamp()[ia]);
      FC_FUNC(g_if_sigma_ccov_ccoo_no6_x0_type0_noeri,G_IF_SIGMA_CCOV_CCOO_NO6_X0_TYPE0_NOERI)
        (sa, ia, sa0, ia0, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
      retval.acc_amp2(ia, S2b);
    } // End ia
    } // End sa
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
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
  // -- Title : sigma_ccov_ccoo

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W25(i,a0,a1,a) += (    1.00000000) V2(a1,a3,a2,a) D2(i,a0,a2,a3) 
  // |-- [    1] --| S2(w,x,i,a) += (    2.00000000) T2(w,x,a0,a1) W25(i,a0,a1,a) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(ia1);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W25aa_sigma_ccov_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa1^sa));
    FC_FUNC(g_if_sigma_ccov_ccoo_no0_x0_type1_eri_o,G_IF_SIGMA_CCOV_CCOO_NO0_X0_TYPE1_ERI_O)
      (sa, ia, sa1, ia1, V2_sym.cptr(), W25aa_sigma_ccov_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccoo_no0_x1_type1_eri_o,G_IF_SIGMA_CCOV_CCOO_NO0_X1_TYPE1_ERI_O)
      (sa, ia, sa1, ia1, T2b.cptr(), W25aa_sigma_ccov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W26(i,a0,a1,a) += (    1.00000000) V2(a1,a3,a2,a) D2(i,a3,a2,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    2.00000000) T2(x,w,a0,a1) W26(i,a0,a1,a) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(ia1);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W26aa_sigma_ccov_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa1^sa));
    FC_FUNC(g_if_sigma_ccov_ccoo_no1_x0_type1_eri_o,G_IF_SIGMA_CCOV_CCOO_NO1_X0_TYPE1_ERI_O)
      (sa, ia, sa1, ia1, V2_sym.cptr(), W26aa_sigma_ccov_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccoo_no1_x1_type1_eri_o,G_IF_SIGMA_CCOV_CCOO_NO1_X1_TYPE1_ERI_O)
      (sa, ia, sa1, ia1, T2b.cptr(), W26aa_sigma_ccov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W31(i,a0,a2,a) += (    1.00000000) V2(a2,a,a3,a1) D2(i,a0,a3,a1) 
  // |-- [    1] --| S2(w,x,i,a) += (   -4.00000000) T2(w,x,a0,a2) W31(i,a0,a2,a) 
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
    orz::DTensor W31aa_sigma_ccov_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa2^sa));
    FC_FUNC(g_if_sigma_ccov_ccoo_no2_x0_type1_eri_o,G_IF_SIGMA_CCOV_CCOO_NO2_X0_TYPE1_ERI_O)
      (sa, ia, sa2, ia2, V2_sym.cptr(), W31aa_sigma_ccov_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccoo_no2_x1_type1_eri_o,G_IF_SIGMA_CCOV_CCOO_NO2_X1_TYPE1_ERI_O)
      (sa, ia, sa2, ia2, T2b.cptr(), W31aa_sigma_ccov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W32(i,a0,a2,a) += (    1.00000000) V2(a2,a,a3,a1) D2(i,a0,a3,a1) 
  // |-- [    1] --| S2(w,x,i,a) += (    2.00000000) T2(x,w,a0,a2) W32(i,a0,a2,a) 
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
    orz::DTensor W32aa_sigma_ccov_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa2^sa));
    FC_FUNC(g_if_sigma_ccov_ccoo_no3_x0_type1_eri_o,G_IF_SIGMA_CCOV_CCOO_NO3_X0_TYPE1_ERI_O)
      (sa, ia, sa2, ia2, V2_sym.cptr(), W32aa_sigma_ccov_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccoo_no3_x1_type1_eri_o,G_IF_SIGMA_CCOV_CCOO_NO3_X1_TYPE1_ERI_O)
      (sa, ia, sa2, ia2, T2b.cptr(), W32aa_sigma_ccov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W33(a0,i,a2,a) += (    1.00000000) V2(a2,a,i,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -4.00000000) T2(w,x,a0,a2) W33(a0,i,a2,a) 
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
    orz::DTensor W33aa_sigma_ccov_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa2^sa));
    FC_FUNC(g_if_sigma_ccov_ccoo_no4_x0_type1_eri_o,G_IF_SIGMA_CCOV_CCOO_NO4_X0_TYPE1_ERI_O)
      (sa, ia, sa2, ia2, V2_sym.cptr(), W33aa_sigma_ccov_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccoo_no4_x1_type1_eri_o,G_IF_SIGMA_CCOV_CCOO_NO4_X1_TYPE1_ERI_O)
      (sa, ia, sa2, ia2, T2b.cptr(), W33aa_sigma_ccov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W34(a0,i,a2,a) += (    1.00000000) V2(a2,a,i,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    2.00000000) T2(x,w,a0,a2) W34(a0,i,a2,a) 
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
    orz::DTensor W34aa_sigma_ccov_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa2^sa));
    FC_FUNC(g_if_sigma_ccov_ccoo_no5_x0_type1_eri_o,G_IF_SIGMA_CCOV_CCOO_NO5_X0_TYPE1_ERI_O)
      (sa, ia, sa2, ia2, V2_sym.cptr(), W34aa_sigma_ccov_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccoo_no5_x1_type1_eri_o,G_IF_SIGMA_CCOV_CCOO_NO5_X1_TYPE1_ERI_O)
      (sa, ia, sa2, ia2, T2b.cptr(), W34aa_sigma_ccov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W35(a0,a) += (    1.00000000) V2(a0,a2,a1,a) D1(a2,a1) 
  // |-- [    1] --| S2(w,x,i,a) += (    2.00000000) T2(x,w,i,a0) W35(a0,a) 
  int sa0(s_eri);
  int ia0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(ia0);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    double W35_sigma_ccov_ccoo(0);
    FC_FUNC(g_if_sigma_ccov_ccoo_no6_x0_type1_eri_o,G_IF_SIGMA_CCOV_CCOO_NO6_X0_TYPE1_ERI_O)
      (sa, ia, sa0, ia0, V2_sym.cptr(), &W35_sigma_ccov_ccoo, nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccoo_no6_x1_type1_eri_o,G_IF_SIGMA_CCOV_CCOO_NO6_X1_TYPE1_ERI_O)
      (sa, ia, sa0, ia0, T2b.cptr(), &W35_sigma_ccov_ccoo, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W36(a0,a) += (    1.00000000) V2(a0,a2,a1,a) D1(a2,a1) 
  // |-- [    1] --| S2(w,x,i,a) += (   -4.00000000) T2(w,x,i,a0) W36(a0,a) 
  int sa0(s_eri);
  int ia0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(ia0);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    double W36_sigma_ccov_ccoo(0);
    FC_FUNC(g_if_sigma_ccov_ccoo_no7_x0_type1_eri_o,G_IF_SIGMA_CCOV_CCOO_NO7_X0_TYPE1_ERI_O)
      (sa, ia, sa0, ia0, V2_sym.cptr(), &W36_sigma_ccov_ccoo, nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccoo_no7_x1_type1_eri_o,G_IF_SIGMA_CCOV_CCOO_NO7_X1_TYPE1_ERI_O)
      (sa, ia, sa0, ia0, T2b.cptr(), &W36_sigma_ccov_ccoo, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| S2(w,x,i,a) += (   -8.00000000) V2(a0,a,x,c0) T2(w,c0,i,a0) 
  int sa0(s_eri);
  int ia0(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(ia0);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccoo_no8_x0_type1_eri_o,G_IF_SIGMA_CCOV_CCOO_NO8_X0_TYPE1_ERI_O)
      (sa, ia, sa0, ia0, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| S2(w,x,i,a) += (    4.00000000) V2(a0,a,w,c0) T2(x,c0,i,a0) 
  int sa0(s_eri);
  int ia0(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(ia0);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccoo_no9_x0_type1_eri_o,G_IF_SIGMA_CCOV_CCOO_NO9_X0_TYPE1_ERI_O)
      (sa, ia, sa0, ia0, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| S2(w,x,i,a) += (    4.00000000) V2(a0,a,x,c0) T2(c0,w,i,a0) 
  int sa0(s_eri);
  int ia0(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(ia0);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccoo_no10_x0_type1_eri_o,G_IF_SIGMA_CCOV_CCOO_NO10_X0_TYPE1_ERI_O)
      (sa, ia, sa0, ia0, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| S2(w,x,i,a) += (   -8.00000000) V2(a0,a,w,c0) T2(c0,x,i,a0) 
  int sa0(s_eri);
  int ia0(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(ia0);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccoo_no11_x0_type1_eri_o,G_IF_SIGMA_CCOV_CCOO_NO11_X0_TYPE1_ERI_O)
      (sa, ia, sa0, ia0, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| W37(a1,a) += (    1.00000000) V2(a1,a,a2,a0) D1(a2,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    8.00000000) T2(w,x,i,a1) W37(a1,a) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(ia1);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    double W37_sigma_ccov_ccoo(0);
    FC_FUNC(g_if_sigma_ccov_ccoo_no12_x0_type1_eri_o,G_IF_SIGMA_CCOV_CCOO_NO12_X0_TYPE1_ERI_O)
      (sa, ia, sa1, ia1, V2_sym.cptr(), &W37_sigma_ccov_ccoo, nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccoo_no12_x1_type1_eri_o,G_IF_SIGMA_CCOV_CCOO_NO12_X1_TYPE1_ERI_O)
      (sa, ia, sa1, ia1, T2b.cptr(), &W37_sigma_ccov_ccoo, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| W38(a1,a) += (    1.00000000) V2(a1,a,a2,a0) D1(a2,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -4.00000000) T2(x,w,i,a1) W38(a1,a) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(ia1);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    double W38_sigma_ccov_ccoo(0);
    FC_FUNC(g_if_sigma_ccov_ccoo_no13_x0_type1_eri_o,G_IF_SIGMA_CCOV_CCOO_NO13_X0_TYPE1_ERI_O)
      (sa, ia, sa1, ia1, V2_sym.cptr(), &W38_sigma_ccov_ccoo, nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccoo_no13_x1_type1_eri_o,G_IF_SIGMA_CCOV_CCOO_NO13_X1_TYPE1_ERI_O)
      (sa, ia, sa1, ia1, T2b.cptr(), &W38_sigma_ccov_ccoo, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   14] -- 
  // |-- [    0] --| W39(a0,i,a1,a) += (    1.00000000) V2(a1,i,a2,a) D1(a2,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    2.00000000) T2(w,x,a0,a1) W39(a0,i,a1,a) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(ia1);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W39aa_sigma_ccov_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa1^sa));
    FC_FUNC(g_if_sigma_ccov_ccoo_no14_x0_type1_eri_o,G_IF_SIGMA_CCOV_CCOO_NO14_X0_TYPE1_ERI_O)
      (sa, ia, sa1, ia1, V2_sym.cptr(), W39aa_sigma_ccov_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccoo_no14_x1_type1_eri_o,G_IF_SIGMA_CCOV_CCOO_NO14_X1_TYPE1_ERI_O)
      (sa, ia, sa1, ia1, T2b.cptr(), W39aa_sigma_ccov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   15] -- 
  // |-- [    0] --| W40(a0,i,a1,a) += (    1.00000000) V2(a1,i,a2,a) D1(a2,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -4.00000000) T2(x,w,a0,a1) W40(a0,i,a1,a) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(ia1);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W40aa_sigma_ccov_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa1^sa));
    FC_FUNC(g_if_sigma_ccov_ccoo_no15_x0_type1_eri_o,G_IF_SIGMA_CCOV_CCOO_NO15_X0_TYPE1_ERI_O)
      (sa, ia, sa1, ia1, V2_sym.cptr(), W40aa_sigma_ccov_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccoo_no15_x1_type1_eri_o,G_IF_SIGMA_CCOV_CCOO_NO15_X1_TYPE1_ERI_O)
      (sa, ia, sa1, ia1, T2b.cptr(), W40aa_sigma_ccov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   16] -- 
  // |-- [    0] --| W41(i,a1,a0,a) += (    1.00000000) V2(a0,a2,a1,a) D1(i,a2) 
  // |-- [    1] --| S2(w,x,i,a) += (   -4.00000000) T2(x,w,a1,a0) W41(i,a1,a0,a) 
  int sa0(s_eri);
  int ia0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(ia0);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W41aa_sigma_ccov_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_ccov_ccoo_no16_x0_type1_eri_o,G_IF_SIGMA_CCOV_CCOO_NO16_X0_TYPE1_ERI_O)
      (sa, ia, sa0, ia0, V2_sym.cptr(), W41aa_sigma_ccov_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccoo_no16_x1_type1_eri_o,G_IF_SIGMA_CCOV_CCOO_NO16_X1_TYPE1_ERI_O)
      (sa, ia, sa0, ia0, T2b.cptr(), W41aa_sigma_ccov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   17] -- 
  // |-- [    0] --| W42(i,a1,a0,a) += (    1.00000000) V2(a0,a2,a1,a) D1(i,a2) 
  // |-- [    1] --| S2(w,x,i,a) += (    2.00000000) T2(w,x,a1,a0) W42(i,a1,a0,a) 
  int sa0(s_eri);
  int ia0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(ia0);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W42aa_sigma_ccov_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_ccov_ccoo_no17_x0_type1_eri_o,G_IF_SIGMA_CCOV_CCOO_NO17_X0_TYPE1_ERI_O)
      (sa, ia, sa0, ia0, V2_sym.cptr(), W42aa_sigma_ccov_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccoo_no17_x1_type1_eri_o,G_IF_SIGMA_CCOV_CCOO_NO17_X1_TYPE1_ERI_O)
      (sa, ia, sa0, ia0, T2b.cptr(), W42aa_sigma_ccov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   18] -- 
  // |-- [    0] --| S2(w,x,i,a) += (    8.00000000) V2(a0,i,a1,a) T2(x,w,a1,a0) 
  int sa0(s_eri);
  int ia0(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(ia0);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccoo_no18_x0_type1_eri_o,G_IF_SIGMA_CCOV_CCOO_NO18_X0_TYPE1_ERI_O)
      (sa, ia, sa0, ia0, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   19] -- 
  // |-- [    0] --| S2(w,x,i,a) += (   -4.00000000) V2(a1,i,a0,a) T2(w,x,a0,a1) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(ia1);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccoo_no19_x0_type1_eri_o,G_IF_SIGMA_CCOV_CCOO_NO19_X0_TYPE1_ERI_O)
      (sa, ia, sa1, ia1, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // --  Title : sigma_ccov_ccoo
  //  >> Intermediates for the external contractions are defined here << 
  orz::DTensor W5ccaa_sigma_ccov_ccoo(orz::mr::sizeof_sympack_Xccaa(symblockinfo, 0));
  orz::DTensor W6ccaa_sigma_ccov_ccoo(orz::mr::sizeof_sympack_Xccaa(symblockinfo, 0));
  orz::DTensor W7ccaa_sigma_ccov_ccoo(orz::mr::sizeof_sympack_Xccaa(symblockinfo, 0));
  orz::DTensor W8ccaa_sigma_ccov_ccoo(orz::mr::sizeof_sympack_Xccaa(symblockinfo, 0));
//-@type(2).declaration(end)

  //*-- Entering to take the type 0 contractions --*//
//-@type(0).contraction(begin)

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W5(w,c0,i,a2) += (    1.00000000) T2(w,c0,a1,a0) D2(i,a1,a2,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccov_ccoo_no0_x0_type0_eri_v,G_IF_SIGMA_CCOV_CCOO_NO0_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W5ccaa_sigma_ccov_ccoo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W6(x,c0,i,a2) += (    1.00000000) T2(x,c0,a1,a0) D2(i,a1,a2,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccov_ccoo_no1_x0_type0_eri_v,G_IF_SIGMA_CCOV_CCOO_NO1_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W6ccaa_sigma_ccov_ccoo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W7(w,c0,i,a2) += (    1.00000000) T2(w,c0,a1,a0) D2(i,a1,a2,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccov_ccoo_no2_x0_type0_eri_v,G_IF_SIGMA_CCOV_CCOO_NO2_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W7ccaa_sigma_ccov_ccoo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W8(x,c0,i,a2) += (    1.00000000) T2(x,c0,a0,a1) D2(i,a1,a2,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_ccov_ccoo_no3_x0_type0_eri_v,G_IF_SIGMA_CCOV_CCOO_NO3_X0_TYPE0_ERI_V)
      (sa1, ia1, T2b.cptr(), W8ccaa_sigma_ccov_ccoo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
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
  // -- Title : sigma_ccov_ccoo

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| S2(w,x,i,a) += (    4.00000000) V2(a,x,c0,a2) W5(w,c0,i,a2) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  FC_FUNC(g_if_sigma_ccov_ccoo_no0_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO0_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W5ccaa_sigma_ccov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| S2(w,x,i,a) += (   -2.00000000) V2(a,w,c0,a2) W6(x,c0,i,a2) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  FC_FUNC(g_if_sigma_ccov_ccoo_no1_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO1_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W6ccaa_sigma_ccov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| S2(w,x,i,a) += (   -2.00000000) V2(a,a2,x,c0) W7(w,c0,i,a2) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  FC_FUNC(g_if_sigma_ccov_ccoo_no2_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO2_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W7ccaa_sigma_ccov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| S2(w,x,i,a) += (   -2.00000000) V2(a,a2,w,c0) W8(x,c0,i,a2) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  FC_FUNC(g_if_sigma_ccov_ccoo_no3_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO3_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W8ccaa_sigma_ccov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W9(i,a1,a0,a) += (    1.00000000) V2(a,a3,a4,a2) D3(i,a1,a4,a2,a3,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    2.00000000) T2(x,w,a0,a1) W9(i,a1,a0,a) 
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
    orz::DTensor W9aa_sigma_ccov_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa1^sa));
    FC_FUNC(g_if_sigma_ccov_ccoo_no4_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO4_X0_TYPE1_ERI_V)
      (sa, ia, sa1, ia1, V2_sym.cptr(), W9aa_sigma_ccov_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccoo_no4_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO4_X1_TYPE1_ERI_V)
      (sa, ia, sa1, ia1, T2b.cptr(), W9aa_sigma_ccov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W10(a1,a0,i,a) += (    1.00000000) V2(a,a3,i,a2) D2(a3,a1,a2,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    2.00000000) T2(x,w,a1,a0) W10(a1,a0,i,a) 
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
    orz::DTensor W10aa_sigma_ccov_ccoo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_ccov_ccoo_no5_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO5_X0_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, V2_sym.cptr(), W10aa_sigma_ccov_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccoo_no5_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO5_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), W10aa_sigma_ccov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W11(x,c0,a0,a) += (    1.00000000) V2(a,x,c0,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -8.00000000) T2(w,c0,i,a0) W11(x,c0,a0,a) 
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
    orz::DTensor W11cc_sigma_ccov_ccoo(orz::mr::sizeof_sympack_Xcc(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_ccov_ccoo_no6_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO6_X0_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, V2_sym.cptr(), W11cc_sigma_ccov_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccoo_no6_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO6_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), W11cc_sigma_ccov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W12(w,c0,a0,a) += (    1.00000000) V2(a,w,c0,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    4.00000000) T2(x,c0,i,a0) W12(w,c0,a0,a) 
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
    orz::DTensor W12cc_sigma_ccov_ccoo(orz::mr::sizeof_sympack_Xcc(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_ccov_ccoo_no7_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO7_X0_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, V2_sym.cptr(), W12cc_sigma_ccov_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccoo_no7_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO7_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), W12cc_sigma_ccov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| W13(w,c0,a0,a) += (    1.00000000) V2(a,w,c0,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) T2(x,c0,a0,i) W13(w,c0,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W13cca_sigma_ccov_ccoo(orz::mr::sizeof_sympack_Xcca(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ccov_ccoo_no8_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO8_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W13cca_sigma_ccov_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    FC_FUNC(g_if_sigma_ccov_ccoo_no8_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO8_X1_TYPE1_ERI_V)
      (sa, ia, si, ii, T2b.cptr(), W13cca_sigma_ccov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W14(x,c0,a0,a) += (    1.00000000) V2(a,x,c0,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    4.00000000) T2(w,c0,a0,i) W14(x,c0,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W14cca_sigma_ccov_ccoo(orz::mr::sizeof_sympack_Xcca(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ccov_ccoo_no9_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO9_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W14cca_sigma_ccov_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    FC_FUNC(g_if_sigma_ccov_ccoo_no9_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO9_X1_TYPE1_ERI_V)
      (sa, ia, si, ii, T2b.cptr(), W14cca_sigma_ccov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| W15(w,c0,a0,a) += (    1.00000000) V2(a,a1,w,c0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) T2(x,c0,i,a0) W15(w,c0,a0,a) 
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
    orz::DTensor W15cc_sigma_ccov_ccoo(orz::mr::sizeof_sympack_Xcc(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_ccov_ccoo_no10_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO10_X0_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, V2_sym.cptr(), W15cc_sigma_ccov_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccoo_no10_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO10_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), W15cc_sigma_ccov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| W16(x,c0,a0,a) += (    1.00000000) V2(a,a1,x,c0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    4.00000000) T2(w,c0,i,a0) W16(x,c0,a0,a) 
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
    orz::DTensor W16cc_sigma_ccov_ccoo(orz::mr::sizeof_sympack_Xcc(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_ccov_ccoo_no11_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO11_X0_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, V2_sym.cptr(), W16cc_sigma_ccov_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccoo_no11_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO11_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), W16cc_sigma_ccov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| W17(w,c0,a0,a) += (    1.00000000) V2(a,a1,w,c0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    4.00000000) T2(x,c0,a0,i) W17(w,c0,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W17cca_sigma_ccov_ccoo(orz::mr::sizeof_sympack_Xcca(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ccov_ccoo_no12_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO12_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W17cca_sigma_ccov_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    FC_FUNC(g_if_sigma_ccov_ccoo_no12_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO12_X1_TYPE1_ERI_V)
      (sa, ia, si, ii, T2b.cptr(), W17cca_sigma_ccov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| W18(x,c0,a0,a) += (    1.00000000) V2(a,a1,x,c0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) T2(w,c0,a0,i) W18(x,c0,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W18cca_sigma_ccov_ccoo(orz::mr::sizeof_sympack_Xcca(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ccov_ccoo_no13_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO13_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W18cca_sigma_ccov_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    FC_FUNC(g_if_sigma_ccov_ccoo_no13_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO13_X1_TYPE1_ERI_V)
      (sa, ia, si, ii, T2b.cptr(), W18cca_sigma_ccov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   14] -- 
  // |-- [    0] --| W19(a0,a) += (    1.00000000) V2(a,a2,a3,a1) D2(a3,a1,a2,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    2.00000000) T2(x,w,i,a0) W19(a0,a) 
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
    double W19_sigma_ccov_ccoo(0);
    FC_FUNC(g_if_sigma_ccov_ccoo_no14_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO14_X0_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, V2_sym.cptr(), &W19_sigma_ccov_ccoo, nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccoo_no14_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO14_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), &W19_sigma_ccov_ccoo, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   15] -- 
  // |-- [    0] --| W20(a0,a) += (    1.00000000) V2(a,a2,a3,a1) D2(a3,a1,a2,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -4.00000000) T2(x,w,a0,i) W20(a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W20a_sigma_ccov_ccoo(orz::mr::sizeof_sympack_Xa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ccov_ccoo_no15_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO15_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W20a_sigma_ccov_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    FC_FUNC(g_if_sigma_ccov_ccoo_no15_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO15_X1_TYPE1_ERI_V)
      (sa, ia, si, ii, T2b.cptr(), W20a_sigma_ccov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   16] -- 
  // |-- [    0] --| W21(w,x,a0,a) += (    1.00000000) V2(a,x,c0,a1) T2(w,c0,a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    4.00000000) D1(i,a0) W21(w,x,a0,a) 
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
    orz::DTensor W21cc_sigma_ccov_ccoo(orz::mr::sizeof_sympack_Xcc(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_ccov_ccoo_no16_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO16_X0_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), V2_sym.cptr(), W21cc_sigma_ccov_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccoo_no16_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO16_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, W21cc_sigma_ccov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   17] -- 
  // |-- [    0] --| W22(x,w,a0,a) += (    1.00000000) V2(a,w,c0,a1) T2(x,c0,a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) D1(i,a0) W22(x,w,a0,a) 
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
    orz::DTensor W22cc_sigma_ccov_ccoo(orz::mr::sizeof_sympack_Xcc(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_ccov_ccoo_no17_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO17_X0_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), V2_sym.cptr(), W22cc_sigma_ccov_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccoo_no17_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO17_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, W22cc_sigma_ccov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   18] -- 
  // |-- [    0] --| W23(x,w,a0,a) += (    1.00000000) V2(a,w,c0,a1) T2(x,c0,a0,a1) 
  // |-- [    1] --| S2(w,x,i,a) += (    4.00000000) D1(i,a0) W23(x,w,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W23cca_sigma_ccov_ccoo(orz::mr::sizeof_sympack_Xcca(symblockinfo, sa));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_ccov_ccoo_no18_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO18_X0_TYPE1_ERI_V)
      (sa, ia, sa1, ia1, T2b.cptr(), V2_sym.cptr(), W23cca_sigma_ccov_ccoo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  FC_FUNC(g_if_sigma_ccov_ccoo_no18_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO18_X1_TYPE1_ERI_V)
    (sa, ia, W23cca_sigma_ccov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   19] -- 
  // |-- [    0] --| W24(w,x,a0,a) += (    1.00000000) V2(a,x,c0,a1) T2(w,c0,a0,a1) 
  // |-- [    1] --| S2(w,x,i,a) += (   -8.00000000) D1(i,a0) W24(w,x,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W24cca_sigma_ccov_ccoo(orz::mr::sizeof_sympack_Xcca(symblockinfo, sa));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_ccov_ccoo_no19_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO19_X0_TYPE1_ERI_V)
      (sa, ia, sa1, ia1, T2b.cptr(), V2_sym.cptr(), W24cca_sigma_ccov_ccoo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  FC_FUNC(g_if_sigma_ccov_ccoo_no19_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO19_X1_TYPE1_ERI_V)
    (sa, ia, W24cca_sigma_ccov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   20] -- 
  // |-- [    0] --| W27(x,w,a0,a) += (    1.00000000) V2(a,a1,w,c0) T2(x,c0,a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    4.00000000) D1(i,a0) W27(x,w,a0,a) 
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
    orz::DTensor W27cc_sigma_ccov_ccoo(orz::mr::sizeof_sympack_Xcc(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_ccov_ccoo_no20_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO20_X0_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), V2_sym.cptr(), W27cc_sigma_ccov_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccoo_no20_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO20_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, W27cc_sigma_ccov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   21] -- 
  // |-- [    0] --| W28(w,x,a0,a) += (    1.00000000) V2(a,a1,x,c0) T2(w,c0,a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) D1(i,a0) W28(w,x,a0,a) 
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
    orz::DTensor W28cc_sigma_ccov_ccoo(orz::mr::sizeof_sympack_Xcc(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_ccov_ccoo_no21_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO21_X0_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), V2_sym.cptr(), W28cc_sigma_ccov_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccoo_no21_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO21_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, W28cc_sigma_ccov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   22] -- 
  // |-- [    0] --| W29(x,w,a0,a) += (    1.00000000) V2(a,a1,w,c0) T2(x,c0,a0,a1) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) D1(i,a0) W29(x,w,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W29cca_sigma_ccov_ccoo(orz::mr::sizeof_sympack_Xcca(symblockinfo, sa));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_ccov_ccoo_no22_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO22_X0_TYPE1_ERI_V)
      (sa, ia, sa1, ia1, T2b.cptr(), V2_sym.cptr(), W29cca_sigma_ccov_ccoo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  FC_FUNC(g_if_sigma_ccov_ccoo_no22_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO22_X1_TYPE1_ERI_V)
    (sa, ia, W29cca_sigma_ccov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   23] -- 
  // |-- [    0] --| W30(w,x,a0,a) += (    1.00000000) V2(a,a1,x,c0) T2(w,c0,a0,a1) 
  // |-- [    1] --| S2(w,x,i,a) += (    4.00000000) D1(i,a0) W30(w,x,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W30cca_sigma_ccov_ccoo(orz::mr::sizeof_sympack_Xcca(symblockinfo, sa));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_ccov_ccoo_no23_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO23_X0_TYPE1_ERI_V)
      (sa, ia, sa1, ia1, T2b.cptr(), V2_sym.cptr(), W30cca_sigma_ccov_ccoo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  FC_FUNC(g_if_sigma_ccov_ccoo_no23_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO23_X1_TYPE1_ERI_V)
    (sa, ia, W30cca_sigma_ccov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   24] -- 
  // |-- [    0] --| S2(w,x,i,a) += (   16.00000000) V2(a,x,c0,a0) T2(w,c0,i,a0) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccov_ccoo_no24_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO24_X0_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   25] -- 
  // |-- [    0] --| S2(w,x,i,a) += (   -8.00000000) V2(a,w,c0,a0) T2(x,c0,i,a0) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccov_ccoo_no25_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO25_X0_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   26] -- 
  // |-- [    0] --| S2(w,x,i,a) += (    4.00000000) V2(a,w,c0,a0) T2(x,c0,a0,i) 
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
    FC_FUNC(g_if_sigma_ccov_ccoo_no26_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO26_X0_TYPE1_ERI_V)
      (sa, ia, si, ii, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   27] -- 
  // |-- [    0] --| S2(w,x,i,a) += (   -8.00000000) V2(a,x,c0,a0) T2(w,c0,a0,i) 
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
    FC_FUNC(g_if_sigma_ccov_ccoo_no27_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCOO_NO27_X0_TYPE1_ERI_V)
      (sa, ia, si, ii, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
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
  orz::mr::Femto::my_timer.push_back(boost::make_tuple("sigma_ccov_ccoo", time_sigma.elapsed_cputime(), time_sigma.elapsed_wallclocktime()));
  orz::mr::Femto::file_timing << "* " << boost::format("%20s : %10.7f %10.7f ") % orz::mr::Femto::my_timer.back().get<0>() % orz::mr::Femto::my_timer.back().get<1>() % orz::mr::Femto::my_timer.back().get<2>() << endl;
  flush(orz::mr::Femto::file_timing);

  return retval; 
} 
