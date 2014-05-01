                                                                                
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
#include <sci/icmr/Femto/elems/c_sigma_ccov_covv.h>                                  
                                                                                
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
                                                                                
orz::mr::BareAmpPack orz::mr::Femto::sigma_ccov_covv(const orz::mr::Input &ctinp,                                    
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
  // -- Title : sigma_ccov_covv

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W0(x,i,a0,v0) += (    1.00000000) V2(x,a2,v0,a1) D2(i,a1,a0,a2) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) T2(w,a0,v0,a) W0(x,i,a0,v0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W0aav_sigma_ccov_covv(orz::mr::sizeof_sympack_Xaav(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccov_covv_no0_x0_type1_eri_c,G_IF_SIGMA_CCOV_COVV_NO0_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W0aav_sigma_ccov_covv.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_covv_no0_x1_type1_eri_c,G_IF_SIGMA_CCOV_COVV_NO0_X1_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), W0aav_sigma_ccov_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W1(w,i,a0,v0) += (    1.00000000) V2(w,a2,v0,a1) D2(i,a2,a0,a1) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) T2(x,a0,v0,a) W1(w,i,a0,v0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W1aav_sigma_ccov_covv(orz::mr::sizeof_sympack_Xaav(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccov_covv_no1_x0_type1_eri_c,G_IF_SIGMA_CCOV_COVV_NO1_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W1aav_sigma_ccov_covv.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_covv_no1_x1_type1_eri_c,G_IF_SIGMA_CCOV_COVV_NO1_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), W1aav_sigma_ccov_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W2(x,a0,i,v0) += (    1.00000000) V2(x,i,v0,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) T2(w,a0,v0,a) W2(x,a0,i,v0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W2aav_sigma_ccov_covv(orz::mr::sizeof_sympack_Xaav(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccov_covv_no2_x0_type1_eri_c,G_IF_SIGMA_CCOV_COVV_NO2_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W2aav_sigma_ccov_covv.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_covv_no2_x1_type1_eri_c,G_IF_SIGMA_CCOV_COVV_NO2_X1_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), W2aav_sigma_ccov_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W3(w,a0,i,v0) += (    1.00000000) V2(w,i,v0,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) T2(x,a0,v0,a) W3(w,a0,i,v0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W3aav_sigma_ccov_covv(orz::mr::sizeof_sympack_Xaav(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccov_covv_no3_x0_type1_eri_c,G_IF_SIGMA_CCOV_COVV_NO3_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W3aav_sigma_ccov_covv.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_covv_no3_x1_type1_eri_c,G_IF_SIGMA_CCOV_COVV_NO3_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), W3aav_sigma_ccov_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W4(w,a0,i,v0) += (    1.00000000) V2(w,a1,v0,i) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) T2(x,a0,v0,a) W4(w,a0,i,v0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W4aav_sigma_ccov_covv(orz::mr::sizeof_sympack_Xaav(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccov_covv_no4_x0_type1_eri_c,G_IF_SIGMA_CCOV_COVV_NO4_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W4aav_sigma_ccov_covv.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_covv_no4_x1_type1_eri_c,G_IF_SIGMA_CCOV_COVV_NO4_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), W4aav_sigma_ccov_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W5(x,a0,i,v0) += (    1.00000000) V2(x,a1,v0,i) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) T2(w,a0,v0,a) W5(x,a0,i,v0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W5aav_sigma_ccov_covv(orz::mr::sizeof_sympack_Xaav(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccov_covv_no5_x0_type1_eri_c,G_IF_SIGMA_CCOV_COVV_NO5_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W5aav_sigma_ccov_covv.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_covv_no5_x1_type1_eri_c,G_IF_SIGMA_CCOV_COVV_NO5_X1_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), W5aav_sigma_ccov_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W6(w,i,a0,v0) += (    1.00000000) V2(w,a2,v0,a1) D2(i,a2,a0,a1) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) T2(a0,x,v0,a) W6(w,i,a0,v0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W6aav_sigma_ccov_covv(orz::mr::sizeof_sympack_Xaav(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccov_covv_no6_x0_type1_eri_c,G_IF_SIGMA_CCOV_COVV_NO6_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W6aav_sigma_ccov_covv.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_covv_no6_x1_type1_eri_c,G_IF_SIGMA_CCOV_COVV_NO6_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), W6aav_sigma_ccov_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W7(x,i,a0,v0) += (    1.00000000) V2(x,a2,v0,a1) D2(i,a2,a0,a1) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) T2(a0,w,v0,a) W7(x,i,a0,v0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W7aav_sigma_ccov_covv(orz::mr::sizeof_sympack_Xaav(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccov_covv_no7_x0_type1_eri_c,G_IF_SIGMA_CCOV_COVV_NO7_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W7aav_sigma_ccov_covv.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_covv_no7_x1_type1_eri_c,G_IF_SIGMA_CCOV_COVV_NO7_X1_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), W7aav_sigma_ccov_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| W8(w,a0,i,v0) += (    1.00000000) V2(w,i,v0,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    4.00000000) T2(a0,x,v0,a) W8(w,a0,i,v0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W8aav_sigma_ccov_covv(orz::mr::sizeof_sympack_Xaav(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccov_covv_no8_x0_type1_eri_c,G_IF_SIGMA_CCOV_COVV_NO8_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W8aav_sigma_ccov_covv.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_covv_no8_x1_type1_eri_c,G_IF_SIGMA_CCOV_COVV_NO8_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), W8aav_sigma_ccov_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W9(x,a0,i,v0) += (    1.00000000) V2(x,i,v0,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) T2(a0,w,v0,a) W9(x,a0,i,v0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W9aav_sigma_ccov_covv(orz::mr::sizeof_sympack_Xaav(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccov_covv_no9_x0_type1_eri_c,G_IF_SIGMA_CCOV_COVV_NO9_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W9aav_sigma_ccov_covv.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_covv_no9_x1_type1_eri_c,G_IF_SIGMA_CCOV_COVV_NO9_X1_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), W9aav_sigma_ccov_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| W10(w,a0,i,v0) += (    1.00000000) V2(w,a1,v0,i) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) T2(a0,x,v0,a) W10(w,a0,i,v0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W10aav_sigma_ccov_covv(orz::mr::sizeof_sympack_Xaav(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccov_covv_no10_x0_type1_eri_c,G_IF_SIGMA_CCOV_COVV_NO10_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W10aav_sigma_ccov_covv.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_covv_no10_x1_type1_eri_c,G_IF_SIGMA_CCOV_COVV_NO10_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), W10aav_sigma_ccov_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| W11(x,a0,i,v0) += (    1.00000000) V2(x,a1,v0,i) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    1.00000000) T2(a0,w,v0,a) W11(x,a0,i,v0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W11aav_sigma_ccov_covv(orz::mr::sizeof_sympack_Xaav(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccov_covv_no11_x0_type1_eri_c,G_IF_SIGMA_CCOV_COVV_NO11_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W11aav_sigma_ccov_covv.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_covv_no11_x1_type1_eri_c,G_IF_SIGMA_CCOV_COVV_NO11_X1_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), W11aav_sigma_ccov_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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

//-@ERI.contractions(end)

//-@D4C.contractions(begin)

//-@D4C.contractions(end)

  // Do timing!
  orz::mr::Femto::my_timer.push_back(boost::make_tuple("sigma_ccov_covv", time_sigma.elapsed_cputime(), time_sigma.elapsed_wallclocktime()));
  orz::mr::Femto::file_timing << "* " << boost::format("%20s : %10.7f %10.7f ") % orz::mr::Femto::my_timer.back().get<0>() % orz::mr::Femto::my_timer.back().get<1>() % orz::mr::Femto::my_timer.back().get<2>() << endl;
  flush(orz::mr::Femto::file_timing);

  return retval; 
} 
