                                                                                
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
#include <sci/icmr/Femto/elems/c_sigma_ccov_ccvv.h>                                  
                                                                                
using std::cout;                                                                
using std::endl;                                                                
                                                                                
#define FLOPCOUNT                                                               
                                                                                
using orz::mr::Femto::Fc0;                                                      
using orz::mr::Femto::h1_int;                                                   
using orz::mr::Femto::h6_int;                                                   
                                                                                
//  `7MM"""YMM                         mm               
//    MM    `7                         MM                  
//    MM   d  .gP"Ya `7MMpMMMb.pMMMb.mmMMmm ,pW"Wq.      
//    MM""MM ,M'   Yb  MM    MM    MM  MM  6W'   `Wb   
//    MM   Y 8M""""""  MM    MM    MM  MM  8M     M8 
//    MM     YM.    ,  MM    MM    MM  MM  YA.   ,A9       
//  .JMML.    `Mbmmd'.JMML  JMML  JMML.`Mbmo`Ybmd9'        

//                                   Generated date : Sun Apr 20 10:26:25 2014

                                                                                
// ***************************************************************************  
// orz::mr::mrci                                                                
// ***************************************************************************  
									    /*!        
   @brief CT input                                                              
									     */        
                                                                                
orz::mr::BareAmpPack orz::mr::Femto::sigma_ccov_ccvv(const orz::mr::Input &ctinp,                                    
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
  // |-- [    0] --| W0(i,v0) += (    1.00000000) D1(i,a0) Fc1(v0,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    2.00000000) T2(x,w,v0,a) W0(i,v0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W0av_sigma_ccov_ccvv(orz::mr::sizeof_sympack_Xav(symblockinfo, 0));
  FC_FUNC(g_if_sigma_ccov_ccvv_no0_x0_type0_noeri,G_IF_SIGMA_CCOV_CCVV_NO0_X0_TYPE0_NOERI)
    (W0av_sigma_ccov_ccvv.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccvv_no0_x1_type0_noeri,G_IF_SIGMA_CCOV_CCVV_NO0_X1_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W0av_sigma_ccov_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W1(i,v0) += (    1.00000000) D1(i,a0) Fc1(v0,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -4.00000000) T2(w,x,v0,a) W1(i,v0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W1av_sigma_ccov_ccvv(orz::mr::sizeof_sympack_Xav(symblockinfo, 0));
  FC_FUNC(g_if_sigma_ccov_ccvv_no1_x0_type0_noeri,G_IF_SIGMA_CCOV_CCVV_NO1_X0_TYPE0_NOERI)
    (W1av_sigma_ccov_ccvv.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccvv_no1_x1_type0_noeri,G_IF_SIGMA_CCOV_CCVV_NO1_X1_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W1av_sigma_ccov_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| S2(w,x,i,a) += (    8.00000000) T2(w,x,v0,a) Fc1(v0,i) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccvv_no2_x0_type0_noeri,G_IF_SIGMA_CCOV_CCVV_NO2_X0_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| S2(w,x,i,a) += (   -4.00000000) T2(x,w,v0,a) Fc1(v0,i) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccvv_no3_x0_type0_noeri,G_IF_SIGMA_CCOV_CCVV_NO3_X0_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- Title : sigma_ccov_ccvv

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W2(w,c0,i,v0) += (    1.00000000) V2(w,a0,v0,c0) D1(i,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -8.00000000) T2(c0,x,v0,a) W2(w,c0,i,v0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W2cav_sigma_ccov_ccvv(orz::mr::sizeof_sympack_Xcav(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccov_ccvv_no0_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCVV_NO0_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W2cav_sigma_ccov_ccvv.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccvv_no0_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCVV_NO0_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), W2cav_sigma_ccov_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W3(x,c0,i,v0) += (    1.00000000) V2(x,a0,v0,c0) D1(i,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    4.00000000) T2(c0,w,v0,a) W3(x,c0,i,v0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W3cav_sigma_ccov_ccvv(orz::mr::sizeof_sympack_Xcav(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccov_ccvv_no1_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCVV_NO1_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W3cav_sigma_ccov_ccvv.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccvv_no1_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCVV_NO1_X1_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), W3cav_sigma_ccov_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W4(w,c0,i,v0) += (    1.00000000) V2(w,a0,v0,c0) D1(i,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    4.00000000) T2(x,c0,v0,a) W4(w,c0,i,v0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W4cav_sigma_ccov_ccvv(orz::mr::sizeof_sympack_Xcav(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccov_ccvv_no2_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCVV_NO2_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W4cav_sigma_ccov_ccvv.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccvv_no2_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCVV_NO2_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), W4cav_sigma_ccov_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W5(x,c0,i,v0) += (    1.00000000) V2(x,a0,v0,c0) D1(i,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) T2(w,c0,v0,a) W5(x,c0,i,v0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W5cav_sigma_ccov_ccvv(orz::mr::sizeof_sympack_Xcav(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccov_ccvv_no3_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCVV_NO3_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W5cav_sigma_ccov_ccvv.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccvv_no3_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCVV_NO3_X1_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), W5cav_sigma_ccov_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W6(w,c0,i,v0) += (    1.00000000) V2(w,c0,v0,a0) D1(i,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    4.00000000) T2(c0,x,v0,a) W6(w,c0,i,v0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W6cav_sigma_ccov_ccvv(orz::mr::sizeof_sympack_Xcav(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccov_ccvv_no4_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCVV_NO4_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W6cav_sigma_ccov_ccvv.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccvv_no4_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCVV_NO4_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), W6cav_sigma_ccov_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W7(x,c0,i,v0) += (    1.00000000) V2(x,c0,v0,a0) D1(i,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) T2(c0,w,v0,a) W7(x,c0,i,v0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W7cav_sigma_ccov_ccvv(orz::mr::sizeof_sympack_Xcav(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccov_ccvv_no5_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCVV_NO5_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W7cav_sigma_ccov_ccvv.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccvv_no5_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCVV_NO5_X1_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), W7cav_sigma_ccov_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W8(w,c0,i,v0) += (    1.00000000) V2(w,c0,v0,a0) D1(i,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -2.00000000) T2(x,c0,v0,a) W8(w,c0,i,v0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W8cav_sigma_ccov_ccvv(orz::mr::sizeof_sympack_Xcav(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccov_ccvv_no6_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCVV_NO6_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W8cav_sigma_ccov_ccvv.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccvv_no6_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCVV_NO6_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), W8cav_sigma_ccov_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W9(x,c0,i,v0) += (    1.00000000) V2(x,c0,v0,a0) D1(i,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    4.00000000) T2(w,c0,v0,a) W9(x,c0,i,v0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W9cav_sigma_ccov_ccvv(orz::mr::sizeof_sympack_Xcav(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccov_ccvv_no7_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCVV_NO7_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W9cav_sigma_ccov_ccvv.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccvv_no7_x1_type1_eri_c,G_IF_SIGMA_CCOV_CCVV_NO7_X1_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), W9cav_sigma_ccov_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| S2(w,x,i,a) += (   16.00000000) V2(w,i,v0,c0) T2(c0,x,v0,a) 
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
    FC_FUNC(g_if_sigma_ccov_ccvv_no8_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCVV_NO8_X0_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| S2(w,x,i,a) += (   -8.00000000) V2(x,i,v0,c0) T2(c0,w,v0,a) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccvv_no9_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCVV_NO9_X0_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| S2(w,x,i,a) += (   -8.00000000) V2(w,i,v0,c0) T2(x,c0,v0,a) 
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
    FC_FUNC(g_if_sigma_ccov_ccvv_no10_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCVV_NO10_X0_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| S2(w,x,i,a) += (    4.00000000) V2(x,i,v0,c0) T2(w,c0,v0,a) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccvv_no11_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCVV_NO11_X0_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| S2(w,x,i,a) += (   -8.00000000) V2(x,c0,v0,i) T2(w,c0,v0,a) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccvv_no12_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCVV_NO12_X0_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| S2(w,x,i,a) += (    4.00000000) V2(w,c0,v0,i) T2(x,c0,v0,a) 
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
    FC_FUNC(g_if_sigma_ccov_ccvv_no13_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCVV_NO13_X0_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   14] -- 
  // |-- [    0] --| S2(w,x,i,a) += (    4.00000000) V2(x,c0,v0,i) T2(c0,w,v0,a) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccvv_no14_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCVV_NO14_X0_TYPE1_ERI_C)
      (sa, ia, sx, ix, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   15] -- 
  // |-- [    0] --| S2(w,x,i,a) += (   -8.00000000) V2(w,c0,v0,i) T2(c0,x,v0,a) 
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
    FC_FUNC(g_if_sigma_ccov_ccvv_no15_x0_type1_eri_c,G_IF_SIGMA_CCOV_CCVV_NO15_X0_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- Title : sigma_ccov_ccvv

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W10(i,v0) += (    1.00000000) V2(v0,a1,a2,a0) D2(i,a1,a2,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -4.00000000) T2(x,w,a,v0) W10(i,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W10a_sigma_ccov_ccvv(orz::mr::sizeof_sympack_Xa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ccov_ccvv_no0_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCVV_NO0_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W10a_sigma_ccov_ccvv.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccvv_no0_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCVV_NO0_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), W10a_sigma_ccov_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W11(i,v0) += (    1.00000000) V2(v0,a1,a2,a0) D2(i,a1,a2,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    2.00000000) T2(x,w,v0,a) W11(i,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W11a_sigma_ccov_ccvv(orz::mr::sizeof_sympack_Xa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ccov_ccvv_no1_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCVV_NO1_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W11a_sigma_ccov_ccvv.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccvv_no1_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCVV_NO1_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), W11a_sigma_ccov_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W12(i,v0) += (    1.00000000) V2(v0,a1,i,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -4.00000000) T2(x,w,a,v0) W12(i,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W12a_sigma_ccov_ccvv(orz::mr::sizeof_sympack_Xa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ccov_ccvv_no2_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCVV_NO2_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W12a_sigma_ccov_ccvv.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccvv_no2_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCVV_NO2_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), W12a_sigma_ccov_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W13(i,v0) += (    1.00000000) V2(v0,a1,i,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    2.00000000) T2(x,w,v0,a) W13(i,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W13a_sigma_ccov_ccvv(orz::mr::sizeof_sympack_Xa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ccov_ccvv_no3_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCVV_NO3_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W13a_sigma_ccov_ccvv.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccvv_no3_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCVV_NO3_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), W13a_sigma_ccov_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W14(i,v0) += (    1.00000000) V2(v0,i,a1,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    8.00000000) T2(x,w,a,v0) W14(i,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W14a_sigma_ccov_ccvv(orz::mr::sizeof_sympack_Xa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ccov_ccvv_no4_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCVV_NO4_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W14a_sigma_ccov_ccvv.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccvv_no4_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCVV_NO4_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), W14a_sigma_ccov_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W15(i,v0) += (    1.00000000) V2(v0,i,a1,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -4.00000000) T2(x,w,v0,a) W15(i,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W15a_sigma_ccov_ccvv(orz::mr::sizeof_sympack_Xa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ccov_ccvv_no5_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCVV_NO5_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W15a_sigma_ccov_ccvv.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccvv_no5_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCVV_NO5_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), W15a_sigma_ccov_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W16(i,v1,v0,a) += (    1.00000000) V2(v1,a0,v0,a) D1(i,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (   -4.00000000) T2(x,w,v0,v1) W16(i,v1,v0,a) 
  int sv1(s_eri);
  int iv1(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv1);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W16av_sigma_ccov_ccvv(orz::mr::sizeof_sympack_Xav(symblockinfo, sv1^sa));
    FC_FUNC(g_if_sigma_ccov_ccvv_no6_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCVV_NO6_X0_TYPE1_ERI_V)
      (sa, ia, sv1, iv1, V2_sym.cptr(), W16av_sigma_ccov_ccvv.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccvv_no6_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCVV_NO6_X1_TYPE1_ERI_V)
      (sa, ia, sv1, iv1, T2b.cptr(), W16av_sigma_ccov_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W17(i,v1,v0,a) += (    1.00000000) V2(v1,a0,v0,a) D1(i,a0) 
  // |-- [    1] --| S2(w,x,i,a) += (    2.00000000) T2(w,x,v0,v1) W17(i,v1,v0,a) 
  int sv1(s_eri);
  int iv1(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv1);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W17av_sigma_ccov_ccvv(orz::mr::sizeof_sympack_Xav(symblockinfo, sv1^sa));
    FC_FUNC(g_if_sigma_ccov_ccvv_no7_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCVV_NO7_X0_TYPE1_ERI_V)
      (sa, ia, sv1, iv1, V2_sym.cptr(), W17av_sigma_ccov_ccvv.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccov_ccvv_no7_x1_type1_eri_v,G_IF_SIGMA_CCOV_CCVV_NO7_X1_TYPE1_ERI_V)
      (sa, ia, sv1, iv1, T2b.cptr(), W17av_sigma_ccov_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| S2(w,x,i,a) += (    8.00000000) V2(v1,a,v0,i) T2(w,x,v0,v1) 
  int sv1(s_eri);
  int iv1(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv1);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccvv_no8_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCVV_NO8_X0_TYPE1_ERI_V)
      (sa, ia, sv1, iv1, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| S2(w,x,i,a) += (   -4.00000000) V2(v1,i,v0,a) T2(w,x,v0,v1) 
  int sv1(s_eri);
  int iv1(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv1);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccov_ccvv_no9_x0_type1_eri_v,G_IF_SIGMA_CCOV_CCVV_NO9_X0_TYPE1_ERI_V)
      (sa, ia, sv1, iv1, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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

//-@D4C.contractions(end)

  // Do timing!
  orz::mr::Femto::my_timer.push_back(boost::make_tuple("sigma_ccov_ccvv", time_sigma.elapsed_cputime(), time_sigma.elapsed_wallclocktime()));
  orz::mr::Femto::file_timing << "* " << boost::format("%20s : %10.7f %10.7f ") % orz::mr::Femto::my_timer.back().get<0>() % orz::mr::Femto::my_timer.back().get<1>() % orz::mr::Femto::my_timer.back().get<2>() << endl;
  flush(orz::mr::Femto::file_timing);

  return retval; 
} 
