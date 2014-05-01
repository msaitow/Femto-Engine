                                                                                
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
#include <sci/icmr/Femto/elems/c_sigma_ccvv_ccov.h>                                  
                                                                                
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
                                                                                
orz::mr::BareAmpPack orz::mr::Femto::sigma_ccvv_ccov(const orz::mr::Input &ctinp,                                    
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
  // |-- [    0] --| W0(a0,a) += (    1.00000000) D1(a1,a0) Fc1(a1,a) 
  // |-- [    1] --| S2(w,x,a,b) += (    1.00000000) T2(x,w,a0,b) W0(a0,a) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W0av_sigma_ccvv_ccov(orz::mr::sizeof_sympack_Xav(symblockinfo, 0));
  FC_FUNC(g_if_sigma_ccvv_ccov_no0_x0_type0_noeri,G_IF_SIGMA_CCVV_CCOV_NO0_X0_TYPE0_NOERI)
    (W0av_sigma_ccvv_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccov_no0_x1_type0_noeri,G_IF_SIGMA_CCVV_CCOV_NO0_X1_TYPE0_NOERI)
      (sb, ib, T2b.cptr(), W0av_sigma_ccvv_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W1(a0,a) += (    1.00000000) D1(a1,a0) Fc1(a1,a) 
  // |-- [    1] --| S2(w,x,a,b) += (   -2.00000000) T2(w,x,a0,b) W1(a0,a) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W1av_sigma_ccvv_ccov(orz::mr::sizeof_sympack_Xav(symblockinfo, 0));
  FC_FUNC(g_if_sigma_ccvv_ccov_no1_x0_type0_noeri,G_IF_SIGMA_CCVV_CCOV_NO1_X0_TYPE0_NOERI)
    (W1av_sigma_ccvv_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccov_no1_x1_type0_noeri,G_IF_SIGMA_CCVV_CCOV_NO1_X1_TYPE0_NOERI)
      (sb, ib, T2b.cptr(), W1av_sigma_ccvv_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W2(a0,b) += (    1.00000000) D1(a1,a0) Fc1(b,a1) 
  // |-- [    1] --| S2(w,x,a,b) += (   -2.00000000) T2(x,w,a0,a) W2(a0,b) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    orz::DTensor W2a_sigma_ccvv_ccov(orz::mr::sizeof_sympack_Xa(symblockinfo, sb));
    FC_FUNC(g_if_sigma_ccvv_ccov_no2_x0_type0_noeri,G_IF_SIGMA_CCVV_CCOV_NO2_X0_TYPE0_NOERI)
      (sb, ib, W2a_sigma_ccvv_ccov.cptr(), nir, nsym, psym, &flops);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      T2b = T2.get_amp2(ia);
      FC_FUNC(g_if_sigma_ccvv_ccov_no2_x1_type0_noeri,G_IF_SIGMA_CCVV_CCOV_NO2_X1_TYPE0_NOERI)
        (sa, ia, sb, ib, T2b.cptr(), W2a_sigma_ccvv_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    } // End ia
    } // End sa
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W3(a0,b) += (    1.00000000) D1(a1,a0) Fc1(b,a1) 
  // |-- [    1] --| S2(w,x,a,b) += (    1.00000000) T2(x,w,a,a0) W3(a0,b) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    for(int sb = 0;sb < nir;++sb){ 
    for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
      S2b = orz::DTensor(retval.namps_iamp()[ib]);
      double W3_sigma_ccvv_ccov(0);
      FC_FUNC(g_if_sigma_ccvv_ccov_no3_x0_type0_noeri,G_IF_SIGMA_CCVV_CCOV_NO3_X0_TYPE0_NOERI)
        (sa0, ia0, sb, ib, &W3_sigma_ccvv_ccov, nir, nsym, psym, &flops);
      FC_FUNC(g_if_sigma_ccvv_ccov_no3_x1_type0_noeri,G_IF_SIGMA_CCVV_CCOV_NO3_X1_TYPE0_NOERI)
        (sa0, ia0, sb, ib, T2b.cptr(), &W3_sigma_ccvv_ccov, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
      retval.acc_amp2(ib, S2b);
    } // End ib
    } // End sb
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| S2(w,x,a,b) += (    4.00000000) T2(x,w,a0,a) Fc1(b,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    for(int sb = 0;sb < nir;++sb){ 
    for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
      S2b = orz::DTensor(retval.namps_iamp()[ib]);
      FC_FUNC(g_if_sigma_ccvv_ccov_no4_x0_type0_noeri,G_IF_SIGMA_CCVV_CCOV_NO4_X0_TYPE0_NOERI)
        (sa, ia, sb, ib, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
      retval.acc_amp2(ib, S2b);
    } // End ib
    } // End sb
  // --> @[a, "virtual"] [notNeeded]
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| S2(w,x,a,b) += (   -2.00000000) T2(x,w,a,a0) Fc1(b,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    for(int sb = 0;sb < nir;++sb){ 
    for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
      S2b = orz::DTensor(retval.namps_iamp()[ib]);
      FC_FUNC(g_if_sigma_ccvv_ccov_no5_x0_type0_noeri,G_IF_SIGMA_CCVV_CCOV_NO5_X0_TYPE0_NOERI)
        (sa0, ia0, sb, ib, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
      retval.acc_amp2(ib, S2b);
    } // End ib
    } // End sb
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| S2(w,x,a,b) += (    4.00000000) T2(w,x,a0,b) Fc1(a0,a) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccov_no6_x0_type0_noeri,G_IF_SIGMA_CCVV_CCOV_NO6_X0_TYPE0_NOERI)
      (sb, ib, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| S2(w,x,a,b) += (   -2.00000000) T2(x,w,a0,b) Fc1(a0,a) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccov_no7_x0_type0_noeri,G_IF_SIGMA_CCVV_CCOV_NO7_X0_TYPE0_NOERI)
      (sb, ib, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
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
  // -- Title : sigma_ccvv_ccov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W4(w,c0,a0,a) += (    1.00000000) V2(w,a,c0,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (   -4.00000000) T2(c0,x,a0,b) W4(w,c0,a0,a) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W4cav_sigma_ccvv_ccov(orz::mr::sizeof_sympack_Xcav(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccvv_ccov_no0_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO0_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W4cav_sigma_ccvv_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccov_no0_x1_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO0_X1_TYPE1_ERI_C)
      (sb, ib, sw, iw, T2b.cptr(), W4cav_sigma_ccvv_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W5(x,c0,a0,a) += (    1.00000000) V2(x,a,c0,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (    2.00000000) T2(c0,w,a0,b) W5(x,c0,a0,a) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W5cav_sigma_ccvv_ccov(orz::mr::sizeof_sympack_Xcav(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccvv_ccov_no1_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO1_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W5cav_sigma_ccvv_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccov_no1_x1_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO1_X1_TYPE1_ERI_C)
      (sb, ib, sx, ix, T2b.cptr(), W5cav_sigma_ccvv_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W6(w,c0,a0,a) += (    1.00000000) V2(w,a,c0,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (    2.00000000) T2(x,c0,a0,b) W6(w,c0,a0,a) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W6cav_sigma_ccvv_ccov(orz::mr::sizeof_sympack_Xcav(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccvv_ccov_no2_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO2_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W6cav_sigma_ccvv_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccov_no2_x1_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO2_X1_TYPE1_ERI_C)
      (sb, ib, sw, iw, T2b.cptr(), W6cav_sigma_ccvv_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W7(x,c0,a0,a) += (    1.00000000) V2(x,a,c0,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (   -1.00000000) T2(w,c0,a0,b) W7(x,c0,a0,a) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W7cav_sigma_ccvv_ccov(orz::mr::sizeof_sympack_Xcav(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccvv_ccov_no3_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO3_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W7cav_sigma_ccvv_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccov_no3_x1_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO3_X1_TYPE1_ERI_C)
      (sb, ib, sx, ix, T2b.cptr(), W7cav_sigma_ccvv_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W8(w,c0,a0,a) += (    1.00000000) V2(w,c0,a1,a) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (    2.00000000) T2(c0,x,a0,b) W8(w,c0,a0,a) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W8cav_sigma_ccvv_ccov(orz::mr::sizeof_sympack_Xcav(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccvv_ccov_no4_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO4_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W8cav_sigma_ccvv_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccov_no4_x1_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO4_X1_TYPE1_ERI_C)
      (sb, ib, sw, iw, T2b.cptr(), W8cav_sigma_ccvv_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W9(x,c0,a0,a) += (    1.00000000) V2(x,c0,a1,a) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (   -1.00000000) T2(c0,w,a0,b) W9(x,c0,a0,a) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W9cav_sigma_ccvv_ccov(orz::mr::sizeof_sympack_Xcav(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccvv_ccov_no5_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO5_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W9cav_sigma_ccvv_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccov_no5_x1_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO5_X1_TYPE1_ERI_C)
      (sb, ib, sx, ix, T2b.cptr(), W9cav_sigma_ccvv_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W10(w,c0,a0,a) += (    1.00000000) V2(w,c0,a1,a) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (   -1.00000000) T2(x,c0,a0,b) W10(w,c0,a0,a) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W10cav_sigma_ccvv_ccov(orz::mr::sizeof_sympack_Xcav(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccvv_ccov_no6_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO6_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W10cav_sigma_ccvv_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccov_no6_x1_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO6_X1_TYPE1_ERI_C)
      (sb, ib, sw, iw, T2b.cptr(), W10cav_sigma_ccvv_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W11(x,c0,a0,a) += (    1.00000000) V2(x,c0,a1,a) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (    2.00000000) T2(w,c0,a0,b) W11(x,c0,a0,a) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W11cav_sigma_ccvv_ccov(orz::mr::sizeof_sympack_Xcav(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccvv_ccov_no7_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO7_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W11cav_sigma_ccvv_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccov_no7_x1_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO7_X1_TYPE1_ERI_C)
      (sb, ib, sx, ix, T2b.cptr(), W11cav_sigma_ccvv_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| S2(w,x,a,b) += (    8.00000000) V2(w,a,c0,a0) T2(c0,x,a0,b) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccov_no8_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO8_X0_TYPE1_ERI_C)
      (sb, ib, sw, iw, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| S2(w,x,a,b) += (   -4.00000000) V2(x,a,c0,a0) T2(c0,w,a0,b) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccov_no9_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO9_X0_TYPE1_ERI_C)
      (sb, ib, sx, ix, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| S2(w,x,a,b) += (   -4.00000000) V2(w,a,c0,a0) T2(x,c0,a0,b) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccov_no10_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO10_X0_TYPE1_ERI_C)
      (sb, ib, sw, iw, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| S2(w,x,a,b) += (    2.00000000) V2(x,a,c0,a0) T2(w,c0,a0,b) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccov_no11_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO11_X0_TYPE1_ERI_C)
      (sb, ib, sx, ix, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| S2(w,x,a,b) += (   -4.00000000) V2(x,c0,a0,a) T2(w,c0,a0,b) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccov_no12_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO12_X0_TYPE1_ERI_C)
      (sb, ib, sx, ix, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| S2(w,x,a,b) += (    2.00000000) V2(w,c0,a0,a) T2(x,c0,a0,b) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccov_no13_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO13_X0_TYPE1_ERI_C)
      (sb, ib, sw, iw, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   14] -- 
  // |-- [    0] --| S2(w,x,a,b) += (    2.00000000) V2(x,c0,a0,a) T2(c0,w,a0,b) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccov_no14_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO14_X0_TYPE1_ERI_C)
      (sb, ib, sx, ix, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   15] -- 
  // |-- [    0] --| S2(w,x,a,b) += (   -4.00000000) V2(w,c0,a0,a) T2(c0,x,a0,b) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccov_no15_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCOV_NO15_X0_TYPE1_ERI_C)
      (sb, ib, sw, iw, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
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
  // -- Title : sigma_ccvv_ccov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W24(a0,a) += (    1.00000000) V2(a0,a2,a1,a) D1(a2,a1) 
  // |-- [    1] --| S2(w,x,a,b) += (   -2.00000000) T2(x,w,b,a0) W24(a0,a) 
  int sa0(s_eri);
  int ia0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia0);
  orz::DTensor W24v_sigma_ccvv_ccov(orz::mr::sizeof_sympack_Xv(symblockinfo, sa0));
  FC_FUNC(g_if_sigma_ccvv_ccov_no0_x0_type1_eri_o,G_IF_SIGMA_CCVV_CCOV_NO0_X0_TYPE1_ERI_O)
    (sa0, ia0, V2_sym.cptr(), W24v_sigma_ccvv_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccov_no0_x1_type1_eri_o,G_IF_SIGMA_CCVV_CCOV_NO0_X1_TYPE1_ERI_O)
      (sa0, ia0, sb, ib, T2b.cptr(), W24v_sigma_ccvv_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W25(a0,a) += (    1.00000000) V2(a0,a2,a1,a) D1(a2,a1) 
  // |-- [    1] --| S2(w,x,a,b) += (    1.00000000) T2(x,w,a0,b) W25(a0,a) 
  int sa0(s_eri);
  int ia0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W25v_sigma_ccvv_ccov(orz::mr::sizeof_sympack_Xv(symblockinfo, sa0));
  FC_FUNC(g_if_sigma_ccvv_ccov_no1_x0_type1_eri_o,G_IF_SIGMA_CCVV_CCOV_NO1_X0_TYPE1_ERI_O)
    (sa0, ia0, V2_sym.cptr(), W25v_sigma_ccvv_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccov_no1_x1_type1_eri_o,G_IF_SIGMA_CCVV_CCOV_NO1_X1_TYPE1_ERI_O)
      (sa0, ia0, sb, ib, T2b.cptr(), W25v_sigma_ccvv_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W30(a1,a) += (    1.00000000) V2(a1,a,a2,a0) D1(a2,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (    4.00000000) T2(x,w,b,a1) W30(a1,a) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia1);
  orz::DTensor W30v_sigma_ccvv_ccov(orz::mr::sizeof_sympack_Xv(symblockinfo, sa1));
  FC_FUNC(g_if_sigma_ccvv_ccov_no2_x0_type1_eri_o,G_IF_SIGMA_CCVV_CCOV_NO2_X0_TYPE1_ERI_O)
    (sa1, ia1, V2_sym.cptr(), W30v_sigma_ccvv_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccov_no2_x1_type1_eri_o,G_IF_SIGMA_CCVV_CCOV_NO2_X1_TYPE1_ERI_O)
      (sa1, ia1, sb, ib, T2b.cptr(), W30v_sigma_ccvv_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W31(a1,a) += (    1.00000000) V2(a1,a,a2,a0) D1(a2,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (   -2.00000000) T2(x,w,a1,b) W31(a1,a) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W31v_sigma_ccvv_ccov(orz::mr::sizeof_sympack_Xv(symblockinfo, sa1));
  FC_FUNC(g_if_sigma_ccvv_ccov_no3_x0_type1_eri_o,G_IF_SIGMA_CCVV_CCOV_NO3_X0_TYPE1_ERI_O)
    (sa1, ia1, V2_sym.cptr(), W31v_sigma_ccvv_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccov_no3_x1_type1_eri_o,G_IF_SIGMA_CCVV_CCOV_NO3_X1_TYPE1_ERI_O)
      (sa1, ia1, sb, ib, T2b.cptr(), W31v_sigma_ccvv_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
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
  // -- Title : sigma_ccvv_ccov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W12(a0,a) += (    1.00000000) V2(a,a2,a3,a1) D2(a3,a1,a2,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (   -2.00000000) T2(w,x,a0,b) W12(a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W12a_sigma_ccvv_ccov(orz::mr::sizeof_sympack_Xa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ccvv_ccov_no0_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO0_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W12a_sigma_ccvv_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccov_no0_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO0_X1_TYPE1_ERI_V)
      (sa, ia, sb, ib, T2b.cptr(), W12a_sigma_ccvv_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W13(a0,a) += (    1.00000000) V2(a,a2,a3,a1) D2(a3,a1,a2,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (    1.00000000) T2(x,w,a0,b) W13(a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W13a_sigma_ccvv_ccov(orz::mr::sizeof_sympack_Xa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ccvv_ccov_no1_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO1_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W13a_sigma_ccvv_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccov_no1_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO1_X1_TYPE1_ERI_V)
      (sa, ia, sb, ib, T2b.cptr(), W13a_sigma_ccvv_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W14(x,c0,a0,b) += (    1.00000000) V2(b,x,c0,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (   -4.00000000) T2(w,c0,a,a0) W14(x,c0,a0,b) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W14cc_sigma_ccvv_ccov(orz::mr::sizeof_sympack_Xcc(symblockinfo, sa0^sb));
    FC_FUNC(g_if_sigma_ccvv_ccov_no2_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO2_X0_TYPE1_ERI_V)
      (sa0, ia0, sb, ib, V2_sym.cptr(), W14cc_sigma_ccvv_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccvv_ccov_no2_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO2_X1_TYPE1_ERI_V)
      (sa0, ia0, sb, ib, T2b.cptr(), W14cc_sigma_ccvv_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W15(w,c0,a0,b) += (    1.00000000) V2(b,w,c0,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (    2.00000000) T2(x,c0,a,a0) W15(w,c0,a0,b) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W15cc_sigma_ccvv_ccov(orz::mr::sizeof_sympack_Xcc(symblockinfo, sa0^sb));
    FC_FUNC(g_if_sigma_ccvv_ccov_no3_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO3_X0_TYPE1_ERI_V)
      (sa0, ia0, sb, ib, V2_sym.cptr(), W15cc_sigma_ccvv_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccvv_ccov_no3_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO3_X1_TYPE1_ERI_V)
      (sa0, ia0, sb, ib, T2b.cptr(), W15cc_sigma_ccvv_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W16(w,c0,a0,b) += (    1.00000000) V2(b,w,c0,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (   -1.00000000) T2(x,c0,a0,a) W16(w,c0,a0,b) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  orz::DTensor W16cca_sigma_ccvv_ccov(orz::mr::sizeof_sympack_Xcca(symblockinfo, sb));
  FC_FUNC(g_if_sigma_ccvv_ccov_no4_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO4_X0_TYPE1_ERI_V)
    (sb, ib, V2_sym.cptr(), W16cca_sigma_ccvv_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    FC_FUNC(g_if_sigma_ccvv_ccov_no4_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO4_X1_TYPE1_ERI_V)
      (sa, ia, sb, ib, T2b.cptr(), W16cca_sigma_ccvv_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
  } // End ia
  } // End sa
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W17(x,c0,a0,b) += (    1.00000000) V2(b,x,c0,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (    2.00000000) T2(w,c0,a0,a) W17(x,c0,a0,b) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  orz::DTensor W17cca_sigma_ccvv_ccov(orz::mr::sizeof_sympack_Xcca(symblockinfo, sb));
  FC_FUNC(g_if_sigma_ccvv_ccov_no5_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO5_X0_TYPE1_ERI_V)
    (sb, ib, V2_sym.cptr(), W17cca_sigma_ccvv_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    FC_FUNC(g_if_sigma_ccvv_ccov_no5_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO5_X1_TYPE1_ERI_V)
      (sa, ia, sb, ib, T2b.cptr(), W17cca_sigma_ccvv_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
  } // End ia
  } // End sa
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W18(w,c0,a0,b) += (    1.00000000) V2(b,a1,w,c0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (   -1.00000000) T2(x,c0,a,a0) W18(w,c0,a0,b) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W18cc_sigma_ccvv_ccov(orz::mr::sizeof_sympack_Xcc(symblockinfo, sa0^sb));
    FC_FUNC(g_if_sigma_ccvv_ccov_no6_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO6_X0_TYPE1_ERI_V)
      (sa0, ia0, sb, ib, V2_sym.cptr(), W18cc_sigma_ccvv_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccvv_ccov_no6_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO6_X1_TYPE1_ERI_V)
      (sa0, ia0, sb, ib, T2b.cptr(), W18cc_sigma_ccvv_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W19(x,c0,a0,b) += (    1.00000000) V2(b,a1,x,c0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (    2.00000000) T2(w,c0,a,a0) W19(x,c0,a0,b) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W19cc_sigma_ccvv_ccov(orz::mr::sizeof_sympack_Xcc(symblockinfo, sa0^sb));
    FC_FUNC(g_if_sigma_ccvv_ccov_no7_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO7_X0_TYPE1_ERI_V)
      (sa0, ia0, sb, ib, V2_sym.cptr(), W19cc_sigma_ccvv_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccvv_ccov_no7_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO7_X1_TYPE1_ERI_V)
      (sa0, ia0, sb, ib, T2b.cptr(), W19cc_sigma_ccvv_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| W20(w,c0,a0,b) += (    1.00000000) V2(b,a1,w,c0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (    2.00000000) T2(x,c0,a0,a) W20(w,c0,a0,b) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  orz::DTensor W20cca_sigma_ccvv_ccov(orz::mr::sizeof_sympack_Xcca(symblockinfo, sb));
  FC_FUNC(g_if_sigma_ccvv_ccov_no8_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO8_X0_TYPE1_ERI_V)
    (sb, ib, V2_sym.cptr(), W20cca_sigma_ccvv_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    FC_FUNC(g_if_sigma_ccvv_ccov_no8_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO8_X1_TYPE1_ERI_V)
      (sa, ia, sb, ib, T2b.cptr(), W20cca_sigma_ccvv_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
  } // End ia
  } // End sa
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W21(x,c0,a0,b) += (    1.00000000) V2(b,a1,x,c0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (   -1.00000000) T2(w,c0,a0,a) W21(x,c0,a0,b) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  orz::DTensor W21cca_sigma_ccvv_ccov(orz::mr::sizeof_sympack_Xcca(symblockinfo, sb));
  FC_FUNC(g_if_sigma_ccvv_ccov_no9_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO9_X0_TYPE1_ERI_V)
    (sb, ib, V2_sym.cptr(), W21cca_sigma_ccvv_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    FC_FUNC(g_if_sigma_ccvv_ccov_no9_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO9_X1_TYPE1_ERI_V)
      (sa, ia, sb, ib, T2b.cptr(), W21cca_sigma_ccvv_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
  } // End ia
  } // End sa
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| W22(a0,b) += (    1.00000000) V2(b,a2,a3,a1) D2(a3,a1,a2,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (    1.00000000) T2(x,w,a,a0) W22(a0,b) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    double W22_sigma_ccvv_ccov(0);
    FC_FUNC(g_if_sigma_ccvv_ccov_no10_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO10_X0_TYPE1_ERI_V)
      (sa0, ia0, sb, ib, V2_sym.cptr(), &W22_sigma_ccvv_ccov, nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccvv_ccov_no10_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO10_X1_TYPE1_ERI_V)
      (sa0, ia0, sb, ib, T2b.cptr(), &W22_sigma_ccvv_ccov, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| W23(a0,b) += (    1.00000000) V2(b,a2,a3,a1) D2(a3,a1,a2,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (   -2.00000000) T2(x,w,a0,a) W23(a0,b) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  orz::DTensor W23a_sigma_ccvv_ccov(orz::mr::sizeof_sympack_Xa(symblockinfo, sb));
  FC_FUNC(g_if_sigma_ccvv_ccov_no11_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO11_X0_TYPE1_ERI_V)
    (sb, ib, V2_sym.cptr(), W23a_sigma_ccvv_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    FC_FUNC(g_if_sigma_ccvv_ccov_no11_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO11_X1_TYPE1_ERI_V)
      (sa, ia, sb, ib, T2b.cptr(), W23a_sigma_ccvv_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
  } // End ia
  } // End sa
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| S2(w,x,a,b) += (    8.00000000) V2(b,x,c0,a0) T2(w,c0,a,a0) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccvv_ccov_no12_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO12_X0_TYPE1_ERI_V)
      (sa0, ia0, sb, ib, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| S2(w,x,a,b) += (   -4.00000000) V2(b,w,c0,a0) T2(x,c0,a,a0) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccvv_ccov_no13_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO13_X0_TYPE1_ERI_V)
      (sa0, ia0, sb, ib, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   14] -- 
  // |-- [    0] --| S2(w,x,a,b) += (    2.00000000) V2(b,w,c0,a0) T2(x,c0,a0,a) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    FC_FUNC(g_if_sigma_ccvv_ccov_no14_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO14_X0_TYPE1_ERI_V)
      (sa, ia, sb, ib, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
  } // End ia
  } // End sa
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   15] -- 
  // |-- [    0] --| S2(w,x,a,b) += (   -4.00000000) V2(b,x,c0,a0) T2(w,c0,a0,a) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    FC_FUNC(g_if_sigma_ccvv_ccov_no15_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO15_X0_TYPE1_ERI_V)
      (sa, ia, sb, ib, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
  } // End ia
  } // End sa
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   16] -- 
  // |-- [    0] --| W26(a0,b) += (    1.00000000) V2(b,a1,a2,a0) D1(a2,a1) 
  // |-- [    1] --| S2(w,x,a,b) += (    1.00000000) T2(x,w,a,a0) W26(a0,b) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    double W26_sigma_ccvv_ccov(0);
    FC_FUNC(g_if_sigma_ccvv_ccov_no16_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO16_X0_TYPE1_ERI_V)
      (sa0, ia0, sb, ib, V2_sym.cptr(), &W26_sigma_ccvv_ccov, nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccvv_ccov_no16_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO16_X1_TYPE1_ERI_V)
      (sa0, ia0, sb, ib, T2b.cptr(), &W26_sigma_ccvv_ccov, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   17] -- 
  // |-- [    0] --| W27(a0,b) += (    1.00000000) V2(b,a1,a2,a0) D1(a2,a1) 
  // |-- [    1] --| S2(w,x,a,b) += (   -2.00000000) T2(x,w,a0,a) W27(a0,b) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  orz::DTensor W27a_sigma_ccvv_ccov(orz::mr::sizeof_sympack_Xa(symblockinfo, sb));
  FC_FUNC(g_if_sigma_ccvv_ccov_no17_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO17_X0_TYPE1_ERI_V)
    (sb, ib, V2_sym.cptr(), W27a_sigma_ccvv_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    FC_FUNC(g_if_sigma_ccvv_ccov_no17_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO17_X1_TYPE1_ERI_V)
      (sa, ia, sb, ib, T2b.cptr(), W27a_sigma_ccvv_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
  } // End ia
  } // End sa
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   18] -- 
  // |-- [    0] --| S2(w,x,a,b) += (   -4.00000000) V2(b,a0,x,c0) T2(w,c0,a,a0) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccvv_ccov_no18_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO18_X0_TYPE1_ERI_V)
      (sa0, ia0, sb, ib, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   19] -- 
  // |-- [    0] --| S2(w,x,a,b) += (    2.00000000) V2(b,a0,w,c0) T2(x,c0,a,a0) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccvv_ccov_no19_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO19_X0_TYPE1_ERI_V)
      (sa0, ia0, sb, ib, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   20] -- 
  // |-- [    0] --| S2(w,x,a,b) += (    2.00000000) V2(b,a0,x,c0) T2(w,c0,a0,a) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    FC_FUNC(g_if_sigma_ccvv_ccov_no20_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO20_X0_TYPE1_ERI_V)
      (sa, ia, sb, ib, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
  } // End ia
  } // End sa
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   21] -- 
  // |-- [    0] --| S2(w,x,a,b) += (   -4.00000000) V2(b,a0,w,c0) T2(x,c0,a0,a) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    FC_FUNC(g_if_sigma_ccvv_ccov_no21_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO21_X0_TYPE1_ERI_V)
      (sa, ia, sb, ib, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
  } // End ia
  } // End sa
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   22] -- 
  // |-- [    0] --| W28(a1,b) += (    1.00000000) V2(b,a1,a2,a0) D1(a2,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (    4.00000000) T2(x,w,a1,a) W28(a1,b) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  orz::DTensor W28a_sigma_ccvv_ccov(orz::mr::sizeof_sympack_Xa(symblockinfo, sb));
  FC_FUNC(g_if_sigma_ccvv_ccov_no22_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO22_X0_TYPE1_ERI_V)
    (sb, ib, V2_sym.cptr(), W28a_sigma_ccvv_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    FC_FUNC(g_if_sigma_ccvv_ccov_no22_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO22_X1_TYPE1_ERI_V)
      (sa, ia, sb, ib, T2b.cptr(), W28a_sigma_ccvv_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
  } // End ia
  } // End sa
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   23] -- 
  // |-- [    0] --| W29(a1,b) += (    1.00000000) V2(b,a1,a2,a0) D1(a2,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (   -2.00000000) T2(x,w,a,a1) W29(a1,b) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    double W29_sigma_ccvv_ccov(0);
    FC_FUNC(g_if_sigma_ccvv_ccov_no23_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO23_X0_TYPE1_ERI_V)
      (sa1, ia1, sb, ib, V2_sym.cptr(), &W29_sigma_ccvv_ccov, nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccvv_ccov_no23_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO23_X1_TYPE1_ERI_V)
      (sa1, ia1, sb, ib, T2b.cptr(), &W29_sigma_ccvv_ccov, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   24] -- 
  // |-- [    0] --| W32(a0,v0,b,a) += (    1.00000000) V2(v0,b,a1,a) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (   -2.00000000) T2(w,x,a0,v0) W32(a0,v0,b,a) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    orz::DTensor W32av_sigma_ccvv_ccov(orz::mr::sizeof_sympack_Xav(symblockinfo, sv0^sb));
    FC_FUNC(g_if_sigma_ccvv_ccov_no24_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO24_X0_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, V2_sym.cptr(), W32av_sigma_ccvv_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccvv_ccov_no24_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO24_X1_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), W32av_sigma_ccvv_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   25] -- 
  // |-- [    0] --| W33(a0,v0,b,a) += (    1.00000000) V2(v0,b,a1,a) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (    1.00000000) T2(x,w,a0,v0) W33(a0,v0,b,a) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    orz::DTensor W33av_sigma_ccvv_ccov(orz::mr::sizeof_sympack_Xav(symblockinfo, sv0^sb));
    FC_FUNC(g_if_sigma_ccvv_ccov_no25_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO25_X0_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, V2_sym.cptr(), W33av_sigma_ccvv_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccvv_ccov_no25_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO25_X1_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), W33av_sigma_ccvv_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   26] -- 
  // |-- [    0] --| W34(a0,v0,b,a) += (    1.00000000) V2(v0,a,b,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (    1.00000000) T2(w,x,a0,v0) W34(a0,v0,b,a) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    orz::DTensor W34av_sigma_ccvv_ccov(orz::mr::sizeof_sympack_Xav(symblockinfo, sv0^sb));
    FC_FUNC(g_if_sigma_ccvv_ccov_no26_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO26_X0_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, V2_sym.cptr(), W34av_sigma_ccvv_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccvv_ccov_no26_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO26_X1_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), W34av_sigma_ccvv_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   27] -- 
  // |-- [    0] --| W35(a0,v0,b,a) += (    1.00000000) V2(v0,a,b,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (   -2.00000000) T2(x,w,a0,v0) W35(a0,v0,b,a) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    orz::DTensor W35av_sigma_ccvv_ccov(orz::mr::sizeof_sympack_Xav(symblockinfo, sv0^sb));
    FC_FUNC(g_if_sigma_ccvv_ccov_no27_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO27_X0_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, V2_sym.cptr(), W35av_sigma_ccvv_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccvv_ccov_no27_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO27_X1_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), W35av_sigma_ccvv_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   28] -- 
  // |-- [    0] --| S2(w,x,a,b) += (    4.00000000) V2(v0,b,a0,a) T2(w,x,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccov_no28_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO28_X0_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   29] -- 
  // |-- [    0] --| S2(w,x,a,b) += (   -2.00000000) V2(v0,b,a0,a) T2(x,w,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccov_no29_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO29_X0_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   30] -- 
  // |-- [    0] --| S2(w,x,a,b) += (    4.00000000) V2(v0,a,b,a0) T2(x,w,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccov_no30_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO30_X0_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   31] -- 
  // |-- [    0] --| S2(w,x,a,b) += (   -2.00000000) V2(v0,a,b,a0) T2(w,x,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccov_no31_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCOV_NO31_X0_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
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
  orz::mr::Femto::my_timer.push_back(boost::make_tuple("sigma_ccvv_ccov", time_sigma.elapsed_cputime(), time_sigma.elapsed_wallclocktime()));
  orz::mr::Femto::file_timing << "* " << boost::format("%20s : %10.7f %10.7f ") % orz::mr::Femto::my_timer.back().get<0>() % orz::mr::Femto::my_timer.back().get<1>() % orz::mr::Femto::my_timer.back().get<2>() << endl;
  flush(orz::mr::Femto::file_timing);

  return retval; 
} 
