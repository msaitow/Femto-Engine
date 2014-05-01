                                                                                
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
#include <sci/icmr/Femto/elems/c_sigma_ccvv_ccvv.h>                                  
                                                                                
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
                                                                                
orz::mr::BareAmpPack orz::mr::Femto::sigma_ccvv_ccvv(const orz::mr::Input &ctinp,                                    
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
  // |-- [    0] --| S2(w,x,a,b) += (    8.00000000) Fc0 T2(w,x,a,b) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x0_type0_noeri,G_IF_SIGMA_CCVV_CCVV_NO0_X0_TYPE0_NOERI)
      (sb, ib, &Fc0, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| S2(w,x,a,b) += (   -4.00000000) Fc0 T2(x,w,a,b) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x0_type0_noeri,G_IF_SIGMA_CCVV_CCVV_NO1_X0_TYPE0_NOERI)
      (sb, ib, &Fc0, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| S2(w,x,a,b) += (   -8.00000000) T2(w,c0,a,b) Fc1(x,c0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no2_x0_type0_noeri,G_IF_SIGMA_CCVV_CCVV_NO2_X0_TYPE0_NOERI)
      (sb, ib, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| S2(w,x,a,b) += (    4.00000000) T2(x,c0,a,b) Fc1(w,c0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no3_x0_type0_noeri,G_IF_SIGMA_CCVV_CCVV_NO3_X0_TYPE0_NOERI)
      (sb, ib, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| S2(w,x,a,b) += (    4.00000000) T2(c0,w,a,b) Fc1(x,c0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no4_x0_type0_noeri,G_IF_SIGMA_CCVV_CCVV_NO4_X0_TYPE0_NOERI)
      (sb, ib, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| S2(w,x,a,b) += (   -8.00000000) T2(c0,x,a,b) Fc1(w,c0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no5_x0_type0_noeri,G_IF_SIGMA_CCVV_CCVV_NO5_X0_TYPE0_NOERI)
      (sb, ib, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W0() += (    1.00000000) D1(a1,a0) Fc1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (    8.00000000) T2(w,x,a,b) W0() 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  double W0_sigma_ccvv_ccvv(0);
  FC_FUNC(g_if_sigma_ccvv_ccvv_no6_x0_type0_noeri,G_IF_SIGMA_CCVV_CCVV_NO6_X0_TYPE0_NOERI)
    (&W0_sigma_ccvv_ccvv, nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no6_x1_type0_noeri,G_IF_SIGMA_CCVV_CCVV_NO6_X1_TYPE0_NOERI)
      (sb, ib, T2b.cptr(), &W0_sigma_ccvv_ccvv, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W1() += (    1.00000000) D1(a1,a0) Fc1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (   -4.00000000) T2(x,w,a,b) W1() 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  double W1_sigma_ccvv_ccvv(0);
  FC_FUNC(g_if_sigma_ccvv_ccvv_no7_x0_type0_noeri,G_IF_SIGMA_CCVV_CCVV_NO7_X0_TYPE0_NOERI)
    (&W1_sigma_ccvv_ccvv, nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no7_x1_type0_noeri,G_IF_SIGMA_CCVV_CCVV_NO7_X1_TYPE0_NOERI)
      (sb, ib, T2b.cptr(), &W1_sigma_ccvv_ccvv, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| S2(x,w,b,a) += (    8.00000000) T2(x,w,v0,a) Fc1(v0,b) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no8_x0_type0_noeri,G_IF_SIGMA_CCVV_CCVV_NO8_X0_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| S2(x,w,b,a) += (   -4.00000000) T2(w,x,v0,a) Fc1(v0,b) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no9_x0_type0_noeri,G_IF_SIGMA_CCVV_CCVV_NO9_X0_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| S2(w,x,a,b) += (    8.00000000) T2(w,x,v0,b) Fc1(v0,a) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no10_x0_type0_noeri,G_IF_SIGMA_CCVV_CCVV_NO10_X0_TYPE0_NOERI)
      (sb, ib, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| S2(w,x,a,b) += (   -4.00000000) T2(x,w,v0,b) Fc1(v0,a) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no11_x0_type0_noeri,G_IF_SIGMA_CCVV_CCVV_NO11_X0_TYPE0_NOERI)
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
  // -- Title : sigma_ccvv_ccvv

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| S2(w,x,a,b) += (    8.00000000) V2(x,c1,w,c0) T2(c0,c1,a,b) 
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCVV_NO0_X0_TYPE1_ERI_C)
      (sb, ib, sx, ix, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| S2(w,x,a,b) += (   -4.00000000) V2(x,c0,w,c1) T2(c0,c1,a,b) 
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCVV_NO1_X0_TYPE1_ERI_C)
      (sb, ib, sx, ix, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W2(x,c0) += (    1.00000000) V2(x,a1,c0,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (    4.00000000) T2(w,c0,a,b) W2(x,c0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W2c_sigma_ccvv_ccvv(orz::mr::sizeof_sympack_Xc(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccvv_ccvv_no2_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCVV_NO2_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W2c_sigma_ccvv_ccvv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no2_x1_type1_eri_c,G_IF_SIGMA_CCVV_CCVV_NO2_X1_TYPE1_ERI_C)
      (sb, ib, sx, ix, T2b.cptr(), W2c_sigma_ccvv_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W3(w,c0) += (    1.00000000) V2(w,a1,c0,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (   -2.00000000) T2(x,c0,a,b) W3(w,c0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W3c_sigma_ccvv_ccvv(orz::mr::sizeof_sympack_Xc(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccvv_ccvv_no3_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCVV_NO3_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W3c_sigma_ccvv_ccvv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no3_x1_type1_eri_c,G_IF_SIGMA_CCVV_CCVV_NO3_X1_TYPE1_ERI_C)
      (sb, ib, sw, iw, T2b.cptr(), W3c_sigma_ccvv_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W4(x,c0) += (    1.00000000) V2(x,a1,c0,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (   -2.00000000) T2(c0,w,a,b) W4(x,c0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W4c_sigma_ccvv_ccvv(orz::mr::sizeof_sympack_Xc(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccvv_ccvv_no4_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCVV_NO4_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W4c_sigma_ccvv_ccvv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no4_x1_type1_eri_c,G_IF_SIGMA_CCVV_CCVV_NO4_X1_TYPE1_ERI_C)
      (sb, ib, sx, ix, T2b.cptr(), W4c_sigma_ccvv_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W5(w,c0) += (    1.00000000) V2(w,a1,c0,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (    4.00000000) T2(c0,x,a,b) W5(w,c0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W5c_sigma_ccvv_ccvv(orz::mr::sizeof_sympack_Xc(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccvv_ccvv_no5_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCVV_NO5_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W5c_sigma_ccvv_ccvv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no5_x1_type1_eri_c,G_IF_SIGMA_CCVV_CCVV_NO5_X1_TYPE1_ERI_C)
      (sb, ib, sw, iw, T2b.cptr(), W5c_sigma_ccvv_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W6(x,c0) += (    1.00000000) V2(x,c0,a1,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (   -8.00000000) T2(w,c0,a,b) W6(x,c0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W6c_sigma_ccvv_ccvv(orz::mr::sizeof_sympack_Xc(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccvv_ccvv_no6_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCVV_NO6_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W6c_sigma_ccvv_ccvv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no6_x1_type1_eri_c,G_IF_SIGMA_CCVV_CCVV_NO6_X1_TYPE1_ERI_C)
      (sb, ib, sx, ix, T2b.cptr(), W6c_sigma_ccvv_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W7(w,c0) += (    1.00000000) V2(w,c0,a1,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (    4.00000000) T2(x,c0,a,b) W7(w,c0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W7c_sigma_ccvv_ccvv(orz::mr::sizeof_sympack_Xc(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccvv_ccvv_no7_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCVV_NO7_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W7c_sigma_ccvv_ccvv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no7_x1_type1_eri_c,G_IF_SIGMA_CCVV_CCVV_NO7_X1_TYPE1_ERI_C)
      (sb, ib, sw, iw, T2b.cptr(), W7c_sigma_ccvv_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| W8(x,c0) += (    1.00000000) V2(x,c0,a1,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (    4.00000000) T2(c0,w,a,b) W8(x,c0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W8c_sigma_ccvv_ccvv(orz::mr::sizeof_sympack_Xc(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccvv_ccvv_no8_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCVV_NO8_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W8c_sigma_ccvv_ccvv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no8_x1_type1_eri_c,G_IF_SIGMA_CCVV_CCVV_NO8_X1_TYPE1_ERI_C)
      (sb, ib, sx, ix, T2b.cptr(), W8c_sigma_ccvv_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W9(w,c0) += (    1.00000000) V2(w,c0,a1,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (   -8.00000000) T2(c0,x,a,b) W9(w,c0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W9c_sigma_ccvv_ccvv(orz::mr::sizeof_sympack_Xc(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccvv_ccvv_no9_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCVV_NO9_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W9c_sigma_ccvv_ccvv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no9_x1_type1_eri_c,G_IF_SIGMA_CCVV_CCVV_NO9_X1_TYPE1_ERI_C)
      (sb, ib, sw, iw, T2b.cptr(), W9c_sigma_ccvv_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| S2(w,x,a,b) += (   16.00000000) V2(w,a,v0,c0) T2(c0,x,v0,b) 
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_no10_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCVV_NO10_X0_TYPE1_ERI_C)
      (sb, ib, sw, iw, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| S2(w,x,a,b) += (   -8.00000000) V2(x,a,v0,c0) T2(c0,w,v0,b) 
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_no11_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCVV_NO11_X0_TYPE1_ERI_C)
      (sb, ib, sx, ix, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| S2(w,x,a,b) += (   -8.00000000) V2(w,a,v0,c0) T2(x,c0,v0,b) 
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_no12_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCVV_NO12_X0_TYPE1_ERI_C)
      (sb, ib, sw, iw, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| S2(w,x,a,b) += (    4.00000000) V2(x,a,v0,c0) T2(w,c0,v0,b) 
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_no13_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCVV_NO13_X0_TYPE1_ERI_C)
      (sb, ib, sx, ix, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   14] -- 
  // |-- [    0] --| S2(w,x,a,b) += (   -8.00000000) V2(x,c0,v0,a) T2(w,c0,v0,b) 
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_no14_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCVV_NO14_X0_TYPE1_ERI_C)
      (sb, ib, sx, ix, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   15] -- 
  // |-- [    0] --| S2(w,x,a,b) += (    4.00000000) V2(w,c0,v0,a) T2(x,c0,v0,b) 
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_no15_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCVV_NO15_X0_TYPE1_ERI_C)
      (sb, ib, sw, iw, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   16] -- 
  // |-- [    0] --| S2(w,x,a,b) += (    4.00000000) V2(x,c0,v0,a) T2(c0,w,v0,b) 
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_no16_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCVV_NO16_X0_TYPE1_ERI_C)
      (sb, ib, sx, ix, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   17] -- 
  // |-- [    0] --| S2(w,x,a,b) += (   -8.00000000) V2(w,c0,v0,a) T2(c0,x,v0,b) 
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_no17_x0_type1_eri_c,G_IF_SIGMA_CCVV_CCVV_NO17_X0_TYPE1_ERI_C)
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

//-@type(2).declaration(begin)
  // --  Title : sigma_ccvv_ccvv
  //  >> Intermediates for the external contractions are defined here << 
  double W10_sigma_ccvv_ccvv(0);
  double W11_sigma_ccvv_ccvv(0);
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
  // -- Title : sigma_ccvv_ccvv

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W10() += (    1.00000000) V2(a3,a1,a2,a0) D2(a3,a1,a2,a0) 
  int sa3(s_eri);
  int ia3(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x0_type1_eri_o,G_IF_SIGMA_CCVV_CCVV_NO0_X0_TYPE1_ERI_O)
    (sa3, ia3, V2_sym.cptr(), &W10_sigma_ccvv_ccvv, nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W11() += (    1.00000000) V2(a3,a1,a2,a0) D2(a3,a1,a2,a0) 
  int sa3(s_eri);
  int ia3(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x0_type1_eri_o,G_IF_SIGMA_CCVV_CCVV_NO1_X0_TYPE1_ERI_O)
    (sa3, ia3, V2_sym.cptr(), &W11_sigma_ccvv_ccvv, nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

//-@type(1).contraction(end)

  } // End myrank
  }
  }
  orz::world().barrier();

//-@type(2).contraction(begin)
  // -- Title : sigma_ccvv_ccvv
  //*-- Entering to take the type 2 contractions --*//

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| S2(w,x,a,b) += (    4.00000000) T2(w,x,a,b) W10() 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x0_type2_eri_o,G_IF_SIGMA_CCVV_CCVV_NO0_X0_TYPE2_ERI_O)
      (sb, ib, T2b.cptr(), &W10_sigma_ccvv_ccvv, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| S2(w,x,a,b) += (   -2.00000000) T2(x,w,a,b) W11() 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x0_type2_eri_o,G_IF_SIGMA_CCVV_CCVV_NO1_X0_TYPE2_ERI_O)
      (sb, ib, T2b.cptr(), &W11_sigma_ccvv_ccvv, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
//-@type(2).contraction(end)
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
  // -- Title : sigma_ccvv_ccvv

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W12(v0,a) += (    1.00000000) V2(v0,a0,a1,a) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (   -4.00000000) T2(x,w,b,v0) W12(v0,a) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W12v_sigma_ccvv_ccvv(orz::mr::sizeof_sympack_Xv(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCVV_NO0_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W12v_sigma_ccvv_ccvv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCVV_NO0_X1_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), W12v_sigma_ccvv_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W13(v0,a) += (    1.00000000) V2(v0,a0,a1,a) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (    2.00000000) T2(x,w,v0,b) W13(v0,a) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W13v_sigma_ccvv_ccvv(orz::mr::sizeof_sympack_Xv(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCVV_NO1_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W13v_sigma_ccvv_ccvv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCVV_NO1_X1_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), W13v_sigma_ccvv_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| S2(w,x,a,b) += (   16.00000000) V2(b,x,v0,c0) T2(w,c0,a,v0) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no2_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCVV_NO2_X0_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| S2(w,x,a,b) += (   -8.00000000) V2(b,w,v0,c0) T2(x,c0,a,v0) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no3_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCVV_NO3_X0_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| S2(w,x,a,b) += (    4.00000000) V2(b,w,v0,c0) T2(x,c0,v0,a) 
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_no4_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCVV_NO4_X0_TYPE1_ERI_V)
      (sa, ia, sb, ib, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
  } // End ia
  } // End sa
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| S2(w,x,a,b) += (   -8.00000000) V2(b,x,v0,c0) T2(w,c0,v0,a) 
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_no5_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCVV_NO5_X0_TYPE1_ERI_V)
      (sa, ia, sb, ib, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
  } // End ia
  } // End sa
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W14(v0,b) += (    1.00000000) V2(v0,a1,b,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (    2.00000000) T2(x,w,a,v0) W14(v0,b) 
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
    double W14_sigma_ccvv_ccvv(0);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no6_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCVV_NO6_X0_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, V2_sym.cptr(), &W14_sigma_ccvv_ccvv, nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no6_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCVV_NO6_X1_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), &W14_sigma_ccvv_ccvv, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W15(v0,b) += (    1.00000000) V2(v0,a1,b,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (   -4.00000000) T2(w,x,a,v0) W15(v0,b) 
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
    double W15_sigma_ccvv_ccvv(0);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no7_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCVV_NO7_X0_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, V2_sym.cptr(), &W15_sigma_ccvv_ccvv, nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no7_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCVV_NO7_X1_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), &W15_sigma_ccvv_ccvv, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| S2(w,x,a,b) += (   -8.00000000) V2(v0,b,x,c0) T2(w,c0,a,v0) 
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_no8_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCVV_NO8_X0_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| S2(w,x,a,b) += (    4.00000000) V2(v0,b,w,c0) T2(x,c0,a,v0) 
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_no9_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCVV_NO9_X0_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| S2(w,x,a,b) += (    4.00000000) V2(v0,b,x,c0) T2(c0,w,a,v0) 
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_no10_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCVV_NO10_X0_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| S2(w,x,a,b) += (   -8.00000000) V2(v0,b,w,c0) T2(c0,x,a,v0) 
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_no11_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCVV_NO11_X0_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| W16(v0,b) += (    1.00000000) V2(v0,b,a1,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (    8.00000000) T2(w,x,a,v0) W16(v0,b) 
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
    double W16_sigma_ccvv_ccvv(0);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no12_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCVV_NO12_X0_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, V2_sym.cptr(), &W16_sigma_ccvv_ccvv, nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no12_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCVV_NO12_X1_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), &W16_sigma_ccvv_ccvv, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| W17(v0,b) += (    1.00000000) V2(v0,b,a1,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (   -4.00000000) T2(x,w,a,v0) W17(v0,b) 
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
    double W17_sigma_ccvv_ccvv(0);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no13_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCVV_NO13_X0_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, V2_sym.cptr(), &W17_sigma_ccvv_ccvv, nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no13_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCVV_NO13_X1_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), &W17_sigma_ccvv_ccvv, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   14] -- 
  // |-- [    0] --| W18(v0,a) += (    1.00000000) V2(v0,a,a1,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (    8.00000000) T2(x,w,b,v0) W18(v0,a) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W18v_sigma_ccvv_ccvv(orz::mr::sizeof_sympack_Xv(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ccvv_ccvv_no14_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCVV_NO14_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W18v_sigma_ccvv_ccvv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no14_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCVV_NO14_X1_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), W18v_sigma_ccvv_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   15] -- 
  // |-- [    0] --| W19(v0,a) += (    1.00000000) V2(v0,a,a1,a0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (   -4.00000000) T2(x,w,v0,b) W19(v0,a) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W19v_sigma_ccvv_ccvv(orz::mr::sizeof_sympack_Xv(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_ccvv_ccvv_no15_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCVV_NO15_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W19v_sigma_ccvv_ccvv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no15_x1_type1_eri_v,G_IF_SIGMA_CCVV_CCVV_NO15_X1_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), W19v_sigma_ccvv_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   16] -- 
  // |-- [    0] --| S2(w,x,a,b) += (    8.00000000) V2(v1,b,v0,a) T2(w,x,v0,v1) 
  int sv1(s_eri);
  int iv1(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv1);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no16_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCVV_NO16_X0_TYPE1_ERI_V)
      (sb, ib, sv1, iv1, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   17] -- 
  // |-- [    0] --| S2(w,x,a,b) += (   -4.00000000) V2(v1,a,v0,b) T2(w,x,v0,v1) 
  int sv1(s_eri);
  int iv1(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv1);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no17_x0_type1_eri_v,G_IF_SIGMA_CCVV_CCVV_NO17_X0_TYPE1_ERI_V)
      (sb, ib, sv1, iv1, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  orz::mr::Femto::my_timer.push_back(boost::make_tuple("sigma_ccvv_ccvv", time_sigma.elapsed_cputime(), time_sigma.elapsed_wallclocktime()));
  orz::mr::Femto::file_timing << "* " << boost::format("%20s : %10.7f %10.7f ") % orz::mr::Femto::my_timer.back().get<0>() % orz::mr::Femto::my_timer.back().get<1>() % orz::mr::Femto::my_timer.back().get<2>() << endl;
  flush(orz::mr::Femto::file_timing);

  return retval; 
} 
