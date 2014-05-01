                                                                                
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
#include <sci/icmr/Femto/elems/c_sigma_covv_covv.h>                                  
                                                                                
using std::cout;                                                                
using std::endl;                                                                
                                                                                
#define FLOPCOUNT                                                               
                                                                                
using orz::mr::Femto::Fc0;                                                      
using orz::mr::Femto::h1_int;                                                   
using orz::mr::Femto::h6_int;                                                   
                                                                                
//  8 8888888888   8 8888888888            ,8.       ,8.    8888888 8888888888 ,o888888o.     
//  8 8888         8 8888                 ,888.     ,888.         8 8888    . 8888     `88.   
//  8 8888         8 8888                .`8888.   .`8888.        8 8888   ,8 8888       `8b  
//  8 8888         8 8888               ,8.`8888. ,8.`8888.       8 8888   88 8888        `8b 
//  8 888888888888 8 888888888888      ,8'8.`8888,8^8.`8888.      8 8888   88 8888         88 
//  8 8888         8 8888             ,8' `8.`8888' `8.`8888.     8 8888   88 8888         88 
//  8 8888         8 8888            ,8'   `8.`88'   `8.`8888.    8 8888   88 8888        ,8P 
//  8 8888         8 8888           ,8'     `8.`'     `8.`8888.   8 8888   `8 8888       ,8P  
//  8 8888         8 8888          ,8'       `8        `8.`8888.  8 8888    ` 8888     ,88'   
//  8 8888         8 888888888888 ,8'         `         `8.`8888. 8 8888       `8888888P'     

//                                   Generated date : Sun Apr 20 10:26:20 2014

                                                                                
// ***************************************************************************  
// orz::mr::mrci                                                                
// ***************************************************************************  
									    /*!        
   @brief CT input                                                              
									     */        
                                                                                
orz::mr::BareAmpPack orz::mr::Femto::sigma_covv_covv(const orz::mr::Input &ctinp,                                    
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
  // |-- [    0] --| S2(w,i,a,b) += (    2.00000000) Fc0 T2(w,a0,a,b) D1(i,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_covv_no0_x0_type0_noeri,G_IF_SIGMA_COVV_COVV_NO0_X0_TYPE0_NOERI)
      (sb, ib, &Fc0, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| S2(w,i,a,b) += (   -1.00000000) Fc0 T2(a0,w,a,b) D1(i,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_covv_no1_x0_type0_noeri,G_IF_SIGMA_COVV_COVV_NO1_X0_TYPE0_NOERI)
      (sb, ib, &Fc0, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W0(c0,i,a,b) += (    1.00000000) T2(c0,a0,a,b) D1(i,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (   -2.00000000) Fc1(w,c0) W0(c0,i,a,b) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    orz::DTensor W0cav_sigma_covv_covv(orz::mr::sizeof_sympack_Xcav(symblockinfo, sb));
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_covv_no2_x0_type0_noeri,G_IF_SIGMA_COVV_COVV_NO2_X0_TYPE0_NOERI)
      (sb, ib, T2b.cptr(), W0cav_sigma_covv_covv.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_covv_covv_no2_x1_type0_noeri,G_IF_SIGMA_COVV_COVV_NO2_X1_TYPE0_NOERI)
      (sb, ib, W0cav_sigma_covv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W1(i,a0) += (    1.00000000) D2(i,a0,a2,a1) Fc1(a2,a1) 
  // |-- [    1] --| S2(w,i,a,b) += (    2.00000000) T2(w,a0,a,b) W1(i,a0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W1aa_sigma_covv_covv(orz::mr::sizeof_sympack_Xaa(symblockinfo, 0));
  FC_FUNC(g_if_sigma_covv_covv_no3_x0_type0_noeri,G_IF_SIGMA_COVV_COVV_NO3_X0_TYPE0_NOERI)
    (W1aa_sigma_covv_covv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_covv_no3_x1_type0_noeri,G_IF_SIGMA_COVV_COVV_NO3_X1_TYPE0_NOERI)
      (sb, ib, T2b.cptr(), W1aa_sigma_covv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W2(c0,i,b,a) += (    1.00000000) T2(a0,c0,a,b) D1(i,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (    1.00000000) Fc1(w,c0) W2(c0,i,b,a) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    orz::DTensor W2cav_sigma_covv_covv(orz::mr::sizeof_sympack_Xcav(symblockinfo, sb));
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_covv_no4_x0_type0_noeri,G_IF_SIGMA_COVV_COVV_NO4_X0_TYPE0_NOERI)
      (sb, ib, T2b.cptr(), W2cav_sigma_covv_covv.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_covv_covv_no4_x1_type0_noeri,G_IF_SIGMA_COVV_COVV_NO4_X1_TYPE0_NOERI)
      (sb, ib, W2cav_sigma_covv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W3(i,a0) += (    1.00000000) D2(i,a0,a2,a1) Fc1(a2,a1) 
  // |-- [    1] --| S2(w,i,a,b) += (   -1.00000000) T2(a0,w,a,b) W3(i,a0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W3aa_sigma_covv_covv(orz::mr::sizeof_sympack_Xaa(symblockinfo, 0));
  FC_FUNC(g_if_sigma_covv_covv_no5_x0_type0_noeri,G_IF_SIGMA_COVV_COVV_NO5_X0_TYPE0_NOERI)
    (W3aa_sigma_covv_covv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_covv_no5_x1_type0_noeri,G_IF_SIGMA_COVV_COVV_NO5_X1_TYPE0_NOERI)
      (sb, ib, T2b.cptr(), W3aa_sigma_covv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W4(w,i,v0,a) += (    1.00000000) T2(w,a0,v0,a) D1(i,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (   -1.00000000) Fc1(v0,b) W4(w,i,v0,a) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W4cav_sigma_covv_covv(orz::mr::sizeof_sympack_Xcav(symblockinfo, sa));
    FC_FUNC(g_if_sigma_covv_covv_no6_x0_type0_noeri,G_IF_SIGMA_COVV_COVV_NO6_X0_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W4cav_sigma_covv_covv.cptr(), nir, nsym, psym, &flops);
    for(int sb = 0;sb < nir;++sb){ 
    for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
      S2b = orz::DTensor(retval.namps_iamp()[ib]);
      FC_FUNC(g_if_sigma_covv_covv_no6_x1_type0_noeri,G_IF_SIGMA_COVV_COVV_NO6_X1_TYPE0_NOERI)
        (sa, ia, sb, ib, W4cav_sigma_covv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- No. [    7] -- 
  // |-- [    0] --| W5(w,i,v0,b) += (    1.00000000) T2(w,a0,v0,b) D1(i,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (    2.00000000) Fc1(v0,a) W5(w,i,v0,b) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    orz::DTensor W5cav_sigma_covv_covv(orz::mr::sizeof_sympack_Xcav(symblockinfo, sb));
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_covv_no7_x0_type0_noeri,G_IF_SIGMA_COVV_COVV_NO7_X0_TYPE0_NOERI)
      (sb, ib, T2b.cptr(), W5cav_sigma_covv_covv.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_covv_covv_no7_x1_type0_noeri,G_IF_SIGMA_COVV_COVV_NO7_X1_TYPE0_NOERI)
      (sb, ib, W5cav_sigma_covv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| W6(w,i,a,v0) += (    1.00000000) T2(w,a0,a,v0) D1(i,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (    2.00000000) Fc1(v0,b) W6(w,i,a,v0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    orz::DTensor W6cav_sigma_covv_covv(orz::mr::sizeof_sympack_Xcav(symblockinfo, sv0));
    FC_FUNC(g_if_sigma_covv_covv_no8_x0_type0_noeri,G_IF_SIGMA_COVV_COVV_NO8_X0_TYPE0_NOERI)
      (sv0, iv0, T2b.cptr(), W6cav_sigma_covv_covv.cptr(), nir, nsym, psym, &flops);
    for(int sb = 0;sb < nir;++sb){ 
    for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
      S2b = orz::DTensor(retval.namps_iamp()[ib]);
      FC_FUNC(g_if_sigma_covv_covv_no8_x1_type0_noeri,G_IF_SIGMA_COVV_COVV_NO8_X1_TYPE0_NOERI)
        (sb, ib, sv0, iv0, W6cav_sigma_covv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
      retval.acc_amp2(ib, S2b);
    } // End ib
    } // End sb
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W7(w,i,b,v0) += (    1.00000000) T2(a0,w,v0,b) D1(i,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (   -1.00000000) Fc1(v0,a) W7(w,i,b,v0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    orz::DTensor W7cav_sigma_covv_covv(orz::mr::sizeof_sympack_Xcav(symblockinfo, sb));
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_covv_no9_x0_type0_noeri,G_IF_SIGMA_COVV_COVV_NO9_X0_TYPE0_NOERI)
      (sb, ib, T2b.cptr(), W7cav_sigma_covv_covv.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_covv_covv_no9_x1_type0_noeri,G_IF_SIGMA_COVV_COVV_NO9_X1_TYPE0_NOERI)
      (sb, ib, W7cav_sigma_covv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- Title : sigma_covv_covv

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W8(w,c0,i,a0) += (    1.00000000) V2(w,a2,c0,a1) D2(i,a0,a2,a1) 
  // |-- [    1] --| S2(w,i,a,b) += (    1.00000000) T2(c0,a0,a,b) W8(w,c0,i,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W8caa_sigma_covv_covv(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_covv_covv_no0_x0_type1_eri_c,G_IF_SIGMA_COVV_COVV_NO0_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W8caa_sigma_covv_covv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_covv_no0_x1_type1_eri_c,G_IF_SIGMA_COVV_COVV_NO0_X1_TYPE1_ERI_C)
      (sb, ib, sw, iw, T2b.cptr(), W8caa_sigma_covv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W9(w,c0,i,a0) += (    1.00000000) V2(w,c0,a2,a1) D2(i,a0,a2,a1) 
  // |-- [    1] --| S2(w,i,a,b) += (   -2.00000000) T2(c0,a0,a,b) W9(w,c0,i,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W9caa_sigma_covv_covv(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_covv_covv_no1_x0_type1_eri_c,G_IF_SIGMA_COVV_COVV_NO1_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W9caa_sigma_covv_covv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_covv_no1_x1_type1_eri_c,G_IF_SIGMA_COVV_COVV_NO1_X1_TYPE1_ERI_C)
      (sb, ib, sw, iw, T2b.cptr(), W9caa_sigma_covv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W11(w,c0,i,a0) += (    1.00000000) V2(w,a2,c0,a1) D2(i,a1,a2,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (    1.00000000) T2(a0,c0,a,b) W11(w,c0,i,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W11caa_sigma_covv_covv(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_covv_covv_no2_x0_type1_eri_c,G_IF_SIGMA_COVV_COVV_NO2_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W11caa_sigma_covv_covv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_covv_no2_x1_type1_eri_c,G_IF_SIGMA_COVV_COVV_NO2_X1_TYPE1_ERI_C)
      (sb, ib, sw, iw, T2b.cptr(), W11caa_sigma_covv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W12(w,c0,i,a0) += (    1.00000000) V2(w,c0,a2,a1) D2(i,a0,a2,a1) 
  // |-- [    1] --| S2(w,i,a,b) += (    1.00000000) T2(a0,c0,a,b) W12(w,c0,i,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W12caa_sigma_covv_covv(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_covv_covv_no3_x0_type1_eri_c,G_IF_SIGMA_COVV_COVV_NO3_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W12caa_sigma_covv_covv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_covv_no3_x1_type1_eri_c,G_IF_SIGMA_COVV_COVV_NO3_X1_TYPE1_ERI_C)
      (sb, ib, sw, iw, T2b.cptr(), W12caa_sigma_covv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W14(w,a0,b,a) += (    1.00000000) V2(w,a,v0,c0) T2(c0,a0,v0,b) 
  // |-- [    1] --| S2(w,i,a,b) += (    4.00000000) D1(i,a0) W14(w,a0,b,a) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    orz::DTensor W14av_sigma_covv_covv(orz::mr::sizeof_sympack_Xav(symblockinfo, sw^sb));
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_covv_no4_x0_type1_eri_c,G_IF_SIGMA_COVV_COVV_NO4_X0_TYPE1_ERI_C)
      (sb, ib, sw, iw, T2b.cptr(), V2_sym.cptr(), W14av_sigma_covv_covv.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_covv_covv_no4_x1_type1_eri_c,G_IF_SIGMA_COVV_COVV_NO4_X1_TYPE1_ERI_C)
      (sb, ib, sw, iw, W14av_sigma_covv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W20(w,a0,b,a) += (    1.00000000) V2(w,c0,v0,a) T2(c0,a0,v0,b) 
  // |-- [    1] --| S2(w,i,a,b) += (   -2.00000000) D1(i,a0) W20(w,a0,b,a) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    orz::DTensor W20av_sigma_covv_covv(orz::mr::sizeof_sympack_Xav(symblockinfo, sw^sb));
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_covv_no5_x0_type1_eri_c,G_IF_SIGMA_COVV_COVV_NO5_X0_TYPE1_ERI_C)
      (sb, ib, sw, iw, T2b.cptr(), V2_sym.cptr(), W20av_sigma_covv_covv.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_covv_covv_no5_x1_type1_eri_c,G_IF_SIGMA_COVV_COVV_NO5_X1_TYPE1_ERI_C)
      (sb, ib, sw, iw, W20av_sigma_covv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W22(w,a0,b,a) += (    1.00000000) V2(w,a,v0,c0) T2(a0,c0,v0,b) 
  // |-- [    1] --| S2(w,i,a,b) += (   -2.00000000) D1(i,a0) W22(w,a0,b,a) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    orz::DTensor W22av_sigma_covv_covv(orz::mr::sizeof_sympack_Xav(symblockinfo, sw^sb));
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_covv_no6_x0_type1_eri_c,G_IF_SIGMA_COVV_COVV_NO6_X0_TYPE1_ERI_C)
      (sb, ib, sw, iw, T2b.cptr(), V2_sym.cptr(), W22av_sigma_covv_covv.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_covv_covv_no6_x1_type1_eri_c,G_IF_SIGMA_COVV_COVV_NO6_X1_TYPE1_ERI_C)
      (sb, ib, sw, iw, W22av_sigma_covv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W28(w,a0,b,a) += (    1.00000000) V2(w,c0,v0,a) T2(a0,c0,v0,b) 
  // |-- [    1] --| S2(w,i,a,b) += (    1.00000000) D1(i,a0) W28(w,a0,b,a) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    orz::DTensor W28av_sigma_covv_covv(orz::mr::sizeof_sympack_Xav(symblockinfo, sw^sb));
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_covv_no7_x0_type1_eri_c,G_IF_SIGMA_COVV_COVV_NO7_X0_TYPE1_ERI_C)
      (sb, ib, sw, iw, T2b.cptr(), V2_sym.cptr(), W28av_sigma_covv_covv.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_covv_covv_no7_x1_type1_eri_c,G_IF_SIGMA_COVV_COVV_NO7_X1_TYPE1_ERI_C)
      (sb, ib, sw, iw, W28av_sigma_covv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // --  Title : sigma_covv_covv
  //  >> Intermediates for the external contractions are defined here << 
  orz::DTensor W10aa_sigma_covv_covv(orz::mr::sizeof_sympack_Xaa(symblockinfo, 0));
  orz::DTensor W13aa_sigma_covv_covv(orz::mr::sizeof_sympack_Xaa(symblockinfo, 0));
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
  // -- Title : sigma_covv_covv

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W10(i,a0) += (    1.00000000) V2(a4,a2,a3,a1) D3(i,a0,a4,a2,a3,a1) 
  int sa4(s_eri);
  int ia4(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  FC_FUNC(g_if_sigma_covv_covv_no0_x0_type1_eri_o,G_IF_SIGMA_COVV_COVV_NO0_X0_TYPE1_ERI_O)
    (sa4, ia4, V2_sym.cptr(), W10aa_sigma_covv_covv.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W13(i,a0) += (    1.00000000) V2(a4,a2,a3,a1) D3(i,a0,a4,a2,a3,a1) 
  int sa4(s_eri);
  int ia4(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  FC_FUNC(g_if_sigma_covv_covv_no1_x0_type1_eri_o,G_IF_SIGMA_COVV_COVV_NO1_X0_TYPE1_ERI_O)
    (sa4, ia4, V2_sym.cptr(), W13aa_sigma_covv_covv.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

//-@type(1).contraction(end)

  } // End myrank
  }
  }
  orz::world().barrier();

//-@type(2).contraction(begin)
  // -- Title : sigma_covv_covv
  //*-- Entering to take the type 2 contractions --*//

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| S2(w,i,a,b) += (    1.00000000) T2(w,a0,a,b) W10(i,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_covv_no0_x0_type2_eri_o,G_IF_SIGMA_COVV_COVV_NO0_X0_TYPE2_ERI_O)
      (sb, ib, T2b.cptr(), W10aa_sigma_covv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| S2(w,i,a,b) += (   -0.50000000) T2(a0,w,a,b) W13(i,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_covv_no1_x0_type2_eri_o,G_IF_SIGMA_COVV_COVV_NO1_X0_TYPE2_ERI_O)
      (sb, ib, T2b.cptr(), W13aa_sigma_covv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- Title : sigma_covv_covv

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W15(i,a0,v0,a) += (    1.00000000) V2(v0,a1,a2,a) D2(i,a0,a1,a2) 
  // |-- [    1] --| S2(w,i,a,b) += (   -1.00000000) T2(w,a0,v0,b) W15(i,a0,v0,a) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W15aav_sigma_covv_covv(orz::mr::sizeof_sympack_Xaav(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_covv_covv_no0_x0_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO0_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W15aav_sigma_covv_covv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_covv_no0_x1_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO0_X1_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), W15aav_sigma_covv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W16(w,a0,a,b) += (    1.00000000) V2(b,w,v0,c0) T2(c0,a0,v0,a) 
  // |-- [    1] --| S2(w,i,a,b) += (   -2.00000000) D1(i,a0) W16(w,a0,a,b) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W16ca_sigma_covv_covv(orz::mr::sizeof_sympack_Xca(symblockinfo, sa^sb));
    FC_FUNC(g_if_sigma_covv_covv_no1_x0_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO1_X0_TYPE1_ERI_V)
      (sa, ia, sb, ib, T2b.cptr(), V2_sym.cptr(), W16ca_sigma_covv_covv.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_covv_covv_no1_x1_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO1_X1_TYPE1_ERI_V)
      (sa, ia, sb, ib, W16ca_sigma_covv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
  } // End ia
  } // End sa
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W17(i,a0,v0,b) += (    1.00000000) V2(v0,a2,b,a1) D2(i,a1,a2,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (   -1.00000000) T2(a0,w,a,v0) W17(i,a0,v0,b) 
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
    orz::DTensor W17aa_sigma_covv_covv(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sb));
    FC_FUNC(g_if_sigma_covv_covv_no2_x0_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO2_X0_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, V2_sym.cptr(), W17aa_sigma_covv_covv.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_covv_covv_no2_x1_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO2_X1_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), W17aa_sigma_covv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W18(w,a0,a,b) += (    1.00000000) V2(b,v0,w,c0) T2(c0,a0,v0,a) 
  // |-- [    1] --| S2(w,i,a,b) += (    1.00000000) D1(i,a0) W18(w,a0,a,b) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W18ca_sigma_covv_covv(orz::mr::sizeof_sympack_Xca(symblockinfo, sa^sb));
    FC_FUNC(g_if_sigma_covv_covv_no3_x0_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO3_X0_TYPE1_ERI_V)
      (sa, ia, sb, ib, T2b.cptr(), V2_sym.cptr(), W18ca_sigma_covv_covv.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_covv_covv_no3_x1_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO3_X1_TYPE1_ERI_V)
      (sa, ia, sb, ib, W18ca_sigma_covv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
  } // End ia
  } // End sa
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W19(i,a0,v0,b) += (    1.00000000) V2(v0,b,a2,a1) D2(i,a0,a2,a1) 
  // |-- [    1] --| S2(w,i,a,b) += (   -1.00000000) T2(a0,w,a,v0) W19(i,a0,v0,b) 
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
    orz::DTensor W19aa_sigma_covv_covv(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sb));
    FC_FUNC(g_if_sigma_covv_covv_no4_x0_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO4_X0_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, V2_sym.cptr(), W19aa_sigma_covv_covv.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_covv_covv_no4_x1_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO4_X1_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), W19aa_sigma_covv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W21(i,a0,v0,a) += (    1.00000000) V2(v0,a,a2,a1) D2(i,a0,a2,a1) 
  // |-- [    1] --| S2(w,i,a,b) += (    2.00000000) T2(w,a0,v0,b) W21(i,a0,v0,a) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W21aav_sigma_covv_covv(orz::mr::sizeof_sympack_Xaav(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_covv_covv_no5_x0_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO5_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W21aav_sigma_covv_covv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_covv_no5_x1_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO5_X1_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), W21aav_sigma_covv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W23(i,a0,v0,a) += (    1.00000000) V2(v0,a1,a2,a) D2(i,a2,a1,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (   -1.00000000) T2(w,a0,b,v0) W23(i,a0,v0,a) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W23aav_sigma_covv_covv(orz::mr::sizeof_sympack_Xaav(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_covv_covv_no6_x0_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO6_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W23aav_sigma_covv_covv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_covv_no6_x1_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO6_X1_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), W23aav_sigma_covv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W24(w,a0,a,b) += (    1.00000000) V2(b,w,v0,c0) T2(c0,a0,a,v0) 
  // |-- [    1] --| S2(w,i,a,b) += (    1.00000000) D1(i,a0) W24(w,a0,a,b) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  orz::DTensor W24cav_sigma_covv_covv(orz::mr::sizeof_sympack_Xcav(symblockinfo, sb));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_covv_covv_no7_x0_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO7_X0_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W24cav_sigma_covv_covv.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_covv_covv_no7_x1_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO7_X1_TYPE1_ERI_V)
    (sb, ib, W24cav_sigma_covv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| W25(i,a0,v0,b) += (    1.00000000) V2(v0,a2,b,a1) D2(i,a1,a2,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (    2.00000000) T2(w,a0,a,v0) W25(i,a0,v0,b) 
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
    orz::DTensor W25aa_sigma_covv_covv(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sb));
    FC_FUNC(g_if_sigma_covv_covv_no8_x0_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO8_X0_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, V2_sym.cptr(), W25aa_sigma_covv_covv.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_covv_covv_no8_x1_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO8_X1_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), W25aa_sigma_covv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W26(w,a0,a,b) += (    1.00000000) V2(b,v0,w,c0) T2(c0,a0,a,v0) 
  // |-- [    1] --| S2(w,i,a,b) += (   -2.00000000) D1(i,a0) W26(w,a0,a,b) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  orz::DTensor W26cav_sigma_covv_covv(orz::mr::sizeof_sympack_Xcav(symblockinfo, sb));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_covv_covv_no9_x0_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO9_X0_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W26cav_sigma_covv_covv.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_covv_covv_no9_x1_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO9_X1_TYPE1_ERI_V)
    (sb, ib, W26cav_sigma_covv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| W27(i,a0,v0,b) += (    1.00000000) V2(v0,b,a2,a1) D2(i,a0,a2,a1) 
  // |-- [    1] --| S2(w,i,a,b) += (    2.00000000) T2(w,a0,a,v0) W27(i,a0,v0,b) 
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
    orz::DTensor W27aa_sigma_covv_covv(orz::mr::sizeof_sympack_Xaa(symblockinfo, sv0^sb));
    FC_FUNC(g_if_sigma_covv_covv_no10_x0_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO10_X0_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, V2_sym.cptr(), W27aa_sigma_covv_covv.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_covv_covv_no10_x1_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO10_X1_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), W27aa_sigma_covv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| W29(i,a0,v0,a) += (    1.00000000) V2(v0,a,a2,a1) D2(i,a0,a2,a1) 
  // |-- [    1] --| S2(w,i,a,b) += (   -1.00000000) T2(w,a0,b,v0) W29(i,a0,v0,a) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W29aav_sigma_covv_covv(orz::mr::sizeof_sympack_Xaav(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_covv_covv_no11_x0_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO11_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W29aav_sigma_covv_covv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_covv_no11_x1_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO11_X1_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), W29aav_sigma_covv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| W30(w,a0,b,a) += (    1.00000000) V2(b,v1,v0,a) T2(w,a0,v0,v1) 
  // |-- [    1] --| S2(w,i,a,b) += (    2.00000000) D1(i,a0) W30(w,a0,b,a) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  orz::DTensor W30cav_sigma_covv_covv(orz::mr::sizeof_sympack_Xcav(symblockinfo, sb));
  for(int sv1 = 0;sv1 < nir;++sv1){ 
  for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
    T2b = T2.get_amp2(iv1);
    FC_FUNC(g_if_sigma_covv_covv_no12_x0_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO12_X0_TYPE1_ERI_V)
      (sb, ib, sv1, iv1, T2b.cptr(), V2_sym.cptr(), W30cav_sigma_covv_covv.cptr(), nir, nsym, psym, &flops);
  // --> @[v1, "virtual"] [notNeeded]
  } // End iv1
  } // End sv1
  FC_FUNC(g_if_sigma_covv_covv_no12_x1_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO12_X1_TYPE1_ERI_V)
    (sb, ib, W30cav_sigma_covv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| W31(w,a0,b,a) += (    1.00000000) V2(b,v1,v0,a) T2(w,a0,v1,v0) 
  // |-- [    1] --| S2(w,i,a,b) += (   -1.00000000) D1(i,a0) W31(w,a0,b,a) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  orz::DTensor W31cav_sigma_covv_covv(orz::mr::sizeof_sympack_Xcav(symblockinfo, sb));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_covv_covv_no13_x0_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO13_X0_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W31cav_sigma_covv_covv.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_covv_covv_no13_x1_type1_eri_v,G_IF_SIGMA_COVV_COVV_NO13_X1_TYPE1_ERI_V)
    (sb, ib, W31cav_sigma_covv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ib, S2b);
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
  orz::mr::Femto::my_timer.push_back(boost::make_tuple("sigma_covv_covv", time_sigma.elapsed_cputime(), time_sigma.elapsed_wallclocktime()));
  orz::mr::Femto::file_timing << "* " << boost::format("%20s : %10.7f %10.7f ") % orz::mr::Femto::my_timer.back().get<0>() % orz::mr::Femto::my_timer.back().get<1>() % orz::mr::Femto::my_timer.back().get<2>() << endl;
  flush(orz::mr::Femto::file_timing);

  return retval; 
} 
