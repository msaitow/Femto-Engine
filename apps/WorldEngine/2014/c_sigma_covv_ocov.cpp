                                                                                
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
#include <sci/icmr/Femto/elems/c_sigma_covv_ocov.h>                                  
                                                                                
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
                                                                                
orz::mr::BareAmpPack orz::mr::Femto::sigma_covv_ocov(const orz::mr::Input &ctinp,                                    
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
  // |-- [    0] --| W0(w,i,a2,b) += (    1.00000000) T2(a0,w,a1,b) D2(i,a2,a1,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (   -1.00000000) Fc1(a2,a) W0(w,i,a2,b) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    orz::DTensor W0caa_sigma_covv_ocov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sb));
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_ocov_no0_x0_type0_noeri,G_IF_SIGMA_COVV_OCOV_NO0_X0_TYPE0_NOERI)
      (sb, ib, T2b.cptr(), W0caa_sigma_covv_ocov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_covv_ocov_no0_x1_type0_noeri,G_IF_SIGMA_COVV_OCOV_NO0_X1_TYPE0_NOERI)
      (sb, ib, W0caa_sigma_covv_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W1(w,i,a2,a) += (    1.00000000) T2(w,a0,a,a1) D2(i,a2,a1,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (    2.00000000) Fc1(b,a2) W1(w,i,a2,a) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W1caav_sigma_covv_ocov(orz::mr::sizeof_sympack_Xcaav(symblockinfo, 0));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_covv_ocov_no1_x0_type0_noeri,G_IF_SIGMA_COVV_OCOV_NO1_X0_TYPE0_NOERI)
      (sa1, ia1, T2b.cptr(), W1caav_sigma_covv_ocov.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_ocov_no1_x1_type0_noeri,G_IF_SIGMA_COVV_OCOV_NO1_X1_TYPE0_NOERI)
      (sb, ib, W1caav_sigma_covv_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W2(w,i,a1,a) += (    1.00000000) T2(w,a0,a,a1) D1(i,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (    2.00000000) Fc1(b,a1) W2(w,i,a1,a) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W2cav_sigma_covv_ocov(orz::mr::sizeof_sympack_Xcav(symblockinfo, sa1));
    FC_FUNC(g_if_sigma_covv_ocov_no2_x0_type0_noeri,G_IF_SIGMA_COVV_OCOV_NO2_X0_TYPE0_NOERI)
      (sa1, ia1, T2b.cptr(), W2cav_sigma_covv_ocov.cptr(), nir, nsym, psym, &flops);
    for(int sb = 0;sb < nir;++sb){ 
    for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
      S2b = orz::DTensor(retval.namps_iamp()[ib]);
      FC_FUNC(g_if_sigma_covv_ocov_no2_x1_type0_noeri,G_IF_SIGMA_COVV_OCOV_NO2_X1_TYPE0_NOERI)
        (sa1, ia1, sb, ib, W2cav_sigma_covv_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
      retval.acc_amp2(ib, S2b);
    } // End ib
    } // End sb
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W3(w,i,a1,b) += (    1.00000000) T2(a0,w,a1,b) D1(i,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (   -1.00000000) Fc1(a1,a) W3(w,i,a1,b) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    orz::DTensor W3caa_sigma_covv_ocov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sb));
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_ocov_no3_x0_type0_noeri,G_IF_SIGMA_COVV_OCOV_NO3_X0_TYPE0_NOERI)
      (sb, ib, T2b.cptr(), W3caa_sigma_covv_ocov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_covv_ocov_no3_x1_type0_noeri,G_IF_SIGMA_COVV_OCOV_NO3_X1_TYPE0_NOERI)
      (sb, ib, W3caa_sigma_covv_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- Title : sigma_covv_ocov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W4(c0,i,a2,b) += (    1.00000000) T2(a0,c0,a1,b) D2(i,a2,a1,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (   -2.00000000) V2(c0,a2,w,a) W4(c0,i,a2,b) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    orz::DTensor W4aa_sigma_covv_ocov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sb));
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_ocov_no0_x0_type1_eri_c,G_IF_SIGMA_COVV_OCOV_NO0_X0_TYPE1_ERI_C)
      (sb, ib, sc0, ic0, T2b.cptr(), W4aa_sigma_covv_ocov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_covv_ocov_no0_x1_type1_eri_c,G_IF_SIGMA_COVV_OCOV_NO0_X1_TYPE1_ERI_C)
      (sb, ib, sc0, ic0, V2_sym.cptr(), W4aa_sigma_covv_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W5(c0,i,a2,b) += (    1.00000000) T2(a0,c0,a1,b) D2(i,a2,a1,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (    1.00000000) V2(c0,w,a2,a) W5(c0,i,a2,b) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    orz::DTensor W5aa_sigma_covv_ocov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sb));
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_ocov_no1_x0_type1_eri_c,G_IF_SIGMA_COVV_OCOV_NO1_X0_TYPE1_ERI_C)
      (sb, ib, sc0, ic0, T2b.cptr(), W5aa_sigma_covv_ocov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_covv_ocov_no1_x1_type1_eri_c,G_IF_SIGMA_COVV_OCOV_NO1_X1_TYPE1_ERI_C)
      (sb, ib, sc0, ic0, V2_sym.cptr(), W5aa_sigma_covv_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W10(w,a0,b,a) += (    1.00000000) V2(w,a,c0,a1) T2(a0,c0,a1,b) 
  // |-- [    1] --| S2(w,i,a,b) += (   -2.00000000) D1(i,a0) W10(w,a0,b,a) 
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
    orz::DTensor W10av_sigma_covv_ocov(orz::mr::sizeof_sympack_Xav(symblockinfo, sw^sb));
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_ocov_no2_x0_type1_eri_c,G_IF_SIGMA_COVV_OCOV_NO2_X0_TYPE1_ERI_C)
      (sb, ib, sw, iw, T2b.cptr(), V2_sym.cptr(), W10av_sigma_covv_ocov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_covv_ocov_no2_x1_type1_eri_c,G_IF_SIGMA_COVV_OCOV_NO2_X1_TYPE1_ERI_C)
      (sb, ib, sw, iw, W10av_sigma_covv_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W16(w,a0,b,a) += (    1.00000000) V2(w,c0,a1,a) T2(a0,c0,a1,b) 
  // |-- [    1] --| S2(w,i,a,b) += (    1.00000000) D1(i,a0) W16(w,a0,b,a) 
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
    orz::DTensor W16av_sigma_covv_ocov(orz::mr::sizeof_sympack_Xav(symblockinfo, sw^sb));
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_ocov_no3_x0_type1_eri_c,G_IF_SIGMA_COVV_OCOV_NO3_X0_TYPE1_ERI_C)
      (sb, ib, sw, iw, T2b.cptr(), V2_sym.cptr(), W16av_sigma_covv_ocov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_covv_ocov_no3_x1_type1_eri_c,G_IF_SIGMA_COVV_OCOV_NO3_X1_TYPE1_ERI_C)
      (sb, ib, sw, iw, W16av_sigma_covv_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- Title : sigma_covv_ocov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W11(i,a0,a1,a) += (    1.00000000) V2(a1,a3,a2,a) D2(i,a2,a3,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (   -1.00000000) T2(w,a0,b,a1) W11(i,a0,a1,a) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia1);
  orz::DTensor W11aav_sigma_covv_ocov(orz::mr::sizeof_sympack_Xaav(symblockinfo, sa1));
  FC_FUNC(g_if_sigma_covv_ocov_no0_x0_type1_eri_o,G_IF_SIGMA_COVV_OCOV_NO0_X0_TYPE1_ERI_O)
    (sa1, ia1, V2_sym.cptr(), W11aav_sigma_covv_ocov.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_ocov_no0_x1_type1_eri_o,G_IF_SIGMA_COVV_OCOV_NO0_X1_TYPE1_ERI_O)
      (sa1, ia1, sb, ib, T2b.cptr(), W11aav_sigma_covv_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W17(i,a0,a2,a) += (    1.00000000) V2(a2,a,a3,a1) D2(i,a0,a3,a1) 
  // |-- [    1] --| S2(w,i,a,b) += (   -1.00000000) T2(w,a0,b,a2) W17(i,a0,a2,a) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia2);
  orz::DTensor W17aav_sigma_covv_ocov(orz::mr::sizeof_sympack_Xaav(symblockinfo, sa2));
  FC_FUNC(g_if_sigma_covv_ocov_no1_x0_type1_eri_o,G_IF_SIGMA_COVV_OCOV_NO1_X0_TYPE1_ERI_O)
    (sa2, ia2, V2_sym.cptr(), W17aav_sigma_covv_ocov.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_ocov_no1_x1_type1_eri_o,G_IF_SIGMA_COVV_OCOV_NO1_X1_TYPE1_ERI_O)
      (sa2, ia2, sb, ib, T2b.cptr(), W17aav_sigma_covv_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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

//-@type(2).declaration(begin)
  // --  Title : sigma_covv_ocov
  //  >> Intermediates for the external contractions are defined here << 
  orz::DTensor W7caav_sigma_covv_ocov(orz::mr::sizeof_sympack_Xcaav(symblockinfo, 0));
  orz::DTensor W8caav_sigma_covv_ocov(orz::mr::sizeof_sympack_Xcaav(symblockinfo, 0));
//-@type(2).declaration(end)

  //*-- Entering to take the type 0 contractions --*//
//-@type(0).contraction(begin)

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W7(c0,i,a2,a) += (    1.00000000) T2(c0,a0,a,a1) D2(i,a2,a1,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_covv_ocov_no0_x0_type0_eri_v,G_IF_SIGMA_COVV_OCOV_NO0_X0_TYPE0_ERI_V)
      (sa1, ia1, T2b.cptr(), W7caav_sigma_covv_ocov.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W8(c0,i,a2,a) += (    1.00000000) T2(c0,a0,a,a1) D2(i,a2,a1,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_covv_ocov_no1_x0_type0_eri_v,G_IF_SIGMA_COVV_OCOV_NO1_X0_TYPE0_ERI_V)
      (sa1, ia1, T2b.cptr(), W8caav_sigma_covv_ocov.cptr(), nir, nsym, psym, &flops);
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
  // -- Title : sigma_covv_ocov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W6(i,a1,a0,a) += (    1.00000000) V2(a,a3,a4,a2) D3(i,a3,a4,a2,a1,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (   -1.00000000) T2(a0,w,a1,b) W6(i,a1,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W6aaa_sigma_covv_ocov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_covv_ocov_no0_x0_type1_eri_v,G_IF_SIGMA_COVV_OCOV_NO0_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W6aaa_sigma_covv_ocov.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_ocov_no0_x1_type1_eri_v,G_IF_SIGMA_COVV_OCOV_NO0_X1_TYPE1_ERI_V)
      (sa, ia, sb, ib, T2b.cptr(), W6aaa_sigma_covv_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| S2(w,i,a,b) += (    1.00000000) V2(b,w,c0,a2) W7(c0,i,a2,a) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  FC_FUNC(g_if_sigma_covv_ocov_no1_x0_type1_eri_v,G_IF_SIGMA_COVV_OCOV_NO1_X0_TYPE1_ERI_V)
    (sb, ib, V2_sym.cptr(), W7caav_sigma_covv_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| S2(w,i,a,b) += (   -2.00000000) V2(b,a2,w,c0) W8(c0,i,a2,a) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  FC_FUNC(g_if_sigma_covv_ocov_no2_x0_type1_eri_v,G_IF_SIGMA_COVV_OCOV_NO2_X0_TYPE1_ERI_V)
    (sb, ib, V2_sym.cptr(), W8caav_sigma_covv_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W9(i,a1,a0,b) += (    1.00000000) V2(b,a3,a4,a2) D3(i,a3,a4,a2,a1,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (    2.00000000) T2(w,a0,a,a1) W9(i,a1,a0,b) 
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
    orz::DTensor W9aa_sigma_covv_ocov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa1^sb));
    FC_FUNC(g_if_sigma_covv_ocov_no3_x0_type1_eri_v,G_IF_SIGMA_COVV_OCOV_NO3_X0_TYPE1_ERI_V)
      (sa1, ia1, sb, ib, V2_sym.cptr(), W9aa_sigma_covv_ocov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_covv_ocov_no3_x1_type1_eri_v,G_IF_SIGMA_COVV_OCOV_NO3_X1_TYPE1_ERI_V)
      (sa1, ia1, sb, ib, T2b.cptr(), W9aa_sigma_covv_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W12(w,a0,a,b) += (    1.00000000) V2(b,w,c0,a1) T2(c0,a0,a,a1) 
  // |-- [    1] --| S2(w,i,a,b) += (    1.00000000) D1(i,a0) W12(w,a0,a,b) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  orz::DTensor W12cav_sigma_covv_ocov(orz::mr::sizeof_sympack_Xcav(symblockinfo, sb));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_covv_ocov_no4_x0_type1_eri_v,G_IF_SIGMA_COVV_OCOV_NO4_X0_TYPE1_ERI_V)
      (sa1, ia1, sb, ib, T2b.cptr(), V2_sym.cptr(), W12cav_sigma_covv_ocov.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  FC_FUNC(g_if_sigma_covv_ocov_no4_x1_type1_eri_v,G_IF_SIGMA_COVV_OCOV_NO4_X1_TYPE1_ERI_V)
    (sb, ib, W12cav_sigma_covv_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W13(i,a0,a1,b) += (    1.00000000) V2(b,a2,a3,a1) D2(i,a2,a3,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (    2.00000000) T2(w,a0,a,a1) W13(i,a0,a1,b) 
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
    orz::DTensor W13aa_sigma_covv_ocov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa1^sb));
    FC_FUNC(g_if_sigma_covv_ocov_no5_x0_type1_eri_v,G_IF_SIGMA_COVV_OCOV_NO5_X0_TYPE1_ERI_V)
      (sa1, ia1, sb, ib, V2_sym.cptr(), W13aa_sigma_covv_ocov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_covv_ocov_no5_x1_type1_eri_v,G_IF_SIGMA_COVV_OCOV_NO5_X1_TYPE1_ERI_V)
      (sa1, ia1, sb, ib, T2b.cptr(), W13aa_sigma_covv_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W14(w,a0,a,b) += (    1.00000000) V2(b,a1,w,c0) T2(c0,a0,a,a1) 
  // |-- [    1] --| S2(w,i,a,b) += (   -2.00000000) D1(i,a0) W14(w,a0,a,b) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  orz::DTensor W14cav_sigma_covv_ocov(orz::mr::sizeof_sympack_Xcav(symblockinfo, sb));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_covv_ocov_no6_x0_type1_eri_v,G_IF_SIGMA_COVV_OCOV_NO6_X0_TYPE1_ERI_V)
      (sa1, ia1, sb, ib, T2b.cptr(), V2_sym.cptr(), W14cav_sigma_covv_ocov.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  FC_FUNC(g_if_sigma_covv_ocov_no6_x1_type1_eri_v,G_IF_SIGMA_COVV_OCOV_NO6_X1_TYPE1_ERI_V)
    (sb, ib, W14cav_sigma_covv_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W15(i,a0,a2,b) += (    1.00000000) V2(b,a2,a3,a1) D2(i,a0,a3,a1) 
  // |-- [    1] --| S2(w,i,a,b) += (    2.00000000) T2(w,a0,a,a2) W15(i,a0,a2,b) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  for(int sa2 = 0;sa2 < nir;++sa2){ 
  for(int ia2 = symblockinfo.psym()(sa2,I_O,I_BEGIN);ia2 <= symblockinfo.psym()(sa2,I_O,I_END);++ia2){ 
    T2b = T2.get_amp2(ia2);
    orz::DTensor W15aa_sigma_covv_ocov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa2^sb));
    FC_FUNC(g_if_sigma_covv_ocov_no7_x0_type1_eri_v,G_IF_SIGMA_COVV_OCOV_NO7_X0_TYPE1_ERI_V)
      (sa2, ia2, sb, ib, V2_sym.cptr(), W15aa_sigma_covv_ocov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_covv_ocov_no7_x1_type1_eri_v,G_IF_SIGMA_COVV_OCOV_NO7_X1_TYPE1_ERI_V)
      (sa2, ia2, sb, ib, T2b.cptr(), W15aa_sigma_covv_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a2, "active"] [notNeeded]
  } // End ia2
  } // End sa2
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| W18(w,i,a2,v0) += (    1.00000000) T2(a0,w,a1,v0) D2(i,a2,a1,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (   -1.00000000) V2(v0,b,a2,a) W18(w,i,a2,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W18caa_sigma_covv_ocov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_covv_ocov_no8_x0_type1_eri_v,G_IF_SIGMA_COVV_OCOV_NO8_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W18caa_sigma_covv_ocov.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_ocov_no8_x1_type1_eri_v,G_IF_SIGMA_COVV_OCOV_NO8_X1_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, V2_sym.cptr(), W18caa_sigma_covv_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W19(w,i,a2,v0) += (    1.00000000) T2(a0,w,a1,v0) D2(i,a2,a1,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (    2.00000000) V2(v0,a,b,a2) W19(w,i,a2,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W19caa_sigma_covv_ocov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_covv_ocov_no9_x0_type1_eri_v,G_IF_SIGMA_COVV_OCOV_NO9_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W19caa_sigma_covv_ocov.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_ocov_no9_x1_type1_eri_v,G_IF_SIGMA_COVV_OCOV_NO9_X1_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, V2_sym.cptr(), W19caa_sigma_covv_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| W20(w,a0,b,a) += (    1.00000000) V2(b,v0,a1,a) T2(w,a0,v0,a1) 
  // |-- [    1] --| S2(w,i,a,b) += (   -1.00000000) D1(i,a0) W20(w,a0,b,a) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  orz::DTensor W20cav_sigma_covv_ocov(orz::mr::sizeof_sympack_Xcav(symblockinfo, sb));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_covv_ocov_no10_x0_type1_eri_v,G_IF_SIGMA_COVV_OCOV_NO10_X0_TYPE1_ERI_V)
      (sa1, ia1, sb, ib, T2b.cptr(), V2_sym.cptr(), W20cav_sigma_covv_ocov.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  FC_FUNC(g_if_sigma_covv_ocov_no10_x1_type1_eri_v,G_IF_SIGMA_COVV_OCOV_NO10_X1_TYPE1_ERI_V)
    (sb, ib, W20cav_sigma_covv_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| W21(w,a0,b,a) += (    1.00000000) V2(b,a1,v0,a) T2(w,a0,v0,a1) 
  // |-- [    1] --| S2(w,i,a,b) += (    2.00000000) D1(i,a0) W21(w,a0,b,a) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  orz::DTensor W21cav_sigma_covv_ocov(orz::mr::sizeof_sympack_Xcav(symblockinfo, sb));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_covv_ocov_no11_x0_type1_eri_v,G_IF_SIGMA_COVV_OCOV_NO11_X0_TYPE1_ERI_V)
      (sa1, ia1, sb, ib, T2b.cptr(), V2_sym.cptr(), W21cav_sigma_covv_ocov.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  FC_FUNC(g_if_sigma_covv_ocov_no11_x1_type1_eri_v,G_IF_SIGMA_COVV_OCOV_NO11_X1_TYPE1_ERI_V)
    (sb, ib, W21cav_sigma_covv_ocov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  orz::mr::Femto::my_timer.push_back(boost::make_tuple("sigma_covv_ocov", time_sigma.elapsed_cputime(), time_sigma.elapsed_wallclocktime()));
  orz::mr::Femto::file_timing << "* " << boost::format("%20s : %10.7f %10.7f ") % orz::mr::Femto::my_timer.back().get<0>() % orz::mr::Femto::my_timer.back().get<1>() % orz::mr::Femto::my_timer.back().get<2>() << endl;
  flush(orz::mr::Femto::file_timing);

  return retval; 
} 
