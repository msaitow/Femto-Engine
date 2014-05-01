                                                                                
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
#include <sci/icmr/Femto/elems/c_sigma_covv_ccvv.h>                                  
                                                                                
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
                                                                                
orz::mr::BareAmpPack orz::mr::Femto::sigma_covv_ccvv(const orz::mr::Input &ctinp,                                    
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
  // |-- [    0] --| W0(c0,i) += (    1.00000000) D1(i,a0) Fc1(c0,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (    2.00000000) T2(c0,w,a,b) W0(c0,i) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W0ca_sigma_covv_ccvv(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  FC_FUNC(g_if_sigma_covv_ccvv_no0_x0_type0_noeri,G_IF_SIGMA_COVV_CCVV_NO0_X0_TYPE0_NOERI)
    (W0ca_sigma_covv_ccvv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_ccvv_no0_x1_type0_noeri,G_IF_SIGMA_COVV_CCVV_NO0_X1_TYPE0_NOERI)
      (sb, ib, T2b.cptr(), W0ca_sigma_covv_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W1(c0,i) += (    1.00000000) D1(i,a0) Fc1(c0,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (   -4.00000000) T2(w,c0,a,b) W1(c0,i) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W1ca_sigma_covv_ccvv(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  FC_FUNC(g_if_sigma_covv_ccvv_no1_x0_type0_noeri,G_IF_SIGMA_COVV_CCVV_NO1_X0_TYPE0_NOERI)
    (W1ca_sigma_covv_ccvv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_ccvv_no1_x1_type0_noeri,G_IF_SIGMA_COVV_CCVV_NO1_X1_TYPE0_NOERI)
      (sb, ib, T2b.cptr(), W1ca_sigma_covv_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- Title : sigma_covv_ccvv

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W2(w,c1,c0,i) += (    1.00000000) V2(w,c0,c1,a0) D1(i,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (   -2.00000000) T2(c1,c0,a,b) W2(w,c1,c0,i) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W2cca_sigma_covv_ccvv(orz::mr::sizeof_sympack_Xcca(symblockinfo, sw));
  FC_FUNC(g_if_sigma_covv_ccvv_no0_x0_type1_eri_c,G_IF_SIGMA_COVV_CCVV_NO0_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W2cca_sigma_covv_ccvv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_ccvv_no0_x1_type1_eri_c,G_IF_SIGMA_COVV_CCVV_NO0_X1_TYPE1_ERI_C)
      (sb, ib, sw, iw, T2b.cptr(), W2cca_sigma_covv_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W3(w,c1,c0,i) += (    1.00000000) V2(w,c0,c1,a0) D1(i,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (    4.00000000) T2(c0,c1,a,b) W3(w,c1,c0,i) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W3cca_sigma_covv_ccvv(orz::mr::sizeof_sympack_Xcca(symblockinfo, sw));
  FC_FUNC(g_if_sigma_covv_ccvv_no1_x0_type1_eri_c,G_IF_SIGMA_COVV_CCVV_NO1_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W3cca_sigma_covv_ccvv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_ccvv_no1_x1_type1_eri_c,G_IF_SIGMA_COVV_CCVV_NO1_X1_TYPE1_ERI_C)
      (sb, ib, sw, iw, T2b.cptr(), W3cca_sigma_covv_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W4(c0,i) += (    1.00000000) V2(c0,a1,a2,a0) D2(i,a1,a2,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (   -4.00000000) T2(w,c0,a,b) W4(c0,i) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W4a_sigma_covv_ccvv(orz::mr::sizeof_sympack_Xa(symblockinfo, sc0));
  FC_FUNC(g_if_sigma_covv_ccvv_no2_x0_type1_eri_c,G_IF_SIGMA_COVV_CCVV_NO2_X0_TYPE1_ERI_C)
    (sc0, ic0, V2_sym.cptr(), W4a_sigma_covv_ccvv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_ccvv_no2_x1_type1_eri_c,G_IF_SIGMA_COVV_CCVV_NO2_X1_TYPE1_ERI_C)
      (sb, ib, sc0, ic0, T2b.cptr(), W4a_sigma_covv_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W5(c0,i) += (    1.00000000) V2(c0,a1,a2,a0) D2(i,a1,a2,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (    2.00000000) T2(c0,w,a,b) W5(c0,i) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W5a_sigma_covv_ccvv(orz::mr::sizeof_sympack_Xa(symblockinfo, sc0));
  FC_FUNC(g_if_sigma_covv_ccvv_no3_x0_type1_eri_c,G_IF_SIGMA_COVV_CCVV_NO3_X0_TYPE1_ERI_C)
    (sc0, ic0, V2_sym.cptr(), W5a_sigma_covv_ccvv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_ccvv_no3_x1_type1_eri_c,G_IF_SIGMA_COVV_CCVV_NO3_X1_TYPE1_ERI_C)
      (sb, ib, sc0, ic0, T2b.cptr(), W5a_sigma_covv_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- Title : sigma_covv_ccvv

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W6(c0,i,v0,a) += (    1.00000000) V2(v0,c0,a0,a) D1(i,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (   -4.00000000) T2(w,c0,b,v0) W6(c0,i,v0,a) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W6cav_sigma_covv_ccvv(orz::mr::sizeof_sympack_Xcav(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_covv_ccvv_no0_x0_type1_eri_v,G_IF_SIGMA_COVV_CCVV_NO0_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W6cav_sigma_covv_ccvv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_ccvv_no0_x1_type1_eri_v,G_IF_SIGMA_COVV_CCVV_NO0_X1_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), W6cav_sigma_covv_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W7(c0,i,v0,a) += (    1.00000000) V2(v0,c0,a0,a) D1(i,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (    2.00000000) T2(w,c0,v0,b) W7(c0,i,v0,a) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W7cav_sigma_covv_ccvv(orz::mr::sizeof_sympack_Xcav(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_covv_ccvv_no1_x0_type1_eri_v,G_IF_SIGMA_COVV_CCVV_NO1_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W7cav_sigma_covv_ccvv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_ccvv_no1_x1_type1_eri_v,G_IF_SIGMA_COVV_CCVV_NO1_X1_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), W7cav_sigma_covv_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W8(c0,i,v0,b) += (    1.00000000) V2(v0,c0,b,a0) D1(i,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (    8.00000000) T2(w,c0,a,v0) W8(c0,i,v0,b) 
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
    orz::DTensor W8ca_sigma_covv_ccvv(orz::mr::sizeof_sympack_Xca(symblockinfo, sv0^sb));
    FC_FUNC(g_if_sigma_covv_ccvv_no2_x0_type1_eri_v,G_IF_SIGMA_COVV_CCVV_NO2_X0_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, V2_sym.cptr(), W8ca_sigma_covv_ccvv.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_covv_ccvv_no2_x1_type1_eri_v,G_IF_SIGMA_COVV_CCVV_NO2_X1_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), W8ca_sigma_covv_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W9(c0,i,v0,b) += (    1.00000000) V2(v0,c0,b,a0) D1(i,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (   -4.00000000) T2(c0,w,a,v0) W9(c0,i,v0,b) 
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
    orz::DTensor W9ca_sigma_covv_ccvv(orz::mr::sizeof_sympack_Xca(symblockinfo, sv0^sb));
    FC_FUNC(g_if_sigma_covv_ccvv_no3_x0_type1_eri_v,G_IF_SIGMA_COVV_CCVV_NO3_X0_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, V2_sym.cptr(), W9ca_sigma_covv_ccvv.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_covv_ccvv_no3_x1_type1_eri_v,G_IF_SIGMA_COVV_CCVV_NO3_X1_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), W9ca_sigma_covv_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W10(c0,i,v0,b) += (    1.00000000) V2(v0,b,c0,a0) D1(i,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (   -4.00000000) T2(w,c0,a,v0) W10(c0,i,v0,b) 
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
    orz::DTensor W10ca_sigma_covv_ccvv(orz::mr::sizeof_sympack_Xca(symblockinfo, sv0^sb));
    FC_FUNC(g_if_sigma_covv_ccvv_no4_x0_type1_eri_v,G_IF_SIGMA_COVV_CCVV_NO4_X0_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, V2_sym.cptr(), W10ca_sigma_covv_ccvv.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_covv_ccvv_no4_x1_type1_eri_v,G_IF_SIGMA_COVV_CCVV_NO4_X1_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), W10ca_sigma_covv_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W11(c0,i,v0,b) += (    1.00000000) V2(v0,b,c0,a0) D1(i,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (    2.00000000) T2(c0,w,a,v0) W11(c0,i,v0,b) 
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
    orz::DTensor W11ca_sigma_covv_ccvv(orz::mr::sizeof_sympack_Xca(symblockinfo, sv0^sb));
    FC_FUNC(g_if_sigma_covv_ccvv_no5_x0_type1_eri_v,G_IF_SIGMA_COVV_CCVV_NO5_X0_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, V2_sym.cptr(), W11ca_sigma_covv_ccvv.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_covv_ccvv_no5_x1_type1_eri_v,G_IF_SIGMA_COVV_CCVV_NO5_X1_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), W11ca_sigma_covv_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W12(c0,i,v0,a) += (    1.00000000) V2(v0,a,c0,a0) D1(i,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (   -4.00000000) T2(w,c0,v0,b) W12(c0,i,v0,a) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W12cav_sigma_covv_ccvv(orz::mr::sizeof_sympack_Xcav(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_covv_ccvv_no6_x0_type1_eri_v,G_IF_SIGMA_COVV_CCVV_NO6_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W12cav_sigma_covv_ccvv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_ccvv_no6_x1_type1_eri_v,G_IF_SIGMA_COVV_CCVV_NO6_X1_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), W12cav_sigma_covv_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W13(c0,i,v0,a) += (    1.00000000) V2(v0,a,c0,a0) D1(i,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (    2.00000000) T2(w,c0,b,v0) W13(c0,i,v0,a) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W13cav_sigma_covv_ccvv(orz::mr::sizeof_sympack_Xcav(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_covv_ccvv_no7_x0_type1_eri_v,G_IF_SIGMA_COVV_CCVV_NO7_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W13cav_sigma_covv_ccvv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_ccvv_no7_x1_type1_eri_v,G_IF_SIGMA_COVV_CCVV_NO7_X1_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), W13cav_sigma_covv_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  orz::mr::Femto::my_timer.push_back(boost::make_tuple("sigma_covv_ccvv", time_sigma.elapsed_cputime(), time_sigma.elapsed_wallclocktime()));
  orz::mr::Femto::file_timing << "* " << boost::format("%20s : %10.7f %10.7f ") % orz::mr::Femto::my_timer.back().get<0>() % orz::mr::Femto::my_timer.back().get<1>() % orz::mr::Femto::my_timer.back().get<2>() << endl;
  flush(orz::mr::Femto::file_timing);

  return retval; 
} 
