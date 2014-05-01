                                                                                
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
#include <sci/icmr/Femto/elems/c_sigma_covv_ooov.h>                                  
                                                                                
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
                                                                                
orz::mr::BareAmpPack orz::mr::Femto::sigma_covv_ooov(const orz::mr::Input &ctinp,                                    
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
  // |-- [    0] --| W0(i,b) += (    1.00000000) T2(a0,a1,a2,b) D2(i,a1,a2,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (    2.00000000) Fc1(w,a) W0(i,b) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    orz::DTensor W0a_sigma_covv_ooov(orz::mr::sizeof_sympack_Xa(symblockinfo, sb));
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_ooov_no0_x0_type0_noeri,G_IF_SIGMA_COVV_OOOV_NO0_X0_TYPE0_NOERI)
      (sb, ib, T2b.cptr(), W0a_sigma_covv_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_covv_ooov_no0_x1_type0_noeri,G_IF_SIGMA_COVV_OOOV_NO0_X1_TYPE0_NOERI)
      (sb, ib, W0a_sigma_covv_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W1(i,a) += (    1.00000000) T2(a1,a0,a,a2) D2(i,a1,a2,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (   -1.00000000) Fc1(w,b) W1(i,a) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W1av_sigma_covv_ooov(orz::mr::sizeof_sympack_Xav(symblockinfo, 0));
  for(int sa2 = 0;sa2 < nir;++sa2){ 
  for(int ia2 = symblockinfo.psym()(sa2,I_O,I_BEGIN);ia2 <= symblockinfo.psym()(sa2,I_O,I_END);++ia2){ 
    T2b = T2.get_amp2(ia2);
    FC_FUNC(g_if_sigma_covv_ooov_no1_x0_type0_noeri,G_IF_SIGMA_COVV_OOOV_NO1_X0_TYPE0_NOERI)
      (sa2, ia2, T2b.cptr(), W1av_sigma_covv_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a2, "active"] [notNeeded]
  } // End ia2
  } // End sa2
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_ooov_no1_x1_type0_noeri,G_IF_SIGMA_COVV_OOOV_NO1_X1_TYPE0_NOERI)
      (sb, ib, W1av_sigma_covv_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
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
  // -- Title : sigma_covv_ooov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W2(i,a4,a3,b) += (    1.00000000) T2(a0,a1,a2,b) D3(i,a1,a4,a3,a2,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (    2.00000000) V2(a4,a3,w,a) W2(i,a4,a3,b) 
  int sa4(s_eri);
  int ia4(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    orz::DTensor W2aa_sigma_covv_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa4^sb));
    FC_FUNC(g_if_sigma_covv_ooov_no0_x0_type1_eri_o,G_IF_SIGMA_COVV_OOOV_NO0_X0_TYPE1_ERI_O)
      (sa4, ia4, sb, ib, T2b.cptr(), W2aa_sigma_covv_ooov.cptr(), nir, nsym, psym, &flops);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_ooov_no0_x1_type1_eri_o,G_IF_SIGMA_COVV_OOOV_NO0_X1_TYPE1_ERI_O)
      (sa4, ia4, sb, ib, V2_sym.cptr(), W2aa_sigma_covv_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W3(i,a3,a4,b) += (    1.00000000) T2(a0,a1,a2,b) D3(i,a1,a3,a4,a2,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (   -1.00000000) V2(a3,w,a4,a) W3(i,a3,a4,b) 
  int sa3(s_eri);
  int ia3(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    orz::DTensor W3aa_sigma_covv_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa3^sb));
    FC_FUNC(g_if_sigma_covv_ooov_no1_x0_type1_eri_o,G_IF_SIGMA_COVV_OOOV_NO1_X0_TYPE1_ERI_O)
      (sa3, ia3, sb, ib, T2b.cptr(), W3aa_sigma_covv_ooov.cptr(), nir, nsym, psym, &flops);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_ooov_no1_x1_type1_eri_o,G_IF_SIGMA_COVV_OOOV_NO1_X1_TYPE1_ERI_O)
      (sa3, ia3, sb, ib, V2_sym.cptr(), W3aa_sigma_covv_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W6(i,a3,a2,b) += (    1.00000000) T2(a1,a0,b,a2) D2(i,a1,a3,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (    2.00000000) V2(a2,a3,w,a) W6(i,a3,a2,b) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  T2b = T2.get_amp2(ia2);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    orz::DTensor W6aa_sigma_covv_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa2^sb));
    FC_FUNC(g_if_sigma_covv_ooov_no2_x0_type1_eri_o,G_IF_SIGMA_COVV_OOOV_NO2_X0_TYPE1_ERI_O)
      (sa2, ia2, sb, ib, T2b.cptr(), W6aa_sigma_covv_ooov.cptr(), nir, nsym, psym, &flops);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_ooov_no2_x1_type1_eri_o,G_IF_SIGMA_COVV_OOOV_NO2_X1_TYPE1_ERI_O)
      (sa2, ia2, sb, ib, V2_sym.cptr(), W6aa_sigma_covv_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W7(i,a3,a2,a) += (    1.00000000) T2(a1,a0,a,a2) D2(i,a1,a3,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (   -1.00000000) V2(a2,a3,w,b) W7(i,a3,a2,a) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia2);
  orz::DTensor W7aav_sigma_covv_ooov(orz::mr::sizeof_sympack_Xaav(symblockinfo, sa2));
  FC_FUNC(g_if_sigma_covv_ooov_no3_x0_type1_eri_o,G_IF_SIGMA_COVV_OOOV_NO3_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), W7aav_sigma_covv_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_ooov_no3_x1_type1_eri_o,G_IF_SIGMA_COVV_OOOV_NO3_X1_TYPE1_ERI_O)
      (sa2, ia2, sb, ib, V2_sym.cptr(), W7aav_sigma_covv_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W8(i,a3,a2,a) += (    1.00000000) T2(a0,a1,a,a2) D2(i,a1,a3,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (   -1.00000000) V2(a2,b,w,a3) W8(i,a3,a2,a) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia2);
  orz::DTensor W8aav_sigma_covv_ooov(orz::mr::sizeof_sympack_Xaav(symblockinfo, sa2));
  FC_FUNC(g_if_sigma_covv_ooov_no4_x0_type1_eri_o,G_IF_SIGMA_COVV_OOOV_NO4_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), W8aav_sigma_covv_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_ooov_no4_x1_type1_eri_o,G_IF_SIGMA_COVV_OOOV_NO4_X1_TYPE1_ERI_O)
      (sa2, ia2, sb, ib, V2_sym.cptr(), W8aav_sigma_covv_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W9(i,a2,a3,b) += (    1.00000000) T2(a1,a0,b,a3) D2(i,a1,a2,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (   -1.00000000) V2(a3,a,w,a2) W9(i,a2,a3,b) 
  int sa3(s_eri);
  int ia3(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  T2b = T2.get_amp2(ia3);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    orz::DTensor W9aa_sigma_covv_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa3^sb));
    FC_FUNC(g_if_sigma_covv_ooov_no5_x0_type1_eri_o,G_IF_SIGMA_COVV_OOOV_NO5_X0_TYPE1_ERI_O)
      (sa3, ia3, sb, ib, T2b.cptr(), W9aa_sigma_covv_ooov.cptr(), nir, nsym, psym, &flops);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_ooov_no5_x1_type1_eri_o,G_IF_SIGMA_COVV_OOOV_NO5_X1_TYPE1_ERI_O)
      (sa3, ia3, sb, ib, V2_sym.cptr(), W9aa_sigma_covv_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // --  Title : sigma_covv_ooov
  //  >> Intermediates for the external contractions are defined here << 
  orz::DTensor W4aaav_sigma_covv_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
  orz::DTensor W5aaav_sigma_covv_ooov(orz::mr::sizeof_sympack_Xaaav(symblockinfo, 0));
//-@type(2).declaration(end)

  //*-- Entering to take the type 0 contractions --*//
//-@type(0).contraction(begin)

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W4(i,a4,a3,a) += (    1.00000000) T2(a1,a0,a,a2) D3(i,a1,a4,a3,a2,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa2 = 0;sa2 < nir;++sa2){ 
  for(int ia2 = symblockinfo.psym()(sa2,I_O,I_BEGIN);ia2 <= symblockinfo.psym()(sa2,I_O,I_END);++ia2){ 
    T2b = T2.get_amp2(ia2);
    FC_FUNC(g_if_sigma_covv_ooov_no0_x0_type0_eri_v,G_IF_SIGMA_COVV_OOOV_NO0_X0_TYPE0_ERI_V)
      (sa2, ia2, T2b.cptr(), W4aaav_sigma_covv_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a2, "active"] [notNeeded]
  } // End ia2
  } // End sa2
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W5(i,a4,a3,a) += (    1.00000000) T2(a1,a0,a,a2) D3(i,a3,a4,a1,a2,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa2 = 0;sa2 < nir;++sa2){ 
  for(int ia2 = symblockinfo.psym()(sa2,I_O,I_BEGIN);ia2 <= symblockinfo.psym()(sa2,I_O,I_END);++ia2){ 
    T2b = T2.get_amp2(ia2);
    FC_FUNC(g_if_sigma_covv_ooov_no1_x0_type0_eri_v,G_IF_SIGMA_COVV_OOOV_NO1_X0_TYPE0_ERI_V)
      (sa2, ia2, T2b.cptr(), W5aaav_sigma_covv_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a2, "active"] [notNeeded]
  } // End ia2
  } // End sa2
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
  // -- Title : sigma_covv_ooov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| S2(w,i,a,b) += (   -1.00000000) V2(b,w,a4,a3) W4(i,a4,a3,a) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  FC_FUNC(g_if_sigma_covv_ooov_no0_x0_type1_eri_v,G_IF_SIGMA_COVV_OOOV_NO0_X0_TYPE1_ERI_V)
    (sb, ib, V2_sym.cptr(), W4aaav_sigma_covv_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| S2(w,i,a,b) += (   -1.00000000) V2(b,a3,w,a4) W5(i,a4,a3,a) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  FC_FUNC(g_if_sigma_covv_ooov_no1_x0_type1_eri_v,G_IF_SIGMA_COVV_OOOV_NO1_X0_TYPE1_ERI_V)
    (sb, ib, V2_sym.cptr(), W5aaav_sigma_covv_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W10(i,v0) += (    1.00000000) T2(a0,a1,a2,v0) D2(i,a1,a2,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (    2.00000000) V2(v0,b,w,a) W10(i,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W10a_sigma_covv_ooov(orz::mr::sizeof_sympack_Xa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_covv_ooov_no2_x0_type1_eri_v,G_IF_SIGMA_COVV_OOOV_NO2_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W10a_sigma_covv_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_ooov_no2_x1_type1_eri_v,G_IF_SIGMA_COVV_OOOV_NO2_X1_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, V2_sym.cptr(), W10a_sigma_covv_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W11(i,v0) += (    1.00000000) T2(a0,a1,a2,v0) D2(i,a1,a2,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (   -1.00000000) V2(v0,a,w,b) W11(i,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W11a_sigma_covv_ooov(orz::mr::sizeof_sympack_Xa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_covv_ooov_no3_x0_type1_eri_v,G_IF_SIGMA_COVV_OOOV_NO3_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W11a_sigma_covv_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_ooov_no3_x1_type1_eri_v,G_IF_SIGMA_COVV_OOOV_NO3_X1_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, V2_sym.cptr(), W11a_sigma_covv_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  orz::mr::Femto::my_timer.push_back(boost::make_tuple("sigma_covv_ooov", time_sigma.elapsed_cputime(), time_sigma.elapsed_wallclocktime()));
  orz::mr::Femto::file_timing << "* " << boost::format("%20s : %10.7f %10.7f ") % orz::mr::Femto::my_timer.back().get<0>() % orz::mr::Femto::my_timer.back().get<1>() % orz::mr::Femto::my_timer.back().get<2>() << endl;
  flush(orz::mr::Femto::file_timing);

  return retval; 
} 
