                                                                                
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
#include <sci/icmr/Femto/elems/c_sigma_covv_oovv.h>                                  
                                                                                
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
                                                                                
orz::mr::BareAmpPack orz::mr::Femto::sigma_covv_oovv(const orz::mr::Input &ctinp,                                    
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
  // |-- [    0] --| W0(i,a2,b,a) += (    1.00000000) T2(a0,a1,a,b) D2(i,a1,a2,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (   -2.00000000) Fc1(w,a2) W0(i,a2,b,a) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    orz::DTensor W0aav_sigma_covv_oovv(orz::mr::sizeof_sympack_Xaav(symblockinfo, sb));
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_oovv_no0_x0_type0_noeri,G_IF_SIGMA_COVV_OOVV_NO0_X0_TYPE0_NOERI)
      (sb, ib, T2b.cptr(), W0aav_sigma_covv_oovv.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_covv_oovv_no0_x1_type0_noeri,G_IF_SIGMA_COVV_OOVV_NO0_X1_TYPE0_NOERI)
      (sb, ib, W0aav_sigma_covv_oovv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- Title : sigma_covv_oovv

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W1(w,i,a1,a0) += (    1.00000000) V2(w,a3,a4,a2) D3(i,a1,a4,a2,a3,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (   -2.00000000) T2(a0,a1,a,b) W1(w,i,a1,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W1aaa_sigma_covv_oovv(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_covv_oovv_no0_x0_type1_eri_c,G_IF_SIGMA_COVV_OOVV_NO0_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W1aaa_sigma_covv_oovv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_oovv_no0_x1_type1_eri_c,G_IF_SIGMA_COVV_OOVV_NO0_X1_TYPE1_ERI_C)
      (sb, ib, sw, iw, T2b.cptr(), W1aaa_sigma_covv_oovv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- Title : sigma_covv_oovv

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W2(i,a2,b,v0) += (    1.00000000) T2(a1,a0,b,v0) D2(i,a1,a2,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (    4.00000000) V2(v0,a2,w,a) W2(i,a2,b,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    orz::DTensor W2aa_sigma_covv_oovv(orz::mr::sizeof_sympack_Xaa(symblockinfo, sb^sv0));
    FC_FUNC(g_if_sigma_covv_oovv_no0_x0_type1_eri_v,G_IF_SIGMA_COVV_OOVV_NO0_X0_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), W2aa_sigma_covv_oovv.cptr(), nir, nsym, psym, &flops);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_oovv_no0_x1_type1_eri_v,G_IF_SIGMA_COVV_OOVV_NO0_X1_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, V2_sym.cptr(), W2aa_sigma_covv_oovv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W3(i,a2,a,v0) += (    1.00000000) T2(a1,a0,a,v0) D2(i,a1,a2,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (   -2.00000000) V2(v0,a2,w,b) W3(i,a2,a,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W3aav_sigma_covv_oovv(orz::mr::sizeof_sympack_Xaav(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_covv_oovv_no1_x0_type1_eri_v,G_IF_SIGMA_COVV_OOVV_NO1_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W3aav_sigma_covv_oovv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_oovv_no1_x1_type1_eri_v,G_IF_SIGMA_COVV_OOVV_NO1_X1_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, V2_sym.cptr(), W3aav_sigma_covv_oovv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W4(i,a2,v0,a) += (    1.00000000) T2(a0,a1,a,v0) D2(i,a1,a2,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (   -2.00000000) V2(v0,b,w,a2) W4(i,a2,v0,a) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W4aav_sigma_covv_oovv(orz::mr::sizeof_sympack_Xaav(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_covv_oovv_no2_x0_type1_eri_v,G_IF_SIGMA_COVV_OOVV_NO2_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W4aav_sigma_covv_oovv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_oovv_no2_x1_type1_eri_v,G_IF_SIGMA_COVV_OOVV_NO2_X1_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, V2_sym.cptr(), W4aav_sigma_covv_oovv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W5(i,a2,b,v0) += (    1.00000000) T2(a1,a0,b,v0) D2(i,a1,a2,a0) 
  // |-- [    1] --| S2(w,i,a,b) += (   -2.00000000) V2(v0,a,w,a2) W5(i,a2,b,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    orz::DTensor W5aa_sigma_covv_oovv(orz::mr::sizeof_sympack_Xaa(symblockinfo, sb^sv0));
    FC_FUNC(g_if_sigma_covv_oovv_no3_x0_type1_eri_v,G_IF_SIGMA_COVV_OOVV_NO3_X0_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), W5aa_sigma_covv_oovv.cptr(), nir, nsym, psym, &flops);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_covv_oovv_no3_x1_type1_eri_v,G_IF_SIGMA_COVV_OOVV_NO3_X1_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, V2_sym.cptr(), W5aa_sigma_covv_oovv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  orz::mr::Femto::my_timer.push_back(boost::make_tuple("sigma_covv_oovv", time_sigma.elapsed_cputime(), time_sigma.elapsed_wallclocktime()));
  orz::mr::Femto::file_timing << "* " << boost::format("%20s : %10.7f %10.7f ") % orz::mr::Femto::my_timer.back().get<0>() % orz::mr::Femto::my_timer.back().get<1>() % orz::mr::Femto::my_timer.back().get<2>() << endl;
  flush(orz::mr::Femto::file_timing);

  return retval; 
} 
