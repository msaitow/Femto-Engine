                                                                                
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
#include <sci/icmr/Femto/elems/c_sigma_ccvv_covv.h>                                  
                                                                                
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
                                                                                
orz::mr::BareAmpPack orz::mr::Femto::sigma_ccvv_covv(const orz::mr::Input &ctinp,                                    
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
  // |-- [    0] --| W0(w,a0) += (    1.00000000) D1(a1,a0) Fc1(w,a1) 
  // |-- [    1] --| S2(w,x,a,b) += (    1.00000000) T2(x,a0,a,b) W0(w,a0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W0ca_sigma_ccvv_covv(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  FC_FUNC(g_if_sigma_ccvv_covv_no0_x0_type0_noeri,G_IF_SIGMA_CCVV_COVV_NO0_X0_TYPE0_NOERI)
    (W0ca_sigma_ccvv_covv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_covv_no0_x1_type0_noeri,G_IF_SIGMA_CCVV_COVV_NO0_X1_TYPE0_NOERI)
      (sb, ib, T2b.cptr(), W0ca_sigma_ccvv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W1(x,a0) += (    1.00000000) D1(a1,a0) Fc1(x,a1) 
  // |-- [    1] --| S2(w,x,a,b) += (   -2.00000000) T2(w,a0,a,b) W1(x,a0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W1ca_sigma_ccvv_covv(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  FC_FUNC(g_if_sigma_ccvv_covv_no1_x0_type0_noeri,G_IF_SIGMA_CCVV_COVV_NO1_X0_TYPE0_NOERI)
    (W1ca_sigma_ccvv_covv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_covv_no1_x1_type0_noeri,G_IF_SIGMA_CCVV_COVV_NO1_X1_TYPE0_NOERI)
      (sb, ib, T2b.cptr(), W1ca_sigma_ccvv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W2(w,a0) += (    1.00000000) D1(a1,a0) Fc1(w,a1) 
  // |-- [    1] --| S2(w,x,a,b) += (   -2.00000000) T2(a0,x,a,b) W2(w,a0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W2ca_sigma_ccvv_covv(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  FC_FUNC(g_if_sigma_ccvv_covv_no2_x0_type0_noeri,G_IF_SIGMA_CCVV_COVV_NO2_X0_TYPE0_NOERI)
    (W2ca_sigma_ccvv_covv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_covv_no2_x1_type0_noeri,G_IF_SIGMA_CCVV_COVV_NO2_X1_TYPE0_NOERI)
      (sb, ib, T2b.cptr(), W2ca_sigma_ccvv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W3(x,a0) += (    1.00000000) D1(a1,a0) Fc1(x,a1) 
  // |-- [    1] --| S2(w,x,a,b) += (    1.00000000) T2(a0,w,a,b) W3(x,a0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W3ca_sigma_ccvv_covv(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  FC_FUNC(g_if_sigma_ccvv_covv_no3_x0_type0_noeri,G_IF_SIGMA_CCVV_COVV_NO3_X0_TYPE0_NOERI)
    (W3ca_sigma_ccvv_covv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_covv_no3_x1_type0_noeri,G_IF_SIGMA_CCVV_COVV_NO3_X1_TYPE0_NOERI)
      (sb, ib, T2b.cptr(), W3ca_sigma_ccvv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- Title : sigma_ccvv_covv

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W4(x,w,c0,a0) += (    1.00000000) V2(x,c0,w,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (   -1.00000000) T2(c0,a0,a,b) W4(x,w,c0,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W4cca_sigma_ccvv_covv(orz::mr::sizeof_sympack_Xcca(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccvv_covv_no0_x0_type1_eri_c,G_IF_SIGMA_CCVV_COVV_NO0_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W4cca_sigma_ccvv_covv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_covv_no0_x1_type1_eri_c,G_IF_SIGMA_CCVV_COVV_NO0_X1_TYPE1_ERI_C)
      (sb, ib, sx, ix, T2b.cptr(), W4cca_sigma_ccvv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W5(x,w,c0,a0) += (    1.00000000) V2(x,a1,w,c0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (    2.00000000) T2(c0,a0,a,b) W5(x,w,c0,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W5cca_sigma_ccvv_covv(orz::mr::sizeof_sympack_Xcca(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccvv_covv_no1_x0_type1_eri_c,G_IF_SIGMA_CCVV_COVV_NO1_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W5cca_sigma_ccvv_covv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_covv_no1_x1_type1_eri_c,G_IF_SIGMA_CCVV_COVV_NO1_X1_TYPE1_ERI_C)
      (sb, ib, sx, ix, T2b.cptr(), W5cca_sigma_ccvv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W6(w,a0) += (    1.00000000) V2(w,a2,a3,a1) D2(a3,a1,a2,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (    1.00000000) T2(x,a0,a,b) W6(w,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W6a_sigma_ccvv_covv(orz::mr::sizeof_sympack_Xa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccvv_covv_no2_x0_type1_eri_c,G_IF_SIGMA_CCVV_COVV_NO2_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W6a_sigma_ccvv_covv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_covv_no2_x1_type1_eri_c,G_IF_SIGMA_CCVV_COVV_NO2_X1_TYPE1_ERI_C)
      (sb, ib, sw, iw, T2b.cptr(), W6a_sigma_ccvv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W7(x,a0) += (    1.00000000) V2(x,a2,a3,a1) D2(a3,a1,a2,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (   -2.00000000) T2(w,a0,a,b) W7(x,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W7a_sigma_ccvv_covv(orz::mr::sizeof_sympack_Xa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccvv_covv_no3_x0_type1_eri_c,G_IF_SIGMA_CCVV_COVV_NO3_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W7a_sigma_ccvv_covv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_covv_no3_x1_type1_eri_c,G_IF_SIGMA_CCVV_COVV_NO3_X1_TYPE1_ERI_C)
      (sb, ib, sx, ix, T2b.cptr(), W7a_sigma_ccvv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W8(x,w,c0,a0) += (    1.00000000) V2(x,c0,w,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (    2.00000000) T2(a0,c0,a,b) W8(x,w,c0,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W8cca_sigma_ccvv_covv(orz::mr::sizeof_sympack_Xcca(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccvv_covv_no4_x0_type1_eri_c,G_IF_SIGMA_CCVV_COVV_NO4_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W8cca_sigma_ccvv_covv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_covv_no4_x1_type1_eri_c,G_IF_SIGMA_CCVV_COVV_NO4_X1_TYPE1_ERI_C)
      (sb, ib, sx, ix, T2b.cptr(), W8cca_sigma_ccvv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W9(x,w,c0,a0) += (    1.00000000) V2(x,a1,w,c0) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (   -1.00000000) T2(a0,c0,a,b) W9(x,w,c0,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W9cca_sigma_ccvv_covv(orz::mr::sizeof_sympack_Xcca(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccvv_covv_no5_x0_type1_eri_c,G_IF_SIGMA_CCVV_COVV_NO5_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W9cca_sigma_ccvv_covv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_covv_no5_x1_type1_eri_c,G_IF_SIGMA_CCVV_COVV_NO5_X1_TYPE1_ERI_C)
      (sb, ib, sx, ix, T2b.cptr(), W9cca_sigma_ccvv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W10(w,a0) += (    1.00000000) V2(w,a2,a3,a1) D2(a3,a1,a2,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (   -2.00000000) T2(a0,x,a,b) W10(w,a0) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W10a_sigma_ccvv_covv(orz::mr::sizeof_sympack_Xa(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccvv_covv_no6_x0_type1_eri_c,G_IF_SIGMA_CCVV_COVV_NO6_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W10a_sigma_ccvv_covv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_covv_no6_x1_type1_eri_c,G_IF_SIGMA_CCVV_COVV_NO6_X1_TYPE1_ERI_C)
      (sb, ib, sw, iw, T2b.cptr(), W10a_sigma_ccvv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W11(x,a0) += (    1.00000000) V2(x,a2,a3,a1) D2(a3,a1,a2,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (    1.00000000) T2(a0,w,a,b) W11(x,a0) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W11a_sigma_ccvv_covv(orz::mr::sizeof_sympack_Xa(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccvv_covv_no7_x0_type1_eri_c,G_IF_SIGMA_CCVV_COVV_NO7_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W11a_sigma_ccvv_covv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_covv_no7_x1_type1_eri_c,G_IF_SIGMA_CCVV_COVV_NO7_X1_TYPE1_ERI_C)
      (sb, ib, sx, ix, T2b.cptr(), W11a_sigma_ccvv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| W12(x,a0,v0,a) += (    1.00000000) V2(x,a,v0,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (    1.00000000) T2(w,a0,v0,b) W12(x,a0,v0,a) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W12avv_sigma_ccvv_covv(orz::mr::sizeof_sympack_Xavv(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccvv_covv_no8_x0_type1_eri_c,G_IF_SIGMA_CCVV_COVV_NO8_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W12avv_sigma_ccvv_covv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_covv_no8_x1_type1_eri_c,G_IF_SIGMA_CCVV_COVV_NO8_X1_TYPE1_ERI_C)
      (sb, ib, sx, ix, T2b.cptr(), W12avv_sigma_ccvv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W13(w,a0,v0,a) += (    1.00000000) V2(w,a,v0,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (   -2.00000000) T2(x,a0,v0,b) W13(w,a0,v0,a) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W13avv_sigma_ccvv_covv(orz::mr::sizeof_sympack_Xavv(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccvv_covv_no9_x0_type1_eri_c,G_IF_SIGMA_CCVV_COVV_NO9_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W13avv_sigma_ccvv_covv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_covv_no9_x1_type1_eri_c,G_IF_SIGMA_CCVV_COVV_NO9_X1_TYPE1_ERI_C)
      (sb, ib, sw, iw, T2b.cptr(), W13avv_sigma_ccvv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| W18(w,a0,v0,a) += (    1.00000000) V2(w,a1,v0,a) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (    1.00000000) T2(x,a0,v0,b) W18(w,a0,v0,a) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W18avv_sigma_ccvv_covv(orz::mr::sizeof_sympack_Xavv(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccvv_covv_no10_x0_type1_eri_c,G_IF_SIGMA_CCVV_COVV_NO10_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W18avv_sigma_ccvv_covv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_covv_no10_x1_type1_eri_c,G_IF_SIGMA_CCVV_COVV_NO10_X1_TYPE1_ERI_C)
      (sb, ib, sw, iw, T2b.cptr(), W18avv_sigma_ccvv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| W19(x,a0,v0,a) += (    1.00000000) V2(x,a1,v0,a) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (   -2.00000000) T2(w,a0,v0,b) W19(x,a0,v0,a) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W19avv_sigma_ccvv_covv(orz::mr::sizeof_sympack_Xavv(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccvv_covv_no11_x0_type1_eri_c,G_IF_SIGMA_CCVV_COVV_NO11_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W19avv_sigma_ccvv_covv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_covv_no11_x1_type1_eri_c,G_IF_SIGMA_CCVV_COVV_NO11_X1_TYPE1_ERI_C)
      (sb, ib, sx, ix, T2b.cptr(), W19avv_sigma_ccvv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| W20(w,a0,v0,a) += (    1.00000000) V2(w,a,v0,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (    4.00000000) T2(a0,x,v0,b) W20(w,a0,v0,a) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W20avv_sigma_ccvv_covv(orz::mr::sizeof_sympack_Xavv(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccvv_covv_no12_x0_type1_eri_c,G_IF_SIGMA_CCVV_COVV_NO12_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W20avv_sigma_ccvv_covv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_covv_no12_x1_type1_eri_c,G_IF_SIGMA_CCVV_COVV_NO12_X1_TYPE1_ERI_C)
      (sb, ib, sw, iw, T2b.cptr(), W20avv_sigma_ccvv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| W21(x,a0,v0,a) += (    1.00000000) V2(x,a,v0,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (   -2.00000000) T2(a0,w,v0,b) W21(x,a0,v0,a) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W21avv_sigma_ccvv_covv(orz::mr::sizeof_sympack_Xavv(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccvv_covv_no13_x0_type1_eri_c,G_IF_SIGMA_CCVV_COVV_NO13_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W21avv_sigma_ccvv_covv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_covv_no13_x1_type1_eri_c,G_IF_SIGMA_CCVV_COVV_NO13_X1_TYPE1_ERI_C)
      (sb, ib, sx, ix, T2b.cptr(), W21avv_sigma_ccvv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   14] -- 
  // |-- [    0] --| W26(w,a0,v0,a) += (    1.00000000) V2(w,a1,v0,a) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (   -2.00000000) T2(a0,x,v0,b) W26(w,a0,v0,a) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W26avv_sigma_ccvv_covv(orz::mr::sizeof_sympack_Xavv(symblockinfo, sw));
  FC_FUNC(g_if_sigma_ccvv_covv_no14_x0_type1_eri_c,G_IF_SIGMA_CCVV_COVV_NO14_X0_TYPE1_ERI_C)
    (sw, iw, V2_sym.cptr(), W26avv_sigma_ccvv_covv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_covv_no14_x1_type1_eri_c,G_IF_SIGMA_CCVV_COVV_NO14_X1_TYPE1_ERI_C)
      (sb, ib, sw, iw, T2b.cptr(), W26avv_sigma_ccvv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   15] -- 
  // |-- [    0] --| W27(x,a0,v0,a) += (    1.00000000) V2(x,a1,v0,a) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (    1.00000000) T2(a0,w,v0,b) W27(x,a0,v0,a) 
  int sx(s_eri);
  int ix(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W27avv_sigma_ccvv_covv(orz::mr::sizeof_sympack_Xavv(symblockinfo, sx));
  FC_FUNC(g_if_sigma_ccvv_covv_no15_x0_type1_eri_c,G_IF_SIGMA_CCVV_COVV_NO15_X0_TYPE1_ERI_C)
    (sx, ix, V2_sym.cptr(), W27avv_sigma_ccvv_covv.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_covv_no15_x1_type1_eri_c,G_IF_SIGMA_CCVV_COVV_NO15_X1_TYPE1_ERI_C)
      (sb, ib, sx, ix, T2b.cptr(), W27avv_sigma_ccvv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- Title : sigma_ccvv_covv

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W14(x,a0,v0,b) += (    1.00000000) V2(b,x,v0,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (   -2.00000000) T2(w,a0,v0,a) W14(x,a0,v0,b) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  orz::DTensor W14cav_sigma_ccvv_covv(orz::mr::sizeof_sympack_Xcav(symblockinfo, sb));
  FC_FUNC(g_if_sigma_ccvv_covv_no0_x0_type1_eri_v,G_IF_SIGMA_CCVV_COVV_NO0_X0_TYPE1_ERI_V)
    (sb, ib, V2_sym.cptr(), W14cav_sigma_ccvv_covv.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    FC_FUNC(g_if_sigma_ccvv_covv_no0_x1_type1_eri_v,G_IF_SIGMA_CCVV_COVV_NO0_X1_TYPE1_ERI_V)
      (sa, ia, sb, ib, T2b.cptr(), W14cav_sigma_ccvv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
  } // End ia
  } // End sa
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W15(w,a0,v0,b) += (    1.00000000) V2(b,w,v0,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (    1.00000000) T2(x,a0,v0,a) W15(w,a0,v0,b) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  orz::DTensor W15cav_sigma_ccvv_covv(orz::mr::sizeof_sympack_Xcav(symblockinfo, sb));
  FC_FUNC(g_if_sigma_ccvv_covv_no1_x0_type1_eri_v,G_IF_SIGMA_CCVV_COVV_NO1_X0_TYPE1_ERI_V)
    (sb, ib, V2_sym.cptr(), W15cav_sigma_ccvv_covv.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    FC_FUNC(g_if_sigma_ccvv_covv_no1_x1_type1_eri_v,G_IF_SIGMA_CCVV_COVV_NO1_X1_TYPE1_ERI_V)
      (sa, ia, sb, ib, T2b.cptr(), W15cav_sigma_ccvv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
  } // End ia
  } // End sa
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W16(w,a0,v0,b) += (    1.00000000) V2(v0,b,w,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (   -2.00000000) T2(a0,x,a,v0) W16(w,a0,v0,b) 
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
    orz::DTensor W16ca_sigma_ccvv_covv(orz::mr::sizeof_sympack_Xca(symblockinfo, sv0^sb));
    FC_FUNC(g_if_sigma_ccvv_covv_no2_x0_type1_eri_v,G_IF_SIGMA_CCVV_COVV_NO2_X0_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, V2_sym.cptr(), W16ca_sigma_ccvv_covv.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccvv_covv_no2_x1_type1_eri_v,G_IF_SIGMA_CCVV_COVV_NO2_X1_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), W16ca_sigma_ccvv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W17(x,a0,v0,b) += (    1.00000000) V2(v0,b,x,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (    1.00000000) T2(a0,w,a,v0) W17(x,a0,v0,b) 
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
    orz::DTensor W17ca_sigma_ccvv_covv(orz::mr::sizeof_sympack_Xca(symblockinfo, sv0^sb));
    FC_FUNC(g_if_sigma_ccvv_covv_no3_x0_type1_eri_v,G_IF_SIGMA_CCVV_COVV_NO3_X0_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, V2_sym.cptr(), W17ca_sigma_ccvv_covv.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccvv_covv_no3_x1_type1_eri_v,G_IF_SIGMA_CCVV_COVV_NO3_X1_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), W17ca_sigma_ccvv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W22(x,a0,v0,b) += (    1.00000000) V2(b,x,v0,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (    4.00000000) T2(w,a0,a,v0) W22(x,a0,v0,b) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    orz::DTensor W22ca_sigma_ccvv_covv(orz::mr::sizeof_sympack_Xca(symblockinfo, sv0^sb));
    FC_FUNC(g_if_sigma_ccvv_covv_no4_x0_type1_eri_v,G_IF_SIGMA_CCVV_COVV_NO4_X0_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, V2_sym.cptr(), W22ca_sigma_ccvv_covv.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccvv_covv_no4_x1_type1_eri_v,G_IF_SIGMA_CCVV_COVV_NO4_X1_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), W22ca_sigma_ccvv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W23(w,a0,v0,b) += (    1.00000000) V2(b,w,v0,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (   -2.00000000) T2(x,a0,a,v0) W23(w,a0,v0,b) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    orz::DTensor W23ca_sigma_ccvv_covv(orz::mr::sizeof_sympack_Xca(symblockinfo, sv0^sb));
    FC_FUNC(g_if_sigma_ccvv_covv_no5_x0_type1_eri_v,G_IF_SIGMA_CCVV_COVV_NO5_X0_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, V2_sym.cptr(), W23ca_sigma_ccvv_covv.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccvv_covv_no5_x1_type1_eri_v,G_IF_SIGMA_CCVV_COVV_NO5_X1_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), W23ca_sigma_ccvv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W24(w,a0,v0,b) += (    1.00000000) V2(v0,b,w,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (    1.00000000) T2(x,a0,a,v0) W24(w,a0,v0,b) 
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
    orz::DTensor W24ca_sigma_ccvv_covv(orz::mr::sizeof_sympack_Xca(symblockinfo, sv0^sb));
    FC_FUNC(g_if_sigma_ccvv_covv_no6_x0_type1_eri_v,G_IF_SIGMA_CCVV_COVV_NO6_X0_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, V2_sym.cptr(), W24ca_sigma_ccvv_covv.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccvv_covv_no6_x1_type1_eri_v,G_IF_SIGMA_CCVV_COVV_NO6_X1_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), W24ca_sigma_ccvv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W25(x,a0,v0,b) += (    1.00000000) V2(v0,b,x,a1) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (   -2.00000000) T2(w,a0,a,v0) W25(x,a0,v0,b) 
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
    orz::DTensor W25ca_sigma_ccvv_covv(orz::mr::sizeof_sympack_Xca(symblockinfo, sv0^sb));
    FC_FUNC(g_if_sigma_ccvv_covv_no7_x0_type1_eri_v,G_IF_SIGMA_CCVV_COVV_NO7_X0_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, V2_sym.cptr(), W25ca_sigma_ccvv_covv.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccvv_covv_no7_x1_type1_eri_v,G_IF_SIGMA_CCVV_COVV_NO7_X1_TYPE1_ERI_V)
      (sb, ib, sv0, iv0, T2b.cptr(), W25ca_sigma_ccvv_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  orz::mr::Femto::my_timer.push_back(boost::make_tuple("sigma_ccvv_covv", time_sigma.elapsed_cputime(), time_sigma.elapsed_wallclocktime()));
  orz::mr::Femto::file_timing << "* " << boost::format("%20s : %10.7f %10.7f ") % orz::mr::Femto::my_timer.back().get<0>() % orz::mr::Femto::my_timer.back().get<1>() % orz::mr::Femto::my_timer.back().get<2>() << endl;
  flush(orz::mr::Femto::file_timing);

  return retval; 
} 
