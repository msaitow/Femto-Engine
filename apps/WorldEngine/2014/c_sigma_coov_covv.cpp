                                                                                
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
#include <sci/icmr/Femto/elems/c_sigma_coov_covv.h>                                  
                                                                                
using std::cout;                                                                
using std::endl;                                                                
                                                                                
#define FLOPCOUNT                                                               
                                                                                
using orz::mr::Femto::Fc0;                                                      
using orz::mr::Femto::h1_int;                                                   
using orz::mr::Femto::h6_int;                                                   
                                                                                
//  `MMMMMMM                                         
//   MM    \                         /              
//   MM       ____  ___  __    __   /M      _____    
//   MM   ,  6MMMMb `MM 6MMb  6MMb /MMMMM  6MMMMMb   
//   MMMMMM 6M'  `Mb MM69 `MM69 `Mb MM    6M'   `Mb  
//   MM   ` MM    MM MM'   MM'   MM MM    MM     MM  
//   MM     MMMMMMMM MM    MM    MM MM    MM     MM  
//   MM     MM       MM    MM    MM MM    MM     MM  
//   MM     YM    d9 MM    MM    MM YM.  ,YM.   ,M9  
//  _MM_     YMMMM9 _MM_  _MM_  _MM_ YMMM9 YMMMMM9   

//                                   Generated date : Sun Apr 20 10:26:19 2014

                                                                                
// ***************************************************************************  
// orz::mr::mrci                                                                
// ***************************************************************************  
									    /*!        
   @brief CT input                                                              
									     */        
                                                                                
orz::mr::BareAmpPack orz::mr::Femto::sigma_coov_covv(const orz::mr::Input &ctinp,                                    
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
  // |-- [    0] --| W0(w,a1,a0,a) += (    1.00000000) T2(w,a0,v0,a) Fc1(v0,a1) 
  // |-- [    1] --| S2(w,i,j,a) += (   -1.00000000) D2(j,a1,a0,i) W0(w,a1,a0,a) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W0caa_sigma_coov_covv(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_covv_no0_x0_type0_noeri,G_IF_SIGMA_COOV_COVV_NO0_X0_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W0caa_sigma_coov_covv.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_coov_covv_no0_x1_type0_noeri,G_IF_SIGMA_COOV_COVV_NO0_X1_TYPE0_NOERI)
      (sa, ia, W0caa_sigma_coov_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W1(w,i,v0,a) += (    1.00000000) T2(w,a0,v0,a) D1(i,a0) 
  // |-- [    1] --| S2(w,i,j,a) += (    2.00000000) Fc1(v0,j) W1(w,i,v0,a) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W1cav_sigma_coov_covv(orz::mr::sizeof_sympack_Xcav(symblockinfo, sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_covv_no1_x0_type0_noeri,G_IF_SIGMA_COOV_COVV_NO1_X0_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W1cav_sigma_coov_covv.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_coov_covv_no1_x1_type0_noeri,G_IF_SIGMA_COOV_COVV_NO1_X1_TYPE0_NOERI)
      (sa, ia, W1cav_sigma_coov_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W2(w,a1,a0,a) += (    1.00000000) T2(a0,w,v0,a) Fc1(v0,a1) 
  // |-- [    1] --| S2(w,i,j,a) += (   -1.00000000) D2(j,i,a0,a1) W2(w,a1,a0,a) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W2caa_sigma_coov_covv(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_covv_no2_x0_type0_noeri,G_IF_SIGMA_COOV_COVV_NO2_X0_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W2caa_sigma_coov_covv.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_coov_covv_no2_x1_type0_noeri,G_IF_SIGMA_COOV_COVV_NO2_X1_TYPE0_NOERI)
      (sa, ia, W2caa_sigma_coov_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W3(w,i,a,v0) += (    1.00000000) T2(a0,w,v0,a) D1(i,a0) 
  // |-- [    1] --| S2(w,i,j,a) += (   -1.00000000) Fc1(v0,j) W3(w,i,a,v0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W3cav_sigma_coov_covv(orz::mr::sizeof_sympack_Xcav(symblockinfo, sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_covv_no3_x0_type0_noeri,G_IF_SIGMA_COOV_COVV_NO3_X0_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W3cav_sigma_coov_covv.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_coov_covv_no3_x1_type0_noeri,G_IF_SIGMA_COOV_COVV_NO3_X1_TYPE0_NOERI)
      (sa, ia, W3cav_sigma_coov_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- Title : sigma_coov_covv

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W4(w,a0,a1,a) += (    1.00000000) V2(w,a1,v0,c0) T2(c0,a0,v0,a) 
  // |-- [    1] --| S2(w,i,j,a) += (   -2.00000000) D2(j,a1,a0,i) W4(w,a0,a1,a) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W4aa_sigma_coov_covv(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_covv_no0_x0_type1_eri_c,G_IF_SIGMA_COOV_COVV_NO0_X0_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), V2_sym.cptr(), W4aa_sigma_coov_covv.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_coov_covv_no0_x1_type1_eri_c,G_IF_SIGMA_COOV_COVV_NO0_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, W4aa_sigma_coov_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W5(w,a0,a1,a) += (    1.00000000) V2(w,c0,v0,a1) T2(c0,a0,v0,a) 
  // |-- [    1] --| S2(w,i,j,a) += (    1.00000000) D2(j,a1,a0,i) W5(w,a0,a1,a) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W5aa_sigma_coov_covv(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_covv_no1_x0_type1_eri_c,G_IF_SIGMA_COOV_COVV_NO1_X0_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), V2_sym.cptr(), W5aa_sigma_coov_covv.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_coov_covv_no1_x1_type1_eri_c,G_IF_SIGMA_COOV_COVV_NO1_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, W5aa_sigma_coov_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W7(w,a0,j,a) += (    1.00000000) V2(w,j,v0,c0) T2(c0,a0,v0,a) 
  // |-- [    1] --| S2(w,i,j,a) += (    4.00000000) D1(i,a0) W7(w,a0,j,a) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W7aa_sigma_coov_covv(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_covv_no2_x0_type1_eri_c,G_IF_SIGMA_COOV_COVV_NO2_X0_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), V2_sym.cptr(), W7aa_sigma_coov_covv.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_coov_covv_no2_x1_type1_eri_c,G_IF_SIGMA_COOV_COVV_NO2_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, W7aa_sigma_coov_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W9(w,a0,j,a) += (    1.00000000) V2(w,c0,v0,j) T2(c0,a0,v0,a) 
  // |-- [    1] --| S2(w,i,j,a) += (   -2.00000000) D1(i,a0) W9(w,a0,j,a) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W9aa_sigma_coov_covv(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_covv_no3_x0_type1_eri_c,G_IF_SIGMA_COOV_COVV_NO3_X0_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), V2_sym.cptr(), W9aa_sigma_coov_covv.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_coov_covv_no3_x1_type1_eri_c,G_IF_SIGMA_COOV_COVV_NO3_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, W9aa_sigma_coov_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W11(w,a0,a1,a) += (    1.00000000) V2(w,a1,v0,c0) T2(a0,c0,v0,a) 
  // |-- [    1] --| S2(w,i,j,a) += (    1.00000000) D2(j,a1,a0,i) W11(w,a0,a1,a) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W11aa_sigma_coov_covv(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_covv_no4_x0_type1_eri_c,G_IF_SIGMA_COOV_COVV_NO4_X0_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), V2_sym.cptr(), W11aa_sigma_coov_covv.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_coov_covv_no4_x1_type1_eri_c,G_IF_SIGMA_COOV_COVV_NO4_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, W11aa_sigma_coov_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W12(w,a0,a1,a) += (    1.00000000) V2(w,c0,v0,a1) T2(a0,c0,v0,a) 
  // |-- [    1] --| S2(w,i,j,a) += (    1.00000000) D2(j,i,a0,a1) W12(w,a0,a1,a) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W12aa_sigma_coov_covv(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_covv_no5_x0_type1_eri_c,G_IF_SIGMA_COOV_COVV_NO5_X0_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), V2_sym.cptr(), W12aa_sigma_coov_covv.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_coov_covv_no5_x1_type1_eri_c,G_IF_SIGMA_COOV_COVV_NO5_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, W12aa_sigma_coov_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W14(w,a0,j,a) += (    1.00000000) V2(w,j,v0,c0) T2(a0,c0,v0,a) 
  // |-- [    1] --| S2(w,i,j,a) += (   -2.00000000) D1(i,a0) W14(w,a0,j,a) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W14aa_sigma_coov_covv(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_covv_no6_x0_type1_eri_c,G_IF_SIGMA_COOV_COVV_NO6_X0_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), V2_sym.cptr(), W14aa_sigma_coov_covv.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_coov_covv_no6_x1_type1_eri_c,G_IF_SIGMA_COOV_COVV_NO6_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, W14aa_sigma_coov_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W16(w,a0,j,a) += (    1.00000000) V2(w,c0,v0,j) T2(a0,c0,v0,a) 
  // |-- [    1] --| S2(w,i,j,a) += (    1.00000000) D1(i,a0) W16(w,a0,j,a) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W16aa_sigma_coov_covv(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_covv_no7_x0_type1_eri_c,G_IF_SIGMA_COOV_COVV_NO7_X0_TYPE1_ERI_C)
      (sa, ia, sw, iw, T2b.cptr(), V2_sym.cptr(), W16aa_sigma_coov_covv.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_coov_covv_no7_x1_type1_eri_c,G_IF_SIGMA_COOV_COVV_NO7_X1_TYPE1_ERI_C)
      (sa, ia, sw, iw, W16aa_sigma_coov_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- Title : sigma_coov_covv

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W6(j,a0,i,v0) += (    1.00000000) V2(v0,a2,a3,a1) D3(j,a2,a3,a1,a0,i) 
  // |-- [    1] --| S2(w,i,j,a) += (   -1.00000000) T2(w,a0,v0,a) W6(j,a0,i,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W6aaa_sigma_coov_covv(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_coov_covv_no0_x0_type1_eri_v,G_IF_SIGMA_COOV_COVV_NO0_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W6aaa_sigma_coov_covv.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_covv_no0_x1_type1_eri_v,G_IF_SIGMA_COOV_COVV_NO0_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), W6aaa_sigma_coov_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W8(i,a0,j,v0) += (    1.00000000) V2(v0,a2,j,a1) D2(i,a0,a2,a1) 
  // |-- [    1] --| S2(w,i,j,a) += (   -1.00000000) T2(w,a0,v0,a) W8(i,a0,j,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W8aaa_sigma_coov_covv(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_coov_covv_no1_x0_type1_eri_v,G_IF_SIGMA_COOV_COVV_NO1_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W8aaa_sigma_coov_covv.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_covv_no1_x1_type1_eri_v,G_IF_SIGMA_COOV_COVV_NO1_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), W8aaa_sigma_coov_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W10(i,a0,j,v0) += (    1.00000000) V2(v0,j,a2,a1) D2(i,a0,a2,a1) 
  // |-- [    1] --| S2(w,i,j,a) += (    2.00000000) T2(w,a0,v0,a) W10(i,a0,j,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W10aaa_sigma_coov_covv(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_coov_covv_no2_x0_type1_eri_v,G_IF_SIGMA_COOV_COVV_NO2_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W10aaa_sigma_coov_covv.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_covv_no2_x1_type1_eri_v,G_IF_SIGMA_COOV_COVV_NO2_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), W10aaa_sigma_coov_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W13(j,a0,i,v0) += (    1.00000000) V2(v0,a2,a3,a1) D3(j,i,a3,a1,a0,a2) 
  // |-- [    1] --| S2(w,i,j,a) += (   -1.00000000) T2(w,a0,a,v0) W13(j,a0,i,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W13aaa_sigma_coov_covv(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_coov_covv_no3_x0_type1_eri_v,G_IF_SIGMA_COOV_COVV_NO3_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W13aaa_sigma_coov_covv.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_covv_no3_x1_type1_eri_v,G_IF_SIGMA_COOV_COVV_NO3_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), W13aaa_sigma_coov_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W15(i,a0,j,v0) += (    1.00000000) V2(v0,a2,j,a1) D2(i,a1,a2,a0) 
  // |-- [    1] --| S2(w,i,j,a) += (   -1.00000000) T2(w,a0,a,v0) W15(i,a0,j,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W15aaa_sigma_coov_covv(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_coov_covv_no4_x0_type1_eri_v,G_IF_SIGMA_COOV_COVV_NO4_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W15aaa_sigma_coov_covv.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_covv_no4_x1_type1_eri_v,G_IF_SIGMA_COOV_COVV_NO4_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), W15aaa_sigma_coov_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W17(i,a0,j,v0) += (    1.00000000) V2(v0,j,a2,a1) D2(i,a0,a2,a1) 
  // |-- [    1] --| S2(w,i,j,a) += (   -1.00000000) T2(w,a0,a,v0) W17(i,a0,j,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W17aaa_sigma_coov_covv(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_coov_covv_no5_x0_type1_eri_v,G_IF_SIGMA_COOV_COVV_NO5_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W17aaa_sigma_coov_covv.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_covv_no5_x1_type1_eri_v,G_IF_SIGMA_COOV_COVV_NO5_X1_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), W17aaa_sigma_coov_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W18(w,a0,a1,a) += (    1.00000000) V2(a,v0,v1,a1) T2(w,a0,v0,v1) 
  // |-- [    1] --| S2(w,i,j,a) += (   -1.00000000) D2(j,i,a0,a1) W18(w,a0,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W18caa_sigma_coov_covv(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  for(int sv1 = 0;sv1 < nir;++sv1){ 
  for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
    T2b = T2.get_amp2(iv1);
    FC_FUNC(g_if_sigma_coov_covv_no6_x0_type1_eri_v,G_IF_SIGMA_COOV_COVV_NO6_X0_TYPE1_ERI_V)
      (sa, ia, sv1, iv1, T2b.cptr(), V2_sym.cptr(), W18caa_sigma_coov_covv.cptr(), nir, nsym, psym, &flops);
  // --> @[v1, "virtual"] [notNeeded]
  } // End iv1
  } // End sv1
  FC_FUNC(g_if_sigma_coov_covv_no6_x1_type1_eri_v,G_IF_SIGMA_COOV_COVV_NO6_X1_TYPE1_ERI_V)
    (sa, ia, W18caa_sigma_coov_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W19(w,a0,a1,a) += (    1.00000000) V2(a,v0,v1,a1) T2(w,a0,v1,v0) 
  // |-- [    1] --| S2(w,i,j,a) += (   -1.00000000) D2(j,a1,a0,i) W19(w,a0,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W19caa_sigma_coov_covv(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_coov_covv_no7_x0_type1_eri_v,G_IF_SIGMA_COOV_COVV_NO7_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W19caa_sigma_coov_covv.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_coov_covv_no7_x1_type1_eri_v,G_IF_SIGMA_COOV_COVV_NO7_X1_TYPE1_ERI_V)
    (sa, ia, W19caa_sigma_coov_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| W20(w,a0,j,a) += (    1.00000000) V2(a,v0,v1,j) T2(w,a0,v1,v0) 
  // |-- [    1] --| S2(w,i,j,a) += (    2.00000000) D1(i,a0) W20(w,a0,j,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W20caa_sigma_coov_covv(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_coov_covv_no8_x0_type1_eri_v,G_IF_SIGMA_COOV_COVV_NO8_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W20caa_sigma_coov_covv.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_coov_covv_no8_x1_type1_eri_v,G_IF_SIGMA_COOV_COVV_NO8_X1_TYPE1_ERI_V)
    (sa, ia, W20caa_sigma_coov_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W21(w,a0,j,a) += (    1.00000000) V2(a,v0,v1,j) T2(w,a0,v0,v1) 
  // |-- [    1] --| S2(w,i,j,a) += (   -1.00000000) D1(i,a0) W21(w,a0,j,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W21caa_sigma_coov_covv(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  for(int sv1 = 0;sv1 < nir;++sv1){ 
  for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
    T2b = T2.get_amp2(iv1);
    FC_FUNC(g_if_sigma_coov_covv_no9_x0_type1_eri_v,G_IF_SIGMA_COOV_COVV_NO9_X0_TYPE1_ERI_V)
      (sa, ia, sv1, iv1, T2b.cptr(), V2_sym.cptr(), W21caa_sigma_coov_covv.cptr(), nir, nsym, psym, &flops);
  // --> @[v1, "virtual"] [notNeeded]
  } // End iv1
  } // End sv1
  FC_FUNC(g_if_sigma_coov_covv_no9_x1_type1_eri_v,G_IF_SIGMA_COOV_COVV_NO9_X1_TYPE1_ERI_V)
    (sa, ia, W21caa_sigma_coov_covv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
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
  orz::mr::Femto::my_timer.push_back(boost::make_tuple("sigma_coov_covv", time_sigma.elapsed_cputime(), time_sigma.elapsed_wallclocktime()));
  orz::mr::Femto::file_timing << "* " << boost::format("%20s : %10.7f %10.7f ") % orz::mr::Femto::my_timer.back().get<0>() % orz::mr::Femto::my_timer.back().get<1>() % orz::mr::Femto::my_timer.back().get<2>() << endl;
  flush(orz::mr::Femto::file_timing);

  return retval; 
} 
