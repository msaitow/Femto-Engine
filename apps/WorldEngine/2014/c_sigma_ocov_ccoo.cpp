                                                                                
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
#include <sci/icmr/Femto/elems/c_sigma_ocov_ccoo.h>                                  
                                                                                
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

//                                   Generated date : Sun Apr 20 10:26:12 2014

                                                                                
// ***************************************************************************  
// orz::mr::mrci                                                                
// ***************************************************************************  
									    /*!        
   @brief CT input                                                              
									     */        
                                                                                
orz::mr::BareAmpPack orz::mr::Femto::sigma_ocov_ccoo(const orz::mr::Input &ctinp,                                    
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
  // -- Title : sigma_ocov_ccoo

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W0(w,c0,j,a3,a2,i) += (    1.00000000) T2(w,c0,a1,a0) D3(j,i,a3,a1,a2,a0) 
  // |-- [    1] --| S2(i,w,j,a) += (    2.00000000) V2(c0,a2,a3,a) W0(w,c0,j,a3,a2,i) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W0caaaa_sigma_ocov_ccoo(orz::mr::sizeof_sympack_Xcaaaa(symblockinfo, sc0));
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ocov_ccoo_no0_x0_type1_eri_c,G_IF_SIGMA_OCOV_CCOO_NO0_X0_TYPE1_ERI_C)
      (sa0, ia0, sc0, ic0, T2b.cptr(), W0caaaa_sigma_ocov_ccoo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ocov_ccoo_no0_x1_type1_eri_c,G_IF_SIGMA_OCOV_CCOO_NO0_X1_TYPE1_ERI_C)
      (sa, ia, sc0, ic0, V2_sym.cptr(), W0caaaa_sigma_ocov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- Title : sigma_ocov_ccoo

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W7(c0,i,a0,a) += (    1.00000000) V2(a0,c0,a1,a) D1(i,a1) 
  // |-- [    1] --| S2(i,w,j,a) += (   -4.00000000) T2(w,c0,j,a0) W7(c0,i,a0,a) 
  int sa0(s_eri);
  int ia0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(ia0);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W7ca_sigma_ocov_ccoo(orz::mr::sizeof_sympack_Xca(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_ocov_ccoo_no0_x0_type1_eri_o,G_IF_SIGMA_OCOV_CCOO_NO0_X0_TYPE1_ERI_O)
      (sa, ia, sa0, ia0, V2_sym.cptr(), W7ca_sigma_ocov_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ocov_ccoo_no0_x1_type1_eri_o,G_IF_SIGMA_OCOV_CCOO_NO0_X1_TYPE1_ERI_O)
      (sa, ia, sa0, ia0, T2b.cptr(), W7ca_sigma_ocov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W8(c0,i,a0,a) += (    1.00000000) V2(a0,c0,a1,a) D1(i,a1) 
  // |-- [    1] --| S2(i,w,j,a) += (    2.00000000) T2(c0,w,j,a0) W8(c0,i,a0,a) 
  int sa0(s_eri);
  int ia0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(ia0);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W8ca_sigma_ocov_ccoo(orz::mr::sizeof_sympack_Xca(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_ocov_ccoo_no1_x0_type1_eri_o,G_IF_SIGMA_OCOV_CCOO_NO1_X0_TYPE1_ERI_O)
      (sa, ia, sa0, ia0, V2_sym.cptr(), W8ca_sigma_ocov_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ocov_ccoo_no1_x1_type1_eri_o,G_IF_SIGMA_OCOV_CCOO_NO1_X1_TYPE1_ERI_O)
      (sa, ia, sa0, ia0, T2b.cptr(), W8ca_sigma_ocov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W9(c0,i,a1,a) += (    1.00000000) V2(a1,a,c0,a0) D1(i,a0) 
  // |-- [    1] --| S2(i,w,j,a) += (   -4.00000000) T2(c0,w,j,a1) W9(c0,i,a1,a) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(ia1);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W9ca_sigma_ocov_ccoo(orz::mr::sizeof_sympack_Xca(symblockinfo, sa1^sa));
    FC_FUNC(g_if_sigma_ocov_ccoo_no2_x0_type1_eri_o,G_IF_SIGMA_OCOV_CCOO_NO2_X0_TYPE1_ERI_O)
      (sa, ia, sa1, ia1, V2_sym.cptr(), W9ca_sigma_ocov_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ocov_ccoo_no2_x1_type1_eri_o,G_IF_SIGMA_OCOV_CCOO_NO2_X1_TYPE1_ERI_O)
      (sa, ia, sa1, ia1, T2b.cptr(), W9ca_sigma_ocov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W10(c0,i,a1,a) += (    1.00000000) V2(a1,a,c0,a0) D1(i,a0) 
  // |-- [    1] --| S2(i,w,j,a) += (    2.00000000) T2(w,c0,j,a1) W10(c0,i,a1,a) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(ia1);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W10ca_sigma_ocov_ccoo(orz::mr::sizeof_sympack_Xca(symblockinfo, sa1^sa));
    FC_FUNC(g_if_sigma_ocov_ccoo_no3_x0_type1_eri_o,G_IF_SIGMA_OCOV_CCOO_NO3_X0_TYPE1_ERI_O)
      (sa, ia, sa1, ia1, V2_sym.cptr(), W10ca_sigma_ocov_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ocov_ccoo_no3_x1_type1_eri_o,G_IF_SIGMA_OCOV_CCOO_NO3_X1_TYPE1_ERI_O)
      (sa, ia, sa1, ia1, T2b.cptr(), W10ca_sigma_ocov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- Title : sigma_ocov_ccoo

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W1(c0,i,a0,a) += (    1.00000000) V2(a,a2,c0,a1) D2(i,a2,a0,a1) 
  // |-- [    1] --| S2(i,w,j,a) += (    2.00000000) T2(w,c0,j,a0) W1(c0,i,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W1ca_sigma_ocov_ccoo(orz::mr::sizeof_sympack_Xca(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_ocov_ccoo_no0_x0_type1_eri_v,G_IF_SIGMA_OCOV_CCOO_NO0_X0_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, V2_sym.cptr(), W1ca_sigma_ocov_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ocov_ccoo_no0_x1_type1_eri_v,G_IF_SIGMA_OCOV_CCOO_NO0_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), W1ca_sigma_ocov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W2(c0,i,a0,a) += (    1.00000000) V2(a,a2,c0,a1) D2(i,a1,a0,a2) 
  // |-- [    1] --| S2(i,w,j,a) += (    2.00000000) T2(w,c0,a0,j) W2(c0,i,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W2caa_sigma_ocov_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ocov_ccoo_no1_x0_type1_eri_v,G_IF_SIGMA_OCOV_CCOO_NO1_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W2caa_sigma_ocov_ccoo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    FC_FUNC(g_if_sigma_ocov_ccoo_no1_x1_type1_eri_v,G_IF_SIGMA_OCOV_CCOO_NO1_X1_TYPE1_ERI_V)
      (sa, ia, sj, ij, T2b.cptr(), W2caa_sigma_ocov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
  } // End ij
  } // End sj
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W3(w,a0,a2,a) += (    1.00000000) V2(a,a2,c0,a1) T2(w,c0,a1,a0) 
  // |-- [    1] --| S2(i,w,j,a) += (    2.00000000) D2(j,i,a2,a0) W3(w,a0,a2,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W3ca_sigma_ocov_ccoo(orz::mr::sizeof_sympack_Xca(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_ocov_ccoo_no2_x0_type1_eri_v,G_IF_SIGMA_OCOV_CCOO_NO2_X0_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), V2_sym.cptr(), W3ca_sigma_ocov_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ocov_ccoo_no2_x1_type1_eri_v,G_IF_SIGMA_OCOV_CCOO_NO2_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, W3ca_sigma_ocov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W4(w,a0,a2,a) += (    1.00000000) V2(a,a2,c0,a1) T2(w,c0,a0,a1) 
  // |-- [    1] --| S2(i,w,j,a) += (   -4.00000000) D2(j,i,a2,a0) W4(w,a0,a2,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W4caa_sigma_ocov_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_ocov_ccoo_no3_x0_type1_eri_v,G_IF_SIGMA_OCOV_CCOO_NO3_X0_TYPE1_ERI_V)
      (sa, ia, sa1, ia1, T2b.cptr(), V2_sym.cptr(), W4caa_sigma_ocov_ccoo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  FC_FUNC(g_if_sigma_ocov_ccoo_no3_x1_type1_eri_v,G_IF_SIGMA_OCOV_CCOO_NO3_X1_TYPE1_ERI_V)
    (sa, ia, W4caa_sigma_ocov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W5(w,a0,a1,a) += (    1.00000000) V2(a,a2,c0,a1) T2(w,c0,a2,a0) 
  // |-- [    1] --| S2(i,w,j,a) += (   -4.00000000) D2(j,i,a1,a0) W5(w,a0,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W5ca_sigma_ocov_ccoo(orz::mr::sizeof_sympack_Xca(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_ocov_ccoo_no4_x0_type1_eri_v,G_IF_SIGMA_OCOV_CCOO_NO4_X0_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), V2_sym.cptr(), W5ca_sigma_ocov_ccoo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ocov_ccoo_no4_x1_type1_eri_v,G_IF_SIGMA_OCOV_CCOO_NO4_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, W5ca_sigma_ocov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W6(w,a0,a1,a) += (    1.00000000) V2(a,a2,c0,a1) T2(w,c0,a0,a2) 
  // |-- [    1] --| S2(i,w,j,a) += (    2.00000000) D2(j,i,a1,a0) W6(w,a0,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W6caa_sigma_ocov_ccoo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  for(int sa2 = 0;sa2 < nir;++sa2){ 
  for(int ia2 = symblockinfo.psym()(sa2,I_O,I_BEGIN);ia2 <= symblockinfo.psym()(sa2,I_O,I_END);++ia2){ 
    T2b = T2.get_amp2(ia2);
    FC_FUNC(g_if_sigma_ocov_ccoo_no5_x0_type1_eri_v,G_IF_SIGMA_OCOV_CCOO_NO5_X0_TYPE1_ERI_V)
      (sa, ia, sa2, ia2, T2b.cptr(), V2_sym.cptr(), W6caa_sigma_ocov_ccoo.cptr(), nir, nsym, psym, &flops);
  // --> @[a2, "active"] [notNeeded]
  } // End ia2
  } // End sa2
  FC_FUNC(g_if_sigma_ocov_ccoo_no5_x1_type1_eri_v,G_IF_SIGMA_OCOV_CCOO_NO5_X1_TYPE1_ERI_V)
    (sa, ia, W6caa_sigma_ocov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W11(w,a) += (    1.00000000) V2(a,a1,c0,a0) T2(w,c0,a1,a0) 
  // |-- [    1] --| S2(i,w,j,a) += (    8.00000000) D1(j,i) W11(w,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W11c_sigma_ocov_ccoo(orz::mr::sizeof_sympack_Xc(symblockinfo, sa));
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ocov_ccoo_no6_x0_type1_eri_v,G_IF_SIGMA_OCOV_CCOO_NO6_X0_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), V2_sym.cptr(), W11c_sigma_ocov_ccoo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  FC_FUNC(g_if_sigma_ocov_ccoo_no6_x1_type1_eri_v,G_IF_SIGMA_OCOV_CCOO_NO6_X1_TYPE1_ERI_V)
    (sa, ia, W11c_sigma_ocov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W12(w,a) += (    1.00000000) V2(a,a1,c0,a0) T2(w,c0,a0,a1) 
  // |-- [    1] --| S2(i,w,j,a) += (   -4.00000000) D1(j,i) W12(w,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W12c_sigma_ocov_ccoo(orz::mr::sizeof_sympack_Xc(symblockinfo, sa));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_ocov_ccoo_no7_x0_type1_eri_v,G_IF_SIGMA_OCOV_CCOO_NO7_X0_TYPE1_ERI_V)
      (sa, ia, sa1, ia1, T2b.cptr(), V2_sym.cptr(), W12c_sigma_ocov_ccoo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  FC_FUNC(g_if_sigma_ocov_ccoo_no7_x1_type1_eri_v,G_IF_SIGMA_OCOV_CCOO_NO7_X1_TYPE1_ERI_V)
    (sa, ia, W12c_sigma_ocov_ccoo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  orz::mr::Femto::my_timer.push_back(boost::make_tuple("sigma_ocov_ccoo", time_sigma.elapsed_cputime(), time_sigma.elapsed_wallclocktime()));
  orz::mr::Femto::file_timing << "* " << boost::format("%20s : %10.7f %10.7f ") % orz::mr::Femto::my_timer.back().get<0>() % orz::mr::Femto::my_timer.back().get<1>() % orz::mr::Femto::my_timer.back().get<2>() << endl;
  flush(orz::mr::Femto::file_timing);

  return retval; 
} 
