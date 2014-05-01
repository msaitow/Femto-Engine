                                                                                
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
#include <sci/icmr/Femto/elems/c_sigma_ooov_coov.h>                                  
                                                                                
using std::cout;                                                                
using std::endl;                                                                
                                                                                
#define FLOPCOUNT                                                               
                                                                                
using orz::mr::Femto::Fc0;                                                      
using orz::mr::Femto::h1_int;                                                   
using orz::mr::Femto::h6_int;                                                   
                                                                                
//    o__ __o__/_                            o                         
//   <|    v                                <|>                        
//   < >                                    < >                        
//    |         o__  __o   \o__ __o__ __o    |        o__ __o         
//    o__/_    /v      |>   |     |     |>   o__/_   /v     v\        
//    |       />      //   / \   / \   / \   |      />       <\    
//   <o>      \o    o/     \o/   \o/   \o/   |      \         /   
//    |        v\  /v __o   |     |     |    o       o       o        
//   / \        <\/> __/>  / \   / \   / \   <\__    <\__ __/>  

//                                   Generated date : Sun Apr 20 10:26:09 2014

                                                                                
// ***************************************************************************  
// orz::mr::mrci                                                                
// ***************************************************************************  
									    /*!        
   @brief CT input                                                              
									     */        
                                                                                
orz::mr::BareAmpPack orz::mr::Femto::sigma_ooov_coov(const orz::mr::Input &ctinp,                                    
                                  const orz::mr::SymBlockInfo &symblockinfo,                                 
                                  const orz::mr::HintMO &hintmo,                                             
                                  const int alloc_type,                                                      
                                  const orz::mr::RdmPack &rdmPack,                                           
                                  const orz::DTensor &rdm4,                                                  
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
  orz::DTensor rdm4_sym;                                                                                         
  orz::DTensor rdm4_ij_sliced(ctinp.use_d4cum_of() ? nocc*nocc*nocc*nocc*nocc*nocc : 0);                         
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
  // |-- [    0] --| W0(a2,a1,a0,a) += (    1.00000000) T2(c0,a1,a0,a) Fc1(c0,a2) 
  // |-- [    1] --| S2(i,j,k,a) += (   -1.00000000) D3(k,i,a2,a0,a1,j) W0(a2,a1,a0,a) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W0aaa_sigma_ooov_coov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_coov_no0_x0_type0_noeri,G_IF_SIGMA_OOOV_COOV_NO0_X0_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W0aaa_sigma_ooov_coov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_coov_no0_x1_type0_noeri,G_IF_SIGMA_OOOV_COOV_NO0_X1_TYPE0_NOERI)
      (sa, ia, W0aaa_sigma_ooov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W1(a1,a0,k,a) += (    1.00000000) T2(c0,a0,k,a) Fc1(c0,a1) 
  // |-- [    1] --| S2(i,j,k,a) += (   -1.00000000) D2(j,a0,i,a1) W1(a1,a0,k,a) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W1aaa_sigma_ooov_coov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_coov_no1_x0_type0_noeri,G_IF_SIGMA_OOOV_COOV_NO1_X0_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W1aaa_sigma_ooov_coov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_coov_no1_x1_type0_noeri,G_IF_SIGMA_OOOV_COOV_NO1_X1_TYPE0_NOERI)
      (sa, ia, W1aaa_sigma_ooov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W2(a0,a) += (    1.00000000) T2(c0,a0,a1,a) Fc1(c0,a1) 
  // |-- [    1] --| S2(i,j,k,a) += (    2.00000000) D2(k,i,a0,j) W2(a0,a) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W2a_sigma_ooov_coov(orz::mr::sizeof_sympack_Xa(symblockinfo, sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_coov_no2_x0_type0_noeri,G_IF_SIGMA_OOOV_COOV_NO2_X0_TYPE0_NOERI)
      (sa, ia, T2b.cptr(), W2a_sigma_ooov_coov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_coov_no2_x1_type0_noeri,G_IF_SIGMA_OOOV_COOV_NO2_X1_TYPE0_NOERI)
      (sa, ia, W2a_sigma_ooov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- Title : sigma_ooov_coov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W4(c0,j,i,a0) += (    1.00000000) V2(c0,a2,a3,a1) D3(j,a0,i,a2,a3,a1) 
  // |-- [    1] --| S2(i,j,k,a) += (   -1.00000000) T2(c0,a0,k,a) W4(c0,j,i,a0) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W4aaa_sigma_ooov_coov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sc0));
  FC_FUNC(g_if_sigma_ooov_coov_no0_x0_type1_eri_c,G_IF_SIGMA_OOOV_COOV_NO0_X0_TYPE1_ERI_C)
    (sc0, ic0, V2_sym.cptr(), W4aaa_sigma_ooov_coov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_coov_no0_x1_type1_eri_c,G_IF_SIGMA_OOOV_COOV_NO0_X1_TYPE1_ERI_C)
      (sa, ia, sc0, ic0, T2b.cptr(), W4aaa_sigma_ooov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- Title : sigma_ooov_coov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W3(c0,j,i,a1,a0,k) += (    1.00000000) V2(k,a3,c0,a2) D3(j,a0,i,a3,a1,a2) 
  // |-- [    1] --| S2(i,j,k,a) += (   -1.00000000) T2(c0,a0,a1,a) W3(c0,j,i,a1,a0,k) 
  int sk(s_eri);
  int ik(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W3caaaa_sigma_ooov_coov(orz::mr::sizeof_sympack_Xcaaaa(symblockinfo, sk));
  FC_FUNC(g_if_sigma_ooov_coov_no0_x0_type1_eri_o,G_IF_SIGMA_OOOV_COOV_NO0_X0_TYPE1_ERI_O)
    (sk, ik, V2_sym.cptr(), W3caaaa_sigma_ooov_coov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_coov_no0_x1_type1_eri_o,G_IF_SIGMA_OOOV_COOV_NO0_X1_TYPE1_ERI_O)
      (sa, ia, sk, ik, T2b.cptr(), W3caaaa_sigma_ooov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W5(a0,a3,a1,a) += (    1.00000000) V2(a3,a1,c0,a2) T2(c0,a0,a2,a) 
  // |-- [    1] --| S2(i,j,k,a) += (    2.00000000) D3(k,i,a3,a1,a0,j) W5(a0,a3,a1,a) 
  int sa3(s_eri);
  int ia3(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W5aa_sigma_ooov_coov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa3^sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_coov_no1_x0_type1_eri_o,G_IF_SIGMA_OOOV_COOV_NO1_X0_TYPE1_ERI_O)
      (sa, ia, sa3, ia3, T2b.cptr(), V2_sym.cptr(), W5aa_sigma_ooov_coov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_coov_no1_x1_type1_eri_o,G_IF_SIGMA_OOOV_COOV_NO1_X1_TYPE1_ERI_O)
      (sa, ia, sa3, ia3, W5aa_sigma_ooov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W6(a0,a3,a2,a) += (    1.00000000) V2(a2,c0,a3,a1) T2(c0,a0,a1,a) 
  // |-- [    1] --| S2(i,j,k,a) += (   -1.00000000) D3(k,i,a2,a3,a0,j) W6(a0,a3,a2,a) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W6aa_sigma_ooov_coov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa2^sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_coov_no2_x0_type1_eri_o,G_IF_SIGMA_OOOV_COOV_NO2_X0_TYPE1_ERI_O)
      (sa, ia, sa2, ia2, T2b.cptr(), V2_sym.cptr(), W6aa_sigma_ooov_coov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_coov_no2_x1_type1_eri_o,G_IF_SIGMA_OOOV_COOV_NO2_X1_TYPE1_ERI_O)
      (sa, ia, sa2, ia2, W6aa_sigma_ooov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W7(a0,k,a2,a) += (    1.00000000) V2(k,a2,c0,a1) T2(c0,a0,a1,a) 
  // |-- [    1] --| S2(i,j,k,a) += (    2.00000000) D2(j,a0,i,a2) W7(a0,k,a2,a) 
  int sk(s_eri);
  int ik(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W7aa_sigma_ooov_coov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sk^sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_coov_no3_x0_type1_eri_o,G_IF_SIGMA_OOOV_COOV_NO3_X0_TYPE1_ERI_O)
      (sa, ia, sk, ik, T2b.cptr(), V2_sym.cptr(), W7aa_sigma_ooov_coov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_coov_no3_x1_type1_eri_o,G_IF_SIGMA_OOOV_COOV_NO3_X1_TYPE1_ERI_O)
      (sa, ia, sk, ik, W7aa_sigma_ooov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W8(a0,k,a1,a) += (    1.00000000) V2(k,a2,c0,a1) T2(c0,a0,a2,a) 
  // |-- [    1] --| S2(i,j,k,a) += (   -1.00000000) D2(j,a0,i,a1) W8(a0,k,a1,a) 
  int sk(s_eri);
  int ik(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W8aa_sigma_ooov_coov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sk^sa));
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_coov_no4_x0_type1_eri_o,G_IF_SIGMA_OOOV_COOV_NO4_X0_TYPE1_ERI_O)
      (sa, ia, sk, ik, T2b.cptr(), V2_sym.cptr(), W8aa_sigma_ooov_coov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_coov_no4_x1_type1_eri_o,G_IF_SIGMA_OOOV_COOV_NO4_X1_TYPE1_ERI_O)
      (sa, ia, sk, ik, W8aa_sigma_ooov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- Title : sigma_ooov_coov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W9(a1,a0,a2,a) += (    1.00000000) V2(a,a2,v0,c0) T2(c0,a1,a0,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -1.00000000) D3(k,i,a2,j,a1,a0) W9(a1,a0,a2,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W9aaa_sigma_ooov_coov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_ooov_coov_no0_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOV_NO0_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W9aaa_sigma_ooov_coov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_ooov_coov_no0_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOV_NO0_X1_TYPE1_ERI_V)
    (sa, ia, W9aaa_sigma_ooov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W10(a1,a0,a2,a) += (    1.00000000) V2(a,v0,c0,a2) T2(c0,a1,a0,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -1.00000000) D3(k,i,a2,a0,a1,j) W10(a1,a0,a2,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W10aaa_sigma_ooov_coov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_ooov_coov_no1_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOV_NO1_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W10aaa_sigma_ooov_coov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_ooov_coov_no1_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOV_NO1_X1_TYPE1_ERI_V)
    (sa, ia, W10aaa_sigma_ooov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W11(a0,k,a1,a) += (    1.00000000) V2(a,a1,v0,c0) T2(c0,a0,k,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -1.00000000) D2(j,a1,i,a0) W11(a0,k,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W11aaa_sigma_ooov_coov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_ooov_coov_no2_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOV_NO2_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W11aaa_sigma_ooov_coov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_ooov_coov_no2_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOV_NO2_X1_TYPE1_ERI_V)
    (sa, ia, W11aaa_sigma_ooov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W12(a0,k,a1,a) += (    1.00000000) V2(a,v0,c0,a1) T2(c0,a0,k,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -1.00000000) D2(j,a0,i,a1) W12(a0,k,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W12aaa_sigma_ooov_coov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_ooov_coov_no3_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOV_NO3_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W12aaa_sigma_ooov_coov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_ooov_coov_no3_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOV_NO3_X1_TYPE1_ERI_V)
    (sa, ia, W12aaa_sigma_ooov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W13(a0,a) += (    1.00000000) V2(a,a1,v0,c0) T2(c0,a0,a1,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -1.00000000) D2(k,i,a0,j) W13(a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W13a_sigma_ooov_coov(orz::mr::sizeof_sympack_Xa(symblockinfo, sa));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_ooov_coov_no4_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOV_NO4_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W13a_sigma_ooov_coov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_ooov_coov_no4_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOV_NO4_X1_TYPE1_ERI_V)
    (sa, ia, W13a_sigma_ooov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W14(a0,a) += (    1.00000000) V2(a,v0,c0,a1) T2(c0,a0,a1,v0) 
  // |-- [    1] --| S2(i,j,k,a) += (    2.00000000) D2(k,i,a0,j) W14(a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W14a_sigma_ooov_coov(orz::mr::sizeof_sympack_Xa(symblockinfo, sa));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_ooov_coov_no5_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOV_NO5_X0_TYPE1_ERI_V)
      (sa, ia, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W14a_sigma_ooov_coov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_ooov_coov_no5_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOV_NO5_X1_TYPE1_ERI_V)
    (sa, ia, W14a_sigma_ooov_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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

//-@loadD4C(c,begin)
  //*-- FEMTO begins --//*
  // Label : d4c_c
  {


  for(int s_eri = 0;s_eri < nir;++s_eri){ 
  for(int i_eri = symblockinfo.psym()(s_eri,I_C,I_BEGIN);i_eri <= symblockinfo.psym()(s_eri,I_C,I_END);++i_eri){ 
  if(hintmo.iproc_havingimo()[i_eri] == myrank) {           
  orz::DTensor C5;
  orz::LoadBin(ctinp.dir()/(format("D4C_g[%d]")%i_eri).str()) >> C5;

  //*-- Entering to take the type 1 contractions --*//
//-@type(1).contraction(begin)
  // -- Title : sigma_ooov_coov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| S2(i,j,k,a) += (   -1.00000000) T2(c0,a1,a0,a) C5(j,a1,i,k,a0,c0) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_coov_no0_x0_type1_d4c_c,G_IF_SIGMA_OOOV_COOV_NO0_X0_TYPE1_D4C_C)
      (sa, ia, sc0, ic0, C5.cptr(), T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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

//-@loadD4C(c,end)

//-@D4C.contractions(end)

  // Do timing!
  orz::mr::Femto::my_timer.push_back(boost::make_tuple("sigma_ooov_coov", time_sigma.elapsed_cputime(), time_sigma.elapsed_wallclocktime()));
  orz::mr::Femto::file_timing << "* " << boost::format("%20s : %10.7f %10.7f ") % orz::mr::Femto::my_timer.back().get<0>() % orz::mr::Femto::my_timer.back().get<1>() % orz::mr::Femto::my_timer.back().get<2>() << endl;
  flush(orz::mr::Femto::file_timing);

  return retval; 
} 
