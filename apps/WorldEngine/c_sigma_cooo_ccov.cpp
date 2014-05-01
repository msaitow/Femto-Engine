                                                                                
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
#include <sci/icmr/Femto/elems/c_sigma_cooo_ccov.h>                                  
                                                                                
using std::cout;                                                                
using std::endl;                                                                
                                                                                
#define FLOPCOUNT                                                               
                                                                                
using orz::mr::Femto::Fc0;                                                      
using orz::mr::Femto::h1_int;                                                   
using orz::mr::Femto::h6_int;                                                   
                                                                                
//  ___________                __               
//  \_   _____/____    _____ _/  |_  ____      
//   |    __)_/ __ \  /     \\   __\/  _ \ 
//   |     \ \  ___/ |  Y Y  \|  | (  <_> )  
//   \___  /  \___  >|__|_|  /|__|  \____/   
//       \/       \/       \/                

//                                   Generated date : Sun Apr 20 10:26:17 2014

                                                                                
// ***************************************************************************  
// orz::mr::mrci                                                                
// ***************************************************************************  
									    /*!        
   @brief CT input                                                              
									     */        
                                                                                
orz::mr::BareAmpPack orz::mr::Femto::sigma_cooo_ccov(const orz::mr::Input &ctinp,                                    
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
  // |-- [    0] --| W0(w,a0) += (    1.00000000) T2(w,c0,v0,a0) Fc1(v0,c0) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) D2(k,j,a0,i) W0(w,a0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W0c_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xc(symblockinfo, sa0));
    FC_FUNC(g_if_sigma_cooo_ccov_no0_x0_type0_noeri,G_IF_SIGMA_COOO_CCOV_NO0_X0_TYPE0_NOERI)
      (sa0, ia0, T2b.cptr(), W0c_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_ccov_no0_x1_type0_noeri,G_IF_SIGMA_COOO_CCOV_NO0_X1_TYPE0_NOERI)
        (sa0, ia0, sj, ij, W0c_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W1(w,a0) += (    1.00000000) T2(w,c0,a0,v0) Fc1(v0,c0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -2.00000000) D2(k,j,a0,i) W1(w,a0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W1ca_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_cooo_ccov_no1_x0_type0_noeri,G_IF_SIGMA_COOO_CCOV_NO1_X0_TYPE0_NOERI)
      (sv0, iv0, T2b.cptr(), W1ca_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccov_no1_x1_type0_noeri,G_IF_SIGMA_COOO_CCOV_NO1_X1_TYPE0_NOERI)
      (sj, ij, W1ca_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W2(w,i) += (    1.00000000) T2(w,c0,i,v0) Fc1(v0,c0) 
  // |-- [    1] --| S2(w,k,i,j) += (    4.00000000) D1(k,j) W2(w,i) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W2ca_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_cooo_ccov_no2_x0_type0_noeri,G_IF_SIGMA_COOO_CCOV_NO2_X0_TYPE0_NOERI)
      (sv0, iv0, T2b.cptr(), W2ca_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccov_no2_x1_type0_noeri,G_IF_SIGMA_COOO_CCOV_NO2_X1_TYPE0_NOERI)
      (sj, ij, W2ca_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W3(w,i) += (    1.00000000) T2(w,c0,v0,i) Fc1(v0,c0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -2.00000000) D1(k,j) W3(w,i) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    T2b = T2.get_amp2(ii);
    orz::DTensor W3c_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xc(symblockinfo, si));
    FC_FUNC(g_if_sigma_cooo_ccov_no3_x0_type0_noeri,G_IF_SIGMA_COOO_CCOV_NO3_X0_TYPE0_NOERI)
      (si, ii, T2b.cptr(), W3c_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_ccov_no3_x1_type0_noeri,G_IF_SIGMA_COOO_CCOV_NO3_X1_TYPE0_NOERI)
        (si, ii, sj, ij, W3c_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[i, "active"] [notNeeded]
  } // End ii
  } // End si
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W4(w,j) += (    1.00000000) T2(c0,w,v0,j) Fc1(v0,c0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -2.00000000) D1(k,i) W4(w,j) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W4c_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xc(symblockinfo, sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccov_no4_x0_type0_noeri,G_IF_SIGMA_COOO_CCOV_NO4_X0_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W4c_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_ccov_no4_x1_type0_noeri,G_IF_SIGMA_COOO_CCOV_NO4_X1_TYPE0_NOERI)
      (sj, ij, W4c_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W5(w,j) += (    1.00000000) T2(w,c0,v0,j) Fc1(v0,c0) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) D1(k,i) W5(w,j) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W5c_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xc(symblockinfo, sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccov_no5_x0_type0_noeri,G_IF_SIGMA_COOO_CCOV_NO5_X0_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W5c_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_ccov_no5_x1_type0_noeri,G_IF_SIGMA_COOO_CCOV_NO5_X1_TYPE0_NOERI)
      (sj, ij, W5c_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
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
  // -- Title : sigma_cooo_ccov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W24(w,j) += (    1.00000000) V2(w,c1,v0,c0) T2(c0,c1,v0,j) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) D1(k,i) W24(w,j) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    double W24_sigma_cooo_ccov(0);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccov_no0_x0_type1_eri_c,G_IF_SIGMA_COOO_CCOV_NO0_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), V2_sym.cptr(), &W24_sigma_cooo_ccov, nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_ccov_no0_x1_type1_eri_c,G_IF_SIGMA_COOO_CCOV_NO0_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, &W24_sigma_cooo_ccov, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W25(w,j) += (    1.00000000) V2(w,c1,v0,c0) T2(c1,c0,v0,j) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) D1(k,i) W25(w,j) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    double W25_sigma_cooo_ccov(0);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccov_no1_x0_type1_eri_c,G_IF_SIGMA_COOO_CCOV_NO1_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), V2_sym.cptr(), &W25_sigma_cooo_ccov, nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_ccov_no1_x1_type1_eri_c,G_IF_SIGMA_COOO_CCOV_NO1_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, &W25_sigma_cooo_ccov, S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
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
  // -- Title : sigma_cooo_ccov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W22(w,a0,j,a1) += (    1.00000000) V2(j,a1,v0,c0) T2(w,c0,v0,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) D2(k,a1,a0,i) W22(w,a0,j,a1) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W22ca_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xca(symblockinfo, sa0^sj));
    FC_FUNC(g_if_sigma_cooo_ccov_no0_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOV_NO0_X0_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, T2b.cptr(), V2_sym.cptr(), W22ca_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_ccov_no0_x1_type1_eri_o,G_IF_SIGMA_COOO_CCOV_NO0_X1_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, W22ca_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W23(w,a0,j,a1) += (    1.00000000) V2(j,a1,v0,c0) T2(w,c0,a0,v0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -2.00000000) D2(k,a1,a0,i) W23(w,a0,j,a1) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W23caa_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_cooo_ccov_no1_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOV_NO1_X0_TYPE1_ERI_O)
      (sj, ij, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W23caa_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_cooo_ccov_no1_x1_type1_eri_o,G_IF_SIGMA_COOO_CCOV_NO1_X1_TYPE1_ERI_O)
    (sj, ij, W23caa_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W30(w,a0,j,a1) += (    1.00000000) V2(j,v0,c0,a1) T2(w,c0,v0,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) D2(k,i,a0,a1) W30(w,a0,j,a1) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W30ca_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xca(symblockinfo, sa0^sj));
    FC_FUNC(g_if_sigma_cooo_ccov_no2_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOV_NO2_X0_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, T2b.cptr(), V2_sym.cptr(), W30ca_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_ccov_no2_x1_type1_eri_o,G_IF_SIGMA_COOO_CCOV_NO2_X1_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, W30ca_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W31(w,a0,j,a1) += (    1.00000000) V2(j,v0,c0,a1) T2(w,c0,a0,v0) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) D2(k,a1,a0,i) W31(w,a0,j,a1) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W31caa_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_cooo_ccov_no3_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOV_NO3_X0_TYPE1_ERI_O)
      (sj, ij, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W31caa_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_cooo_ccov_no3_x1_type1_eri_o,G_IF_SIGMA_COOO_CCOV_NO3_X1_TYPE1_ERI_O)
    (sj, ij, W31caa_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W48(w,j) += (    1.00000000) V2(j,a0,v0,c0) T2(w,c0,a0,v0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -2.00000000) D1(k,i) W48(w,j) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W48c_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xc(symblockinfo, sj));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_cooo_ccov_no4_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOV_NO4_X0_TYPE1_ERI_O)
      (sj, ij, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W48c_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_cooo_ccov_no4_x1_type1_eri_o,G_IF_SIGMA_COOO_CCOV_NO4_X1_TYPE1_ERI_O)
    (sj, ij, W48c_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W49(w,j) += (    1.00000000) V2(j,a0,v0,c0) T2(w,c0,v0,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) D1(k,i) W49(w,j) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W49c_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xc(symblockinfo, sj));
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_cooo_ccov_no5_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOV_NO5_X0_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, T2b.cptr(), V2_sym.cptr(), W49c_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  FC_FUNC(g_if_sigma_cooo_ccov_no5_x1_type1_eri_o,G_IF_SIGMA_COOO_CCOV_NO5_X1_TYPE1_ERI_O)
    (sj, ij, W49c_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W50(w,j) += (    1.00000000) V2(j,v0,c0,a0) T2(w,c0,v0,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -2.00000000) D1(k,i) W50(w,j) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W50c_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xc(symblockinfo, sj));
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_cooo_ccov_no6_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOV_NO6_X0_TYPE1_ERI_O)
      (sa0, ia0, sj, ij, T2b.cptr(), V2_sym.cptr(), W50c_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  FC_FUNC(g_if_sigma_cooo_ccov_no6_x1_type1_eri_o,G_IF_SIGMA_COOO_CCOV_NO6_X1_TYPE1_ERI_O)
    (sj, ij, W50c_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W51(w,j) += (    1.00000000) V2(j,v0,c0,a0) T2(w,c0,a0,v0) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) D1(k,i) W51(w,j) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W51c_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xc(symblockinfo, sj));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_cooo_ccov_no7_x0_type1_eri_o,G_IF_SIGMA_COOO_CCOV_NO7_X0_TYPE1_ERI_O)
      (sj, ij, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W51c_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_cooo_ccov_no7_x1_type1_eri_o,G_IF_SIGMA_COOO_CCOV_NO7_X1_TYPE1_ERI_O)
    (sj, ij, W51c_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
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
  // --  Title : sigma_cooo_ccov
  //  >> Intermediates for the external contractions are defined here << 
  orz::DTensor W6ca_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W7ca_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W8caaa_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W9caaa_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W10caaa_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W11caaa_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W12caaa_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W13caaa_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W14ca_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W15ca_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W20caaa_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W21caaa_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W40ca_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W41ca_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W42ca_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W43ca_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W44ca_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W45ca_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W46ca_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W47ca_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
//-@type(2).declaration(end)

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
  // -- Title : sigma_cooo_ccov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W6(w,a0) += (    1.00000000) V2(v0,c0,w,c1) T2(c0,c1,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  FC_FUNC(g_if_sigma_cooo_ccov_no0_x0_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO0_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), V2_sym.cptr(), W6ca_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W7(w,a0) += (    1.00000000) V2(v0,c0,w,c1) T2(c1,c0,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  FC_FUNC(g_if_sigma_cooo_ccov_no1_x0_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO1_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), V2_sym.cptr(), W7ca_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W8(w,a0,a2,a1) += (    1.00000000) V2(v0,a2,c0,a1) T2(c0,w,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  FC_FUNC(g_if_sigma_cooo_ccov_no2_x0_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO2_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), V2_sym.cptr(), W8caaa_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W9(w,a0,a2,a1) += (    1.00000000) V2(v0,a2,c0,a1) T2(w,c0,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  FC_FUNC(g_if_sigma_cooo_ccov_no3_x0_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO3_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), V2_sym.cptr(), W9caaa_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W10(w,a0,a2,a1) += (    1.00000000) V2(v0,c0,a2,a1) T2(c0,w,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  FC_FUNC(g_if_sigma_cooo_ccov_no4_x0_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO4_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), V2_sym.cptr(), W10caaa_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W11(w,a0,a2,a1) += (    1.00000000) V2(v0,c0,a2,a1) T2(w,c0,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  FC_FUNC(g_if_sigma_cooo_ccov_no5_x0_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO5_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), V2_sym.cptr(), W11caaa_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W12(w,a0,i,a1) += (    1.00000000) V2(v0,c0,i,a1) T2(c0,w,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  FC_FUNC(g_if_sigma_cooo_ccov_no6_x0_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO6_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), V2_sym.cptr(), W12caaa_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W13(w,a0,i,a1) += (    1.00000000) V2(v0,c0,i,a1) T2(w,c0,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  FC_FUNC(g_if_sigma_cooo_ccov_no7_x0_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO7_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), V2_sym.cptr(), W13caaa_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| W14(w,i) += (    1.00000000) V2(v0,c0,w,c1) T2(c1,c0,i,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  FC_FUNC(g_if_sigma_cooo_ccov_no8_x0_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO8_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), V2_sym.cptr(), W14ca_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W15(w,i) += (    1.00000000) V2(v0,c0,w,c1) T2(c0,c1,i,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  FC_FUNC(g_if_sigma_cooo_ccov_no9_x0_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO9_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), V2_sym.cptr(), W15ca_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| W16(c0,k,j,v0) += (    1.00000000) V2(v0,a1,c0,a0) D2(k,j,a1,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -2.00000000) T2(w,c0,i,v0) W16(c0,k,j,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W16ca_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_cooo_ccov_no10_x0_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO10_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W16ca_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_ccov_no10_x1_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO10_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W16ca_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| W17(c0,k,j,v0) += (    1.00000000) V2(v0,a1,c0,a0) D2(k,j,a1,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) T2(c0,w,i,v0) W17(c0,k,j,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W17ca_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_cooo_ccov_no11_x0_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO11_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W17ca_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_ccov_no11_x1_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO11_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W17ca_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| W18(c0,k,j,v0) += (    1.00000000) V2(v0,c0,a1,a0) D2(k,j,a1,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    4.00000000) T2(w,c0,i,v0) W18(c0,k,j,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W18ca_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_cooo_ccov_no12_x0_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO12_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W18ca_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_ccov_no12_x1_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO12_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W18ca_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| W19(c0,k,j,v0) += (    1.00000000) V2(v0,c0,a1,a0) D2(k,j,a1,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -2.00000000) T2(c0,w,i,v0) W19(c0,k,j,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W19ca_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_cooo_ccov_no13_x0_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO13_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W19ca_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_ccov_no13_x1_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO13_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W19ca_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   14] -- 
  // |-- [    0] --| W20(w,a0,i,a1) += (    1.00000000) V2(v0,i,c0,a1) T2(c0,w,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  FC_FUNC(g_if_sigma_cooo_ccov_no14_x0_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO14_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), V2_sym.cptr(), W20caaa_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   15] -- 
  // |-- [    0] --| W21(w,a0,i,a1) += (    1.00000000) V2(v0,i,c0,a1) T2(w,c0,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  FC_FUNC(g_if_sigma_cooo_ccov_no15_x0_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO15_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), V2_sym.cptr(), W21caaa_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   16] -- 
  // |-- [    0] --| W26(c0,k,i,v0) += (    1.00000000) V2(v0,a1,c0,a0) D2(k,i,a1,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) T2(w,c0,j,v0) W26(c0,k,i,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W26caa_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ccov_no16_x0_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO16_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W26caa_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccov_no16_x1_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO16_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W26caa_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   17] -- 
  // |-- [    0] --| W27(c0,k,i,v0) += (    1.00000000) V2(v0,a1,c0,a0) D2(k,a0,a1,i) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) T2(w,c0,v0,j) W27(c0,k,i,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W27caa_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ccov_no17_x0_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO17_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W27caa_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccov_no17_x1_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO17_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W27caa_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   18] -- 
  // |-- [    0] --| W28(c0,k,i,v0) += (    1.00000000) V2(v0,c0,a1,a0) D2(k,i,a1,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -2.00000000) T2(w,c0,j,v0) W28(c0,k,i,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W28caa_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ccov_no18_x0_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO18_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W28caa_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccov_no18_x1_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO18_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W28caa_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   19] -- 
  // |-- [    0] --| W29(c0,k,i,v0) += (    1.00000000) V2(v0,c0,a1,a0) D2(k,i,a1,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) T2(w,c0,v0,j) W29(c0,k,i,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W29caa_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ccov_no19_x0_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO19_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W29caa_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccov_no19_x1_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO19_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W29caa_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   20] -- 
  // |-- [    0] --| W32(c0,k,i,v0) += (    1.00000000) V2(v0,c0,i,a0) D1(k,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -2.00000000) T2(w,c0,j,v0) W32(c0,k,i,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W32caa_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ccov_no20_x0_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO20_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W32caa_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccov_no20_x1_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO20_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W32caa_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   21] -- 
  // |-- [    0] --| W33(c0,k,i,v0) += (    1.00000000) V2(v0,c0,i,a0) D1(k,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) T2(w,c0,v0,j) W33(c0,k,i,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W33caa_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ccov_no21_x0_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO21_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W33caa_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccov_no21_x1_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO21_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W33caa_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   22] -- 
  // |-- [    0] --| W34(c0,k,j,v0) += (    1.00000000) V2(v0,c0,j,a0) D1(k,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    4.00000000) T2(w,c0,i,v0) W34(c0,k,j,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W34ca_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_cooo_ccov_no22_x0_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO22_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W34ca_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_ccov_no22_x1_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO22_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W34ca_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   23] -- 
  // |-- [    0] --| W35(c0,k,j,v0) += (    1.00000000) V2(v0,c0,j,a0) D1(k,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -2.00000000) T2(c0,w,i,v0) W35(c0,k,j,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W35ca_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_cooo_ccov_no23_x0_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO23_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W35ca_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_ccov_no23_x1_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO23_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W35ca_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   24] -- 
  // |-- [    0] --| W36(c0,k,j,v0) += (    1.00000000) V2(v0,j,c0,a0) D1(k,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -2.00000000) T2(w,c0,i,v0) W36(c0,k,j,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W36ca_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_cooo_ccov_no24_x0_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO24_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W36ca_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_ccov_no24_x1_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO24_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W36ca_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   25] -- 
  // |-- [    0] --| W37(c0,k,j,v0) += (    1.00000000) V2(v0,j,c0,a0) D1(k,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) T2(c0,w,i,v0) W37(c0,k,j,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W37ca_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_cooo_ccov_no25_x0_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO25_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W37ca_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_ccov_no25_x1_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO25_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W37ca_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   26] -- 
  // |-- [    0] --| W38(c0,k,i,v0) += (    1.00000000) V2(v0,i,c0,a0) D1(k,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -2.00000000) T2(w,c0,v0,j) W38(c0,k,i,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W38caa_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ccov_no26_x0_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO26_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W38caa_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccov_no26_x1_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO26_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W38caa_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   27] -- 
  // |-- [    0] --| W39(c0,k,i,v0) += (    1.00000000) V2(v0,i,c0,a0) D1(k,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) T2(w,c0,j,v0) W39(c0,k,i,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W39caa_sigma_cooo_ccov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ccov_no27_x0_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO27_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W39caa_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccov_no27_x1_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO27_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W39caa_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   28] -- 
  // |-- [    0] --| W40(w,a1) += (    1.00000000) V2(v0,a1,c0,a0) T2(w,c0,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  FC_FUNC(g_if_sigma_cooo_ccov_no28_x0_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO28_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), V2_sym.cptr(), W40ca_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   29] -- 
  // |-- [    0] --| W41(w,a1) += (    1.00000000) V2(v0,a1,c0,a0) T2(c0,w,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  FC_FUNC(g_if_sigma_cooo_ccov_no29_x0_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO29_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), V2_sym.cptr(), W41ca_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   30] -- 
  // |-- [    0] --| W42(w,a1) += (    1.00000000) V2(v0,c0,a1,a0) T2(w,c0,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  FC_FUNC(g_if_sigma_cooo_ccov_no30_x0_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO30_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), V2_sym.cptr(), W42ca_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   31] -- 
  // |-- [    0] --| W43(w,a1) += (    1.00000000) V2(v0,c0,a1,a0) T2(c0,w,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  FC_FUNC(g_if_sigma_cooo_ccov_no31_x0_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO31_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), V2_sym.cptr(), W43ca_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   32] -- 
  // |-- [    0] --| W44(w,i) += (    1.00000000) V2(v0,c0,i,a0) T2(w,c0,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  FC_FUNC(g_if_sigma_cooo_ccov_no32_x0_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO32_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), V2_sym.cptr(), W44ca_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   33] -- 
  // |-- [    0] --| W45(w,i) += (    1.00000000) V2(v0,c0,i,a0) T2(c0,w,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  FC_FUNC(g_if_sigma_cooo_ccov_no33_x0_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO33_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), V2_sym.cptr(), W45ca_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   34] -- 
  // |-- [    0] --| W46(w,i) += (    1.00000000) V2(v0,i,c0,a0) T2(c0,w,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  FC_FUNC(g_if_sigma_cooo_ccov_no34_x0_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO34_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), V2_sym.cptr(), W46ca_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   35] -- 
  // |-- [    0] --| W47(w,i) += (    1.00000000) V2(v0,i,c0,a0) T2(w,c0,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  FC_FUNC(g_if_sigma_cooo_ccov_no35_x0_type1_eri_v,G_IF_SIGMA_COOO_CCOV_NO35_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), V2_sym.cptr(), W47ca_sigma_cooo_ccov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

//-@type(1).contraction(end)

  } // End myrank
  }
  }
  orz::world().barrier();

//-@type(2).contraction(begin)
  // -- Title : sigma_cooo_ccov
  //*-- Entering to take the type 2 contractions --*//

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -1.00000000) D2(k,j,a0,i) W6(w,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccov_no0_x0_type2_eri_v,G_IF_SIGMA_COOO_CCOV_NO0_X0_TYPE2_ERI_V)
      (sj, ij, W6ca_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    2.00000000) D2(k,j,a0,i) W7(w,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccov_no1_x0_type2_eri_v,G_IF_SIGMA_COOO_CCOV_NO1_X0_TYPE2_ERI_V)
      (sj, ij, W7ca_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    1.00000000) D3(k,j,a2,i,a0,a1) W8(w,a0,a2,a1) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccov_no2_x0_type2_eri_v,G_IF_SIGMA_COOO_CCOV_NO2_X0_TYPE2_ERI_V)
      (sj, ij, W8caaa_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    1.00000000) D3(k,j,a2,a1,a0,i) W9(w,a0,a2,a1) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccov_no3_x0_type2_eri_v,G_IF_SIGMA_COOO_CCOV_NO3_X0_TYPE2_ERI_V)
      (sj, ij, W9caaa_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    1.00000000) D3(k,j,a2,a1,a0,i) W10(w,a0,a2,a1) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccov_no4_x0_type2_eri_v,G_IF_SIGMA_COOO_CCOV_NO4_X0_TYPE2_ERI_V)
      (sj, ij, W10caaa_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -2.00000000) D3(k,j,a2,a1,a0,i) W11(w,a0,a2,a1) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccov_no5_x0_type2_eri_v,G_IF_SIGMA_COOO_CCOV_NO5_X0_TYPE2_ERI_V)
      (sj, ij, W11caaa_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    1.00000000) D2(k,j,a0,a1) W12(w,a0,i,a1) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccov_no6_x0_type2_eri_v,G_IF_SIGMA_COOO_CCOV_NO6_X0_TYPE2_ERI_V)
      (sj, ij, W12caaa_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -2.00000000) D2(k,j,a0,a1) W13(w,a0,i,a1) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccov_no7_x0_type2_eri_v,G_IF_SIGMA_COOO_CCOV_NO7_X0_TYPE2_ERI_V)
      (sj, ij, W13caaa_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -4.00000000) D1(k,j) W14(w,i) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccov_no8_x0_type2_eri_v,G_IF_SIGMA_COOO_CCOV_NO8_X0_TYPE2_ERI_V)
      (sj, ij, W14ca_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    2.00000000) D1(k,j) W15(w,i) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccov_no9_x0_type2_eri_v,G_IF_SIGMA_COOO_CCOV_NO9_X0_TYPE2_ERI_V)
      (sj, ij, W15ca_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -2.00000000) D2(k,j,a0,a1) W20(w,a0,i,a1) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccov_no10_x0_type2_eri_v,G_IF_SIGMA_COOO_CCOV_NO10_X0_TYPE2_ERI_V)
      (sj, ij, W20caaa_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    1.00000000) D2(k,j,a0,a1) W21(w,a0,i,a1) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccov_no11_x0_type2_eri_v,G_IF_SIGMA_COOO_CCOV_NO11_X0_TYPE2_ERI_V)
      (sj, ij, W21caaa_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    1.00000000) D2(k,j,a1,i) W40(w,a1) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccov_no12_x0_type2_eri_v,G_IF_SIGMA_COOO_CCOV_NO12_X0_TYPE2_ERI_V)
      (sj, ij, W40ca_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -2.00000000) D2(k,j,a1,i) W41(w,a1) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccov_no13_x0_type2_eri_v,G_IF_SIGMA_COOO_CCOV_NO13_X0_TYPE2_ERI_V)
      (sj, ij, W41ca_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   14] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -2.00000000) D2(k,j,a1,i) W42(w,a1) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccov_no14_x0_type2_eri_v,G_IF_SIGMA_COOO_CCOV_NO14_X0_TYPE2_ERI_V)
      (sj, ij, W42ca_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   15] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    1.00000000) D2(k,j,a1,i) W43(w,a1) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccov_no15_x0_type2_eri_v,G_IF_SIGMA_COOO_CCOV_NO15_X0_TYPE2_ERI_V)
      (sj, ij, W43ca_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   16] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    4.00000000) D1(k,j) W44(w,i) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccov_no16_x0_type2_eri_v,G_IF_SIGMA_COOO_CCOV_NO16_X0_TYPE2_ERI_V)
      (sj, ij, W44ca_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   17] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -2.00000000) D1(k,j) W45(w,i) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccov_no17_x0_type2_eri_v,G_IF_SIGMA_COOO_CCOV_NO17_X0_TYPE2_ERI_V)
      (sj, ij, W45ca_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   18] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    4.00000000) D1(k,j) W46(w,i) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccov_no18_x0_type2_eri_v,G_IF_SIGMA_COOO_CCOV_NO18_X0_TYPE2_ERI_V)
      (sj, ij, W46ca_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   19] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -2.00000000) D1(k,j) W47(w,i) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ccov_no19_x0_type2_eri_v,G_IF_SIGMA_COOO_CCOV_NO19_X0_TYPE2_ERI_V)
      (sj, ij, W47ca_sigma_cooo_ccov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
//-@type(2).contraction(end)
  orz::world().barrier();

  } // End Femto
  //*-- FEMTO ends --//*

//-@loadERI(v,end)

//-@ERI.contractions(end)

//-@D4C.contractions(begin)

//-@D4C.contractions(end)

  // Do timing!
  orz::mr::Femto::my_timer.push_back(boost::make_tuple("sigma_cooo_ccov", time_sigma.elapsed_cputime(), time_sigma.elapsed_wallclocktime()));
  orz::mr::Femto::file_timing << "* " << boost::format("%20s : %10.7f %10.7f ") % orz::mr::Femto::my_timer.back().get<0>() % orz::mr::Femto::my_timer.back().get<1>() % orz::mr::Femto::my_timer.back().get<2>() << endl;
  flush(orz::mr::Femto::file_timing);

  return retval; 
} 
