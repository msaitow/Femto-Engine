                                                                                
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
#include <sci/icmr/Femto/elems/c_sigma_cooo_ooov.h>                                  
                                                                                
using std::cout;                                                                
using std::endl;                                                                
                                                                                
#define FLOPCOUNT                                                               
                                                                                
using orz::mr::Femto::Fc0;                                                      
using orz::mr::Femto::h1_int;                                                   
using orz::mr::Femto::h6_int;                                                   
                                                                                
//  8888888888                     888                  
//  888                            888                  
//  888                            888                  
//  8888888  .d88b.  88888b.d88b.  888888  .d88b.       
//  888     d8P  Y8b 888 "888 "88b 888    d88""88b  
//  888     88888888 888  888  888 888    888  888      
//  888     Y8b.     888  888  888 Y88b.  Y88..88P      
//  888      "Y8888  888  888  888  "Y888  "Y88P"   

//                                   Generated date : Sun Apr 20 10:26:13 2014

                                                                                
// ***************************************************************************  
// orz::mr::mrci                                                                
// ***************************************************************************  
									    /*!        
   @brief CT input                                                              
									     */        
                                                                                
orz::mr::BareAmpPack orz::mr::Femto::sigma_cooo_ooov(const orz::mr::Input &ctinp,                                    
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
  // -- Title : sigma_cooo_ooov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W1(w,a1,a0,i,a3,a2) += (    1.00000000) V2(w,a3,v0,a2) T2(a1,a0,i,v0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) D3(k,j,a3,a1,a2,a0) W1(w,a1,a0,i,a3,a2) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W1aaaaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaaaa(symblockinfo, sw));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_cooo_ooov_no0_x0_type1_eri_c,G_IF_SIGMA_COOO_OOOV_NO0_X0_TYPE1_ERI_C)
      (sv0, iv0, sw, iw, T2b.cptr(), V2_sym.cptr(), W1aaaaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no0_x1_type1_eri_c,G_IF_SIGMA_COOO_OOOV_NO0_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W1aaaaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W4(w,a1,a0,j,a3,a2) += (    1.00000000) V2(w,a3,v0,a2) T2(a0,a1,v0,j) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) D3(k,a1,a3,i,a2,a0) W4(w,a1,a0,j,a3,a2) 
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
    orz::DTensor W4aaaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sw^sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no1_x0_type1_eri_c,G_IF_SIGMA_COOO_OOOV_NO1_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), V2_sym.cptr(), W4aaaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_ooov_no1_x1_type1_eri_c,G_IF_SIGMA_COOO_OOOV_NO1_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W4aaaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W14(w,a1,a0,a2,a4,a3) += (    1.00000000) V2(w,a4,v0,a3) T2(a1,a0,v0,a2) 
  // |-- [    1] --| W15(w,j,a1,a0,a4,a3) += (    1.00000000) D1(j,a2) W14(w,a1,a0,a2,a4,a3) 
  // |-- [    2] --| S2(w,k,i,j) += (    0.50000000) D3(a0,k,a1,a3,i,a4) W15(w,j,a1,a0,a4,a3) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 1 <3, 2> 
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 2
  // Type1 anormality found for external indices of the contraction >> 0 <<
  orz::DTensor W14aaaaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaaaa(symblockinfo, sw));
  for(int sa2 = 0;sa2 < nir;++sa2){ 
  for(int ia2 = symblockinfo.psym()(sa2,I_O,I_BEGIN);ia2 <= symblockinfo.psym()(sa2,I_O,I_END);++ia2){ 
    T2b = T2.get_amp2(ia2);
    FC_FUNC(g_if_sigma_cooo_ooov_no2_x0_type1_eri_c,G_IF_SIGMA_COOO_OOOV_NO2_X0_TYPE1_ERI_C)
      (sa2, ia2, sw, iw, T2b.cptr(), V2_sym.cptr(), W14aaaaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a2, "active"] [notNeeded]
  } // End ia2
  } // End sa2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W15aaaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sw^sj));
    FC_FUNC(g_if_sigma_cooo_ooov_no2_x1_type1_eri_c,G_IF_SIGMA_COOO_OOOV_NO2_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W14aaaaa_sigma_cooo_ooov.cptr(), W15aaaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no2_x2_type1_eri_c,G_IF_SIGMA_COOO_OOOV_NO2_X2_TYPE1_ERI_C)
      (sj, ij, sw, iw, W15aaaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W22(w,a1,a0,a2,a4,a3) += (    1.00000000) V2(w,a4,v0,a3) T2(a1,a0,v0,a2) 
  // |-- [    1] --| W23(w,i,a1,a0,a4,a3) += (    1.00000000) D1(i,a2) W22(w,a1,a0,a2,a4,a3) 
  // |-- [    2] --| S2(w,k,i,j) += (    0.50000000) D3(a0,a4,a1,a3,j,k) W23(w,i,a1,a0,a4,a3) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 1 <2, 1> 
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 2
  orz::DTensor W23aaaaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaaaa(symblockinfo, sw));
  for(int sa2 = 0;sa2 < nir;++sa2){ 
  for(int ia2 = symblockinfo.psym()(sa2,I_O,I_BEGIN);ia2 <= symblockinfo.psym()(sa2,I_O,I_END);++ia2){ 
    T2b = T2.get_amp2(ia2);
    orz::DTensor W22aaaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sw^sa2));
    FC_FUNC(g_if_sigma_cooo_ooov_no3_x0_type1_eri_c,G_IF_SIGMA_COOO_OOOV_NO3_X0_TYPE1_ERI_C)
      (sa2, ia2, sw, iw, T2b.cptr(), V2_sym.cptr(), W22aaaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_ooov_no3_x1_type1_eri_c,G_IF_SIGMA_COOO_OOOV_NO3_X1_TYPE1_ERI_C)
      (sa2, ia2, sw, iw, W22aaaa_sigma_cooo_ooov.cptr(), W23aaaaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  // --> @[a2, "active"] [notNeeded]
  } // End ia2
  } // End sa2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no3_x2_type1_eri_c,G_IF_SIGMA_COOO_OOOV_NO3_X2_TYPE1_ERI_C)
      (sj, ij, sw, iw, W23aaaaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W34(w,a1,a0,a2,a4,a3) += (    1.00000000) V2(w,a4,v0,a3) T2(a1,a0,v0,a2) 
  // |-- [    1] --| W35(w,k,a0,a2,a4,a3) += (    1.00000000) D1(a1,k) W34(w,a1,a0,a2,a4,a3) 
  // |-- [    2] --| S2(w,k,i,j) += (    0.50000000) D3(a0,a2,i,a4,j,a3) W35(w,k,a0,a2,a4,a3) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 0
  for(int sa2 = 0;sa2 < nir;++sa2){ 
  for(int ia2 = symblockinfo.psym()(sa2,I_O,I_BEGIN);ia2 <= symblockinfo.psym()(sa2,I_O,I_END);++ia2){ 
    T2b = T2.get_amp2(ia2);
    orz::DTensor W34aaaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sw^sa2));
    orz::DTensor W35aaaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sw^sa2));
    FC_FUNC(g_if_sigma_cooo_ooov_no4_x0_type1_eri_c,G_IF_SIGMA_COOO_OOOV_NO4_X0_TYPE1_ERI_C)
      (sa2, ia2, sw, iw, T2b.cptr(), V2_sym.cptr(), W34aaaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_ooov_no4_x1_type1_eri_c,G_IF_SIGMA_COOO_OOOV_NO4_X1_TYPE1_ERI_C)
      (sa2, ia2, sw, iw, W34aaaa_sigma_cooo_ooov.cptr(), W35aaaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_ooov_no4_x2_type1_eri_c,G_IF_SIGMA_COOO_OOOV_NO4_X2_TYPE1_ERI_C)
        (sa2, ia2, sj, ij, sw, iw, W35aaaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a2, "active"] [notNeeded]
  } // End ia2
  } // End sa2
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W42(w,a1,a0,a2,a4,a3) += (    1.00000000) V2(w,a4,v0,a3) T2(a1,a0,v0,a2) 
  // |-- [    1] --| W43(w,k,a1,a2,a4,a3) += (    1.00000000) D1(a0,k) W42(w,a1,a0,a2,a4,a3) 
  // |-- [    2] --| S2(w,k,i,j) += (    0.50000000) D3(a1,a3,i,a4,j,a2) W43(w,k,a1,a2,a4,a3) 
  int sw(s_eri);
  int iw(i_eri);
  double flops(0); // Flop count  
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 0
  for(int sa2 = 0;sa2 < nir;++sa2){ 
  for(int ia2 = symblockinfo.psym()(sa2,I_O,I_BEGIN);ia2 <= symblockinfo.psym()(sa2,I_O,I_END);++ia2){ 
    T2b = T2.get_amp2(ia2);
    orz::DTensor W42aaaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sw^sa2));
    orz::DTensor W43aaaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaaa(symblockinfo, sw^sa2));
    FC_FUNC(g_if_sigma_cooo_ooov_no5_x0_type1_eri_c,G_IF_SIGMA_COOO_OOOV_NO5_X0_TYPE1_ERI_C)
      (sa2, ia2, sw, iw, T2b.cptr(), V2_sym.cptr(), W42aaaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_ooov_no5_x1_type1_eri_c,G_IF_SIGMA_COOO_OOOV_NO5_X1_TYPE1_ERI_C)
      (sa2, ia2, sw, iw, W42aaaa_sigma_cooo_ooov.cptr(), W43aaaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_ooov_no5_x2_type1_eri_c,G_IF_SIGMA_COOO_OOOV_NO5_X2_TYPE1_ERI_C)
        (sa2, ia2, sj, ij, sw, iw, W43aaaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[a2, "active"] [notNeeded]
  } // End ia2
  } // End sa2
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

//-@type(2).declaration(begin)
  // --  Title : sigma_cooo_ooov
  //  >> Intermediates for the external contractions are defined here << 
  orz::DTensor W13ca_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W27ca_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W29caaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W31caaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W33caaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W37caaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W39caaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W41caaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W45ca_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W47ca_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W49caaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W51ca_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W53caaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W55caaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W57caaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W59caaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W61caaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W63caaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W65caaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W71caaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W73caaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W75caaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W77caaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W79caaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W83caaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W85caaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W87caaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W89caaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W99ca_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W101ca_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W103caaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W105ca_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W107caaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W109caaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W111caaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W113caaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W115caaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W117caaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W119caaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W125caaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W127caaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W129caaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W131caaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W133caaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W137caaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W139caaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W141caaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W143caaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
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
  // -- Title : sigma_cooo_ooov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W0(k,a3,j,v0) += (    1.00000000) T2(a0,a1,a2,v0) D3(k,j,a3,a1,a2,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) V2(v0,a3,w,i) W0(k,a3,j,v0) 
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
    orz::DTensor W0aa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_cooo_ooov_no0_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO0_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W0aa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_ooov_no0_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO0_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W0aa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W2(k,a3,j,v0) += (    1.00000000) T2(a0,a1,a2,v0) D3(k,j,a3,a1,a2,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) V2(v0,i,w,a3) W2(k,a3,j,v0) 
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
    orz::DTensor W2aa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_cooo_ooov_no1_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO1_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W2aa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_ooov_no1_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO1_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W2aa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W3(k,a3,i,v0) += (    1.00000000) T2(a0,a1,a2,v0) D3(k,i,a3,a1,a2,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) V2(v0,a3,w,j) W3(k,a3,i,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W3aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no2_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO2_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W3aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no2_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO2_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W3aaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W5(k,a3,i,v0) += (    1.00000000) T2(a0,a1,a2,v0) D3(k,a1,a3,i,a2,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) V2(v0,j,w,a3) W5(k,a3,i,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W5aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no3_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO3_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W5aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no3_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO3_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W5aaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W6(k,a2,j,v0) += (    1.00000000) T2(a1,a0,j,v0) D2(k,a1,a2,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) V2(v0,a2,w,i) W6(k,a2,j,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W6aa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_cooo_ooov_no4_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO4_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W6aa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no4_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO4_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W6aa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W7(k,v0) += (    1.00000000) T2(a0,a1,a2,v0) D2(k,a1,a2,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) V2(v0,j,w,i) W7(k,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W7a_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no5_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO5_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W7a_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no5_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO5_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W7a_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W8(k,a2,i,v0) += (    1.00000000) T2(a1,a0,i,v0) D2(k,a1,a2,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) V2(v0,a2,w,j) W8(k,a2,i,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W8aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no6_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO6_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W8aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no6_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO6_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W8aaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W9(k,a2,i,v0) += (    1.00000000) T2(a0,a1,i,v0) D2(k,a1,a2,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) V2(v0,j,w,a2) W9(k,a2,i,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W9aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no7_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO7_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W9aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no7_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO7_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W9aaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| W10(k,v0) += (    1.00000000) T2(a0,a1,a2,v0) D2(k,a1,a2,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) V2(v0,i,w,j) W10(k,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W10a_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no8_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO8_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W10a_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no8_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO8_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W10a_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W11(k,a2,j,v0) += (    1.00000000) T2(a1,a0,j,v0) D2(k,a1,a2,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) V2(v0,i,w,a2) W11(k,a2,j,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W11aa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_cooo_ooov_no9_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO9_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W11aa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no9_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO9_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W11aa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| W12(i,a3,a4,v0) += (    1.00000000) T2(a0,a1,a2,v0) D3(a0,a2,a1,a3,i,a4) 
  // |-- [    1] --| W13(w,i) += (    1.00000000) V2(v0,a3,w,a4) W12(i,a3,a4,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W12aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no10_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO10_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W12aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no10_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO10_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W12aaa_sigma_cooo_ooov.cptr(), W13ca_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| W17(i,k,a4,v0) += (    1.00000000) T2(a0,a1,a2,v0) D3(a0,a2,a1,k,i,a4) 
  // |-- [    1] --| W16(w,j,a4,v0) += (    1.00000000) V2(v0,a3,w,a4) D1(j,a3) 
  // |-- [    2] --| S2(w,k,i,j) += (    0.50000000) W16(w,j,a4,v0) W17(i,k,a4,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W17aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no11_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO11_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W17aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W16ca_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_cooo_ooov_no11_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO11_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W16ca_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no11_x2_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO11_X2_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, W16ca_sigma_cooo_ooov.cptr(), W17aaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| W19(i,a3,k,v0) += (    1.00000000) T2(a0,a1,a2,v0) D3(a0,a2,a1,a3,i,k) 
  // |-- [    1] --| W18(w,j,a3,v0) += (    1.00000000) V2(v0,a3,w,a4) D1(j,a4) 
  // |-- [    2] --| S2(w,k,i,j) += (    0.50000000) W18(w,j,a3,v0) W19(i,a3,k,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W19aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no12_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO12_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W19aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W18ca_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_cooo_ooov_no12_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO12_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W18ca_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no12_x2_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO12_X2_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, W18ca_sigma_cooo_ooov.cptr(), W19aaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| W20(w,i,a3,v0) += (    1.00000000) V2(v0,a3,w,a4) D1(i,a4) 
  // |-- [    1] --| W21(j,a3,k,v0) += (    1.00000000) T2(a0,a1,a2,v0) D3(a0,a2,a1,a3,j,k) 
  // |-- [    2] --| S2(w,k,i,j) += (   -1.00000000) W20(w,i,a3,v0) W21(j,a3,k,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W20caa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no13_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO13_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W20caa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W21aa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_cooo_ooov_no13_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO13_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W21aa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_ooov_no13_x2_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO13_X2_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, W20caa_sigma_cooo_ooov.cptr(), W21aa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   14] -- 
  // |-- [    0] --| W24(w,i,a4,v0) += (    1.00000000) V2(v0,a3,w,a4) D1(i,a3) 
  // |-- [    1] --| W25(j,a4,k,v0) += (    1.00000000) T2(a0,a1,a2,v0) D3(a0,a2,a1,a4,j,k) 
  // |-- [    2] --| S2(w,k,i,j) += (    0.50000000) W24(w,i,a4,v0) W25(j,a4,k,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W24caa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no14_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO14_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W24caa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W25aa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_cooo_ooov_no14_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO14_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W25aa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_ooov_no14_x2_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO14_X2_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, W24caa_sigma_cooo_ooov.cptr(), W25aa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   15] -- 
  // |-- [    0] --| W26(j,a3,a4,v0) += (    1.00000000) T2(a0,a1,a2,v0) D3(a0,a2,a1,a3,j,a4) 
  // |-- [    1] --| W27(w,j) += (    1.00000000) V2(v0,a3,w,a4) W26(j,a3,a4,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W26aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no15_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO15_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W26aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no15_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO15_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W26aaa_sigma_cooo_ooov.cptr(), W27ca_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   16] -- 
  // |-- [    0] --| W28(w,a1,a4,v0) += (    1.00000000) V2(v0,a3,w,a4) D1(a1,a3) 
  // |-- [    1] --| W29(w,a0,a2,a4) += (    1.00000000) T2(a0,a1,a2,v0) W28(w,a1,a4,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W28caa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no16_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO16_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W28caa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no16_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO16_X1_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W28caa_sigma_cooo_ooov.cptr(), W29caaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   17] -- 
  // |-- [    0] --| W30(a0,v0) += (    1.00000000) T2(a0,a1,a2,v0) D1(a1,a2) 
  // |-- [    1] --| W31(w,a4,a3,a0) += (    1.00000000) V2(v0,a3,w,a4) W30(a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W30a_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no17_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO17_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W30a_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no17_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO17_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W30a_sigma_cooo_ooov.cptr(), W31caaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   18] -- 
  // |-- [    0] --| W32(w,a1,a3,v0) += (    1.00000000) V2(v0,a3,w,a4) D1(a1,a4) 
  // |-- [    1] --| W33(w,a0,a2,a3) += (    1.00000000) T2(a0,a1,a2,v0) W32(w,a1,a3,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W32caa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no18_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO18_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W32caa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no18_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO18_X1_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W32caa_sigma_cooo_ooov.cptr(), W33caaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   19] -- 
  // |-- [    0] --| W36(a1,v0) += (    1.00000000) T2(a0,a1,a2,v0) D1(a0,a2) 
  // |-- [    1] --| W37(w,a4,a3,a1) += (    1.00000000) V2(v0,a3,w,a4) W36(a1,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W36a_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no19_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO19_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W36a_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no19_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO19_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W36a_sigma_cooo_ooov.cptr(), W37caaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   20] -- 
  // |-- [    0] --| W38(w,a0,a4,v0) += (    1.00000000) V2(v0,a3,w,a4) D1(a0,a3) 
  // |-- [    1] --| W39(w,a1,a2,a4) += (    1.00000000) T2(a0,a1,a2,v0) W38(w,a0,a4,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W38caa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no20_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO20_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W38caa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no20_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO20_X1_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W38caa_sigma_cooo_ooov.cptr(), W39caaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   21] -- 
  // |-- [    0] --| W40(w,a0,a3,v0) += (    1.00000000) V2(v0,a3,w,a4) D1(a0,a4) 
  // |-- [    1] --| W41(w,a1,a2,a3) += (    1.00000000) T2(a0,a1,a2,v0) W40(w,a0,a3,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W40caa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no21_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO21_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W40caa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no21_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO21_X1_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W40caa_sigma_cooo_ooov.cptr(), W41caaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   22] -- 
  // |-- [    0] --| W44(a3,v0) += (    1.00000000) T2(a0,a1,a2,v0) D2(a0,a2,a1,a3) 
  // |-- [    1] --| W45(w,a4) += (    1.00000000) V2(v0,a3,w,a4) W44(a3,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W44a_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no22_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO22_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W44a_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no22_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO22_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W44a_sigma_cooo_ooov.cptr(), W45ca_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   23] -- 
  // |-- [    0] --| W46(a4,v0) += (    1.00000000) T2(a0,a1,a2,v0) D2(a0,a2,a1,a4) 
  // |-- [    1] --| W47(w,a3) += (    1.00000000) V2(v0,a3,w,a4) W46(a4,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W46a_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no23_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO23_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W46a_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no23_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO23_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W46a_sigma_cooo_ooov.cptr(), W47ca_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   24] -- 
  // |-- [    0] --| W48(k,v0) += (    1.00000000) T2(a0,a1,a2,v0) D2(a0,a2,a1,k) 
  // |-- [    1] --| W49(w,a4,a3,k) += (    1.00000000) V2(v0,a3,w,a4) W48(k,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W48a_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no24_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO24_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W48a_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no24_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO24_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W48a_sigma_cooo_ooov.cptr(), W49caaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   25] -- 
  // |-- [    0] --| W50(a4,a3,a2,v0) += (    1.00000000) T2(a0,a1,a2,v0) D2(a0,a4,a1,a3) 
  // |-- [    1] --| W51(w,a2) += (    1.00000000) V2(v0,a3,w,a4) W50(a4,a3,a2,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W50aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no25_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO25_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W50aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no25_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO25_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W50aaa_sigma_cooo_ooov.cptr(), W51ca_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   26] -- 
  // |-- [    0] --| W52(k,a3,a2,v0) += (    1.00000000) T2(a0,a1,a2,v0) D2(a0,k,a1,a3) 
  // |-- [    1] --| W53(w,a4,k,a2) += (    1.00000000) V2(v0,a3,w,a4) W52(k,a3,a2,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W52aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no26_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO26_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W52aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no26_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO26_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W52aaa_sigma_cooo_ooov.cptr(), W53caaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   27] -- 
  // |-- [    0] --| W54(a4,k,a2,v0) += (    1.00000000) T2(a0,a1,a2,v0) D2(a0,a4,a1,k) 
  // |-- [    1] --| W55(w,a3,k,a2) += (    1.00000000) V2(v0,a3,w,a4) W54(a4,k,a2,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W54aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no27_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO27_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W54aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no27_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO27_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W54aaa_sigma_cooo_ooov.cptr(), W55caaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   28] -- 
  // |-- [    0] --| W56(a4,k,a2,v0) += (    1.00000000) T2(a0,a1,a2,v0) D2(a0,a4,a1,k) 
  // |-- [    1] --| W57(w,a3,k,a2) += (    1.00000000) V2(v0,a3,w,a4) W56(a4,k,a2,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W56aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no28_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO28_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W56aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no28_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO28_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W56aaa_sigma_cooo_ooov.cptr(), W57caaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   29] -- 
  // |-- [    0] --| W58(k,a4,a2,v0) += (    1.00000000) T2(a0,a1,a2,v0) D2(a0,k,a1,a4) 
  // |-- [    1] --| W59(w,a3,k,a2) += (    1.00000000) V2(v0,a3,w,a4) W58(k,a4,a2,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W58aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no29_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO29_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W58aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no29_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO29_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W58aaa_sigma_cooo_ooov.cptr(), W59caaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   30] -- 
  // |-- [    0] --| W60(k,a4,a2,v0) += (    1.00000000) T2(a0,a1,a2,v0) D2(a0,k,a1,a4) 
  // |-- [    1] --| W61(w,a3,k,a2) += (    1.00000000) V2(v0,a3,w,a4) W60(k,a4,a2,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W60aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no30_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO30_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W60aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no30_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO30_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W60aaa_sigma_cooo_ooov.cptr(), W61caaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   31] -- 
  // |-- [    0] --| W62(i,a4,a1,v0) += (    1.00000000) T2(a0,a1,a2,v0) D2(a0,a2,i,a4) 
  // |-- [    1] --| W63(w,a3,i,a1) += (    1.00000000) V2(v0,a3,w,a4) W62(i,a4,a1,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W62aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no31_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO31_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W62aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no31_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO31_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W62aaa_sigma_cooo_ooov.cptr(), W63caaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   32] -- 
  // |-- [    0] --| W64(i,a3,a1,v0) += (    1.00000000) T2(a0,a1,a2,v0) D2(a0,a2,i,a3) 
  // |-- [    1] --| W65(w,a4,i,a1) += (    1.00000000) V2(v0,a3,w,a4) W64(i,a3,a1,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W64aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no32_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO32_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W64aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no32_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO32_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W64aaa_sigma_cooo_ooov.cptr(), W65caaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   33] -- 
  // |-- [    0] --| W67(i,k,a1,v0) += (    1.00000000) T2(a0,a1,a2,v0) D2(a0,a2,i,k) 
  // |-- [    1] --| W66(w,a1,j,v0) += (    1.00000000) V2(v0,a3,w,a4) D2(a1,a3,j,a4) 
  // |-- [    2] --| S2(w,k,i,j) += (   -0.50000000) W66(w,a1,j,v0) W67(i,k,a1,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W67aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no33_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO33_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W67aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W66ca_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_cooo_ooov_no33_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO33_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W66ca_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no33_x2_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO33_X2_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, W66ca_sigma_cooo_ooov.cptr(), W67aaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   34] -- 
  // |-- [    0] --| W68(w,a0,i,v0) += (    1.00000000) V2(v0,a3,w,a4) D2(a0,a3,i,a4) 
  // |-- [    1] --| W69(j,k,a0,v0) += (    1.00000000) T2(a0,a1,a2,v0) D2(a1,a2,j,k) 
  // |-- [    2] --| S2(w,k,i,j) += (   -0.50000000) W68(w,a0,i,v0) W69(j,k,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W68caa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no34_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO34_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W68caa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W69aa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_cooo_ooov_no34_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO34_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W69aa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_ooov_no34_x2_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO34_X2_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, W68caa_sigma_cooo_ooov.cptr(), W69aa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   35] -- 
  // |-- [    0] --| W70(j,a3,a0,v0) += (    1.00000000) T2(a0,a1,a2,v0) D2(a1,a3,j,a2) 
  // |-- [    1] --| W71(w,a4,j,a0) += (    1.00000000) V2(v0,a3,w,a4) W70(j,a3,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W70aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no35_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO35_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W70aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no35_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO35_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W70aaa_sigma_cooo_ooov.cptr(), W71caaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   36] -- 
  // |-- [    0] --| W72(j,a4,a0,v0) += (    1.00000000) T2(a0,a1,a2,v0) D2(a1,a2,j,a4) 
  // |-- [    1] --| W73(w,a3,j,a0) += (    1.00000000) V2(v0,a3,w,a4) W72(j,a4,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W72aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no36_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO36_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W72aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no36_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO36_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W72aaa_sigma_cooo_ooov.cptr(), W73caaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   37] -- 
  // |-- [    0] --| W74(j,a4,a0,v0) += (    1.00000000) T2(a0,a1,a2,v0) D2(a1,a4,j,a2) 
  // |-- [    1] --| W75(w,a3,j,a0) += (    1.00000000) V2(v0,a3,w,a4) W74(j,a4,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W74aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no37_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO37_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W74aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no37_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO37_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W74aaa_sigma_cooo_ooov.cptr(), W75caaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   38] -- 
  // |-- [    0] --| W76(j,a4,a0,v0) += (    1.00000000) T2(a0,a1,a2,v0) D2(a1,a2,j,a4) 
  // |-- [    1] --| W77(w,a3,j,a0) += (    1.00000000) V2(v0,a3,w,a4) W76(j,a4,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W76aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no38_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO38_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W76aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no38_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO38_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W76aaa_sigma_cooo_ooov.cptr(), W77caaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   39] -- 
  // |-- [    0] --| W78(j,a4,a0,v0) += (    1.00000000) T2(a0,a1,a2,v0) D2(a1,a4,j,a2) 
  // |-- [    1] --| W79(w,a3,j,a0) += (    1.00000000) V2(v0,a3,w,a4) W78(j,a4,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W78aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no39_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO39_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W78aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no39_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO39_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W78aaa_sigma_cooo_ooov.cptr(), W79caaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   40] -- 
  // |-- [    0] --| W80(w,i,a1,v0) += (    1.00000000) V2(v0,a3,w,a4) D2(i,a4,a1,a3) 
  // |-- [    1] --| W81(j,k,a1,v0) += (    1.00000000) T2(a0,a1,a2,v0) D2(a0,a2,j,k) 
  // |-- [    2] --| S2(w,k,i,j) += (    1.00000000) W80(w,i,a1,v0) W81(j,k,a1,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W80caa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no40_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO40_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W80caa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W81aa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_cooo_ooov_no40_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO40_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W81aa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_ooov_no40_x2_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO40_X2_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, W80caa_sigma_cooo_ooov.cptr(), W81aa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   41] -- 
  // |-- [    0] --| W82(j,a4,a1,v0) += (    1.00000000) T2(a0,a1,a2,v0) D2(a0,a2,j,a4) 
  // |-- [    1] --| W83(w,a3,j,a1) += (    1.00000000) V2(v0,a3,w,a4) W82(j,a4,a1,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W82aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no41_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO41_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W82aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no41_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO41_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W82aaa_sigma_cooo_ooov.cptr(), W83caaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   42] -- 
  // |-- [    0] --| W84(j,a3,a1,v0) += (    1.00000000) T2(a0,a1,a2,v0) D2(a0,a2,j,a3) 
  // |-- [    1] --| W85(w,a4,j,a1) += (    1.00000000) V2(v0,a3,w,a4) W84(j,a3,a1,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W84aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no42_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO42_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W84aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no42_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO42_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W84aaa_sigma_cooo_ooov.cptr(), W85caaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   43] -- 
  // |-- [    0] --| W86(i,a3,a0,v0) += (    1.00000000) T2(a0,a1,a2,v0) D2(i,a2,a1,a3) 
  // |-- [    1] --| W87(w,a4,i,a0) += (    1.00000000) V2(v0,a3,w,a4) W86(i,a3,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W86aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no43_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO43_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W86aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no43_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO43_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W86aaa_sigma_cooo_ooov.cptr(), W87caaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   44] -- 
  // |-- [    0] --| W88(i,a4,a0,v0) += (    1.00000000) T2(a0,a1,a2,v0) D2(i,a4,a1,a2) 
  // |-- [    1] --| W89(w,a3,i,a0) += (    1.00000000) V2(v0,a3,w,a4) W88(i,a4,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W88aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no44_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO44_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W88aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no44_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO44_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W88aaa_sigma_cooo_ooov.cptr(), W89caaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   45] -- 
  // |-- [    0] --| W91(i,k,a0,v0) += (    1.00000000) T2(a0,a1,a2,v0) D2(i,a2,a1,k) 
  // |-- [    1] --| W90(w,a0,j,v0) += (    1.00000000) V2(v0,a3,w,a4) D2(a0,a4,j,a3) 
  // |-- [    2] --| S2(w,k,i,j) += (    0.33333333) W90(w,a0,j,v0) W91(i,k,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W91aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no45_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO45_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W91aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W90ca_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_cooo_ooov_no45_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO45_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W90ca_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no45_x2_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO45_X2_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, W90ca_sigma_cooo_ooov.cptr(), W91aaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   46] -- 
  // |-- [    0] --| W93(i,k,a0,v0) += (    1.00000000) T2(a0,a1,a2,v0) D2(i,k,a1,a2) 
  // |-- [    1] --| W92(w,a0,j,v0) += (    1.00000000) V2(v0,a3,w,a4) D2(a0,a4,j,a3) 
  // |-- [    2] --| S2(w,k,i,j) += (    0.16666667) W92(w,a0,j,v0) W93(i,k,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W93aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no46_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO46_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W93aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W92ca_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_cooo_ooov_no46_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO46_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W92ca_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no46_x2_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO46_X2_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, W92ca_sigma_cooo_ooov.cptr(), W93aaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   47] -- 
  // |-- [    0] --| W95(i,k,a0,v0) += (    1.00000000) T2(a0,a1,a2,v0) D2(i,a2,a1,k) 
  // |-- [    1] --| W94(w,a0,j,v0) += (    1.00000000) V2(v0,a3,w,a4) D2(a0,a3,j,a4) 
  // |-- [    2] --| S2(w,k,i,j) += (    0.16666667) W94(w,a0,j,v0) W95(i,k,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W95aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no47_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO47_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W95aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W94ca_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_cooo_ooov_no47_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO47_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W94ca_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no47_x2_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO47_X2_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, W94ca_sigma_cooo_ooov.cptr(), W95aaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   48] -- 
  // |-- [    0] --| W97(i,k,a0,v0) += (    1.00000000) T2(a0,a1,a2,v0) D2(i,k,a1,a2) 
  // |-- [    1] --| W96(w,a0,j,v0) += (    1.00000000) V2(v0,a3,w,a4) D2(a0,a3,j,a4) 
  // |-- [    2] --| S2(w,k,i,j) += (    0.33333333) W96(w,a0,j,v0) W97(i,k,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W97aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no48_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO48_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W97aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W96ca_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_cooo_ooov_no48_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO48_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W96ca_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no48_x2_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO48_X2_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, W96ca_sigma_cooo_ooov.cptr(), W97aaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   49] -- 
  // |-- [    0] --| W98(a3,v0) += (    1.00000000) T2(a0,a1,a2,v0) C2(a0,a2,a1,a3) 
  // |-- [    1] --| W99(w,a4) += (    1.00000000) V2(v0,a3,w,a4) W98(a3,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W98a_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no49_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO49_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W98a_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no49_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO49_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W98a_sigma_cooo_ooov.cptr(), W99ca_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   50] -- 
  // |-- [    0] --| W100(a4,v0) += (    1.00000000) T2(a0,a1,a2,v0) C2(a0,a2,a1,a4) 
  // |-- [    1] --| W101(w,a3) += (    1.00000000) V2(v0,a3,w,a4) W100(a4,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W100a_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no50_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO50_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W100a_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no50_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO50_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W100a_sigma_cooo_ooov.cptr(), W101ca_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   51] -- 
  // |-- [    0] --| W102(k,v0) += (    1.00000000) T2(a0,a1,a2,v0) C2(a0,a2,a1,k) 
  // |-- [    1] --| W103(w,a4,a3,k) += (    1.00000000) V2(v0,a3,w,a4) W102(k,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W102a_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no51_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO51_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W102a_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no51_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO51_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W102a_sigma_cooo_ooov.cptr(), W103caaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   52] -- 
  // |-- [    0] --| W104(a4,a3,a2,v0) += (    1.00000000) T2(a0,a1,a2,v0) C2(a0,a4,a1,a3) 
  // |-- [    1] --| W105(w,a2) += (    1.00000000) V2(v0,a3,w,a4) W104(a4,a3,a2,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W104aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no52_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO52_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W104aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no52_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO52_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W104aaa_sigma_cooo_ooov.cptr(), W105ca_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   53] -- 
  // |-- [    0] --| W106(k,a3,a2,v0) += (    1.00000000) T2(a0,a1,a2,v0) C2(a0,k,a1,a3) 
  // |-- [    1] --| W107(w,a4,k,a2) += (    1.00000000) V2(v0,a3,w,a4) W106(k,a3,a2,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W106aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no53_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO53_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W106aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no53_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO53_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W106aaa_sigma_cooo_ooov.cptr(), W107caaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   54] -- 
  // |-- [    0] --| W108(a4,k,a2,v0) += (    1.00000000) T2(a0,a1,a2,v0) C2(a0,a4,a1,k) 
  // |-- [    1] --| W109(w,a3,k,a2) += (    1.00000000) V2(v0,a3,w,a4) W108(a4,k,a2,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W108aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no54_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO54_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W108aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no54_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO54_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W108aaa_sigma_cooo_ooov.cptr(), W109caaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   55] -- 
  // |-- [    0] --| W110(a4,k,a2,v0) += (    1.00000000) T2(a0,a1,a2,v0) C2(a0,a4,a1,k) 
  // |-- [    1] --| W111(w,a3,k,a2) += (    1.00000000) V2(v0,a3,w,a4) W110(a4,k,a2,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W110aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no55_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO55_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W110aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no55_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO55_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W110aaa_sigma_cooo_ooov.cptr(), W111caaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   56] -- 
  // |-- [    0] --| W112(k,a4,a2,v0) += (    1.00000000) T2(a0,a1,a2,v0) C2(a0,k,a1,a4) 
  // |-- [    1] --| W113(w,a3,k,a2) += (    1.00000000) V2(v0,a3,w,a4) W112(k,a4,a2,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W112aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no56_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO56_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W112aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no56_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO56_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W112aaa_sigma_cooo_ooov.cptr(), W113caaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   57] -- 
  // |-- [    0] --| W114(k,a4,a2,v0) += (    1.00000000) T2(a0,a1,a2,v0) C2(a0,k,a1,a4) 
  // |-- [    1] --| W115(w,a3,k,a2) += (    1.00000000) V2(v0,a3,w,a4) W114(k,a4,a2,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W114aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no57_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO57_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W114aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no57_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO57_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W114aaa_sigma_cooo_ooov.cptr(), W115caaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   58] -- 
  // |-- [    0] --| W116(i,a4,a1,v0) += (    1.00000000) T2(a0,a1,a2,v0) C2(a0,a2,i,a4) 
  // |-- [    1] --| W117(w,a3,i,a1) += (    1.00000000) V2(v0,a3,w,a4) W116(i,a4,a1,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W116aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no58_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO58_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W116aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no58_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO58_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W116aaa_sigma_cooo_ooov.cptr(), W117caaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   59] -- 
  // |-- [    0] --| W118(i,a3,a1,v0) += (    1.00000000) T2(a0,a1,a2,v0) C2(a0,a2,i,a3) 
  // |-- [    1] --| W119(w,a4,i,a1) += (    1.00000000) V2(v0,a3,w,a4) W118(i,a3,a1,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W118aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no59_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO59_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W118aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no59_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO59_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W118aaa_sigma_cooo_ooov.cptr(), W119caaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   60] -- 
  // |-- [    0] --| W121(i,k,a1,v0) += (    1.00000000) T2(a0,a1,a2,v0) C2(a0,a2,i,k) 
  // |-- [    1] --| W120(w,a1,j,v0) += (    1.00000000) V2(v0,a3,w,a4) C2(a1,a3,j,a4) 
  // |-- [    2] --| S2(w,k,i,j) += (    1.00000000) W120(w,a1,j,v0) W121(i,k,a1,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W121aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no60_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO60_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W121aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W120ca_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_cooo_ooov_no60_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO60_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W120ca_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no60_x2_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO60_X2_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, W120ca_sigma_cooo_ooov.cptr(), W121aaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   61] -- 
  // |-- [    0] --| W122(w,a0,i,v0) += (    1.00000000) V2(v0,a3,w,a4) C2(a0,a3,i,a4) 
  // |-- [    1] --| W123(j,k,a0,v0) += (    1.00000000) T2(a0,a1,a2,v0) C2(a1,a2,j,k) 
  // |-- [    2] --| S2(w,k,i,j) += (    1.00000000) W122(w,a0,i,v0) W123(j,k,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W122caa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no61_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO61_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W122caa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W123aa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_cooo_ooov_no61_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO61_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W123aa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_ooov_no61_x2_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO61_X2_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, W122caa_sigma_cooo_ooov.cptr(), W123aa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   62] -- 
  // |-- [    0] --| W124(j,a3,a0,v0) += (    1.00000000) T2(a0,a1,a2,v0) C2(a1,a3,j,a2) 
  // |-- [    1] --| W125(w,a4,j,a0) += (    1.00000000) V2(v0,a3,w,a4) W124(j,a3,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W124aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no62_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO62_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W124aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no62_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO62_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W124aaa_sigma_cooo_ooov.cptr(), W125caaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   63] -- 
  // |-- [    0] --| W126(j,a4,a0,v0) += (    1.00000000) T2(a0,a1,a2,v0) C2(a1,a2,j,a4) 
  // |-- [    1] --| W127(w,a3,j,a0) += (    1.00000000) V2(v0,a3,w,a4) W126(j,a4,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W126aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no63_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO63_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W126aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no63_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO63_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W126aaa_sigma_cooo_ooov.cptr(), W127caaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   64] -- 
  // |-- [    0] --| W128(j,a4,a0,v0) += (    1.00000000) T2(a0,a1,a2,v0) C2(a1,a4,j,a2) 
  // |-- [    1] --| W129(w,a3,j,a0) += (    1.00000000) V2(v0,a3,w,a4) W128(j,a4,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W128aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no64_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO64_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W128aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no64_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO64_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W128aaa_sigma_cooo_ooov.cptr(), W129caaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   65] -- 
  // |-- [    0] --| W130(j,a4,a0,v0) += (    1.00000000) T2(a0,a1,a2,v0) C2(a1,a2,j,a4) 
  // |-- [    1] --| W131(w,a3,j,a0) += (    1.00000000) V2(v0,a3,w,a4) W130(j,a4,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W130aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no65_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO65_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W130aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no65_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO65_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W130aaa_sigma_cooo_ooov.cptr(), W131caaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   66] -- 
  // |-- [    0] --| W132(j,a4,a0,v0) += (    1.00000000) T2(a0,a1,a2,v0) C2(a1,a4,j,a2) 
  // |-- [    1] --| W133(w,a3,j,a0) += (    1.00000000) V2(v0,a3,w,a4) W132(j,a4,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W132aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no66_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO66_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W132aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no66_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO66_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W132aaa_sigma_cooo_ooov.cptr(), W133caaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   67] -- 
  // |-- [    0] --| W134(w,i,a1,v0) += (    1.00000000) V2(v0,a3,w,a4) C2(a1,a3,i,a4) 
  // |-- [    1] --| W135(j,k,a1,v0) += (    1.00000000) T2(a0,a1,a2,v0) C2(a0,a2,j,k) 
  // |-- [    2] --| S2(w,k,i,j) += (   -2.00000000) W134(w,i,a1,v0) W135(j,k,a1,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W134caa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no67_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO67_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W134caa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W135aa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_cooo_ooov_no67_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO67_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W135aa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_ooov_no67_x2_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO67_X2_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, W134caa_sigma_cooo_ooov.cptr(), W135aa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   68] -- 
  // |-- [    0] --| W136(j,a4,a1,v0) += (    1.00000000) T2(a0,a1,a2,v0) C2(a0,a2,j,a4) 
  // |-- [    1] --| W137(w,a3,j,a1) += (    1.00000000) V2(v0,a3,w,a4) W136(j,a4,a1,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W136aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no68_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO68_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W136aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no68_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO68_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W136aaa_sigma_cooo_ooov.cptr(), W137caaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   69] -- 
  // |-- [    0] --| W138(j,a3,a1,v0) += (    1.00000000) T2(a0,a1,a2,v0) C2(a0,a2,j,a3) 
  // |-- [    1] --| W139(w,a4,j,a1) += (    1.00000000) V2(v0,a3,w,a4) W138(j,a3,a1,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W138aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no69_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO69_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W138aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no69_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO69_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W138aaa_sigma_cooo_ooov.cptr(), W139caaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   70] -- 
  // |-- [    0] --| W140(i,a3,a0,v0) += (    1.00000000) T2(a0,a1,a2,v0) C2(a1,a3,i,a2) 
  // |-- [    1] --| W141(w,a4,i,a0) += (    1.00000000) V2(v0,a3,w,a4) W140(i,a3,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W140aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no70_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO70_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W140aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no70_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO70_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W140aaa_sigma_cooo_ooov.cptr(), W141caaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   71] -- 
  // |-- [    0] --| W142(i,a4,a0,v0) += (    1.00000000) T2(a0,a1,a2,v0) C2(a1,a2,i,a4) 
  // |-- [    1] --| W143(w,a3,i,a0) += (    1.00000000) V2(v0,a3,w,a4) W142(i,a4,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W142aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no71_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO71_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W142aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_cooo_ooov_no71_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO71_X1_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W142aaa_sigma_cooo_ooov.cptr(), W143caaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   72] -- 
  // |-- [    0] --| W145(i,k,a0,v0) += (    1.00000000) T2(a0,a1,a2,v0) C2(a1,k,i,a2) 
  // |-- [    1] --| W144(w,a0,j,v0) += (    1.00000000) V2(v0,a3,w,a4) C2(a0,a4,j,a3) 
  // |-- [    2] --| S2(w,k,i,j) += (   -0.66666667) W144(w,a0,j,v0) W145(i,k,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W145aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no72_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO72_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W145aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W144ca_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_cooo_ooov_no72_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO72_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W144ca_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no72_x2_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO72_X2_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, W144ca_sigma_cooo_ooov.cptr(), W145aaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   73] -- 
  // |-- [    0] --| W147(i,k,a0,v0) += (    1.00000000) T2(a0,a1,a2,v0) C2(a1,a2,i,k) 
  // |-- [    1] --| W146(w,a0,j,v0) += (    1.00000000) V2(v0,a3,w,a4) C2(a0,a4,j,a3) 
  // |-- [    2] --| S2(w,k,i,j) += (   -0.33333333) W146(w,a0,j,v0) W147(i,k,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W147aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no73_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO73_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W147aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W146ca_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_cooo_ooov_no73_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO73_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W146ca_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no73_x2_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO73_X2_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, W146ca_sigma_cooo_ooov.cptr(), W147aaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   74] -- 
  // |-- [    0] --| W149(i,k,a0,v0) += (    1.00000000) T2(a0,a1,a2,v0) C2(a1,k,i,a2) 
  // |-- [    1] --| W148(w,a0,j,v0) += (    1.00000000) V2(v0,a3,w,a4) C2(a0,a3,j,a4) 
  // |-- [    2] --| S2(w,k,i,j) += (   -0.33333333) W148(w,a0,j,v0) W149(i,k,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W149aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no74_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO74_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W149aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W148ca_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_cooo_ooov_no74_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO74_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W148ca_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no74_x2_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO74_X2_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, W148ca_sigma_cooo_ooov.cptr(), W149aaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   75] -- 
  // |-- [    0] --| W151(i,k,a0,v0) += (    1.00000000) T2(a0,a1,a2,v0) C2(a1,a2,i,k) 
  // |-- [    1] --| W150(w,a0,j,v0) += (    1.00000000) V2(v0,a3,w,a4) C2(a0,a3,j,a4) 
  // |-- [    2] --| S2(w,k,i,j) += (   -0.66666667) W150(w,a0,j,v0) W151(i,k,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [2] L- and R-contractions are not contiguous >> 0
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W151aaa_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_ooov_no75_x0_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO75_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W151aaa_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W150ca_sigma_cooo_ooov(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_cooo_ooov_no75_x1_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO75_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W150ca_sigma_cooo_ooov.cptr(), nir, nsym, psym, &flops);
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no75_x2_type1_eri_v,G_IF_SIGMA_COOO_OOOV_NO75_X2_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, W150ca_sigma_cooo_ooov.cptr(), W151aaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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

//-@type(2).contraction(begin)
  // -- Title : sigma_cooo_ooov
  //*-- Entering to take the type 2 contractions --*//

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -1.00000000) D1(j,k) W13(w,i) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no0_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO0_X0_TYPE2_ERI_V)
      (sj, ij, W13ca_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    0.50000000) D1(i,k) W27(w,j) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no1_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO1_X0_TYPE2_ERI_V)
      (sj, ij, W27ca_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -1.00000000) D3(a0,a2,i,a4,j,k) W29(w,a0,a2,a4) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no2_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO2_X0_TYPE2_ERI_V)
      (sj, ij, W29caaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    0.50000000) D3(a0,a3,i,a4,j,k) W31(w,a4,a3,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no3_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO3_X0_TYPE2_ERI_V)
      (sj, ij, W31caaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    0.50000000) D3(a0,a2,i,a3,j,k) W33(w,a0,a2,a3) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no4_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO4_X0_TYPE2_ERI_V)
      (sj, ij, W33caaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -1.00000000) D3(a1,a3,i,a4,j,k) W37(w,a4,a3,a1) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no5_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO5_X0_TYPE2_ERI_V)
      (sj, ij, W37caaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    0.50000000) D3(a1,a2,i,a4,j,k) W39(w,a1,a2,a4) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no6_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO6_X0_TYPE2_ERI_V)
      (sj, ij, W39caaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    0.50000000) D3(a1,a3,i,a2,j,k) W41(w,a1,a2,a3) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no7_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO7_X0_TYPE2_ERI_V)
      (sj, ij, W41caaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    1.00000000) D2(i,a4,j,k) W45(w,a4) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no8_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO8_X0_TYPE2_ERI_V)
      (sj, ij, W45ca_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -0.50000000) D2(i,a3,j,k) W47(w,a3) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no9_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO9_X0_TYPE2_ERI_V)
      (sj, ij, W47ca_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -0.50000000) D2(i,a4,j,a3) W49(w,a4,a3,k) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no10_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO10_X0_TYPE2_ERI_V)
      (sj, ij, W49caaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -0.50000000) D2(i,a2,j,k) W51(w,a2) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no11_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO11_X0_TYPE2_ERI_V)
      (sj, ij, W51ca_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -0.50000000) D2(i,a4,j,a2) W53(w,a4,k,a2) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no12_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO12_X0_TYPE2_ERI_V)
      (sj, ij, W53caaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    0.33333333) D2(i,a2,j,a3) W55(w,a3,k,a2) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no13_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO13_X0_TYPE2_ERI_V)
      (sj, ij, W55caaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   14] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    0.16666667) D2(i,a3,j,a2) W57(w,a3,k,a2) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no14_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO14_X0_TYPE2_ERI_V)
      (sj, ij, W57caaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   15] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    0.16666667) D2(i,a2,j,a3) W59(w,a3,k,a2) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no15_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO15_X0_TYPE2_ERI_V)
      (sj, ij, W59caaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   16] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    0.33333333) D2(i,a3,j,a2) W61(w,a3,k,a2) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no16_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO16_X0_TYPE2_ERI_V)
      (sj, ij, W61caaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   17] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    1.00000000) D2(a1,a3,j,k) W63(w,a3,i,a1) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no17_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO17_X0_TYPE2_ERI_V)
      (sj, ij, W63caaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   18] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -0.50000000) D2(a1,a4,j,k) W65(w,a4,i,a1) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no18_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO18_X0_TYPE2_ERI_V)
      (sj, ij, W65caaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   19] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -0.50000000) D2(a0,k,i,a4) W71(w,a4,j,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no19_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO19_X0_TYPE2_ERI_V)
      (sj, ij, W71caaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   20] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    0.33333333) D2(a0,a3,i,k) W73(w,a3,j,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no20_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO20_X0_TYPE2_ERI_V)
      (sj, ij, W73caaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   21] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    0.16666667) D2(a0,a3,i,k) W75(w,a3,j,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no21_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO21_X0_TYPE2_ERI_V)
      (sj, ij, W75caaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   22] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    0.16666667) D2(a0,k,i,a3) W77(w,a3,j,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no22_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO22_X0_TYPE2_ERI_V)
      (sj, ij, W77caaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   23] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    0.33333333) D2(a0,k,i,a3) W79(w,a3,j,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no23_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO23_X0_TYPE2_ERI_V)
      (sj, ij, W79caaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   24] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -0.50000000) D2(i,k,a1,a3) W83(w,a3,j,a1) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no24_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO24_X0_TYPE2_ERI_V)
      (sj, ij, W83caaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   25] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -0.50000000) D2(i,a4,a1,k) W85(w,a4,j,a1) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no25_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO25_X0_TYPE2_ERI_V)
      (sj, ij, W85caaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   26] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -0.50000000) D2(a0,a4,j,k) W87(w,a4,i,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no26_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO26_X0_TYPE2_ERI_V)
      (sj, ij, W87caaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   27] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -0.50000000) D2(a0,a3,j,k) W89(w,a3,i,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no27_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO27_X0_TYPE2_ERI_V)
      (sj, ij, W89caaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   28] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -2.00000000) C2(a4,i,k,j) W99(w,a4) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no28_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO28_X0_TYPE2_ERI_V)
      (sj, ij, W99ca_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   29] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    1.00000000) C2(a3,i,k,j) W101(w,a3) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no29_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO29_X0_TYPE2_ERI_V)
      (sj, ij, W101ca_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   30] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    1.00000000) C2(a3,j,a4,i) W103(w,a4,a3,k) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no30_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO30_X0_TYPE2_ERI_V)
      (sj, ij, W103caaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   31] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    1.00000000) C2(a2,i,k,j) W105(w,a2) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no31_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO31_X0_TYPE2_ERI_V)
      (sj, ij, W105ca_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   32] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    1.00000000) C2(a2,j,a4,i) W107(w,a4,k,a2) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no32_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO32_X0_TYPE2_ERI_V)
      (sj, ij, W107caaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   33] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -0.66666667) C2(a2,i,a3,j) W109(w,a3,k,a2) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no33_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO33_X0_TYPE2_ERI_V)
      (sj, ij, W109caaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   34] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -0.33333333) C2(a2,j,a3,i) W111(w,a3,k,a2) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no34_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO34_X0_TYPE2_ERI_V)
      (sj, ij, W111caaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   35] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -0.33333333) C2(a2,i,a3,j) W113(w,a3,k,a2) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no35_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO35_X0_TYPE2_ERI_V)
      (sj, ij, W113caaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   36] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -0.66666667) C2(a2,j,a3,i) W115(w,a3,k,a2) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no36_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO36_X0_TYPE2_ERI_V)
      (sj, ij, W115caaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   37] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -2.00000000) C2(a1,a3,j,k) W117(w,a3,i,a1) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no37_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO37_X0_TYPE2_ERI_V)
      (sj, ij, W117caaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   38] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    1.00000000) C2(a1,a4,j,k) W119(w,a4,i,a1) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no38_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO38_X0_TYPE2_ERI_V)
      (sj, ij, W119caaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   39] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    1.00000000) C2(a0,k,i,a4) W125(w,a4,j,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no39_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO39_X0_TYPE2_ERI_V)
      (sj, ij, W125caaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   40] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -0.66666667) C2(a0,a3,i,k) W127(w,a3,j,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no40_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO40_X0_TYPE2_ERI_V)
      (sj, ij, W127caaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   41] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -0.33333333) C2(a0,a3,i,k) W129(w,a3,j,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no41_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO41_X0_TYPE2_ERI_V)
      (sj, ij, W129caaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   42] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -0.33333333) C2(a0,k,i,a3) W131(w,a3,j,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no42_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO42_X0_TYPE2_ERI_V)
      (sj, ij, W131caaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   43] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -0.66666667) C2(a0,k,i,a3) W133(w,a3,j,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no43_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO43_X0_TYPE2_ERI_V)
      (sj, ij, W133caaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   44] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    1.00000000) C2(a1,a3,i,k) W137(w,a3,j,a1) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no44_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO44_X0_TYPE2_ERI_V)
      (sj, ij, W137caaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   45] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    1.00000000) C2(a1,k,i,a4) W139(w,a4,j,a1) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no45_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO45_X0_TYPE2_ERI_V)
      (sj, ij, W139caaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   46] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    1.00000000) C2(a0,a4,j,k) W141(w,a4,i,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no46_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO46_X0_TYPE2_ERI_V)
      (sj, ij, W141caaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   47] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    1.00000000) C2(a0,a3,j,k) W143(w,a3,i,a0) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_ooov_no47_x0_type2_eri_v,G_IF_SIGMA_COOO_OOOV_NO47_X0_TYPE2_ERI_V)
      (sj, ij, W143caaa_sigma_cooo_ooov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  orz::mr::Femto::my_timer.push_back(boost::make_tuple("sigma_cooo_ooov", time_sigma.elapsed_cputime(), time_sigma.elapsed_wallclocktime()));
  orz::mr::Femto::file_timing << "* " << boost::format("%20s : %10.7f %10.7f ") % orz::mr::Femto::my_timer.back().get<0>() % orz::mr::Femto::my_timer.back().get<1>() % orz::mr::Femto::my_timer.back().get<2>() << endl;
  flush(orz::mr::Femto::file_timing);

  return retval; 
} 
