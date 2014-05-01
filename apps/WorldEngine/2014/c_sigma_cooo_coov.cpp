                                                                                
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
#include <sci/icmr/Femto/elems/c_sigma_cooo_coov.h>                                  
                                                                                
using std::cout;                                                                
using std::endl;                                                                
                                                                                
#define FLOPCOUNT                                                               
                                                                                
using orz::mr::Femto::Fc0;                                                      
using orz::mr::Femto::h1_int;                                                   
using orz::mr::Femto::h6_int;                                                   
                                                                                
//  __/\\\\\\\\\\\\\\\____________________________________________________________________                                   
//   _\/\\\///////////_____________________________________________________________________                                             
//    _\/\\\_______________________________________________________/\\\_____________________                                         
//     _\/\\\\\\\\\\\__________/\\\\\\\\______/\\\\\__/\\\\\_____/\\\\\\\\\\\______/\\\\\____ 
//      _\/\\\///////_________/\\\/////\\\___/\\\///\\\\\///\\\__\////\\\////_____/\\\///\\\__               
//       _\/\\\_______________/\\\\\\\\\\\___\/\\\_\//\\\__\/\\\_____\/\\\________/\\\__\//\\\_       
//        _\/\\\______________\//\\///////____\/\\\__\/\\\__\/\\\_____\/\\\_/\\___\//\\\__/\\\__            
//         _\/\\\_______________\//\\\\\\\\\\__\/\\\__\/\\\__\/\\\_____\//\\\\\_____\///\\\\\/___    
//          _\///_________________\//////////___\///___\///___\///_______\/////________\/////_____                                   

//                                   Generated date : Sun Apr 20 10:26:16 2014

                                                                                
// ***************************************************************************  
// orz::mr::mrci                                                                
// ***************************************************************************  
									    /*!        
   @brief CT input                                                              
									     */        
                                                                                
orz::mr::BareAmpPack orz::mr::Femto::sigma_cooo_coov(const orz::mr::Input &ctinp,                                    
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
  // |-- [    0] --| W0(w,a2,a0,a1) += (    1.00000000) T2(w,a0,a1,v0) Fc1(v0,a2) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) D3(k,j,a2,a0,a1,i) W0(w,a2,a0,a1) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W0caaa_sigma_cooo_coov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_cooo_coov_no0_x0_type0_noeri,G_IF_SIGMA_COOO_COOV_NO0_X0_TYPE0_NOERI)
      (sv0, iv0, T2b.cptr(), W0caaa_sigma_cooo_coov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_coov_no0_x1_type0_noeri,G_IF_SIGMA_COOO_COOV_NO0_X1_TYPE0_NOERI)
      (sj, ij, W0caaa_sigma_cooo_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W1(w,a1,a0,i) += (    1.00000000) T2(w,a0,i,v0) Fc1(v0,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) D2(k,j,a1,a0) W1(w,a1,a0,i) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W1caaa_sigma_cooo_coov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_cooo_coov_no1_x0_type0_noeri,G_IF_SIGMA_COOO_COOV_NO1_X0_TYPE0_NOERI)
      (sv0, iv0, T2b.cptr(), W1caaa_sigma_cooo_coov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_coov_no1_x1_type0_noeri,G_IF_SIGMA_COOO_COOV_NO1_X1_TYPE0_NOERI)
      (sj, ij, W1caaa_sigma_cooo_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W2(w,i,a0,a1) += (    1.00000000) T2(w,a0,a1,v0) Fc1(v0,i) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) D2(k,j,a1,a0) W2(w,i,a0,a1) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W2caaa_sigma_cooo_coov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_cooo_coov_no2_x0_type0_noeri,G_IF_SIGMA_COOO_COOV_NO2_X0_TYPE0_NOERI)
      (sv0, iv0, T2b.cptr(), W2caaa_sigma_cooo_coov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_coov_no2_x1_type0_noeri,G_IF_SIGMA_COOO_COOV_NO2_X1_TYPE0_NOERI)
      (sj, ij, W2caaa_sigma_cooo_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W3(w,a1,a0,j) += (    1.00000000) T2(a0,w,v0,j) Fc1(v0,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) D2(k,i,a1,a0) W3(w,a1,a0,j) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W3caa_sigma_cooo_coov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_coov_no3_x0_type0_noeri,G_IF_SIGMA_COOO_COOV_NO3_X0_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W3caa_sigma_cooo_coov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_coov_no3_x1_type0_noeri,G_IF_SIGMA_COOO_COOV_NO3_X1_TYPE0_NOERI)
      (sj, ij, W3caa_sigma_cooo_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W4(w,j,a0,a1) += (    1.00000000) T2(w,a0,a1,v0) Fc1(v0,j) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) D2(k,a0,a1,i) W4(w,j,a0,a1) 
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    orz::DTensor W4caa_sigma_cooo_coov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
    for(int sv0 = 0;sv0 < nir;++sv0){ 
    for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
      T2b = T2.get_amp2(iv0);
      FC_FUNC(g_if_sigma_cooo_coov_no4_x0_type0_noeri,G_IF_SIGMA_COOO_COOV_NO4_X0_TYPE0_NOERI)
        (sj, ij, sv0, iv0, T2b.cptr(), W4caa_sigma_cooo_coov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
    } // End iv0
    } // End sv0
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_coov_no4_x1_type0_noeri,G_IF_SIGMA_COOO_COOV_NO4_X1_TYPE0_NOERI)
      (sj, ij, W4caa_sigma_cooo_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W5(w,k,i,v0) += (    1.00000000) T2(w,a0,i,v0) D1(k,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) Fc1(v0,j) W5(w,k,i,v0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    orz::DTensor W5caa_sigma_cooo_coov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sv0));
    FC_FUNC(g_if_sigma_cooo_coov_no5_x0_type0_noeri,G_IF_SIGMA_COOO_COOV_NO5_X0_TYPE0_NOERI)
      (sv0, iv0, T2b.cptr(), W5caa_sigma_cooo_coov.cptr(), nir, nsym, psym, &flops);
    for(int sj = 0;sj < nir;++sj){ 
    for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
      S2b = orz::DTensor(retval.namps_iamp()[ij]);
      FC_FUNC(g_if_sigma_cooo_coov_no5_x1_type0_noeri,G_IF_SIGMA_COOO_COOV_NO5_X1_TYPE0_NOERI)
        (sj, ij, sv0, iv0, W5caa_sigma_cooo_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
      retval.acc_amp2(ij, S2b);
    } // End ij
    } // End sj
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W6(w,k,j,v0) += (    1.00000000) T2(a0,w,v0,j) D1(k,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) Fc1(v0,i) W6(w,k,j,v0) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W6cav_sigma_cooo_coov(orz::mr::sizeof_sympack_Xcav(symblockinfo, sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_coov_no6_x0_type0_noeri,G_IF_SIGMA_COOO_COOV_NO6_X0_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W6cav_sigma_cooo_coov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_coov_no6_x1_type0_noeri,G_IF_SIGMA_COOO_COOV_NO6_X1_TYPE0_NOERI)
      (sj, ij, W6cav_sigma_cooo_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- Title : sigma_cooo_coov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W18(w,a0,j,a1) += (    1.00000000) V2(w,a1,v0,c0) T2(a0,c0,v0,j) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) D2(k,a0,a1,i) W18(w,a0,j,a1) 
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
    orz::DTensor W18aa_sigma_cooo_coov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_coov_no0_x0_type1_eri_c,G_IF_SIGMA_COOO_COOV_NO0_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), V2_sym.cptr(), W18aa_sigma_cooo_coov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_coov_no0_x1_type1_eri_c,G_IF_SIGMA_COOO_COOV_NO0_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W18aa_sigma_cooo_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W19(w,a0,j,a1) += (    1.00000000) V2(w,c0,v0,a1) T2(a0,c0,v0,j) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) D2(k,i,a1,a0) W19(w,a0,j,a1) 
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
    orz::DTensor W19aa_sigma_cooo_coov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_coov_no1_x0_type1_eri_c,G_IF_SIGMA_COOO_COOV_NO1_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), V2_sym.cptr(), W19aa_sigma_cooo_coov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_coov_no1_x1_type1_eri_c,G_IF_SIGMA_COOO_COOV_NO1_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W19aa_sigma_cooo_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W23(w,a0,j,i) += (    1.00000000) V2(w,i,v0,c0) T2(a0,c0,v0,j) 
  // |-- [    1] --| S2(w,k,i,j) += (   -2.00000000) D1(k,a0) W23(w,a0,j,i) 
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
    orz::DTensor W23aa_sigma_cooo_coov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_coov_no2_x0_type1_eri_c,G_IF_SIGMA_COOO_COOV_NO2_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), V2_sym.cptr(), W23aa_sigma_cooo_coov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_coov_no2_x1_type1_eri_c,G_IF_SIGMA_COOO_COOV_NO2_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W23aa_sigma_cooo_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W31(w,a0,j,i) += (    1.00000000) V2(w,c0,v0,i) T2(a0,c0,v0,j) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) D1(k,a0) W31(w,a0,j,i) 
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
    orz::DTensor W31aa_sigma_cooo_coov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sw^sj));
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_coov_no3_x0_type1_eri_c,G_IF_SIGMA_COOO_COOV_NO3_X0_TYPE1_ERI_C)
      (sj, ij, sw, iw, T2b.cptr(), V2_sym.cptr(), W31aa_sigma_cooo_coov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_coov_no3_x1_type1_eri_c,G_IF_SIGMA_COOO_COOV_NO3_X1_TYPE1_ERI_C)
      (sj, ij, sw, iw, W31aa_sigma_cooo_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- Title : sigma_cooo_coov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W10(w,a0,a1,i,a3,a2) += (    1.00000000) V2(a3,v0,i,a2) T2(w,a0,a1,v0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) D3(k,j,a3,a0,a1,a2) W10(w,a0,a1,i,a3,a2) 
  int sa3(s_eri);
  int ia3(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W10caaaa_sigma_cooo_coov(orz::mr::sizeof_sympack_Xcaaaa(symblockinfo, sa3));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_cooo_coov_no0_x0_type1_eri_o,G_IF_SIGMA_COOO_COOV_NO0_X0_TYPE1_ERI_O)
      (sa3, ia3, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W10caaaa_sigma_cooo_coov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_coov_no0_x1_type1_eri_o,G_IF_SIGMA_COOO_COOV_NO0_X1_TYPE1_ERI_O)
      (sa3, ia3, sj, ij, W10caaaa_sigma_cooo_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W15(w,a0,a1,a3,i,a2) += (    1.00000000) V2(i,v0,a3,a2) T2(w,a0,a1,v0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) D3(k,j,a3,a2,a1,a0) W15(w,a0,a1,a3,i,a2) 
  int si(s_eri);
  int ii(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W15caaaa_sigma_cooo_coov(orz::mr::sizeof_sympack_Xcaaaa(symblockinfo, si));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_cooo_coov_no1_x0_type1_eri_o,G_IF_SIGMA_COOO_COOV_NO1_X0_TYPE1_ERI_O)
      (si, ii, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W15caaaa_sigma_cooo_coov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_coov_no1_x1_type1_eri_o,G_IF_SIGMA_COOO_COOV_NO1_X1_TYPE1_ERI_O)
      (si, ii, sj, ij, W15caaaa_sigma_cooo_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W16(w,a0,a1,j) += (    1.00000000) V2(j,w,v0,c0) T2(c0,a0,a1,v0) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) D2(k,i,a1,a0) W16(w,a0,a1,j) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W16caa_sigma_cooo_coov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_cooo_coov_no2_x0_type1_eri_o,G_IF_SIGMA_COOO_COOV_NO2_X0_TYPE1_ERI_O)
      (sj, ij, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W16caa_sigma_cooo_coov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_cooo_coov_no2_x1_type1_eri_o,G_IF_SIGMA_COOO_COOV_NO2_X1_TYPE1_ERI_O)
    (sj, ij, W16caa_sigma_cooo_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W17(w,a0,a1,j,a3,a2) += (    1.00000000) V2(j,a2,v0,a3) T2(w,a0,a1,v0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) D3(k,a2,a3,a0,a1,i) W17(w,a0,a1,j,a3,a2) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W17caaaa_sigma_cooo_coov(orz::mr::sizeof_sympack_Xcaaaa(symblockinfo, sj));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_cooo_coov_no3_x0_type1_eri_o,G_IF_SIGMA_COOO_COOV_NO3_X0_TYPE1_ERI_O)
      (sj, ij, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W17caaaa_sigma_cooo_coov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_cooo_coov_no3_x1_type1_eri_o,G_IF_SIGMA_COOO_COOV_NO3_X1_TYPE1_ERI_O)
    (sj, ij, W17caaaa_sigma_cooo_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W21(w,a0,a1,j) += (    1.00000000) V2(j,v0,w,c0) T2(c0,a0,a1,v0) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) D2(k,a0,a1,i) W21(w,a0,a1,j) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W21caa_sigma_cooo_coov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_cooo_coov_no4_x0_type1_eri_o,G_IF_SIGMA_COOO_COOV_NO4_X0_TYPE1_ERI_O)
      (sj, ij, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W21caa_sigma_cooo_coov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_cooo_coov_no4_x1_type1_eri_o,G_IF_SIGMA_COOO_COOV_NO4_X1_TYPE1_ERI_O)
    (sj, ij, W21caa_sigma_cooo_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W22(w,a0,a1,a3,j,a2) += (    1.00000000) V2(j,v0,a3,a2) T2(w,a0,a1,v0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) D3(k,a0,a3,a2,a1,i) W22(w,a0,a1,a3,j,a2) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W22caaaa_sigma_cooo_coov(orz::mr::sizeof_sympack_Xcaaaa(symblockinfo, sj));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_cooo_coov_no5_x0_type1_eri_o,G_IF_SIGMA_COOO_COOV_NO5_X0_TYPE1_ERI_O)
      (sj, ij, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W22caaaa_sigma_cooo_coov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_cooo_coov_no5_x1_type1_eri_o,G_IF_SIGMA_COOO_COOV_NO5_X1_TYPE1_ERI_O)
    (sj, ij, W22caaaa_sigma_cooo_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W26(w,a0,i,j) += (    1.00000000) V2(j,w,v0,c0) T2(c0,a0,i,v0) 
  // |-- [    1] --| S2(w,k,i,j) += (    1.00000000) D1(k,a0) W26(w,a0,i,j) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W26caa_sigma_cooo_coov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_cooo_coov_no6_x0_type1_eri_o,G_IF_SIGMA_COOO_COOV_NO6_X0_TYPE1_ERI_O)
      (sj, ij, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W26caa_sigma_cooo_coov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_cooo_coov_no6_x1_type1_eri_o,G_IF_SIGMA_COOO_COOV_NO6_X1_TYPE1_ERI_O)
    (sj, ij, W26caa_sigma_cooo_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W28(w,a0,i,j) += (    1.00000000) V2(j,v0,w,c0) T2(c0,a0,i,v0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -2.00000000) D1(k,a0) W28(w,a0,i,j) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W28caa_sigma_cooo_coov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_cooo_coov_no7_x0_type1_eri_o,G_IF_SIGMA_COOO_COOV_NO7_X0_TYPE1_ERI_O)
      (sj, ij, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W28caa_sigma_cooo_coov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_cooo_coov_no7_x1_type1_eri_o,G_IF_SIGMA_COOO_COOV_NO7_X1_TYPE1_ERI_O)
    (sj, ij, W28caa_sigma_cooo_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| W36(w,a0,j,a2) += (    1.00000000) V2(j,a1,v0,a2) T2(w,a0,a1,v0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) D2(k,i,a2,a0) W36(w,a0,j,a2) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W36caa_sigma_cooo_coov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_cooo_coov_no8_x0_type1_eri_o,G_IF_SIGMA_COOO_COOV_NO8_X0_TYPE1_ERI_O)
      (sj, ij, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W36caa_sigma_cooo_coov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_cooo_coov_no8_x1_type1_eri_o,G_IF_SIGMA_COOO_COOV_NO8_X1_TYPE1_ERI_O)
    (sj, ij, W36caa_sigma_cooo_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W37(w,a0,a2,j) += (    1.00000000) V2(j,v0,a2,a1) T2(w,a0,a1,v0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) D2(k,a0,a2,i) W37(w,a0,a2,j) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W37caa_sigma_cooo_coov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_cooo_coov_no9_x0_type1_eri_o,G_IF_SIGMA_COOO_COOV_NO9_X0_TYPE1_ERI_O)
      (sj, ij, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W37caa_sigma_cooo_coov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_cooo_coov_no9_x1_type1_eri_o,G_IF_SIGMA_COOO_COOV_NO9_X1_TYPE1_ERI_O)
    (sj, ij, W37caa_sigma_cooo_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| W38(w,a0,i,j) += (    1.00000000) V2(j,v0,i,a1) T2(w,a0,a1,v0) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) D1(k,a0) W38(w,a0,i,j) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W38caa_sigma_cooo_coov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_cooo_coov_no10_x0_type1_eri_o,G_IF_SIGMA_COOO_COOV_NO10_X0_TYPE1_ERI_O)
      (sj, ij, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W38caa_sigma_cooo_coov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_cooo_coov_no10_x1_type1_eri_o,G_IF_SIGMA_COOO_COOV_NO10_X1_TYPE1_ERI_O)
    (sj, ij, W38caa_sigma_cooo_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ij, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| W39(w,a0,j,i) += (    1.00000000) V2(j,a1,v0,i) T2(w,a0,a1,v0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) D1(k,a0) W39(w,a0,j,i) 
  int sj(s_eri);
  int ij(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ij]);
  orz::DTensor W39caa_sigma_cooo_coov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    FC_FUNC(g_if_sigma_cooo_coov_no11_x0_type1_eri_o,G_IF_SIGMA_COOO_COOV_NO11_X0_TYPE1_ERI_O)
      (sj, ij, sv0, iv0, T2b.cptr(), V2_sym.cptr(), W39caa_sigma_cooo_coov.cptr(), nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  FC_FUNC(g_if_sigma_cooo_coov_no11_x1_type1_eri_o,G_IF_SIGMA_COOO_COOV_NO11_X1_TYPE1_ERI_O)
    (sj, ij, W39caa_sigma_cooo_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // --  Title : sigma_cooo_coov
  //  >> Intermediates for the external contractions are defined here << 
  orz::DTensor W7caaa_sigma_cooo_coov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W8caaa_sigma_cooo_coov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W9caaa_sigma_cooo_coov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W11caaa_sigma_cooo_coov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W12caaa_sigma_cooo_coov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W14caaa_sigma_cooo_coov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W33caaa_sigma_cooo_coov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W34caaa_sigma_cooo_coov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W35caaa_sigma_cooo_coov(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
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
  // -- Title : sigma_cooo_coov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W7(w,a0,a1,a2) += (    1.00000000) V2(v0,c0,w,a2) T2(c0,a0,a1,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  FC_FUNC(g_if_sigma_cooo_coov_no0_x0_type1_eri_v,G_IF_SIGMA_COOO_COOV_NO0_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), V2_sym.cptr(), W7caaa_sigma_cooo_coov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W8(w,a0,a1,a2) += (    1.00000000) V2(v0,a2,w,c0) T2(c0,a0,a1,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  FC_FUNC(g_if_sigma_cooo_coov_no1_x0_type1_eri_v,G_IF_SIGMA_COOO_COOV_NO1_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), V2_sym.cptr(), W8caaa_sigma_cooo_coov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W9(w,a0,a1,i) += (    1.00000000) V2(v0,c0,w,i) T2(c0,a0,a1,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  FC_FUNC(g_if_sigma_cooo_coov_no2_x0_type1_eri_v,G_IF_SIGMA_COOO_COOV_NO2_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), V2_sym.cptr(), W9caaa_sigma_cooo_coov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W11(w,a0,i,a1) += (    1.00000000) V2(v0,c0,w,a1) T2(c0,a0,i,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  FC_FUNC(g_if_sigma_cooo_coov_no3_x0_type1_eri_v,G_IF_SIGMA_COOO_COOV_NO3_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), V2_sym.cptr(), W11caaa_sigma_cooo_coov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W12(w,a0,i,a1) += (    1.00000000) V2(v0,a1,w,c0) T2(c0,a0,i,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  FC_FUNC(g_if_sigma_cooo_coov_no4_x0_type1_eri_v,G_IF_SIGMA_COOO_COOV_NO4_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), V2_sym.cptr(), W12caaa_sigma_cooo_coov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W13(k,j,a0,v0) += (    1.00000000) V2(v0,a2,a3,a1) D3(k,j,a3,a1,a2,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) T2(w,a0,i,v0) W13(k,j,a0,v0) 
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
    orz::DTensor W13aa_sigma_cooo_coov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_cooo_coov_no5_x0_type1_eri_v,G_IF_SIGMA_COOO_COOV_NO5_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W13aa_sigma_cooo_coov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_coov_no5_x1_type1_eri_v,G_IF_SIGMA_COOO_COOV_NO5_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W13aa_sigma_cooo_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W14(w,a0,a1,i) += (    1.00000000) V2(v0,i,w,c0) T2(c0,a0,a1,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  FC_FUNC(g_if_sigma_cooo_coov_no6_x0_type1_eri_v,G_IF_SIGMA_COOO_COOV_NO6_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), V2_sym.cptr(), W14caaa_sigma_cooo_coov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W20(k,i,a0,v0) += (    1.00000000) V2(v0,a2,a3,a1) D3(k,i,a3,a1,a2,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) T2(w,a0,j,v0) W20(k,i,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W20aaa_sigma_cooo_coov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_coov_no7_x0_type1_eri_v,G_IF_SIGMA_COOO_COOV_NO7_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W20aaa_sigma_cooo_coov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_coov_no7_x1_type1_eri_v,G_IF_SIGMA_COOO_COOV_NO7_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W20aaa_sigma_cooo_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| W24(k,a0,i,v0) += (    1.00000000) V2(v0,a2,i,a1) D2(k,a1,a2,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) T2(w,a0,j,v0) W24(k,a0,i,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W24aaa_sigma_cooo_coov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_coov_no8_x0_type1_eri_v,G_IF_SIGMA_COOO_COOV_NO8_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W24aaa_sigma_cooo_coov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_coov_no8_x1_type1_eri_v,G_IF_SIGMA_COOO_COOV_NO8_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W24aaa_sigma_cooo_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W25(w,k,a2,v0) += (    1.00000000) T2(w,a0,a1,v0) D2(k,a0,a1,a2) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) V2(v0,j,i,a2) W25(w,k,a2,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W25caa_sigma_cooo_coov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_coov_no9_x0_type1_eri_v,G_IF_SIGMA_COOO_COOV_NO9_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W25caa_sigma_cooo_coov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_coov_no9_x1_type1_eri_v,G_IF_SIGMA_COOO_COOV_NO9_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W25caa_sigma_cooo_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| W27(k,a0,j,v0) += (    1.00000000) V2(v0,a2,j,a1) D2(k,a1,a2,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) T2(w,a0,i,v0) W27(k,a0,j,v0) 
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
    orz::DTensor W27aa_sigma_cooo_coov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_cooo_coov_no10_x0_type1_eri_v,G_IF_SIGMA_COOO_COOV_NO10_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W27aa_sigma_cooo_coov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_coov_no10_x1_type1_eri_v,G_IF_SIGMA_COOO_COOV_NO10_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W27aa_sigma_cooo_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| W29(k,a0,j,v0) += (    1.00000000) V2(v0,j,a2,a1) D2(k,a0,a2,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (    2.00000000) T2(w,a0,i,v0) W29(k,a0,j,v0) 
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
    orz::DTensor W29aa_sigma_cooo_coov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sj^sv0));
    FC_FUNC(g_if_sigma_cooo_coov_no11_x0_type1_eri_v,G_IF_SIGMA_COOO_COOV_NO11_X0_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W29aa_sigma_cooo_coov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_cooo_coov_no11_x1_type1_eri_v,G_IF_SIGMA_COOO_COOV_NO11_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W29aa_sigma_cooo_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| W30(w,k,a2,v0) += (    1.00000000) T2(w,a0,a1,v0) D2(k,a2,a1,a0) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) V2(v0,i,j,a2) W30(w,k,a2,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W30caa_sigma_cooo_coov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_coov_no12_x0_type1_eri_v,G_IF_SIGMA_COOO_COOV_NO12_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W30caa_sigma_cooo_coov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_coov_no12_x1_type1_eri_v,G_IF_SIGMA_COOO_COOV_NO12_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, V2_sym.cptr(), W30caa_sigma_cooo_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| W32(k,a0,i,v0) += (    1.00000000) V2(v0,i,a2,a1) D2(k,a0,a2,a1) 
  // |-- [    1] --| S2(w,k,i,j) += (   -1.00000000) T2(w,a0,j,v0) W32(k,a0,i,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(iv0);
  orz::DTensor W32aaa_sigma_cooo_coov(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_cooo_coov_no13_x0_type1_eri_v,G_IF_SIGMA_COOO_COOV_NO13_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W32aaa_sigma_cooo_coov.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_coov_no13_x1_type1_eri_v,G_IF_SIGMA_COOO_COOV_NO13_X1_TYPE1_ERI_V)
      (sj, ij, sv0, iv0, T2b.cptr(), W32aaa_sigma_cooo_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   14] -- 
  // |-- [    0] --| W33(w,a0,a3,a2) += (    1.00000000) V2(v0,a2,a3,a1) T2(w,a0,a1,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  FC_FUNC(g_if_sigma_cooo_coov_no14_x0_type1_eri_v,G_IF_SIGMA_COOO_COOV_NO14_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), V2_sym.cptr(), W33caaa_sigma_cooo_coov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   15] -- 
  // |-- [    0] --| W34(w,a0,i,a2) += (    1.00000000) V2(v0,a2,i,a1) T2(w,a0,a1,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  FC_FUNC(g_if_sigma_cooo_coov_no15_x0_type1_eri_v,G_IF_SIGMA_COOO_COOV_NO15_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), V2_sym.cptr(), W34caaa_sigma_cooo_coov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   16] -- 
  // |-- [    0] --| W35(w,a0,a2,i) += (    1.00000000) V2(v0,i,a2,a1) T2(w,a0,a1,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  FC_FUNC(g_if_sigma_cooo_coov_no16_x0_type1_eri_v,G_IF_SIGMA_COOO_COOV_NO16_X0_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), V2_sym.cptr(), W35caaa_sigma_cooo_coov.cptr(), nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

//-@type(1).contraction(end)

  } // End myrank
  }
  }
  orz::world().barrier();

//-@type(2).contraction(begin)
  // -- Title : sigma_cooo_coov
  //*-- Entering to take the type 2 contractions --*//

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    1.00000000) D3(k,j,a2,i,a1,a0) W7(w,a0,a1,a2) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_coov_no0_x0_type2_eri_v,G_IF_SIGMA_COOO_COOV_NO0_X0_TYPE2_ERI_V)
      (sj, ij, W7caaa_sigma_cooo_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    1.00000000) D3(k,j,a2,a0,a1,i) W8(w,a0,a1,a2) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_coov_no1_x0_type2_eri_v,G_IF_SIGMA_COOO_COOV_NO1_X0_TYPE2_ERI_V)
      (sj, ij, W8caaa_sigma_cooo_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -2.00000000) D2(k,j,a1,a0) W9(w,a0,a1,i) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_coov_no2_x0_type2_eri_v,G_IF_SIGMA_COOO_COOV_NO2_X0_TYPE2_ERI_V)
      (sj, ij, W9caaa_sigma_cooo_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    1.00000000) D2(k,j,a1,a0) W11(w,a0,i,a1) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_coov_no3_x0_type2_eri_v,G_IF_SIGMA_COOO_COOV_NO3_X0_TYPE2_ERI_V)
      (sj, ij, W11caaa_sigma_cooo_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -2.00000000) D2(k,j,a1,a0) W12(w,a0,i,a1) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_coov_no4_x0_type2_eri_v,G_IF_SIGMA_COOO_COOV_NO4_X0_TYPE2_ERI_V)
      (sj, ij, W12caaa_sigma_cooo_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    1.00000000) D2(k,j,a1,a0) W14(w,a0,a1,i) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_coov_no5_x0_type2_eri_v,G_IF_SIGMA_COOO_COOV_NO5_X0_TYPE2_ERI_V)
      (sj, ij, W14caaa_sigma_cooo_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -1.00000000) D3(k,j,a3,i,a2,a0) W33(w,a0,a3,a2) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_coov_no6_x0_type2_eri_v,G_IF_SIGMA_COOO_COOV_NO6_X0_TYPE2_ERI_V)
      (sj, ij, W33caaa_sigma_cooo_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| S2(w,k,i,j) += (    2.00000000) D2(k,j,a2,a0) W34(w,a0,i,a2) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_coov_no7_x0_type2_eri_v,G_IF_SIGMA_COOO_COOV_NO7_X0_TYPE2_ERI_V)
      (sj, ij, W34caaa_sigma_cooo_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -1.00000000) D2(k,j,a2,a0) W35(w,a0,a2,i) 
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_coov_no8_x0_type2_eri_v,G_IF_SIGMA_COOO_COOV_NO8_X0_TYPE2_ERI_V)
      (sj, ij, W35caaa_sigma_cooo_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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

//-@loadD4C(v,begin)
  //*-- FEMTO begins --//*
  // Label : d4c_v
  {


  for(int s_eri = 0;s_eri < nir;++s_eri){ 
  for(int i_eri = symblockinfo.psym()(s_eri,I_V,I_BEGIN);i_eri <= symblockinfo.psym()(s_eri,I_V,I_END);++i_eri){ 
  if(hintmo.iproc_havingimo()[i_eri] == myrank) {           
  orz::DTensor C5;
  orz::LoadBin(ctinp.dir()/(format("D4C_g[%d]")%i_eri).str()) >> C5;

  //*-- Entering to take the type 1 contractions --*//
//-@type(1).contraction(begin)
  // -- Title : sigma_cooo_coov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| S2(w,k,i,j) += (   -1.00000000) T2(w,a0,a1,v0) C5(i,a1,j,k,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_cooo_coov_no0_x0_type1_d4c_v,G_IF_SIGMA_COOO_COOV_NO0_X0_TYPE1_D4C_V)
      (sj, ij, sv0, iv0, C5.cptr(), T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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

//-@loadD4C(v,end)

//-@D4C.contractions(end)

  // Do timing!
  orz::mr::Femto::my_timer.push_back(boost::make_tuple("sigma_cooo_coov", time_sigma.elapsed_cputime(), time_sigma.elapsed_wallclocktime()));
  orz::mr::Femto::file_timing << "* " << boost::format("%20s : %10.7f %10.7f ") % orz::mr::Femto::my_timer.back().get<0>() % orz::mr::Femto::my_timer.back().get<1>() % orz::mr::Femto::my_timer.back().get<2>() << endl;
  flush(orz::mr::Femto::file_timing);

  return retval; 
} 
