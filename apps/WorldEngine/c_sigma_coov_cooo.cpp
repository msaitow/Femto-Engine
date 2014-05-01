                                                                                
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
#include <sci/icmr/Femto/elems/c_sigma_coov_cooo.h>                                  
                                                                                
using std::cout;                                                                
using std::endl;                                                                
                                                                                
#define FLOPCOUNT                                                               
                                                                                
using orz::mr::Femto::Fc0;                                                      
using orz::mr::Femto::h1_int;                                                   
using orz::mr::Femto::h6_int;                                                   
                                                                                
//       #                #########      #     #   # 
//  ########## ##########         #   #######  #   # 
//      #    #         #          #    # #     #   # 
//      #    #        #   ########     # #     #   # 
//     #     #     # #           #  ##########    #  
//    #   # #       #            #       #       #   
//   #     #         #    ########       #     ##    

//                                   Generated date : Sun Apr 20 10:26:18 2014

                                                                                
// ***************************************************************************  
// orz::mr::mrci                                                                
// ***************************************************************************  
									    /*!        
   @brief CT input                                                              
									     */        
                                                                                
orz::mr::BareAmpPack orz::mr::Femto::sigma_coov_cooo(const orz::mr::Input &ctinp,                                    
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
  // |-- [    0] --| W0(w,j,a3,i) += (    1.00000000) T2(w,a2,a1,a0) D3(j,a1,a3,i,a2,a0) 
  // |-- [    1] --| S2(w,i,j,a) += (   -1.00000000) Fc1(a3,a) W0(w,j,a3,i) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W0caaa_sigma_coov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_coov_cooo_no0_x0_type0_noeri,G_IF_SIGMA_COOV_COOO_NO0_X0_TYPE0_NOERI)
      (sa0, ia0, T2b.cptr(), W0caaa_sigma_coov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_cooo_no0_x1_type0_noeri,G_IF_SIGMA_COOV_COOO_NO0_X1_TYPE0_NOERI)
      (sa, ia, W0caaa_sigma_coov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W1(w,i,a2,j) += (    1.00000000) T2(w,a0,j,a1) D2(i,a2,a1,a0) 
  // |-- [    1] --| S2(w,i,j,a) += (    2.00000000) Fc1(a2,a) W1(w,i,a2,j) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W1caaa_sigma_coov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_coov_cooo_no1_x0_type0_noeri,G_IF_SIGMA_COOV_COOO_NO1_X0_TYPE0_NOERI)
      (sa1, ia1, T2b.cptr(), W1caaa_sigma_coov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_cooo_no1_x1_type0_noeri,G_IF_SIGMA_COOV_COOO_NO1_X1_TYPE0_NOERI)
      (sa, ia, W1caaa_sigma_coov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W2(w,i,a2,j) += (    1.00000000) T2(w,a0,a1,j) D2(i,a2,a1,a0) 
  // |-- [    1] --| S2(w,i,j,a) += (   -1.00000000) Fc1(a2,a) W2(w,i,a2,j) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W2caa_sigma_coov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
    FC_FUNC(g_if_sigma_coov_cooo_no2_x0_type0_noeri,G_IF_SIGMA_COOV_COOO_NO2_X0_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W2caa_sigma_coov_cooo.cptr(), nir, nsym, psym, &flops);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      S2b = orz::DTensor(retval.namps_iamp()[ia]);
      FC_FUNC(g_if_sigma_coov_cooo_no2_x1_type0_noeri,G_IF_SIGMA_COOV_COOO_NO2_X1_TYPE0_NOERI)
        (sa, ia, sj, ij, W2caa_sigma_coov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
      retval.acc_amp2(ia, S2b);
    } // End ia
    } // End sa
  // --> @[j, "active"] [notNeeded]
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W3(w,j,i,a2) += (    1.00000000) T2(w,a1,a2,a0) D2(j,i,a1,a0) 
  // |-- [    1] --| S2(w,i,j,a) += (   -1.00000000) Fc1(a2,a) W3(w,j,i,a2) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W3caaa_sigma_coov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_coov_cooo_no3_x0_type0_noeri,G_IF_SIGMA_COOV_COOO_NO3_X0_TYPE0_NOERI)
      (sa0, ia0, T2b.cptr(), W3caaa_sigma_coov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_cooo_no3_x1_type0_noeri,G_IF_SIGMA_COOV_COOO_NO3_X1_TYPE0_NOERI)
      (sa, ia, W3caaa_sigma_coov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W4(w,i,a1,j) += (    1.00000000) T2(w,a0,a1,j) D1(i,a0) 
  // |-- [    1] --| S2(w,i,j,a) += (   -1.00000000) Fc1(a1,a) W4(w,i,a1,j) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W4caa_sigma_coov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sj));
    FC_FUNC(g_if_sigma_coov_cooo_no4_x0_type0_noeri,G_IF_SIGMA_COOV_COOO_NO4_X0_TYPE0_NOERI)
      (sj, ij, T2b.cptr(), W4caa_sigma_coov_cooo.cptr(), nir, nsym, psym, &flops);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      S2b = orz::DTensor(retval.namps_iamp()[ia]);
      FC_FUNC(g_if_sigma_coov_cooo_no4_x1_type0_noeri,G_IF_SIGMA_COOV_COOO_NO4_X1_TYPE0_NOERI)
        (sa, ia, sj, ij, W4caa_sigma_coov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
      retval.acc_amp2(ia, S2b);
    } // End ia
    } // End sa
  // --> @[j, "active"] [notNeeded]
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W5(w,j,i,a2) += (    1.00000000) T2(w,a1,a0,a2) D2(j,a0,a1,i) 
  // |-- [    1] --| S2(w,i,j,a) += (   -1.00000000) Fc1(a2,a) W5(w,j,i,a2) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int sa2 = 0;sa2 < nir;++sa2){ 
  for(int ia2 = symblockinfo.psym()(sa2,I_O,I_BEGIN);ia2 <= symblockinfo.psym()(sa2,I_O,I_END);++ia2){ 
    T2b = T2.get_amp2(ia2);
    orz::DTensor W5caa_sigma_coov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa2));
    FC_FUNC(g_if_sigma_coov_cooo_no5_x0_type0_noeri,G_IF_SIGMA_COOV_COOO_NO5_X0_TYPE0_NOERI)
      (sa2, ia2, T2b.cptr(), W5caa_sigma_coov_cooo.cptr(), nir, nsym, psym, &flops);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      S2b = orz::DTensor(retval.namps_iamp()[ia]);
      FC_FUNC(g_if_sigma_coov_cooo_no5_x1_type0_noeri,G_IF_SIGMA_COOV_COOO_NO5_X1_TYPE0_NOERI)
        (sa, ia, sa2, ia2, W5caa_sigma_coov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
      retval.acc_amp2(ia, S2b);
    } // End ia
    } // End sa
  // --> @[a2, "active"] [notNeeded]
  } // End ia2
  } // End sa2
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W6(w,i,j,a1) += (    1.00000000) T2(w,a0,j,a1) D1(i,a0) 
  // |-- [    1] --| S2(w,i,j,a) += (    2.00000000) Fc1(a1,a) W6(w,i,j,a1) 
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W6caa_sigma_coov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa1));
    FC_FUNC(g_if_sigma_coov_cooo_no6_x0_type0_noeri,G_IF_SIGMA_COOV_COOO_NO6_X0_TYPE0_NOERI)
      (sa1, ia1, T2b.cptr(), W6caa_sigma_coov_cooo.cptr(), nir, nsym, psym, &flops);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      S2b = orz::DTensor(retval.namps_iamp()[ia]);
      FC_FUNC(g_if_sigma_coov_cooo_no6_x1_type0_noeri,G_IF_SIGMA_COOV_COOO_NO6_X1_TYPE0_NOERI)
        (sa, ia, sa1, ia1, W6caa_sigma_coov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
      retval.acc_amp2(ia, S2b);
    } // End ia
    } // End sa
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
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
  // -- Title : sigma_coov_cooo

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W17(w,j,a3,a4,i,a2) += (    1.00000000) T2(a1,w,a0,a2) D3(j,a4,a3,i,a1,a0) 
  // |-- [    1] --| S2(w,i,j,a) += (   -1.00000000) V2(a2,a4,a3,a) W17(w,j,a3,a4,i,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia2);
  orz::DTensor W17caaaa_sigma_coov_cooo(orz::mr::sizeof_sympack_Xcaaaa(symblockinfo, sa2));
  FC_FUNC(g_if_sigma_coov_cooo_no0_x0_type1_eri_o,G_IF_SIGMA_COOV_COOO_NO0_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), W17caaaa_sigma_coov_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_cooo_no0_x1_type1_eri_o,G_IF_SIGMA_COOV_COOO_NO0_X1_TYPE1_ERI_O)
      (sa, ia, sa2, ia2, V2_sym.cptr(), W17caaaa_sigma_coov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W19(w,j,a4,i,a2,a3) += (    1.00000000) T2(a1,w,a0,a3) D3(j,i,a4,a2,a1,a0) 
  // |-- [    1] --| S2(w,i,j,a) += (   -1.00000000) V2(a3,a,a4,a2) W19(w,j,a4,i,a2,a3) 
  int sa3(s_eri);
  int ia3(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia3);
  orz::DTensor W19caaaa_sigma_coov_cooo(orz::mr::sizeof_sympack_Xcaaaa(symblockinfo, sa3));
  FC_FUNC(g_if_sigma_coov_cooo_no1_x0_type1_eri_o,G_IF_SIGMA_COOV_COOO_NO1_X0_TYPE1_ERI_O)
    (sa3, ia3, T2b.cptr(), W19caaaa_sigma_coov_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_cooo_no1_x1_type1_eri_o,G_IF_SIGMA_COOV_COOO_NO1_X1_TYPE1_ERI_O)
      (sa, ia, sa3, ia3, V2_sym.cptr(), W19caaaa_sigma_coov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W20(w,i,a2,a3) += (    1.00000000) T2(a0,w,a1,a3) D2(i,a2,a1,a0) 
  // |-- [    1] --| S2(w,i,j,a) += (   -1.00000000) V2(a3,a,j,a2) W20(w,i,a2,a3) 
  int sa3(s_eri);
  int ia3(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia3);
  orz::DTensor W20caa_sigma_coov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa3));
  FC_FUNC(g_if_sigma_coov_cooo_no2_x0_type1_eri_o,G_IF_SIGMA_COOV_COOO_NO2_X0_TYPE1_ERI_O)
    (sa3, ia3, T2b.cptr(), W20caa_sigma_coov_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_cooo_no2_x1_type1_eri_o,G_IF_SIGMA_COOV_COOO_NO2_X1_TYPE1_ERI_O)
      (sa, ia, sa3, ia3, V2_sym.cptr(), W20caa_sigma_coov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W22(i,a0,a1,a) += (    1.00000000) V2(a1,a3,a2,a) D2(i,a2,a3,a0) 
  // |-- [    1] --| S2(w,i,j,a) += (   -1.00000000) T2(a0,w,j,a1) W22(i,a0,a1,a) 
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
    orz::DTensor W22aa_sigma_coov_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa1^sa));
    FC_FUNC(g_if_sigma_coov_cooo_no3_x0_type1_eri_o,G_IF_SIGMA_COOV_COOO_NO3_X0_TYPE1_ERI_O)
      (sa, ia, sa1, ia1, V2_sym.cptr(), W22aa_sigma_coov_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_coov_cooo_no3_x1_type1_eri_o,G_IF_SIGMA_COOV_COOO_NO3_X1_TYPE1_ERI_O)
      (sa, ia, sa1, ia1, T2b.cptr(), W22aa_sigma_coov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W24(i,a0,a2,a) += (    1.00000000) V2(a2,a,a3,a1) D2(i,a0,a3,a1) 
  // |-- [    1] --| S2(w,i,j,a) += (   -1.00000000) T2(a0,w,j,a2) W24(i,a0,a2,a) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(ia2);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W24aa_sigma_coov_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa2^sa));
    FC_FUNC(g_if_sigma_coov_cooo_no4_x0_type1_eri_o,G_IF_SIGMA_COOV_COOO_NO4_X0_TYPE1_ERI_O)
      (sa, ia, sa2, ia2, V2_sym.cptr(), W24aa_sigma_coov_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_coov_cooo_no4_x1_type1_eri_o,G_IF_SIGMA_COOV_COOO_NO4_X1_TYPE1_ERI_O)
      (sa, ia, sa2, ia2, T2b.cptr(), W24aa_sigma_coov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W25(w,i,a3,a2) += (    1.00000000) T2(a0,w,a1,a2) D2(i,a3,a1,a0) 
  // |-- [    1] --| S2(w,i,j,a) += (    2.00000000) V2(a2,j,a3,a) W25(w,i,a3,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia2);
  orz::DTensor W25caa_sigma_coov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa2));
  FC_FUNC(g_if_sigma_coov_cooo_no5_x0_type1_eri_o,G_IF_SIGMA_COOV_COOO_NO5_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), W25caa_sigma_coov_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_cooo_no5_x1_type1_eri_o,G_IF_SIGMA_COOV_COOO_NO5_X1_TYPE1_ERI_O)
      (sa, ia, sa2, ia2, V2_sym.cptr(), W25caa_sigma_coov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W27(w,j,a3,i,a4,a2) += (    1.00000000) T2(w,a1,a0,a2) D3(j,a0,a3,i,a1,a4) 
  // |-- [    1] --| S2(w,i,j,a) += (   -1.00000000) V2(a2,a4,a3,a) W27(w,j,a3,i,a4,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia2);
  orz::DTensor W27caaaa_sigma_coov_cooo(orz::mr::sizeof_sympack_Xcaaaa(symblockinfo, sa2));
  FC_FUNC(g_if_sigma_coov_cooo_no6_x0_type1_eri_o,G_IF_SIGMA_COOV_COOO_NO6_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), W27caaaa_sigma_coov_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_cooo_no6_x1_type1_eri_o,G_IF_SIGMA_COOV_COOO_NO6_X1_TYPE1_ERI_O)
      (sa, ia, sa2, ia2, V2_sym.cptr(), W27caaaa_sigma_coov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W29(w,j,a4,a2,i,a3) += (    1.00000000) T2(w,a1,a0,a3) D3(j,a0,a4,a2,a1,i) 
  // |-- [    1] --| S2(w,i,j,a) += (   -1.00000000) V2(a3,a,a4,a2) W29(w,j,a4,a2,i,a3) 
  int sa3(s_eri);
  int ia3(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia3);
  orz::DTensor W29caaaa_sigma_coov_cooo(orz::mr::sizeof_sympack_Xcaaaa(symblockinfo, sa3));
  FC_FUNC(g_if_sigma_coov_cooo_no7_x0_type1_eri_o,G_IF_SIGMA_COOV_COOO_NO7_X0_TYPE1_ERI_O)
    (sa3, ia3, T2b.cptr(), W29caaaa_sigma_coov_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_cooo_no7_x1_type1_eri_o,G_IF_SIGMA_COOV_COOO_NO7_X1_TYPE1_ERI_O)
      (sa, ia, sa3, ia3, V2_sym.cptr(), W29caaaa_sigma_coov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| W30(w,i,a2,a3) += (    1.00000000) T2(w,a0,a1,a3) D2(i,a0,a1,a2) 
  // |-- [    1] --| S2(w,i,j,a) += (   -1.00000000) V2(a3,a,j,a2) W30(w,i,a2,a3) 
  int sa3(s_eri);
  int ia3(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia3);
  orz::DTensor W30caa_sigma_coov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa3));
  FC_FUNC(g_if_sigma_coov_cooo_no8_x0_type1_eri_o,G_IF_SIGMA_COOV_COOO_NO8_X0_TYPE1_ERI_O)
    (sa3, ia3, T2b.cptr(), W30caa_sigma_coov_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_cooo_no8_x1_type1_eri_o,G_IF_SIGMA_COOV_COOO_NO8_X1_TYPE1_ERI_O)
      (sa, ia, sa3, ia3, V2_sym.cptr(), W30caa_sigma_coov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W32(i,a0,a1,a) += (    1.00000000) V2(a1,a3,a2,a) D2(i,a2,a3,a0) 
  // |-- [    1] --| S2(w,i,j,a) += (    2.00000000) T2(w,a0,j,a1) W32(i,a0,a1,a) 
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
    orz::DTensor W32aa_sigma_coov_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa1^sa));
    FC_FUNC(g_if_sigma_coov_cooo_no9_x0_type1_eri_o,G_IF_SIGMA_COOV_COOO_NO9_X0_TYPE1_ERI_O)
      (sa, ia, sa1, ia1, V2_sym.cptr(), W32aa_sigma_coov_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_coov_cooo_no9_x1_type1_eri_o,G_IF_SIGMA_COOV_COOO_NO9_X1_TYPE1_ERI_O)
      (sa, ia, sa1, ia1, T2b.cptr(), W32aa_sigma_coov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| W34(i,a0,a2,a) += (    1.00000000) V2(a2,a,a3,a1) D2(i,a0,a3,a1) 
  // |-- [    1] --| S2(w,i,j,a) += (    2.00000000) T2(w,a0,j,a2) W34(i,a0,a2,a) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(ia2);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor W34aa_sigma_coov_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa2^sa));
    FC_FUNC(g_if_sigma_coov_cooo_no10_x0_type1_eri_o,G_IF_SIGMA_COOV_COOO_NO10_X0_TYPE1_ERI_O)
      (sa, ia, sa2, ia2, V2_sym.cptr(), W34aa_sigma_coov_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_coov_cooo_no10_x1_type1_eri_o,G_IF_SIGMA_COOV_COOO_NO10_X1_TYPE1_ERI_O)
      (sa, ia, sa2, ia2, T2b.cptr(), W34aa_sigma_coov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| W35(w,i,a3,a2) += (    1.00000000) T2(w,a0,a1,a2) D2(i,a3,a1,a0) 
  // |-- [    1] --| S2(w,i,j,a) += (   -1.00000000) V2(a2,j,a3,a) W35(w,i,a3,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia2);
  orz::DTensor W35caa_sigma_coov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa2));
  FC_FUNC(g_if_sigma_coov_cooo_no11_x0_type1_eri_o,G_IF_SIGMA_COOV_COOO_NO11_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), W35caa_sigma_coov_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_coov_cooo_no11_x1_type1_eri_o,G_IF_SIGMA_COOV_COOO_NO11_X1_TYPE1_ERI_O)
      (sa, ia, sa2, ia2, V2_sym.cptr(), W35caa_sigma_coov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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

//-@type(2).declaration(begin)
  // --  Title : sigma_coov_cooo
  //  >> Intermediates for the external contractions are defined here << 
  orz::DTensor W7caaa_sigma_coov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W8caaa_sigma_coov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W9caaa_sigma_coov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W10caaa_sigma_coov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W11caaa_sigma_coov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W13caaa_sigma_coov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W14caaa_sigma_coov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
//-@type(2).declaration(end)

  //*-- Entering to take the type 0 contractions --*//
//-@type(0).contraction(begin)

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W7(c0,j,a3,i) += (    1.00000000) T2(c0,a2,a1,a0) D3(j,i,a3,a1,a2,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_coov_cooo_no0_x0_type0_eri_v,G_IF_SIGMA_COOV_COOO_NO0_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W7caaa_sigma_coov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W8(c0,j,a3,i) += (    1.00000000) T2(c0,a2,a1,a0) D3(j,a1,a3,i,a2,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_coov_cooo_no1_x0_type0_eri_v,G_IF_SIGMA_COOV_COOO_NO1_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W8caaa_sigma_coov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W9(w,i,a4,a3) += (    1.00000000) T2(w,a0,a2,a1) D3(i,a4,a2,a3,a1,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_coov_cooo_no2_x0_type0_eri_v,G_IF_SIGMA_COOV_COOO_NO2_X0_TYPE0_ERI_V)
      (sa1, ia1, T2b.cptr(), W9caaa_sigma_coov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W10(c0,i,a2,j) += (    1.00000000) T2(c0,a0,j,a1) D2(i,a2,a1,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_coov_cooo_no3_x0_type0_eri_v,G_IF_SIGMA_COOV_COOO_NO3_X0_TYPE0_ERI_V)
      (sa1, ia1, T2b.cptr(), W10caaa_sigma_coov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W11(c0,i,a2,j) += (    1.00000000) T2(c0,a0,j,a1) D2(i,a2,a1,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_coov_cooo_no4_x0_type0_eri_v,G_IF_SIGMA_COOV_COOO_NO4_X0_TYPE0_ERI_V)
      (sa1, ia1, T2b.cptr(), W11caaa_sigma_coov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W13(c0,i,a2,j) += (    1.00000000) T2(c0,a0,a1,j) D2(i,a0,a1,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    FC_FUNC(g_if_sigma_coov_cooo_no5_x0_type0_eri_v,G_IF_SIGMA_COOV_COOO_NO5_X0_TYPE0_ERI_V)
      (sj, ij, T2b.cptr(), W13caaa_sigma_coov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W14(c0,i,a2,j) += (    1.00000000) T2(c0,a0,a1,j) D2(i,a2,a1,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    FC_FUNC(g_if_sigma_coov_cooo_no6_x0_type0_eri_v,G_IF_SIGMA_COOV_COOO_NO6_X0_TYPE0_ERI_V)
      (sj, ij, T2b.cptr(), W14caaa_sigma_coov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
  } // End ij
  } // End sj
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
  // -- Title : sigma_coov_cooo

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| S2(w,i,j,a) += (    1.00000000) V2(a,w,c0,a3) W7(c0,j,a3,i) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  FC_FUNC(g_if_sigma_coov_cooo_no0_x0_type1_eri_v,G_IF_SIGMA_COOV_COOO_NO0_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W7caaa_sigma_coov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| S2(w,i,j,a) += (    1.00000000) V2(a,a3,w,c0) W8(c0,j,a3,i) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  FC_FUNC(g_if_sigma_coov_cooo_no1_x0_type1_eri_v,G_IF_SIGMA_COOV_COOO_NO1_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W8caaa_sigma_coov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| S2(w,i,j,a) += (   -1.00000000) V2(a,a4,j,a3) W9(w,i,a4,a3) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  FC_FUNC(g_if_sigma_coov_cooo_no2_x0_type1_eri_v,G_IF_SIGMA_COOV_COOO_NO2_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W9caaa_sigma_coov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| S2(w,i,j,a) += (    1.00000000) V2(a,w,c0,a2) W10(c0,i,a2,j) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  FC_FUNC(g_if_sigma_coov_cooo_no3_x0_type1_eri_v,G_IF_SIGMA_COOV_COOO_NO3_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W10caaa_sigma_coov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| S2(w,i,j,a) += (   -2.00000000) V2(a,a2,w,c0) W11(c0,i,a2,j) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  FC_FUNC(g_if_sigma_coov_cooo_no4_x0_type1_eri_v,G_IF_SIGMA_COOV_COOO_NO4_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W11caaa_sigma_coov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W12(i,a1,a0,a) += (    1.00000000) V2(a,a3,a4,a2) D3(i,a3,a4,a2,a1,a0) 
  // |-- [    1] --| S2(w,i,j,a) += (    2.00000000) T2(w,a0,j,a1) W12(i,a1,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W12aa_sigma_coov_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa1^sa));
    FC_FUNC(g_if_sigma_coov_cooo_no5_x0_type1_eri_v,G_IF_SIGMA_COOV_COOO_NO5_X0_TYPE1_ERI_V)
      (sa, ia, sa1, ia1, V2_sym.cptr(), W12aa_sigma_coov_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_coov_cooo_no5_x1_type1_eri_v,G_IF_SIGMA_COOV_COOO_NO5_X1_TYPE1_ERI_V)
      (sa, ia, sa1, ia1, T2b.cptr(), W12aa_sigma_coov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| S2(w,i,j,a) += (    1.00000000) V2(a,w,c0,a2) W13(c0,i,a2,j) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  FC_FUNC(g_if_sigma_coov_cooo_no6_x0_type1_eri_v,G_IF_SIGMA_COOV_COOO_NO6_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W13caaa_sigma_coov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| S2(w,i,j,a) += (    1.00000000) V2(a,a2,w,c0) W14(c0,i,a2,j) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  FC_FUNC(g_if_sigma_coov_cooo_no7_x0_type1_eri_v,G_IF_SIGMA_COOV_COOO_NO7_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W14caaa_sigma_coov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| W15(i,a1,a0,a) += (    1.00000000) V2(a,a3,a4,a2) D3(i,a3,a4,a2,a1,a0) 
  // |-- [    1] --| S2(w,i,j,a) += (   -1.00000000) T2(w,a0,a1,j) W15(i,a1,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W15aaa_sigma_coov_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_coov_cooo_no8_x0_type1_eri_v,G_IF_SIGMA_COOV_COOO_NO8_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W15aaa_sigma_coov_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    FC_FUNC(g_if_sigma_coov_cooo_no8_x1_type1_eri_v,G_IF_SIGMA_COOV_COOO_NO8_X1_TYPE1_ERI_V)
      (sa, ia, sj, ij, T2b.cptr(), W15aaa_sigma_coov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
  } // End ij
  } // End sj
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W16(w,a1,a0,a) += (    1.00000000) V2(a,w,c0,a2) T2(c0,a1,a2,a0) 
  // |-- [    1] --| S2(w,i,j,a) += (   -2.00000000) D2(j,i,a1,a0) W16(w,a1,a0,a) 
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
    orz::DTensor W16ca_sigma_coov_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_coov_cooo_no9_x0_type1_eri_v,G_IF_SIGMA_COOV_COOO_NO9_X0_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), V2_sym.cptr(), W16ca_sigma_coov_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_coov_cooo_no9_x1_type1_eri_v,G_IF_SIGMA_COOV_COOO_NO9_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, W16ca_sigma_coov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| W18(w,a1,a0,a) += (    1.00000000) V2(a,a2,w,c0) T2(c0,a1,a2,a0) 
  // |-- [    1] --| S2(w,i,j,a) += (    1.00000000) D2(j,i,a1,a0) W18(w,a1,a0,a) 
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
    orz::DTensor W18ca_sigma_coov_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_coov_cooo_no10_x0_type1_eri_v,G_IF_SIGMA_COOV_COOO_NO10_X0_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), V2_sym.cptr(), W18ca_sigma_coov_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_coov_cooo_no10_x1_type1_eri_v,G_IF_SIGMA_COOV_COOO_NO10_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, W18ca_sigma_coov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| W21(w,a0,j,a) += (    1.00000000) V2(a,w,c0,a1) T2(c0,a0,a1,j) 
  // |-- [    1] --| S2(w,i,j,a) += (   -2.00000000) D1(i,a0) W21(w,a0,j,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W21ca_sigma_coov_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sa));
    FC_FUNC(g_if_sigma_coov_cooo_no11_x0_type1_eri_v,G_IF_SIGMA_COOV_COOO_NO11_X0_TYPE1_ERI_V)
      (sa, ia, sj, ij, T2b.cptr(), V2_sym.cptr(), W21ca_sigma_coov_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_coov_cooo_no11_x1_type1_eri_v,G_IF_SIGMA_COOV_COOO_NO11_X1_TYPE1_ERI_V)
      (sa, ia, sj, ij, W21ca_sigma_coov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
  } // End ij
  } // End sj
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| W23(w,a0,j,a) += (    1.00000000) V2(a,a1,w,c0) T2(c0,a0,a1,j) 
  // |-- [    1] --| S2(w,i,j,a) += (    1.00000000) D1(i,a0) W23(w,a0,j,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    T2b = T2.get_amp2(ij);
    orz::DTensor W23ca_sigma_coov_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sj^sa));
    FC_FUNC(g_if_sigma_coov_cooo_no12_x0_type1_eri_v,G_IF_SIGMA_COOV_COOO_NO12_X0_TYPE1_ERI_V)
      (sa, ia, sj, ij, T2b.cptr(), V2_sym.cptr(), W23ca_sigma_coov_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_coov_cooo_no12_x1_type1_eri_v,G_IF_SIGMA_COOV_COOO_NO12_X1_TYPE1_ERI_V)
      (sa, ia, sj, ij, W23ca_sigma_coov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
  } // End ij
  } // End sj
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| W26(w,a1,a0,a) += (    1.00000000) V2(a,w,c0,a2) T2(c0,a1,a0,a2) 
  // |-- [    1] --| S2(w,i,j,a) += (    1.00000000) D2(j,i,a1,a0) W26(w,a1,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W26caa_sigma_coov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  for(int sa2 = 0;sa2 < nir;++sa2){ 
  for(int ia2 = symblockinfo.psym()(sa2,I_O,I_BEGIN);ia2 <= symblockinfo.psym()(sa2,I_O,I_END);++ia2){ 
    T2b = T2.get_amp2(ia2);
    FC_FUNC(g_if_sigma_coov_cooo_no13_x0_type1_eri_v,G_IF_SIGMA_COOV_COOO_NO13_X0_TYPE1_ERI_V)
      (sa, ia, sa2, ia2, T2b.cptr(), V2_sym.cptr(), W26caa_sigma_coov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a2, "active"] [notNeeded]
  } // End ia2
  } // End sa2
  FC_FUNC(g_if_sigma_coov_cooo_no13_x1_type1_eri_v,G_IF_SIGMA_COOV_COOO_NO13_X1_TYPE1_ERI_V)
    (sa, ia, W26caa_sigma_coov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   14] -- 
  // |-- [    0] --| W28(w,a1,a0,a) += (    1.00000000) V2(a,a2,w,c0) T2(c0,a1,a0,a2) 
  // |-- [    1] --| S2(w,i,j,a) += (    1.00000000) D2(j,a0,a1,i) W28(w,a1,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W28caa_sigma_coov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  for(int sa2 = 0;sa2 < nir;++sa2){ 
  for(int ia2 = symblockinfo.psym()(sa2,I_O,I_BEGIN);ia2 <= symblockinfo.psym()(sa2,I_O,I_END);++ia2){ 
    T2b = T2.get_amp2(ia2);
    FC_FUNC(g_if_sigma_coov_cooo_no14_x0_type1_eri_v,G_IF_SIGMA_COOV_COOO_NO14_X0_TYPE1_ERI_V)
      (sa, ia, sa2, ia2, T2b.cptr(), V2_sym.cptr(), W28caa_sigma_coov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a2, "active"] [notNeeded]
  } // End ia2
  } // End sa2
  FC_FUNC(g_if_sigma_coov_cooo_no14_x1_type1_eri_v,G_IF_SIGMA_COOV_COOO_NO14_X1_TYPE1_ERI_V)
    (sa, ia, W28caa_sigma_coov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   15] -- 
  // |-- [    0] --| W31(w,a0,j,a) += (    1.00000000) V2(a,w,c0,a1) T2(c0,a0,j,a1) 
  // |-- [    1] --| S2(w,i,j,a) += (    1.00000000) D1(i,a0) W31(w,a0,j,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W31caa_sigma_coov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_coov_cooo_no15_x0_type1_eri_v,G_IF_SIGMA_COOV_COOO_NO15_X0_TYPE1_ERI_V)
      (sa, ia, sa1, ia1, T2b.cptr(), V2_sym.cptr(), W31caa_sigma_coov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  FC_FUNC(g_if_sigma_coov_cooo_no15_x1_type1_eri_v,G_IF_SIGMA_COOV_COOO_NO15_X1_TYPE1_ERI_V)
    (sa, ia, W31caa_sigma_coov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   16] -- 
  // |-- [    0] --| W33(w,a0,j,a) += (    1.00000000) V2(a,a1,w,c0) T2(c0,a0,j,a1) 
  // |-- [    1] --| S2(w,i,j,a) += (   -2.00000000) D1(i,a0) W33(w,a0,j,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W33caa_sigma_coov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_coov_cooo_no16_x0_type1_eri_v,G_IF_SIGMA_COOV_COOO_NO16_X0_TYPE1_ERI_V)
      (sa, ia, sa1, ia1, T2b.cptr(), V2_sym.cptr(), W33caa_sigma_coov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  FC_FUNC(g_if_sigma_coov_cooo_no16_x1_type1_eri_v,G_IF_SIGMA_COOV_COOO_NO16_X1_TYPE1_ERI_V)
    (sa, ia, W33caa_sigma_coov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   17] -- 
  // |-- [    0] --| W36(w,a0,a3,a) += (    1.00000000) V2(a,a2,a3,a1) T2(w,a0,a2,a1) 
  // |-- [    1] --| S2(w,i,j,a) += (   -1.00000000) D2(j,i,a0,a3) W36(w,a0,a3,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W36caa_sigma_coov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_coov_cooo_no17_x0_type1_eri_v,G_IF_SIGMA_COOV_COOO_NO17_X0_TYPE1_ERI_V)
      (sa, ia, sa1, ia1, T2b.cptr(), V2_sym.cptr(), W36caa_sigma_coov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  FC_FUNC(g_if_sigma_coov_cooo_no17_x1_type1_eri_v,G_IF_SIGMA_COOV_COOO_NO17_X1_TYPE1_ERI_V)
    (sa, ia, W36caa_sigma_coov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   18] -- 
  // |-- [    0] --| W37(w,a0,a3,a) += (    1.00000000) V2(a,a2,a3,a1) T2(w,a0,a1,a2) 
  // |-- [    1] --| S2(w,i,j,a) += (   -1.00000000) D2(j,a3,a0,i) W37(w,a0,a3,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W37caa_sigma_coov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  for(int sa2 = 0;sa2 < nir;++sa2){ 
  for(int ia2 = symblockinfo.psym()(sa2,I_O,I_BEGIN);ia2 <= symblockinfo.psym()(sa2,I_O,I_END);++ia2){ 
    T2b = T2.get_amp2(ia2);
    FC_FUNC(g_if_sigma_coov_cooo_no18_x0_type1_eri_v,G_IF_SIGMA_COOV_COOO_NO18_X0_TYPE1_ERI_V)
      (sa, ia, sa2, ia2, T2b.cptr(), V2_sym.cptr(), W37caa_sigma_coov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a2, "active"] [notNeeded]
  } // End ia2
  } // End sa2
  FC_FUNC(g_if_sigma_coov_cooo_no18_x1_type1_eri_v,G_IF_SIGMA_COOV_COOO_NO18_X1_TYPE1_ERI_V)
    (sa, ia, W37caa_sigma_coov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   19] -- 
  // |-- [    0] --| W38(w,a0,j,a) += (    1.00000000) V2(a,a2,j,a1) T2(w,a0,a1,a2) 
  // |-- [    1] --| S2(w,i,j,a) += (    2.00000000) D1(i,a0) W38(w,a0,j,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W38caa_sigma_coov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  for(int sa2 = 0;sa2 < nir;++sa2){ 
  for(int ia2 = symblockinfo.psym()(sa2,I_O,I_BEGIN);ia2 <= symblockinfo.psym()(sa2,I_O,I_END);++ia2){ 
    T2b = T2.get_amp2(ia2);
    FC_FUNC(g_if_sigma_coov_cooo_no19_x0_type1_eri_v,G_IF_SIGMA_COOV_COOO_NO19_X0_TYPE1_ERI_V)
      (sa, ia, sa2, ia2, T2b.cptr(), V2_sym.cptr(), W38caa_sigma_coov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a2, "active"] [notNeeded]
  } // End ia2
  } // End sa2
  FC_FUNC(g_if_sigma_coov_cooo_no19_x1_type1_eri_v,G_IF_SIGMA_COOV_COOO_NO19_X1_TYPE1_ERI_V)
    (sa, ia, W38caa_sigma_coov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   20] -- 
  // |-- [    0] --| W39(w,a0,j,a) += (    1.00000000) V2(a,a2,j,a1) T2(w,a0,a2,a1) 
  // |-- [    1] --| S2(w,i,j,a) += (   -1.00000000) D1(i,a0) W39(w,a0,j,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W39caa_sigma_coov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_coov_cooo_no20_x0_type1_eri_v,G_IF_SIGMA_COOV_COOO_NO20_X0_TYPE1_ERI_V)
      (sa, ia, sa1, ia1, T2b.cptr(), V2_sym.cptr(), W39caa_sigma_coov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  FC_FUNC(g_if_sigma_coov_cooo_no20_x1_type1_eri_v,G_IF_SIGMA_COOV_COOO_NO20_X1_TYPE1_ERI_V)
    (sa, ia, W39caa_sigma_coov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- Title : sigma_coov_cooo

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| S2(w,i,j,a) += (   -1.00000000) T2(w,a2,a1,a0) C5(a0,a2,a1,j,i,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_coov_cooo_no0_x0_type1_d4c_v,G_IF_SIGMA_COOV_COOO_NO0_X0_TYPE1_D4C_V)
      (sa, ia, sa0, ia0, C5.cptr(), T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ia, S2b);
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
  orz::mr::Femto::my_timer.push_back(boost::make_tuple("sigma_coov_cooo", time_sigma.elapsed_cputime(), time_sigma.elapsed_wallclocktime()));
  orz::mr::Femto::file_timing << "* " << boost::format("%20s : %10.7f %10.7f ") % orz::mr::Femto::my_timer.back().get<0>() % orz::mr::Femto::my_timer.back().get<1>() % orz::mr::Femto::my_timer.back().get<2>() << endl;
  flush(orz::mr::Femto::file_timing);

  return retval; 
} 
