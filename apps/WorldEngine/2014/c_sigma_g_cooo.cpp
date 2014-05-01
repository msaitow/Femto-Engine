                                                                                
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
#include <sci/icmr/Femto/elems/c_sigma_g_cooo.h>                                  
                                                                                
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

//                                   Generated date : Sun Apr 20 10:26:26 2014

                                                                                
// ***************************************************************************  
// orz::mr::mrci                                                                
// ***************************************************************************  
									    /*!        
   @brief CT input                                                              
									     */        
                                                                                
double orz::mr::Femto::sigma_g_cooo(const orz::mr::Input &ctinp,                                                  
                                  const orz::mr::SymBlockInfo &symblockinfo,                                 
                                  const orz::mr::HintMO &hintmo,                                             
                                  const double init_value,                                                   
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
                                                                                                                 
  double S0 = init_value;                                                               
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
  // |-- [    0] --| W0(c0,a2,a1,a0) += (    1.00000000) D2(a3,a1,a2,a0) Fc1(c0,a3) 
  // |-- [    1] --| S0() += (   -1.00000000) T2(c0,a2,a1,a0) W0(c0,a2,a1,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W0caa_sigma_g_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa0));
    FC_FUNC(g_if_sigma_g_cooo_no0_x0_type0_noeri,G_IF_SIGMA_G_COOO_NO0_X0_TYPE0_NOERI)
      (sa0, ia0, W0caa_sigma_g_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_g_cooo_no0_x1_type0_noeri,G_IF_SIGMA_G_COOO_NO0_X1_TYPE0_NOERI)
      (sa0, ia0, T2b.cptr(), W0caa_sigma_g_cooo.cptr(), &S0, nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W1(c0,a1,a0,a2) += (    1.00000000) D1(a1,a0) Fc1(c0,a2) 
  // |-- [    1] --| S0() += (    2.00000000) T2(c0,a1,a2,a0) W1(c0,a1,a0,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W1caa_sigma_g_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa0));
    FC_FUNC(g_if_sigma_g_cooo_no1_x0_type0_noeri,G_IF_SIGMA_G_COOO_NO1_X0_TYPE0_NOERI)
      (sa0, ia0, W1caa_sigma_g_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_g_cooo_no1_x1_type0_noeri,G_IF_SIGMA_G_COOO_NO1_X1_TYPE0_NOERI)
      (sa0, ia0, T2b.cptr(), W1caa_sigma_g_cooo.cptr(), &S0, nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  if(myrank == 0)
  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W2(c0,a1,a0,a2) += (    1.00000000) D1(a1,a0) Fc1(c0,a2) 
  // |-- [    1] --| S0() += (   -1.00000000) T2(c0,a1,a0,a2) W2(c0,a1,a0,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa2 = 0;sa2 < nir;++sa2){ 
  for(int ia2 = symblockinfo.psym()(sa2,I_O,I_BEGIN);ia2 <= symblockinfo.psym()(sa2,I_O,I_END);++ia2){ 
    T2b = T2.get_amp2(ia2);
    orz::DTensor W2caa_sigma_g_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa2));
    FC_FUNC(g_if_sigma_g_cooo_no2_x0_type0_noeri,G_IF_SIGMA_G_COOO_NO2_X0_TYPE0_NOERI)
      (sa2, ia2, W2caa_sigma_g_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_g_cooo_no2_x1_type0_noeri,G_IF_SIGMA_G_COOO_NO2_X1_TYPE0_NOERI)
      (sa2, ia2, T2b.cptr(), W2caa_sigma_g_cooo.cptr(), &S0, nir, nsym, psym, &flops);
  // --> @[a2, "active"] [notNeeded]
  } // End ia2
  } // End sa2
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
  // -- Title : sigma_g_cooo

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W3(c0,a2,a1,a0) += (    1.00000000) V2(c0,a4,a5,a3) D3(a5,a3,a4,a1,a2,a0) 
  // |-- [    1] --| S0() += (   -1.00000000) T2(c0,a2,a1,a0) W3(c0,a2,a1,a0) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W3aa_sigma_g_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sa0));
    FC_FUNC(g_if_sigma_g_cooo_no0_x0_type1_eri_c,G_IF_SIGMA_G_COOO_NO0_X0_TYPE1_ERI_C)
      (sa0, ia0, sc0, ic0, V2_sym.cptr(), W3aa_sigma_g_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_g_cooo_no0_x1_type1_eri_c,G_IF_SIGMA_G_COOO_NO0_X1_TYPE1_ERI_C)
      (sa0, ia0, sc0, ic0, T2b.cptr(), W3aa_sigma_g_cooo.cptr(), &S0, nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W4(c0,a1,a0,a3) += (    1.00000000) V2(c0,a3,a4,a2) D2(a4,a2,a1,a0) 
  // |-- [    1] --| S0() += (    2.00000000) T2(c0,a1,a3,a0) W4(c0,a1,a0,a3) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W4aa_sigma_g_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sa0));
    FC_FUNC(g_if_sigma_g_cooo_no1_x0_type1_eri_c,G_IF_SIGMA_G_COOO_NO1_X0_TYPE1_ERI_C)
      (sa0, ia0, sc0, ic0, V2_sym.cptr(), W4aa_sigma_g_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_g_cooo_no1_x1_type1_eri_c,G_IF_SIGMA_G_COOO_NO1_X1_TYPE1_ERI_C)
      (sa0, ia0, sc0, ic0, T2b.cptr(), W4aa_sigma_g_cooo.cptr(), &S0, nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W5(c0,a1,a0,a2) += (    1.00000000) V2(c0,a3,a4,a2) D2(a4,a3,a1,a0) 
  // |-- [    1] --| S0() += (   -1.00000000) T2(c0,a0,a2,a1) W5(c0,a1,a0,a2) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W5aa_sigma_g_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sa1));
    FC_FUNC(g_if_sigma_g_cooo_no2_x0_type1_eri_c,G_IF_SIGMA_G_COOO_NO2_X0_TYPE1_ERI_C)
      (sa1, ia1, sc0, ic0, V2_sym.cptr(), W5aa_sigma_g_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_g_cooo_no2_x1_type1_eri_c,G_IF_SIGMA_G_COOO_NO2_X1_TYPE1_ERI_C)
      (sa1, ia1, sc0, ic0, T2b.cptr(), W5aa_sigma_g_cooo.cptr(), &S0, nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W6(c0,a1,a0,a3) += (    1.00000000) V2(c0,a3,a4,a2) D2(a4,a2,a1,a0) 
  // |-- [    1] --| S0() += (   -1.00000000) T2(c0,a1,a0,a3) W6(c0,a1,a0,a3) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa3 = 0;sa3 < nir;++sa3){ 
  for(int ia3 = symblockinfo.psym()(sa3,I_O,I_BEGIN);ia3 <= symblockinfo.psym()(sa3,I_O,I_END);++ia3){ 
    T2b = T2.get_amp2(ia3);
    orz::DTensor W6aa_sigma_g_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sa3));
    FC_FUNC(g_if_sigma_g_cooo_no3_x0_type1_eri_c,G_IF_SIGMA_G_COOO_NO3_X0_TYPE1_ERI_C)
      (sa3, ia3, sc0, ic0, V2_sym.cptr(), W6aa_sigma_g_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_g_cooo_no3_x1_type1_eri_c,G_IF_SIGMA_G_COOO_NO3_X1_TYPE1_ERI_C)
      (sa3, ia3, sc0, ic0, T2b.cptr(), W6aa_sigma_g_cooo.cptr(), &S0, nir, nsym, psym, &flops);
  // --> @[a3, "active"] [notNeeded]
  } // End ia3
  } // End sa3
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W7(c0,a1,a0,a2) += (    1.00000000) V2(c0,a3,a4,a2) D2(a4,a0,a1,a3) 
  // |-- [    1] --| S0() += (   -1.00000000) T2(c0,a0,a1,a2) W7(c0,a1,a0,a2) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa2 = 0;sa2 < nir;++sa2){ 
  for(int ia2 = symblockinfo.psym()(sa2,I_O,I_BEGIN);ia2 <= symblockinfo.psym()(sa2,I_O,I_END);++ia2){ 
    T2b = T2.get_amp2(ia2);
    orz::DTensor W7aa_sigma_g_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sa2));
    FC_FUNC(g_if_sigma_g_cooo_no4_x0_type1_eri_c,G_IF_SIGMA_G_COOO_NO4_X0_TYPE1_ERI_C)
      (sa2, ia2, sc0, ic0, V2_sym.cptr(), W7aa_sigma_g_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_g_cooo_no4_x1_type1_eri_c,G_IF_SIGMA_G_COOO_NO4_X1_TYPE1_ERI_C)
      (sa2, ia2, sc0, ic0, T2b.cptr(), W7aa_sigma_g_cooo.cptr(), &S0, nir, nsym, psym, &flops);
  // --> @[a2, "active"] [notNeeded]
  } // End ia2
  } // End sa2
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W8(c0,a0,a2,a1) += (    1.00000000) V2(c0,a2,a3,a1) D1(a3,a0) 
  // |-- [    1] --| S0() += (    2.00000000) T2(c0,a0,a2,a1) W8(c0,a0,a2,a1) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    orz::DTensor W8aa_sigma_g_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sa1));
    FC_FUNC(g_if_sigma_g_cooo_no5_x0_type1_eri_c,G_IF_SIGMA_G_COOO_NO5_X0_TYPE1_ERI_C)
      (sa1, ia1, sc0, ic0, V2_sym.cptr(), W8aa_sigma_g_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_g_cooo_no5_x1_type1_eri_c,G_IF_SIGMA_G_COOO_NO5_X1_TYPE1_ERI_C)
      (sa1, ia1, sc0, ic0, T2b.cptr(), W8aa_sigma_g_cooo.cptr(), &S0, nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W9(c0,a0,a2,a1) += (    1.00000000) V2(c0,a2,a3,a1) D1(a3,a0) 
  // |-- [    1] --| S0() += (   -1.00000000) T2(c0,a0,a1,a2) W9(c0,a0,a2,a1) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa2 = 0;sa2 < nir;++sa2){ 
  for(int ia2 = symblockinfo.psym()(sa2,I_O,I_BEGIN);ia2 <= symblockinfo.psym()(sa2,I_O,I_END);++ia2){ 
    T2b = T2.get_amp2(ia2);
    orz::DTensor W9aa_sigma_g_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sa2));
    FC_FUNC(g_if_sigma_g_cooo_no6_x0_type1_eri_c,G_IF_SIGMA_G_COOO_NO6_X0_TYPE1_ERI_C)
      (sa2, ia2, sc0, ic0, V2_sym.cptr(), W9aa_sigma_g_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_g_cooo_no6_x1_type1_eri_c,G_IF_SIGMA_G_COOO_NO6_X1_TYPE1_ERI_C)
      (sa2, ia2, sc0, ic0, T2b.cptr(), W9aa_sigma_g_cooo.cptr(), &S0, nir, nsym, psym, &flops);
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

//-@ERI.contractions(end)

//-@D4C.contractions(begin)

//-@D4C.contractions(end)

  // Do timing!
  orz::mr::Femto::my_timer.push_back(boost::make_tuple("sigma_g_cooo", time_sigma.elapsed_cputime(), time_sigma.elapsed_wallclocktime()));
  orz::mr::Femto::file_timing << "* " << boost::format("%20s : %10.7f %10.7f ") % orz::mr::Femto::my_timer.back().get<0>() % orz::mr::Femto::my_timer.back().get<1>() % orz::mr::Femto::my_timer.back().get<2>() << endl;
  flush(orz::mr::Femto::file_timing);

  double sum_S0 = 0.0;
  boost::mpi::reduce(orz::world(),S0, sum_S0, std::plus<double>(), 0);
  S0 = sum_S0;

  return  S0;
} 
