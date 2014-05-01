                                                                                
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
#include <sci/icmr/Femto/elems/c_sigma_oovv_coov.h>                                  
                                                                                
using std::cout;                                                                
using std::endl;                                                                
                                                                                
#define FLOPCOUNT                                                               
                                                                                
using orz::mr::Femto::Fc0;                                                      
using orz::mr::Femto::h1_int;                                                   
using orz::mr::Femto::h6_int;                                                   
                                                                                
//      ______                  __           
//     / ____/___   ____ ___   / /_ ____     
//    / /_   / _ \ / __ `__ \ / __// __ \ 
//   / __/  /  __// / / / / // /_ / /_/ /    
//  /_/     \___//_/ /_/ /_/ \__/ \____/  

//                                   Generated date : Sun Apr 20 10:26:10 2014

                                                                                
// ***************************************************************************  
// orz::mr::mrci                                                                
// ***************************************************************************  
									    /*!        
   @brief CT input                                                              
									     */        
                                                                                
orz::mr::BareAmpPack orz::mr::Femto::sigma_oovv_coov(const orz::mr::Input &ctinp,                                    
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
  // -- Title : sigma_oovv_coov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W0(c0,j,i,a1,a0,a) += (    1.00000000) V2(c0,a2,a3,a) D3(j,a0,i,a3,a1,a2) 
  // |-- [    1] --| S2(i,j,a,b) += (   -1.00000000) T2(c0,a0,a1,b) W0(c0,j,i,a1,a0,a) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  orz::DTensor W0aaaav_sigma_oovv_coov(orz::mr::sizeof_sympack_Xaaaav(symblockinfo, sc0));
  FC_FUNC(g_if_sigma_oovv_coov_no0_x0_type1_eri_c,G_IF_SIGMA_OOVV_COOV_NO0_X0_TYPE1_ERI_C)
    (sc0, ic0, V2_sym.cptr(), W0aaaav_sigma_oovv_coov.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_oovv_coov_no0_x1_type1_eri_c,G_IF_SIGMA_OOVV_COOV_NO0_X1_TYPE1_ERI_C)
      (sb, ib, sc0, ic0, T2b.cptr(), W0aaaav_sigma_oovv_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- Title : sigma_oovv_coov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W2(a0,a2,b,a) += (    1.00000000) V2(a2,a,c0,a1) T2(c0,a0,a1,b) 
  // |-- [    1] --| S2(i,j,a,b) += (    2.00000000) D2(j,a0,i,a2) W2(a0,a2,b,a) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    orz::DTensor W2av_sigma_oovv_coov(orz::mr::sizeof_sympack_Xav(symblockinfo, sa2^sb));
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_oovv_coov_no0_x0_type1_eri_o,G_IF_SIGMA_OOVV_COOV_NO0_X0_TYPE1_ERI_O)
      (sa2, ia2, sb, ib, T2b.cptr(), V2_sym.cptr(), W2av_sigma_oovv_coov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_oovv_coov_no0_x1_type1_eri_o,G_IF_SIGMA_OOVV_COOV_NO0_X1_TYPE1_ERI_O)
      (sa2, ia2, sb, ib, W2av_sigma_oovv_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W5(a0,a1,b,a) += (    1.00000000) V2(a1,c0,a2,a) T2(c0,a0,a2,b) 
  // |-- [    1] --| S2(i,j,a,b) += (   -1.00000000) D2(j,a0,i,a1) W5(a0,a1,b,a) 
  int sa1(s_eri);
  int ia1(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    T2b = T2.get_amp2(ib);
    orz::DTensor W5av_sigma_oovv_coov(orz::mr::sizeof_sympack_Xav(symblockinfo, sa1^sb));
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_oovv_coov_no1_x0_type1_eri_o,G_IF_SIGMA_OOVV_COOV_NO1_X0_TYPE1_ERI_O)
      (sa1, ia1, sb, ib, T2b.cptr(), V2_sym.cptr(), W5av_sigma_oovv_coov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_oovv_coov_no1_x1_type1_eri_o,G_IF_SIGMA_OOVV_COOV_NO1_X1_TYPE1_ERI_O)
      (sa1, ia1, sb, ib, W5av_sigma_oovv_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  // -- Title : sigma_oovv_coov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W1(c0,j,i,a1,a0,b) += (    1.00000000) V2(b,a2,c0,a3) D3(j,a2,i,a0,a1,a3) 
  // |-- [    1] --| S2(i,j,a,b) += (   -1.00000000) T2(c0,a0,a1,a) W1(c0,j,i,a1,a0,b) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  orz::DTensor W1caaaa_sigma_oovv_coov(orz::mr::sizeof_sympack_Xcaaaa(symblockinfo, sb));
  FC_FUNC(g_if_sigma_oovv_coov_no0_x0_type1_eri_v,G_IF_SIGMA_OOVV_COOV_NO0_X0_TYPE1_ERI_V)
    (sb, ib, V2_sym.cptr(), W1caaaa_sigma_oovv_coov.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    FC_FUNC(g_if_sigma_oovv_coov_no0_x1_type1_eri_v,G_IF_SIGMA_OOVV_COOV_NO0_X1_TYPE1_ERI_V)
      (sa, ia, sb, ib, T2b.cptr(), W1caaaa_sigma_oovv_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
  } // End ia
  } // End sa
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W3(a0,a1,a,b) += (    1.00000000) V2(b,a1,c0,a2) T2(c0,a0,a2,a) 
  // |-- [    1] --| S2(i,j,a,b) += (    2.00000000) D2(j,a1,i,a0) W3(a0,a1,a,b) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W3aa_sigma_oovv_coov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa^sb));
    FC_FUNC(g_if_sigma_oovv_coov_no1_x0_type1_eri_v,G_IF_SIGMA_OOVV_COOV_NO1_X0_TYPE1_ERI_V)
      (sa, ia, sb, ib, T2b.cptr(), V2_sym.cptr(), W3aa_sigma_oovv_coov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_oovv_coov_no1_x1_type1_eri_v,G_IF_SIGMA_OOVV_COOV_NO1_X1_TYPE1_ERI_V)
      (sa, ia, sb, ib, W3aa_sigma_oovv_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
  } // End ia
  } // End sa
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W4(a0,a2,a,b) += (    1.00000000) V2(b,a1,c0,a2) T2(c0,a0,a1,a) 
  // |-- [    1] --| S2(i,j,a,b) += (   -1.00000000) D2(j,a2,i,a0) W4(a0,a2,a,b) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor W4aa_sigma_oovv_coov(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa^sb));
    FC_FUNC(g_if_sigma_oovv_coov_no2_x0_type1_eri_v,G_IF_SIGMA_OOVV_COOV_NO2_X0_TYPE1_ERI_V)
      (sa, ia, sb, ib, T2b.cptr(), V2_sym.cptr(), W4aa_sigma_oovv_coov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_oovv_coov_no2_x1_type1_eri_v,G_IF_SIGMA_OOVV_COOV_NO2_X1_TYPE1_ERI_V)
      (sa, ia, sb, ib, W4aa_sigma_oovv_coov.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
  } // End ia
  } // End sa
  retval.acc_amp2(ib, S2b);
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
  orz::mr::Femto::my_timer.push_back(boost::make_tuple("sigma_oovv_coov", time_sigma.elapsed_cputime(), time_sigma.elapsed_wallclocktime()));
  orz::mr::Femto::file_timing << "* " << boost::format("%20s : %10.7f %10.7f ") % orz::mr::Femto::my_timer.back().get<0>() % orz::mr::Femto::my_timer.back().get<1>() % orz::mr::Femto::my_timer.back().get<2>() << endl;
  flush(orz::mr::Femto::file_timing);

  return retval; 
} 
