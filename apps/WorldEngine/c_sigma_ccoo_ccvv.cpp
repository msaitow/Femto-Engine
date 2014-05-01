                                                                                
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
#include <sci/icmr/Femto/elems/c_sigma_ccoo_ccvv.h>                                  
                                                                                
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

//                                   Generated date : Sun Apr 20 10:26:23 2014

                                                                                
// ***************************************************************************  
// orz::mr::mrci                                                                
// ***************************************************************************  
									    /*!        
   @brief CT input                                                              
									     */        
                                                                                
orz::mr::BareAmpPack orz::mr::Femto::sigma_ccoo_ccvv(const orz::mr::Input &ctinp,                                    
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
  // -- Title : sigma_ccoo_ccvv

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W0(j,i,v1,v0) += (    1.00000000) V2(v1,a1,v0,a0) D2(j,a1,i,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) T2(w,x,v0,v1) W0(j,i,v1,v0) 
  int sv1(s_eri);
  int iv1(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv1);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W0av_sigma_ccoo_ccvv(orz::mr::sizeof_sympack_Xav(symblockinfo, sj^sv1));
    FC_FUNC(g_if_sigma_ccoo_ccvv_no0_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCVV_NO0_X0_TYPE1_ERI_V)
      (sj, ij, sv1, iv1, V2_sym.cptr(), W0av_sigma_ccoo_ccvv.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccvv_no0_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCVV_NO0_X1_TYPE1_ERI_V)
      (sj, ij, sv1, iv1, T2b.cptr(), W0av_sigma_ccoo_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W1(j,i,v1,v0) += (    1.00000000) V2(v1,i,v0,a0) D1(j,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) T2(w,x,v0,v1) W1(j,i,v1,v0) 
  int sv1(s_eri);
  int iv1(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv1);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W1av_sigma_ccoo_ccvv(orz::mr::sizeof_sympack_Xav(symblockinfo, sj^sv1));
    FC_FUNC(g_if_sigma_ccoo_ccvv_no1_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCVV_NO1_X0_TYPE1_ERI_V)
      (sj, ij, sv1, iv1, V2_sym.cptr(), W1av_sigma_ccoo_ccvv.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccvv_no1_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCVV_NO1_X1_TYPE1_ERI_V)
      (sj, ij, sv1, iv1, T2b.cptr(), W1av_sigma_ccoo_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W2(j,i,v1,v0) += (    1.00000000) V2(v1,i,v0,a0) D1(j,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -4.00000000) T2(x,w,v0,v1) W2(j,i,v1,v0) 
  int sv1(s_eri);
  int iv1(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv1);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W2av_sigma_ccoo_ccvv(orz::mr::sizeof_sympack_Xav(symblockinfo, sj^sv1));
    FC_FUNC(g_if_sigma_ccoo_ccvv_no2_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCVV_NO2_X0_TYPE1_ERI_V)
      (sj, ij, sv1, iv1, V2_sym.cptr(), W2av_sigma_ccoo_ccvv.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccvv_no2_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCVV_NO2_X1_TYPE1_ERI_V)
      (sj, ij, sv1, iv1, T2b.cptr(), W2av_sigma_ccoo_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W3(i,j,v1,v0) += (    1.00000000) V2(v1,j,v0,a0) D1(i,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (   -4.00000000) T2(w,x,v0,v1) W3(i,j,v1,v0) 
  int sv1(s_eri);
  int iv1(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv1);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W3av_sigma_ccoo_ccvv(orz::mr::sizeof_sympack_Xav(symblockinfo, sj^sv1));
    FC_FUNC(g_if_sigma_ccoo_ccvv_no3_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCVV_NO3_X0_TYPE1_ERI_V)
      (sj, ij, sv1, iv1, V2_sym.cptr(), W3av_sigma_ccoo_ccvv.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccvv_no3_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCVV_NO3_X1_TYPE1_ERI_V)
      (sj, ij, sv1, iv1, T2b.cptr(), W3av_sigma_ccoo_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W4(i,j,v1,v0) += (    1.00000000) V2(v1,j,v0,a0) D1(i,a0) 
  // |-- [    1] --| S2(w,x,i,j) += (    2.00000000) T2(x,w,v0,v1) W4(i,j,v1,v0) 
  int sv1(s_eri);
  int iv1(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv1);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    orz::DTensor W4av_sigma_ccoo_ccvv(orz::mr::sizeof_sympack_Xav(symblockinfo, sj^sv1));
    FC_FUNC(g_if_sigma_ccoo_ccvv_no4_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCVV_NO4_X0_TYPE1_ERI_V)
      (sj, ij, sv1, iv1, V2_sym.cptr(), W4av_sigma_ccoo_ccvv.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ccoo_ccvv_no4_x1_type1_eri_v,G_IF_SIGMA_CCOO_CCVV_NO4_X1_TYPE1_ERI_V)
      (sj, ij, sv1, iv1, T2b.cptr(), W4av_sigma_ccoo_ccvv.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| S2(w,x,i,j) += (    8.00000000) V2(v1,j,v0,i) T2(w,x,v0,v1) 
  int sv1(s_eri);
  int iv1(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv1);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccvv_no5_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCVV_NO5_X0_TYPE1_ERI_V)
      (sj, ij, sv1, iv1, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[j, "active"] [notNeeded]
    retval.acc_amp2(ij, S2b);
  } // End ij
  } // End sj
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| S2(w,x,i,j) += (   -4.00000000) V2(v1,i,v0,j) T2(w,x,v0,v1) 
  int sv1(s_eri);
  int iv1(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  T2b = T2.get_amp2(iv1);
  for(int sj = 0;sj < nir;++sj){ 
  for(int ij = symblockinfo.psym()(sj,I_O,I_BEGIN);ij <= symblockinfo.psym()(sj,I_O,I_END);++ij){ 
    S2b = orz::DTensor(retval.namps_iamp()[ij]);
    FC_FUNC(g_if_sigma_ccoo_ccvv_no6_x0_type1_eri_v,G_IF_SIGMA_CCOO_CCVV_NO6_X0_TYPE1_ERI_V)
      (sj, ij, sv1, iv1, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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

//-@loadERI(v,end)

//-@ERI.contractions(end)

//-@D4C.contractions(begin)

//-@D4C.contractions(end)

  // Do timing!
  orz::mr::Femto::my_timer.push_back(boost::make_tuple("sigma_ccoo_ccvv", time_sigma.elapsed_cputime(), time_sigma.elapsed_wallclocktime()));
  orz::mr::Femto::file_timing << "* " << boost::format("%20s : %10.7f %10.7f ") % orz::mr::Femto::my_timer.back().get<0>() % orz::mr::Femto::my_timer.back().get<1>() % orz::mr::Femto::my_timer.back().get<2>() << endl;
  flush(orz::mr::Femto::file_timing);

  return retval; 
} 
