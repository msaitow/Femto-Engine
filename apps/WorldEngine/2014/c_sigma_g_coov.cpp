                                                                                
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
#include <sci/icmr/Femto/elems/c_sigma_g_coov.h>                                  
                                                                                
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
                                                                                
double orz::mr::Femto::sigma_g_coov(const orz::mr::Input &ctinp,                                                  
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
  // |-- [    0] --| W0(c0,a1,a0,v0) += (    1.00000000) D1(a1,a0) Fc1(v0,c0) 
  // |-- [    1] --| S0() += (   -1.00000000) T2(c0,a1,a0,v0) W0(c0,a1,a0,v0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sv0 = 0;sv0 < nir;++sv0){ 
  for(int iv0 = symblockinfo.psym()(sv0,I_V,I_BEGIN);iv0 <= symblockinfo.psym()(sv0,I_V,I_END);++iv0){ 
    T2b = T2.get_amp2(iv0);
    orz::DTensor W0caa_sigma_g_coov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sv0));
    FC_FUNC(g_if_sigma_g_coov_no0_x0_type0_noeri,G_IF_SIGMA_G_COOV_NO0_X0_TYPE0_NOERI)
      (sv0, iv0, W0caa_sigma_g_coov.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_g_coov_no0_x1_type0_noeri,G_IF_SIGMA_G_COOV_NO0_X1_TYPE0_NOERI)
      (sv0, iv0, T2b.cptr(), W0caa_sigma_g_coov.cptr(), &S0, nir, nsym, psym, &flops);
  // --> @[v0, "virtual"] [notNeeded]
  } // End iv0
  } // End sv0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  } // End Femto
  //*-- FEMTO ends --//*

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
  // -- Title : sigma_g_coov

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W1(c0,a1,a0,v0) += (    1.00000000) V2(v0,c0,a3,a2) D2(a3,a2,a1,a0) 
  // |-- [    1] --| S0() += (   -1.00000000) T2(c0,a1,a0,v0) W1(c0,a1,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W1caa_sigma_g_coov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_g_coov_no0_x0_type1_eri_v,G_IF_SIGMA_G_COOV_NO0_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W1caa_sigma_g_coov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_g_coov_no0_x1_type1_eri_v,G_IF_SIGMA_G_COOV_NO0_X1_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W1caa_sigma_g_coov.cptr(), &S0, nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W2(c0,a1,a0,v0) += (    1.00000000) V2(v0,a3,c0,a2) D2(a3,a0,a1,a2) 
  // |-- [    1] --| S0() += (   -1.00000000) T2(c0,a0,a1,v0) W2(c0,a1,a0,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W2caa_sigma_g_coov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_g_coov_no1_x0_type1_eri_v,G_IF_SIGMA_G_COOV_NO1_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W2caa_sigma_g_coov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_g_coov_no1_x1_type1_eri_v,G_IF_SIGMA_G_COOV_NO1_X1_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W2caa_sigma_g_coov.cptr(), &S0, nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W3(c0,a0,a1,v0) += (    1.00000000) V2(v0,a2,c0,a1) D1(a2,a0) 
  // |-- [    1] --| S0() += (    2.00000000) T2(c0,a0,a1,v0) W3(c0,a0,a1,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W3caa_sigma_g_coov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_g_coov_no2_x0_type1_eri_v,G_IF_SIGMA_G_COOV_NO2_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W3caa_sigma_g_coov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_g_coov_no2_x1_type1_eri_v,G_IF_SIGMA_G_COOV_NO2_X1_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W3caa_sigma_g_coov.cptr(), &S0, nir, nsym, psym, &flops);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W4(c0,a0,a1,v0) += (    1.00000000) V2(v0,c0,a2,a1) D1(a2,a0) 
  // |-- [    1] --| S0() += (   -1.00000000) T2(c0,a0,a1,v0) W4(c0,a0,a1,v0) 
  int sv0(s_eri);
  int iv0(i_eri);
  double flops(0); // Flop count  
  // Pref: 2
  T2b = T2.get_amp2(iv0);
  orz::DTensor W4caa_sigma_g_coov(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sv0));
  FC_FUNC(g_if_sigma_g_coov_no3_x0_type1_eri_v,G_IF_SIGMA_G_COOV_NO3_X0_TYPE1_ERI_V)
    (sv0, iv0, V2_sym.cptr(), W4caa_sigma_g_coov.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_g_coov_no3_x1_type1_eri_v,G_IF_SIGMA_G_COOV_NO3_X1_TYPE1_ERI_V)
    (sv0, iv0, T2b.cptr(), W4caa_sigma_g_coov.cptr(), &S0, nir, nsym, psym, &flops);
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
  orz::mr::Femto::my_timer.push_back(boost::make_tuple("sigma_g_coov", time_sigma.elapsed_cputime(), time_sigma.elapsed_wallclocktime()));
  orz::mr::Femto::file_timing << "* " << boost::format("%20s : %10.7f %10.7f ") % orz::mr::Femto::my_timer.back().get<0>() % orz::mr::Femto::my_timer.back().get<1>() % orz::mr::Femto::my_timer.back().get<2>() << endl;
  flush(orz::mr::Femto::file_timing);

  double sum_S0 = 0.0;
  boost::mpi::reduce(orz::world(),S0, sum_S0, std::plus<double>(), 0);
  S0 = sum_S0;

  return  S0;
} 
