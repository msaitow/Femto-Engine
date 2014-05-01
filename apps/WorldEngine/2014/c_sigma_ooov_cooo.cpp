                                                                                
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
#include <sci/icmr/Femto/elems/c_sigma_ooov_cooo.h>                                  
                                                                                
using std::cout;                                                                
using std::endl;                                                                
                                                                                
#define FLOPCOUNT                                                               
                                                                                
using orz::mr::Femto::Fc0;                                                      
using orz::mr::Femto::h1_int;                                                   
using orz::mr::Femto::h6_int;                                                   
                                                                                
//                                                              
//   _______________                                  ______    
//  |          |                 .'. .`. `````|`````.~      ~.  
//  |______    |______         .'   `   `.    |    |          | 
//  |          |             .'           `.  |    |          | 
//  |          |___________.'               `.|     `.______.'  
//                                                              

//                                   Generated date : Sun Apr 20 10:26:08 2014

                                                                                
// ***************************************************************************  
// orz::mr::mrci                                                                
// ***************************************************************************  
									    /*!        
   @brief CT input                                                              
									     */        
                                                                                
orz::mr::BareAmpPack orz::mr::Femto::sigma_ooov_cooo(const orz::mr::Input &ctinp,                                    
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
  // -- Title : sigma_ooov_cooo

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W12(c0,i,a1,a0) += (    1.00000000) T2(c0,a2,a1,a0) D1(i,a2) 
  // |-- [    1] --| W13(c0,j,k,a3,a4,i) += (    1.00000000) D3(a0,k,a1,a3,j,a4) W12(c0,i,a1,a0) 
  // |-- [    2] --| S2(i,j,k,a) += (    0.50000000) V2(c0,a3,a4,a) W13(c0,j,k,a3,a4,i) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 1 <2, 1> 
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 2
  orz::DTensor W13aaaaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaaaa(symblockinfo, sc0));
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W12aa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sa0));
    FC_FUNC(g_if_sigma_ooov_cooo_no0_x0_type1_eri_c,G_IF_SIGMA_OOOV_COOO_NO0_X0_TYPE1_ERI_C)
      (sa0, ia0, sc0, ic0, T2b.cptr(), W12aa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_cooo_no0_x1_type1_eri_c,G_IF_SIGMA_OOOV_COOO_NO0_X1_TYPE1_ERI_C)
      (sa0, ia0, sc0, ic0, W12aa_sigma_ooov_cooo.cptr(), W13aaaaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_cooo_no0_x2_type1_eri_c,G_IF_SIGMA_OOOV_COOO_NO0_X2_TYPE1_ERI_C)
      (sa, ia, sc0, ic0, V2_sym.cptr(), W13aaaaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W20(c0,j,a1,a0) += (    1.00000000) T2(c0,a2,a1,a0) D1(j,a2) 
  // |-- [    1] --| W21(c0,i,a4,a3,k,j) += (    1.00000000) D3(a0,a4,a1,a3,i,k) W20(c0,j,a1,a0) 
  // |-- [    2] --| S2(i,j,k,a) += (    0.50000000) V2(c0,a3,a4,a) W21(c0,i,a4,a3,k,j) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 1 <2, 1> 
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 2
  orz::DTensor W21aaaaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaaaa(symblockinfo, sc0));
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W20aa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sa0));
    FC_FUNC(g_if_sigma_ooov_cooo_no1_x0_type1_eri_c,G_IF_SIGMA_OOOV_COOO_NO1_X0_TYPE1_ERI_C)
      (sa0, ia0, sc0, ic0, T2b.cptr(), W20aa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_cooo_no1_x1_type1_eri_c,G_IF_SIGMA_OOOV_COOO_NO1_X1_TYPE1_ERI_C)
      (sa0, ia0, sc0, ic0, W20aa_sigma_ooov_cooo.cptr(), W21aaaaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_cooo_no1_x2_type1_eri_c,G_IF_SIGMA_OOOV_COOO_NO1_X2_TYPE1_ERI_C)
      (sa, ia, sc0, ic0, V2_sym.cptr(), W21aaaaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W32(c0,k,a2,a0) += (    1.00000000) T2(c0,a2,a1,a0) D1(a1,k) 
  // |-- [    1] --| W33(c0,j,i,a4,a3,k) += (    1.00000000) D3(a0,a2,j,a4,i,a3) W32(c0,k,a2,a0) 
  // |-- [    2] --| S2(i,j,k,a) += (    0.50000000) V2(c0,a3,a4,a) W33(c0,j,i,a4,a3,k) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 1 <2, 1> 
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 2
  orz::DTensor W33aaaaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaaaa(symblockinfo, sc0));
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W32aa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sc0^sa0));
    FC_FUNC(g_if_sigma_ooov_cooo_no2_x0_type1_eri_c,G_IF_SIGMA_OOOV_COOO_NO2_X0_TYPE1_ERI_C)
      (sa0, ia0, sc0, ic0, T2b.cptr(), W32aa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_cooo_no2_x1_type1_eri_c,G_IF_SIGMA_OOOV_COOO_NO2_X1_TYPE1_ERI_C)
      (sa0, ia0, sc0, ic0, W32aa_sigma_ooov_cooo.cptr(), W33aaaaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_cooo_no2_x2_type1_eri_c,G_IF_SIGMA_OOOV_COOO_NO2_X2_TYPE1_ERI_C)
      (sa, ia, sc0, ic0, V2_sym.cptr(), W33aaaaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a, "virtual"] [notNeeded]
    retval.acc_amp2(ia, S2b);
  } // End ia
  } // End sa
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W40(c0,k,a2,a1) += (    1.00000000) T2(c0,a2,a1,a0) D1(a0,k) 
  // |-- [    1] --| W41(c0,j,i,a3,a4,k) += (    1.00000000) D3(a1,a3,j,a4,i,a2) W40(c0,k,a2,a1) 
  // |-- [    2] --| S2(i,j,k,a) += (    0.50000000) V2(c0,a3,a4,a) W41(c0,j,i,a3,a4,k) 
  int sc0(s_eri);
  int ic0(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 2
  orz::DTensor W40aaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sc0));
  orz::DTensor W41aaaaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaaaa(symblockinfo, sc0));
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no3_x0_type1_eri_c,G_IF_SIGMA_OOOV_COOO_NO3_X0_TYPE1_ERI_C)
      (sa0, ia0, sc0, ic0, T2b.cptr(), W40aaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  FC_FUNC(g_if_sigma_ooov_cooo_no3_x1_type1_eri_c,G_IF_SIGMA_OOOV_COOO_NO3_X1_TYPE1_ERI_C)
    (sc0, ic0, W40aaa_sigma_ooov_cooo.cptr(), W41aaaaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_sigma_ooov_cooo_no3_x2_type1_eri_c,G_IF_SIGMA_OOOV_COOO_NO3_X2_TYPE1_ERI_C)
      (sa, ia, sc0, ic0, V2_sym.cptr(), W41aaaaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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

//-@loadERI(v,begin)
  //*-- FEMTO begins --//*
  // Label : eri_v
  {

//-@type(2).declaration(begin)
  // --  Title : sigma_ooov_cooo
  //  >> Intermediates for the external contractions are defined here << 
  orz::DTensor W10caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W15caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W17caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W19caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W23caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W24caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W28ca_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W34ca_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W42ca_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W44ca_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W46ca_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W48caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W50caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W52caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W54caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W56caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W58caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W60caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W62caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W65caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W67caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W68caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W70caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W72caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W74caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W76caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W79caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W80caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W82caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W84caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W86caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W89caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W91caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W93caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W95caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W96ca_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W98ca_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W100ca_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W102caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W104caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W106caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W108caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W110caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W112caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W114caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W116caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W119caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W121caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W122caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W124caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W126caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W128caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W130caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W133caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W134caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W136caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W138caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W140caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W143caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W145caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W147caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W149caaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaa(symblockinfo, 0));
  orz::DTensor W150caaaaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaaaa(symblockinfo, 0));
  orz::DTensor W151caaaaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaaaaa(symblockinfo, 0));
//-@type(2).declaration(end)

  //*-- Entering to take the type 0 contractions --*//
//-@type(0).contraction(begin)

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W10(c0,j,a3,a4) += (    1.00000000) T2(c0,a2,a1,a0) D3(a0,a2,a1,a3,j,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no0_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO0_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W10caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W15(c0,j,k,a4) += (    1.00000000) T2(c0,a2,a1,a0) D3(a0,a2,a1,k,j,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no1_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO1_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W15caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W17(c0,j,a3,k) += (    1.00000000) T2(c0,a2,a1,a0) D3(a0,a2,a1,a3,j,k) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no2_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO2_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W17caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W19(c0,i,a3,k) += (    1.00000000) T2(c0,a2,a1,a0) D3(a0,a2,a1,a3,i,k) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no3_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO3_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W19caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W23(c0,i,a4,k) += (    1.00000000) T2(c0,a2,a1,a0) D3(a0,a2,a1,a4,i,k) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no4_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO4_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W23caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W24(c0,i,a3,a4) += (    1.00000000) T2(c0,a2,a1,a0) D3(a0,a2,a1,a3,i,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no5_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO5_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W24caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W28(c0,a0) += (    1.00000000) T2(c0,a2,a1,a0) D1(a1,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no6_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO6_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W28ca_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W34(c0,a1) += (    1.00000000) T2(c0,a2,a1,a0) D1(a0,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no7_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO7_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W34ca_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| W42(c0,a3) += (    1.00000000) T2(c0,a2,a1,a0) D2(a0,a2,a1,a3) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no8_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO8_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W42ca_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W44(c0,a4) += (    1.00000000) T2(c0,a2,a1,a0) D2(a0,a2,a1,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no9_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO9_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W44ca_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| W46(c0,k) += (    1.00000000) T2(c0,a2,a1,a0) D2(a0,a2,a1,k) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no10_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO10_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W46ca_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| W48(c0,a4,a3,a2) += (    1.00000000) T2(c0,a2,a1,a0) D2(a0,a4,a1,a3) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no11_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO11_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W48caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| W50(c0,k,a3,a2) += (    1.00000000) T2(c0,a2,a1,a0) D2(a0,k,a1,a3) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no12_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO12_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W50caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| W52(c0,a4,k,a2) += (    1.00000000) T2(c0,a2,a1,a0) D2(a0,a4,a1,k) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no13_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO13_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W52caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   14] -- 
  // |-- [    0] --| W54(c0,a4,k,a2) += (    1.00000000) T2(c0,a2,a1,a0) D2(a0,a4,a1,k) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no14_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO14_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W54caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   15] -- 
  // |-- [    0] --| W56(c0,k,a4,a2) += (    1.00000000) T2(c0,a2,a1,a0) D2(a0,k,a1,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no15_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO15_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W56caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   16] -- 
  // |-- [    0] --| W58(c0,k,a4,a2) += (    1.00000000) T2(c0,a2,a1,a0) D2(a0,k,a1,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no16_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO16_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W58caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   17] -- 
  // |-- [    0] --| W60(c0,j,a4,a1) += (    1.00000000) T2(c0,a2,a1,a0) D2(a0,a2,j,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no17_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO17_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W60caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   18] -- 
  // |-- [    0] --| W62(c0,j,a3,a1) += (    1.00000000) T2(c0,a2,a1,a0) D2(a0,a2,j,a3) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no18_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO18_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W62caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   19] -- 
  // |-- [    0] --| W65(c0,j,k,a1) += (    1.00000000) T2(c0,a2,a1,a0) D2(a0,a2,j,k) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no19_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO19_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W65caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   20] -- 
  // |-- [    0] --| W67(c0,i,k,a0) += (    1.00000000) T2(c0,a2,a1,a0) D2(a1,a2,i,k) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no20_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO20_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W67caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   21] -- 
  // |-- [    0] --| W68(c0,i,a3,a0) += (    1.00000000) T2(c0,a2,a1,a0) D2(a1,a3,i,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no21_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO21_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W68caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   22] -- 
  // |-- [    0] --| W70(c0,i,a4,a0) += (    1.00000000) T2(c0,a2,a1,a0) D2(a1,a2,i,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no22_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO22_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W70caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   23] -- 
  // |-- [    0] --| W72(c0,i,a4,a0) += (    1.00000000) T2(c0,a2,a1,a0) D2(a1,a4,i,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no23_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO23_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W72caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   24] -- 
  // |-- [    0] --| W74(c0,i,a4,a0) += (    1.00000000) T2(c0,a2,a1,a0) D2(a1,a2,i,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no24_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO24_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W74caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   25] -- 
  // |-- [    0] --| W76(c0,i,a4,a0) += (    1.00000000) T2(c0,a2,a1,a0) D2(a1,a4,i,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no25_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO25_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W76caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   26] -- 
  // |-- [    0] --| W79(c0,i,k,a1) += (    1.00000000) T2(c0,a2,a1,a0) D2(a0,a2,i,k) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no26_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO26_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W79caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   27] -- 
  // |-- [    0] --| W80(c0,i,a4,a1) += (    1.00000000) T2(c0,a2,a1,a0) D2(a0,a2,i,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no27_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO27_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W80caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   28] -- 
  // |-- [    0] --| W82(c0,i,a3,a1) += (    1.00000000) T2(c0,a2,a1,a0) D2(a0,a2,i,a3) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no28_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO28_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W82caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   29] -- 
  // |-- [    0] --| W84(c0,j,a3,a0) += (    1.00000000) T2(c0,a2,a1,a0) D2(j,a2,a1,a3) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no29_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO29_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W84caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   30] -- 
  // |-- [    0] --| W86(c0,j,a4,a0) += (    1.00000000) T2(c0,a2,a1,a0) D2(j,a4,a1,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no30_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO30_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W86caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   31] -- 
  // |-- [    0] --| W89(c0,j,k,a0) += (    1.00000000) T2(c0,a2,a1,a0) D2(j,a2,a1,k) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no31_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO31_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W89caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   32] -- 
  // |-- [    0] --| W91(c0,j,k,a0) += (    1.00000000) T2(c0,a2,a1,a0) D2(j,k,a1,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no32_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO32_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W91caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   33] -- 
  // |-- [    0] --| W93(c0,j,k,a0) += (    1.00000000) T2(c0,a2,a1,a0) D2(j,a2,a1,k) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no33_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO33_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W93caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   34] -- 
  // |-- [    0] --| W95(c0,j,k,a0) += (    1.00000000) T2(c0,a2,a1,a0) D2(j,k,a1,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no34_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO34_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W95caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   35] -- 
  // |-- [    0] --| W96(c0,a3) += (    1.00000000) T2(c0,a2,a1,a0) C2(a0,a2,a1,a3) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no35_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO35_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W96ca_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   36] -- 
  // |-- [    0] --| W98(c0,a4) += (    1.00000000) T2(c0,a2,a1,a0) C2(a0,a2,a1,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no36_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO36_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W98ca_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   37] -- 
  // |-- [    0] --| W100(c0,k) += (    1.00000000) T2(c0,a2,a1,a0) C2(a0,a2,a1,k) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no37_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO37_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W100ca_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   38] -- 
  // |-- [    0] --| W102(c0,a4,a3,a2) += (    1.00000000) T2(c0,a2,a1,a0) C2(a0,a4,a1,a3) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no38_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO38_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W102caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   39] -- 
  // |-- [    0] --| W104(c0,k,a3,a2) += (    1.00000000) T2(c0,a2,a1,a0) C2(a0,k,a1,a3) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no39_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO39_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W104caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   40] -- 
  // |-- [    0] --| W106(c0,a4,k,a2) += (    1.00000000) T2(c0,a2,a1,a0) C2(a0,a4,a1,k) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no40_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO40_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W106caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   41] -- 
  // |-- [    0] --| W108(c0,a4,k,a2) += (    1.00000000) T2(c0,a2,a1,a0) C2(a0,a4,a1,k) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no41_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO41_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W108caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   42] -- 
  // |-- [    0] --| W110(c0,k,a4,a2) += (    1.00000000) T2(c0,a2,a1,a0) C2(a0,k,a1,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no42_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO42_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W110caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   43] -- 
  // |-- [    0] --| W112(c0,k,a4,a2) += (    1.00000000) T2(c0,a2,a1,a0) C2(a0,k,a1,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no43_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO43_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W112caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   44] -- 
  // |-- [    0] --| W114(c0,j,a4,a1) += (    1.00000000) T2(c0,a2,a1,a0) C2(a0,a2,j,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no44_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO44_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W114caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   45] -- 
  // |-- [    0] --| W116(c0,j,a3,a1) += (    1.00000000) T2(c0,a2,a1,a0) C2(a0,a2,j,a3) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no45_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO45_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W116caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   46] -- 
  // |-- [    0] --| W119(c0,j,k,a1) += (    1.00000000) T2(c0,a2,a1,a0) C2(a0,a2,j,k) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no46_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO46_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W119caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   47] -- 
  // |-- [    0] --| W121(c0,i,k,a0) += (    1.00000000) T2(c0,a2,a1,a0) C2(a1,a2,i,k) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no47_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO47_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W121caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   48] -- 
  // |-- [    0] --| W122(c0,i,a3,a0) += (    1.00000000) T2(c0,a2,a1,a0) C2(a1,a3,i,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no48_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO48_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W122caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   49] -- 
  // |-- [    0] --| W124(c0,i,a4,a0) += (    1.00000000) T2(c0,a2,a1,a0) C2(a1,a2,i,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no49_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO49_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W124caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   50] -- 
  // |-- [    0] --| W126(c0,i,a4,a0) += (    1.00000000) T2(c0,a2,a1,a0) C2(a1,a4,i,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no50_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO50_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W126caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   51] -- 
  // |-- [    0] --| W128(c0,i,a4,a0) += (    1.00000000) T2(c0,a2,a1,a0) C2(a1,a2,i,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no51_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO51_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W128caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   52] -- 
  // |-- [    0] --| W130(c0,i,a4,a0) += (    1.00000000) T2(c0,a2,a1,a0) C2(a1,a4,i,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no52_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO52_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W130caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   53] -- 
  // |-- [    0] --| W133(c0,i,k,a1) += (    1.00000000) T2(c0,a2,a1,a0) C2(a0,a2,i,k) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no53_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO53_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W133caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   54] -- 
  // |-- [    0] --| W134(c0,i,a4,a1) += (    1.00000000) T2(c0,a2,a1,a0) C2(a0,a2,i,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no54_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO54_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W134caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   55] -- 
  // |-- [    0] --| W136(c0,i,a3,a1) += (    1.00000000) T2(c0,a2,a1,a0) C2(a0,a2,i,a3) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no55_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO55_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W136caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   56] -- 
  // |-- [    0] --| W138(c0,j,a3,a0) += (    1.00000000) T2(c0,a2,a1,a0) C2(a1,a3,j,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no56_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO56_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W138caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   57] -- 
  // |-- [    0] --| W140(c0,j,a4,a0) += (    1.00000000) T2(c0,a2,a1,a0) C2(a1,a2,j,a4) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no57_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO57_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W140caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   58] -- 
  // |-- [    0] --| W143(c0,j,k,a0) += (    1.00000000) T2(c0,a2,a1,a0) C2(a1,k,j,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no58_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO58_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W143caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   59] -- 
  // |-- [    0] --| W145(c0,j,k,a0) += (    1.00000000) T2(c0,a2,a1,a0) C2(a1,a2,j,k) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no59_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO59_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W145caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   60] -- 
  // |-- [    0] --| W147(c0,j,k,a0) += (    1.00000000) T2(c0,a2,a1,a0) C2(a1,k,j,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no60_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO60_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W147caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   61] -- 
  // |-- [    0] --| W149(c0,j,k,a0) += (    1.00000000) T2(c0,a2,a1,a0) C2(a1,a2,j,k) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ooov_cooo_no61_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO61_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W149caaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   62] -- 
  // |-- [    0] --| W150(c0,j,i,a3,a2,k) += (    1.00000000) T2(c0,a0,a1,k) D3(j,a3,i,a0,a1,a2) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sk = 0;sk < nir;++sk){ 
  for(int ik = symblockinfo.psym()(sk,I_O,I_BEGIN);ik <= symblockinfo.psym()(sk,I_O,I_END);++ik){ 
    T2b = T2.get_amp2(ik);
    FC_FUNC(g_if_sigma_ooov_cooo_no62_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO62_X0_TYPE0_ERI_V)
      (sk, ik, T2b.cptr(), W150caaaaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[k, "active"] [notNeeded]
  } // End ik
  } // End sk
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   63] -- 
  // |-- [    0] --| W151(c0,j,i,a3,a2,k) += (    1.00000000) T2(a0,c0,a1,k) D3(j,a3,i,a2,a1,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sk = 0;sk < nir;++sk){ 
  for(int ik = symblockinfo.psym()(sk,I_O,I_BEGIN);ik <= symblockinfo.psym()(sk,I_O,I_END);++ik){ 
    T2b = T2.get_amp2(ik);
    FC_FUNC(g_if_sigma_ooov_cooo_no63_x0_type0_eri_v,G_IF_SIGMA_OOOV_COOO_NO63_X0_TYPE0_ERI_V)
      (sk, ik, T2b.cptr(), W151caaaaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[k, "active"] [notNeeded]
  } // End ik
  } // End sk
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
  // -- Title : sigma_ooov_cooo

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W0(a1,a0,a3,a) += (    1.00000000) V2(a,a3,c0,a2) T2(c0,a1,a2,a0) 
  // |-- [    1] --| S2(i,j,k,a) += (    2.00000000) D3(k,i,a3,j,a1,a0) W0(a1,a0,a3,a) 
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
    orz::DTensor W0aa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_ooov_cooo_no0_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO0_X0_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), V2_sym.cptr(), W0aa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_cooo_no0_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO0_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, W0aa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W1(a1,a0,a2,a) += (    1.00000000) V2(a,a3,c0,a2) T2(c0,a1,a3,a0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -1.00000000) D3(k,i,a2,j,a1,a0) W1(a1,a0,a2,a) 
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
    orz::DTensor W1aa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_ooov_cooo_no1_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO1_X0_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), V2_sym.cptr(), W1aa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_cooo_no1_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO1_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, W1aa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W2(a0,k,a2,a) += (    1.00000000) V2(a,a2,c0,a1) T2(c0,a0,a1,k) 
  // |-- [    1] --| S2(i,j,k,a) += (    2.00000000) D2(j,a2,i,a0) W2(a0,k,a2,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  for(int sk = 0;sk < nir;++sk){ 
  for(int ik = symblockinfo.psym()(sk,I_O,I_BEGIN);ik <= symblockinfo.psym()(sk,I_O,I_END);++ik){ 
    T2b = T2.get_amp2(ik);
    orz::DTensor W2aa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sk^sa));
    FC_FUNC(g_if_sigma_ooov_cooo_no2_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO2_X0_TYPE1_ERI_V)
      (sa, ia, sk, ik, T2b.cptr(), V2_sym.cptr(), W2aa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_cooo_no2_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO2_X1_TYPE1_ERI_V)
      (sa, ia, sk, ik, W2aa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[k, "active"] [notNeeded]
  } // End ik
  } // End sk
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W3(a0,k,a1,a) += (    1.00000000) V2(a,a2,c0,a1) T2(c0,a0,a2,k) 
  // |-- [    1] --| S2(i,j,k,a) += (   -1.00000000) D2(j,a1,i,a0) W3(a0,k,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  for(int sk = 0;sk < nir;++sk){ 
  for(int ik = symblockinfo.psym()(sk,I_O,I_BEGIN);ik <= symblockinfo.psym()(sk,I_O,I_END);++ik){ 
    T2b = T2.get_amp2(ik);
    orz::DTensor W3aa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sk^sa));
    FC_FUNC(g_if_sigma_ooov_cooo_no3_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO3_X0_TYPE1_ERI_V)
      (sa, ia, sk, ik, T2b.cptr(), V2_sym.cptr(), W3aa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_cooo_no3_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO3_X1_TYPE1_ERI_V)
      (sa, ia, sk, ik, W3aa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[k, "active"] [notNeeded]
  } // End ik
  } // End sk
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W4(a1,a0,a3,a) += (    1.00000000) V2(a,a3,c0,a2) T2(c0,a1,a0,a2) 
  // |-- [    1] --| S2(i,j,k,a) += (   -1.00000000) D3(k,i,a3,j,a1,a0) W4(a1,a0,a3,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W4aaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  for(int sa2 = 0;sa2 < nir;++sa2){ 
  for(int ia2 = symblockinfo.psym()(sa2,I_O,I_BEGIN);ia2 <= symblockinfo.psym()(sa2,I_O,I_END);++ia2){ 
    T2b = T2.get_amp2(ia2);
    FC_FUNC(g_if_sigma_ooov_cooo_no4_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO4_X0_TYPE1_ERI_V)
      (sa, ia, sa2, ia2, T2b.cptr(), V2_sym.cptr(), W4aaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a2, "active"] [notNeeded]
  } // End ia2
  } // End sa2
  FC_FUNC(g_if_sigma_ooov_cooo_no4_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO4_X1_TYPE1_ERI_V)
    (sa, ia, W4aaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W5(a1,a0,a2,a) += (    1.00000000) V2(a,a3,c0,a2) T2(c0,a1,a0,a3) 
  // |-- [    1] --| S2(i,j,k,a) += (   -1.00000000) D3(k,i,a2,a0,a1,j) W5(a1,a0,a2,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W5aaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  for(int sa3 = 0;sa3 < nir;++sa3){ 
  for(int ia3 = symblockinfo.psym()(sa3,I_O,I_BEGIN);ia3 <= symblockinfo.psym()(sa3,I_O,I_END);++ia3){ 
    T2b = T2.get_amp2(ia3);
    FC_FUNC(g_if_sigma_ooov_cooo_no5_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO5_X0_TYPE1_ERI_V)
      (sa, ia, sa3, ia3, T2b.cptr(), V2_sym.cptr(), W5aaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a3, "active"] [notNeeded]
  } // End ia3
  } // End sa3
  FC_FUNC(g_if_sigma_ooov_cooo_no5_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO5_X1_TYPE1_ERI_V)
    (sa, ia, W5aaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W6(a0,k,a2,a) += (    1.00000000) V2(a,a2,c0,a1) T2(c0,a0,k,a1) 
  // |-- [    1] --| S2(i,j,k,a) += (   -1.00000000) D2(j,a2,i,a0) W6(a0,k,a2,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W6aaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_ooov_cooo_no6_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO6_X0_TYPE1_ERI_V)
      (sa, ia, sa1, ia1, T2b.cptr(), V2_sym.cptr(), W6aaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  FC_FUNC(g_if_sigma_ooov_cooo_no6_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO6_X1_TYPE1_ERI_V)
    (sa, ia, W6aaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W7(a0,k,a1,a) += (    1.00000000) V2(a,a2,c0,a1) T2(c0,a0,k,a2) 
  // |-- [    1] --| S2(i,j,k,a) += (   -1.00000000) D2(j,a0,i,a1) W7(a0,k,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W7aaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  for(int sa2 = 0;sa2 < nir;++sa2){ 
  for(int ia2 = symblockinfo.psym()(sa2,I_O,I_BEGIN);ia2 <= symblockinfo.psym()(sa2,I_O,I_END);++ia2){ 
    T2b = T2.get_amp2(ia2);
    FC_FUNC(g_if_sigma_ooov_cooo_no7_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO7_X0_TYPE1_ERI_V)
      (sa, ia, sa2, ia2, T2b.cptr(), V2_sym.cptr(), W7aaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a2, "active"] [notNeeded]
  } // End ia2
  } // End sa2
  FC_FUNC(g_if_sigma_ooov_cooo_no7_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO7_X1_TYPE1_ERI_V)
    (sa, ia, W7aaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    8] -- 
  // |-- [    0] --| W8(a0,a) += (    1.00000000) V2(a,a2,c0,a1) T2(c0,a0,a2,a1) 
  // |-- [    1] --| S2(i,j,k,a) += (   -1.00000000) D2(k,i,a0,j) W8(a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W8a_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xa(symblockinfo, sa));
  for(int sa1 = 0;sa1 < nir;++sa1){ 
  for(int ia1 = symblockinfo.psym()(sa1,I_O,I_BEGIN);ia1 <= symblockinfo.psym()(sa1,I_O,I_END);++ia1){ 
    T2b = T2.get_amp2(ia1);
    FC_FUNC(g_if_sigma_ooov_cooo_no8_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO8_X0_TYPE1_ERI_V)
      (sa, ia, sa1, ia1, T2b.cptr(), V2_sym.cptr(), W8a_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a1, "active"] [notNeeded]
  } // End ia1
  } // End sa1
  FC_FUNC(g_if_sigma_ooov_cooo_no8_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO8_X1_TYPE1_ERI_V)
    (sa, ia, W8a_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    9] -- 
  // |-- [    0] --| W9(a0,a) += (    1.00000000) V2(a,a2,c0,a1) T2(c0,a0,a1,a2) 
  // |-- [    1] --| S2(i,j,k,a) += (    2.00000000) D2(k,i,a0,j) W9(a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 0 <2, 1> 
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W9a_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xa(symblockinfo, sa));
  for(int sa2 = 0;sa2 < nir;++sa2){ 
  for(int ia2 = symblockinfo.psym()(sa2,I_O,I_BEGIN);ia2 <= symblockinfo.psym()(sa2,I_O,I_END);++ia2){ 
    T2b = T2.get_amp2(ia2);
    FC_FUNC(g_if_sigma_ooov_cooo_no9_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO9_X0_TYPE1_ERI_V)
      (sa, ia, sa2, ia2, T2b.cptr(), V2_sym.cptr(), W9a_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a2, "active"] [notNeeded]
  } // End ia2
  } // End sa2
  FC_FUNC(g_if_sigma_ooov_cooo_no9_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO9_X1_TYPE1_ERI_V)
    (sa, ia, W9a_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   10] -- 
  // |-- [    0] --| W11(j,a) += (    1.00000000) V2(a,a4,c0,a3) W10(c0,j,a3,a4) 
  // |-- [    1] --| S2(i,j,k,a) += (   -1.00000000) D1(i,k) W11(j,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W11a_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no10_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO10_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W10caaa_sigma_ooov_cooo.cptr(), W11a_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no10_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO10_X1_TYPE1_ERI_V)
    (sa, ia, W11a_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   11] -- 
  // |-- [    0] --| W14(c0,i,a4,a) += (    1.00000000) V2(a,a4,c0,a3) D1(i,a3) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.50000000) W14(c0,i,a4,a) W15(c0,j,k,a4) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W14caa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no11_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO11_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W14caa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no11_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO11_X1_TYPE1_ERI_V)
    (sa, ia, W14caa_sigma_ooov_cooo.cptr(), W15caaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   12] -- 
  // |-- [    0] --| W16(c0,i,a3,a) += (    1.00000000) V2(a,a4,c0,a3) D1(i,a4) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.50000000) W16(c0,i,a3,a) W17(c0,j,a3,k) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W16caa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no12_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO12_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W16caa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no12_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO12_X1_TYPE1_ERI_V)
    (sa, ia, W16caa_sigma_ooov_cooo.cptr(), W17caaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   13] -- 
  // |-- [    0] --| W18(c0,j,a3,a) += (    1.00000000) V2(a,a4,c0,a3) D1(j,a4) 
  // |-- [    1] --| S2(i,j,k,a) += (   -1.00000000) W18(c0,j,a3,a) W19(c0,i,a3,k) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W18caa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no13_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO13_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W18caa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no13_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO13_X1_TYPE1_ERI_V)
    (sa, ia, W18caa_sigma_ooov_cooo.cptr(), W19caaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   14] -- 
  // |-- [    0] --| W22(c0,j,a4,a) += (    1.00000000) V2(a,a4,c0,a3) D1(j,a3) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.50000000) W22(c0,j,a4,a) W23(c0,i,a4,k) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W22caa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no14_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO14_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W22caa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no14_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO14_X1_TYPE1_ERI_V)
    (sa, ia, W22caa_sigma_ooov_cooo.cptr(), W23caaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   15] -- 
  // |-- [    0] --| W25(i,a) += (    1.00000000) V2(a,a4,c0,a3) W24(c0,i,a3,a4) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.50000000) D1(j,k) W25(i,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W25a_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no15_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO15_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W24caaa_sigma_ooov_cooo.cptr(), W25a_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no15_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO15_X1_TYPE1_ERI_V)
    (sa, ia, W25a_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   16] -- 
  // |-- [    0] --| W26(c0,a1,a4,a) += (    1.00000000) V2(a,a4,c0,a3) D1(a1,a3) 
  // |-- [    1] --| W27(a2,a0,a4,a) += (    1.00000000) T2(c0,a2,a1,a0) W26(c0,a1,a4,a) 
  // |-- [    2] --| S2(i,j,k,a) += (   -1.00000000) D3(a0,a2,j,a4,i,k) W27(a2,a0,a4,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W26caa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no16_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO16_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W26caa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W27aa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_ooov_cooo_no16_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO16_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), W26caa_sigma_ooov_cooo.cptr(), W27aa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_cooo_no16_x2_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO16_X2_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, W27aa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   17] -- 
  // |-- [    0] --| W29(a4,a3,a0,a) += (    1.00000000) V2(a,a4,c0,a3) W28(c0,a0) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.50000000) D3(a0,a3,j,a4,i,k) W29(a4,a3,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W29aaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no17_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO17_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W28ca_sigma_ooov_cooo.cptr(), W29aaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no17_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO17_X1_TYPE1_ERI_V)
    (sa, ia, W29aaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   18] -- 
  // |-- [    0] --| W30(c0,a1,a3,a) += (    1.00000000) V2(a,a4,c0,a3) D1(a1,a4) 
  // |-- [    1] --| W31(a2,a0,a3,a) += (    1.00000000) T2(c0,a2,a1,a0) W30(c0,a1,a3,a) 
  // |-- [    2] --| S2(i,j,k,a) += (    0.50000000) D3(a0,a2,j,a3,i,k) W31(a2,a0,a3,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 2 (alloc)
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W30caa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no18_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO18_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W30caa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W31aa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaa(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_ooov_cooo_no18_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO18_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), W30caa_sigma_ooov_cooo.cptr(), W31aa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_cooo_no18_x2_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO18_X2_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, W31aa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   19] -- 
  // |-- [    0] --| W35(a4,a3,a1,a) += (    1.00000000) V2(a,a4,c0,a3) W34(c0,a1) 
  // |-- [    1] --| S2(i,j,k,a) += (   -1.00000000) D3(a1,a3,j,a4,i,k) W35(a4,a3,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W35aaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no19_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO19_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W34ca_sigma_ooov_cooo.cptr(), W35aaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no19_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO19_X1_TYPE1_ERI_V)
    (sa, ia, W35aaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   20] -- 
  // |-- [    0] --| W36(c0,a0,a4,a) += (    1.00000000) V2(a,a4,c0,a3) D1(a0,a3) 
  // |-- [    1] --| W37(a2,a1,a4,a) += (    1.00000000) T2(c0,a2,a1,a0) W36(c0,a0,a4,a) 
  // |-- [    2] --| S2(i,j,k,a) += (    0.50000000) D3(a1,a2,j,a4,i,k) W37(a2,a1,a4,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 1 <2, 1> 
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W37aaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W36ca_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_ooov_cooo_no20_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO20_X0_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, V2_sym.cptr(), W36ca_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_cooo_no20_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO20_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), W36ca_sigma_ooov_cooo.cptr(), W37aaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  FC_FUNC(g_if_sigma_ooov_cooo_no20_x2_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO20_X2_TYPE1_ERI_V)
    (sa, ia, W37aaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   21] -- 
  // |-- [    0] --| W38(c0,a0,a3,a) += (    1.00000000) V2(a,a4,c0,a3) D1(a0,a4) 
  // |-- [    1] --| W39(a2,a1,a3,a) += (    1.00000000) T2(c0,a2,a1,a0) W38(c0,a0,a3,a) 
  // |-- [    2] --| S2(i,j,k,a) += (    0.50000000) D3(a1,a3,j,a2,i,k) W39(a2,a1,a3,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// ==> makeContractions: [1] Special care is taken for the contraction 1 <2, 1> 
// |------> 2 (alloc)
// allocSigma.first = 2
  // Pref: 1
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W39aaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    orz::DTensor W38ca_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, sa0^sa));
    FC_FUNC(g_if_sigma_ooov_cooo_no21_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO21_X0_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, V2_sym.cptr(), W38ca_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_ooov_cooo_no21_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO21_X1_TYPE1_ERI_V)
      (sa, ia, sa0, ia0, T2b.cptr(), W38ca_sigma_ooov_cooo.cptr(), W39aaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  FC_FUNC(g_if_sigma_ooov_cooo_no21_x2_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO21_X2_TYPE1_ERI_V)
    (sa, ia, W39aaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   22] -- 
  // |-- [    0] --| W43(a4,a) += (    1.00000000) V2(a,a4,c0,a3) W42(c0,a3) 
  // |-- [    1] --| S2(i,j,k,a) += (    1.00000000) D2(j,a4,i,k) W43(a4,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W43a_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no22_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO22_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W42ca_sigma_ooov_cooo.cptr(), W43a_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no22_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO22_X1_TYPE1_ERI_V)
    (sa, ia, W43a_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   23] -- 
  // |-- [    0] --| W45(a3,a) += (    1.00000000) V2(a,a4,c0,a3) W44(c0,a4) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.50000000) D2(j,a3,i,k) W45(a3,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W45a_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no23_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO23_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W44ca_sigma_ooov_cooo.cptr(), W45a_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no23_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO23_X1_TYPE1_ERI_V)
    (sa, ia, W45a_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   24] -- 
  // |-- [    0] --| W47(a4,a3,k,a) += (    1.00000000) V2(a,a4,c0,a3) W46(c0,k) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.50000000) D2(j,a4,i,a3) W47(a4,a3,k,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W47aaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no24_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO24_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W46ca_sigma_ooov_cooo.cptr(), W47aaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no24_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO24_X1_TYPE1_ERI_V)
    (sa, ia, W47aaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   25] -- 
  // |-- [    0] --| W49(a2,a) += (    1.00000000) V2(a,a4,c0,a3) W48(c0,a4,a3,a2) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.50000000) D2(j,a2,i,k) W49(a2,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W49a_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no25_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO25_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W48caaa_sigma_ooov_cooo.cptr(), W49a_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no25_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO25_X1_TYPE1_ERI_V)
    (sa, ia, W49a_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   26] -- 
  // |-- [    0] --| W51(a4,k,a2,a) += (    1.00000000) V2(a,a4,c0,a3) W50(c0,k,a3,a2) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.50000000) D2(j,a4,i,a2) W51(a4,k,a2,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W51aaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no26_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO26_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W50caaa_sigma_ooov_cooo.cptr(), W51aaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no26_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO26_X1_TYPE1_ERI_V)
    (sa, ia, W51aaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   27] -- 
  // |-- [    0] --| W53(a3,k,a2,a) += (    1.00000000) V2(a,a4,c0,a3) W52(c0,a4,k,a2) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.33333333) D2(j,a2,i,a3) W53(a3,k,a2,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W53aaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no27_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO27_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W52caaa_sigma_ooov_cooo.cptr(), W53aaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no27_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO27_X1_TYPE1_ERI_V)
    (sa, ia, W53aaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   28] -- 
  // |-- [    0] --| W55(a3,k,a2,a) += (    1.00000000) V2(a,a4,c0,a3) W54(c0,a4,k,a2) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.16666667) D2(j,a3,i,a2) W55(a3,k,a2,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W55aaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no28_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO28_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W54caaa_sigma_ooov_cooo.cptr(), W55aaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no28_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO28_X1_TYPE1_ERI_V)
    (sa, ia, W55aaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   29] -- 
  // |-- [    0] --| W57(a3,k,a2,a) += (    1.00000000) V2(a,a4,c0,a3) W56(c0,k,a4,a2) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.16666667) D2(j,a2,i,a3) W57(a3,k,a2,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W57aaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no29_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO29_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W56caaa_sigma_ooov_cooo.cptr(), W57aaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no29_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO29_X1_TYPE1_ERI_V)
    (sa, ia, W57aaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   30] -- 
  // |-- [    0] --| W59(a3,k,a2,a) += (    1.00000000) V2(a,a4,c0,a3) W58(c0,k,a4,a2) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.33333333) D2(j,a3,i,a2) W59(a3,k,a2,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W59aaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no30_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO30_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W58caaa_sigma_ooov_cooo.cptr(), W59aaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no30_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO30_X1_TYPE1_ERI_V)
    (sa, ia, W59aaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   31] -- 
  // |-- [    0] --| W61(a3,j,a1,a) += (    1.00000000) V2(a,a4,c0,a3) W60(c0,j,a4,a1) 
  // |-- [    1] --| S2(i,j,k,a) += (    1.00000000) D2(a1,a3,i,k) W61(a3,j,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W61aaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no31_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO31_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W60caaa_sigma_ooov_cooo.cptr(), W61aaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no31_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO31_X1_TYPE1_ERI_V)
    (sa, ia, W61aaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   32] -- 
  // |-- [    0] --| W63(a4,j,a1,a) += (    1.00000000) V2(a,a4,c0,a3) W62(c0,j,a3,a1) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.50000000) D2(a1,a4,i,k) W63(a4,j,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W63aaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no32_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO32_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W62caaa_sigma_ooov_cooo.cptr(), W63aaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no32_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO32_X1_TYPE1_ERI_V)
    (sa, ia, W63aaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   33] -- 
  // |-- [    0] --| W64(c0,a1,i,a) += (    1.00000000) V2(a,a4,c0,a3) D2(a1,a3,i,a4) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.50000000) W64(c0,a1,i,a) W65(c0,j,k,a1) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W64caa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no33_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO33_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W64caa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no33_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO33_X1_TYPE1_ERI_V)
    (sa, ia, W64caa_sigma_ooov_cooo.cptr(), W65caaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   34] -- 
  // |-- [    0] --| W66(c0,a0,j,a) += (    1.00000000) V2(a,a4,c0,a3) D2(a0,a3,j,a4) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.50000000) W66(c0,a0,j,a) W67(c0,i,k,a0) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W66caa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no34_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO34_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W66caa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no34_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO34_X1_TYPE1_ERI_V)
    (sa, ia, W66caa_sigma_ooov_cooo.cptr(), W67caaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   35] -- 
  // |-- [    0] --| W69(a4,i,a0,a) += (    1.00000000) V2(a,a4,c0,a3) W68(c0,i,a3,a0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.50000000) D2(a0,k,j,a4) W69(a4,i,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W69aaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no35_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO35_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W68caaa_sigma_ooov_cooo.cptr(), W69aaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no35_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO35_X1_TYPE1_ERI_V)
    (sa, ia, W69aaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   36] -- 
  // |-- [    0] --| W71(a3,i,a0,a) += (    1.00000000) V2(a,a4,c0,a3) W70(c0,i,a4,a0) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.33333333) D2(a0,a3,j,k) W71(a3,i,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W71aaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no36_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO36_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W70caaa_sigma_ooov_cooo.cptr(), W71aaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no36_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO36_X1_TYPE1_ERI_V)
    (sa, ia, W71aaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   37] -- 
  // |-- [    0] --| W73(a3,i,a0,a) += (    1.00000000) V2(a,a4,c0,a3) W72(c0,i,a4,a0) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.16666667) D2(a0,a3,j,k) W73(a3,i,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W73aaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no37_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO37_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W72caaa_sigma_ooov_cooo.cptr(), W73aaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no37_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO37_X1_TYPE1_ERI_V)
    (sa, ia, W73aaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   38] -- 
  // |-- [    0] --| W75(a3,i,a0,a) += (    1.00000000) V2(a,a4,c0,a3) W74(c0,i,a4,a0) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.16666667) D2(a0,k,j,a3) W75(a3,i,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W75aaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no38_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO38_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W74caaa_sigma_ooov_cooo.cptr(), W75aaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no38_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO38_X1_TYPE1_ERI_V)
    (sa, ia, W75aaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   39] -- 
  // |-- [    0] --| W77(a3,i,a0,a) += (    1.00000000) V2(a,a4,c0,a3) W76(c0,i,a4,a0) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.33333333) D2(a0,k,j,a3) W77(a3,i,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W77aaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no39_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO39_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W76caaa_sigma_ooov_cooo.cptr(), W77aaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no39_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO39_X1_TYPE1_ERI_V)
    (sa, ia, W77aaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   40] -- 
  // |-- [    0] --| W78(c0,j,a1,a) += (    1.00000000) V2(a,a4,c0,a3) D2(j,a4,a1,a3) 
  // |-- [    1] --| S2(i,j,k,a) += (    1.00000000) W78(c0,j,a1,a) W79(c0,i,k,a1) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W78caa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no40_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO40_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W78caa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no40_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO40_X1_TYPE1_ERI_V)
    (sa, ia, W78caa_sigma_ooov_cooo.cptr(), W79caaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   41] -- 
  // |-- [    0] --| W81(a3,i,a1,a) += (    1.00000000) V2(a,a4,c0,a3) W80(c0,i,a4,a1) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.50000000) D2(j,k,a1,a3) W81(a3,i,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W81aaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no41_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO41_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W80caaa_sigma_ooov_cooo.cptr(), W81aaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no41_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO41_X1_TYPE1_ERI_V)
    (sa, ia, W81aaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   42] -- 
  // |-- [    0] --| W83(a4,i,a1,a) += (    1.00000000) V2(a,a4,c0,a3) W82(c0,i,a3,a1) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.50000000) D2(j,a4,a1,k) W83(a4,i,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W83aaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no42_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO42_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W82caaa_sigma_ooov_cooo.cptr(), W83aaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no42_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO42_X1_TYPE1_ERI_V)
    (sa, ia, W83aaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   43] -- 
  // |-- [    0] --| W85(a4,j,a0,a) += (    1.00000000) V2(a,a4,c0,a3) W84(c0,j,a3,a0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.50000000) D2(a0,a4,i,k) W85(a4,j,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W85aaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no43_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO43_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W84caaa_sigma_ooov_cooo.cptr(), W85aaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no43_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO43_X1_TYPE1_ERI_V)
    (sa, ia, W85aaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   44] -- 
  // |-- [    0] --| W87(a3,j,a0,a) += (    1.00000000) V2(a,a4,c0,a3) W86(c0,j,a4,a0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.50000000) D2(a0,a3,i,k) W87(a3,j,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W87aaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no44_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO44_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W86caaa_sigma_ooov_cooo.cptr(), W87aaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no44_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO44_X1_TYPE1_ERI_V)
    (sa, ia, W87aaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   45] -- 
  // |-- [    0] --| W88(c0,a0,i,a) += (    1.00000000) V2(a,a4,c0,a3) D2(a0,a4,i,a3) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.33333333) W88(c0,a0,i,a) W89(c0,j,k,a0) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W88caa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no45_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO45_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W88caa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no45_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO45_X1_TYPE1_ERI_V)
    (sa, ia, W88caa_sigma_ooov_cooo.cptr(), W89caaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   46] -- 
  // |-- [    0] --| W90(c0,a0,i,a) += (    1.00000000) V2(a,a4,c0,a3) D2(a0,a4,i,a3) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.16666667) W90(c0,a0,i,a) W91(c0,j,k,a0) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W90caa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no46_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO46_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W90caa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no46_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO46_X1_TYPE1_ERI_V)
    (sa, ia, W90caa_sigma_ooov_cooo.cptr(), W91caaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   47] -- 
  // |-- [    0] --| W92(c0,a0,i,a) += (    1.00000000) V2(a,a4,c0,a3) D2(a0,a3,i,a4) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.16666667) W92(c0,a0,i,a) W93(c0,j,k,a0) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W92caa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no47_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO47_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W92caa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no47_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO47_X1_TYPE1_ERI_V)
    (sa, ia, W92caa_sigma_ooov_cooo.cptr(), W93caaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   48] -- 
  // |-- [    0] --| W94(c0,a0,i,a) += (    1.00000000) V2(a,a4,c0,a3) D2(a0,a3,i,a4) 
  // |-- [    1] --| S2(i,j,k,a) += (    0.33333333) W94(c0,a0,i,a) W95(c0,j,k,a0) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W94caa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no48_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO48_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W94caa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no48_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO48_X1_TYPE1_ERI_V)
    (sa, ia, W94caa_sigma_ooov_cooo.cptr(), W95caaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   49] -- 
  // |-- [    0] --| W97(a4,a) += (    1.00000000) V2(a,a4,c0,a3) W96(c0,a3) 
  // |-- [    1] --| S2(i,j,k,a) += (   -2.00000000) C2(a4,j,k,i) W97(a4,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W97a_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no49_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO49_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W96ca_sigma_ooov_cooo.cptr(), W97a_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no49_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO49_X1_TYPE1_ERI_V)
    (sa, ia, W97a_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   50] -- 
  // |-- [    0] --| W99(a3,a) += (    1.00000000) V2(a,a4,c0,a3) W98(c0,a4) 
  // |-- [    1] --| S2(i,j,k,a) += (    1.00000000) C2(a3,j,k,i) W99(a3,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W99a_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no50_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO50_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W98ca_sigma_ooov_cooo.cptr(), W99a_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no50_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO50_X1_TYPE1_ERI_V)
    (sa, ia, W99a_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   51] -- 
  // |-- [    0] --| W101(a4,a3,k,a) += (    1.00000000) V2(a,a4,c0,a3) W100(c0,k) 
  // |-- [    1] --| S2(i,j,k,a) += (    1.00000000) C2(a3,i,a4,j) W101(a4,a3,k,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W101aaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no51_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO51_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W100ca_sigma_ooov_cooo.cptr(), W101aaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no51_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO51_X1_TYPE1_ERI_V)
    (sa, ia, W101aaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   52] -- 
  // |-- [    0] --| W103(a2,a) += (    1.00000000) V2(a,a4,c0,a3) W102(c0,a4,a3,a2) 
  // |-- [    1] --| S2(i,j,k,a) += (    1.00000000) C2(a2,j,k,i) W103(a2,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W103a_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no52_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO52_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W102caaa_sigma_ooov_cooo.cptr(), W103a_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no52_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO52_X1_TYPE1_ERI_V)
    (sa, ia, W103a_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   53] -- 
  // |-- [    0] --| W105(a4,k,a2,a) += (    1.00000000) V2(a,a4,c0,a3) W104(c0,k,a3,a2) 
  // |-- [    1] --| S2(i,j,k,a) += (    1.00000000) C2(a2,i,a4,j) W105(a4,k,a2,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W105aaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no53_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO53_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W104caaa_sigma_ooov_cooo.cptr(), W105aaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no53_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO53_X1_TYPE1_ERI_V)
    (sa, ia, W105aaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   54] -- 
  // |-- [    0] --| W107(a3,k,a2,a) += (    1.00000000) V2(a,a4,c0,a3) W106(c0,a4,k,a2) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.66666667) C2(a2,j,a3,i) W107(a3,k,a2,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W107aaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no54_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO54_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W106caaa_sigma_ooov_cooo.cptr(), W107aaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no54_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO54_X1_TYPE1_ERI_V)
    (sa, ia, W107aaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   55] -- 
  // |-- [    0] --| W109(a3,k,a2,a) += (    1.00000000) V2(a,a4,c0,a3) W108(c0,a4,k,a2) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.33333333) C2(a2,i,a3,j) W109(a3,k,a2,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W109aaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no55_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO55_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W108caaa_sigma_ooov_cooo.cptr(), W109aaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no55_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO55_X1_TYPE1_ERI_V)
    (sa, ia, W109aaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   56] -- 
  // |-- [    0] --| W111(a3,k,a2,a) += (    1.00000000) V2(a,a4,c0,a3) W110(c0,k,a4,a2) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.33333333) C2(a2,j,a3,i) W111(a3,k,a2,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W111aaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no56_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO56_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W110caaa_sigma_ooov_cooo.cptr(), W111aaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no56_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO56_X1_TYPE1_ERI_V)
    (sa, ia, W111aaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   57] -- 
  // |-- [    0] --| W113(a3,k,a2,a) += (    1.00000000) V2(a,a4,c0,a3) W112(c0,k,a4,a2) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.66666667) C2(a2,i,a3,j) W113(a3,k,a2,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W113aaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no57_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO57_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W112caaa_sigma_ooov_cooo.cptr(), W113aaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no57_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO57_X1_TYPE1_ERI_V)
    (sa, ia, W113aaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   58] -- 
  // |-- [    0] --| W115(a3,j,a1,a) += (    1.00000000) V2(a,a4,c0,a3) W114(c0,j,a4,a1) 
  // |-- [    1] --| S2(i,j,k,a) += (   -2.00000000) C2(a1,a3,i,k) W115(a3,j,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W115aaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no58_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO58_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W114caaa_sigma_ooov_cooo.cptr(), W115aaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no58_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO58_X1_TYPE1_ERI_V)
    (sa, ia, W115aaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   59] -- 
  // |-- [    0] --| W117(a4,j,a1,a) += (    1.00000000) V2(a,a4,c0,a3) W116(c0,j,a3,a1) 
  // |-- [    1] --| S2(i,j,k,a) += (    1.00000000) C2(a1,a4,i,k) W117(a4,j,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W117aaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no59_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO59_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W116caaa_sigma_ooov_cooo.cptr(), W117aaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no59_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO59_X1_TYPE1_ERI_V)
    (sa, ia, W117aaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   60] -- 
  // |-- [    0] --| W118(c0,a1,i,a) += (    1.00000000) V2(a,a4,c0,a3) C2(a1,a3,i,a4) 
  // |-- [    1] --| S2(i,j,k,a) += (    1.00000000) W118(c0,a1,i,a) W119(c0,j,k,a1) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W118caa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no60_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO60_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W118caa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no60_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO60_X1_TYPE1_ERI_V)
    (sa, ia, W118caa_sigma_ooov_cooo.cptr(), W119caaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   61] -- 
  // |-- [    0] --| W120(c0,a0,j,a) += (    1.00000000) V2(a,a4,c0,a3) C2(a0,a3,j,a4) 
  // |-- [    1] --| S2(i,j,k,a) += (    1.00000000) W120(c0,a0,j,a) W121(c0,i,k,a0) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W120caa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no61_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO61_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W120caa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no61_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO61_X1_TYPE1_ERI_V)
    (sa, ia, W120caa_sigma_ooov_cooo.cptr(), W121caaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   62] -- 
  // |-- [    0] --| W123(a4,i,a0,a) += (    1.00000000) V2(a,a4,c0,a3) W122(c0,i,a3,a0) 
  // |-- [    1] --| S2(i,j,k,a) += (    1.00000000) C2(a0,k,j,a4) W123(a4,i,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W123aaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no62_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO62_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W122caaa_sigma_ooov_cooo.cptr(), W123aaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no62_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO62_X1_TYPE1_ERI_V)
    (sa, ia, W123aaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   63] -- 
  // |-- [    0] --| W125(a3,i,a0,a) += (    1.00000000) V2(a,a4,c0,a3) W124(c0,i,a4,a0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.66666667) C2(a0,a3,j,k) W125(a3,i,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W125aaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no63_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO63_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W124caaa_sigma_ooov_cooo.cptr(), W125aaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no63_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO63_X1_TYPE1_ERI_V)
    (sa, ia, W125aaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   64] -- 
  // |-- [    0] --| W127(a3,i,a0,a) += (    1.00000000) V2(a,a4,c0,a3) W126(c0,i,a4,a0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.33333333) C2(a0,a3,j,k) W127(a3,i,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W127aaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no64_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO64_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W126caaa_sigma_ooov_cooo.cptr(), W127aaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no64_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO64_X1_TYPE1_ERI_V)
    (sa, ia, W127aaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   65] -- 
  // |-- [    0] --| W129(a3,i,a0,a) += (    1.00000000) V2(a,a4,c0,a3) W128(c0,i,a4,a0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.33333333) C2(a0,k,j,a3) W129(a3,i,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W129aaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no65_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO65_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W128caaa_sigma_ooov_cooo.cptr(), W129aaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no65_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO65_X1_TYPE1_ERI_V)
    (sa, ia, W129aaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   66] -- 
  // |-- [    0] --| W131(a3,i,a0,a) += (    1.00000000) V2(a,a4,c0,a3) W130(c0,i,a4,a0) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.66666667) C2(a0,k,j,a3) W131(a3,i,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W131aaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no66_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO66_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W130caaa_sigma_ooov_cooo.cptr(), W131aaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no66_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO66_X1_TYPE1_ERI_V)
    (sa, ia, W131aaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   67] -- 
  // |-- [    0] --| W132(c0,j,a1,a) += (    1.00000000) V2(a,a4,c0,a3) C2(a1,a3,j,a4) 
  // |-- [    1] --| S2(i,j,k,a) += (   -2.00000000) W132(c0,j,a1,a) W133(c0,i,k,a1) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W132caa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no67_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO67_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W132caa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no67_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO67_X1_TYPE1_ERI_V)
    (sa, ia, W132caa_sigma_ooov_cooo.cptr(), W133caaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   68] -- 
  // |-- [    0] --| W135(a3,i,a1,a) += (    1.00000000) V2(a,a4,c0,a3) W134(c0,i,a4,a1) 
  // |-- [    1] --| S2(i,j,k,a) += (    1.00000000) C2(a1,a3,j,k) W135(a3,i,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W135aaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no68_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO68_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W134caaa_sigma_ooov_cooo.cptr(), W135aaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no68_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO68_X1_TYPE1_ERI_V)
    (sa, ia, W135aaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   69] -- 
  // |-- [    0] --| W137(a4,i,a1,a) += (    1.00000000) V2(a,a4,c0,a3) W136(c0,i,a3,a1) 
  // |-- [    1] --| S2(i,j,k,a) += (    1.00000000) C2(a1,k,j,a4) W137(a4,i,a1,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W137aaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no69_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO69_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W136caaa_sigma_ooov_cooo.cptr(), W137aaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no69_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO69_X1_TYPE1_ERI_V)
    (sa, ia, W137aaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   70] -- 
  // |-- [    0] --| W139(a4,j,a0,a) += (    1.00000000) V2(a,a4,c0,a3) W138(c0,j,a3,a0) 
  // |-- [    1] --| S2(i,j,k,a) += (    1.00000000) C2(a0,a4,i,k) W139(a4,j,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W139aaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no70_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO70_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W138caaa_sigma_ooov_cooo.cptr(), W139aaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no70_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO70_X1_TYPE1_ERI_V)
    (sa, ia, W139aaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   71] -- 
  // |-- [    0] --| W141(a3,j,a0,a) += (    1.00000000) V2(a,a4,c0,a3) W140(c0,j,a4,a0) 
  // |-- [    1] --| S2(i,j,k,a) += (    1.00000000) C2(a0,a3,i,k) W141(a3,j,a0,a) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W141aaa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xaaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no71_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO71_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W140caaa_sigma_ooov_cooo.cptr(), W141aaa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no71_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO71_X1_TYPE1_ERI_V)
    (sa, ia, W141aaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   72] -- 
  // |-- [    0] --| W142(c0,a0,i,a) += (    1.00000000) V2(a,a4,c0,a3) C2(a0,a4,i,a3) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.66666667) W142(c0,a0,i,a) W143(c0,j,k,a0) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W142caa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no72_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO72_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W142caa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no72_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO72_X1_TYPE1_ERI_V)
    (sa, ia, W142caa_sigma_ooov_cooo.cptr(), W143caaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   73] -- 
  // |-- [    0] --| W144(c0,a0,i,a) += (    1.00000000) V2(a,a4,c0,a3) C2(a0,a4,i,a3) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.33333333) W144(c0,a0,i,a) W145(c0,j,k,a0) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W144caa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no73_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO73_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W144caa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no73_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO73_X1_TYPE1_ERI_V)
    (sa, ia, W144caa_sigma_ooov_cooo.cptr(), W145caaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   74] -- 
  // |-- [    0] --| W146(c0,a0,i,a) += (    1.00000000) V2(a,a4,c0,a3) C2(a0,a3,i,a4) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.33333333) W146(c0,a0,i,a) W147(c0,j,k,a0) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W146caa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no74_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO74_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W146caa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no74_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO74_X1_TYPE1_ERI_V)
    (sa, ia, W146caa_sigma_ooov_cooo.cptr(), W147caaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   75] -- 
  // |-- [    0] --| W148(c0,a0,i,a) += (    1.00000000) V2(a,a4,c0,a3) C2(a0,a3,i,a4) 
  // |-- [    1] --| S2(i,j,k,a) += (   -0.66666667) W148(c0,a0,i,a) W149(c0,j,k,a0) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  orz::DTensor W148caa_sigma_ooov_cooo(orz::mr::sizeof_sympack_Xcaa(symblockinfo, sa));
  FC_FUNC(g_if_sigma_ooov_cooo_no75_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO75_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W148caa_sigma_ooov_cooo.cptr(), nir, nsym, psym, &flops);
  FC_FUNC(g_if_sigma_ooov_cooo_no75_x1_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO75_X1_TYPE1_ERI_V)
    (sa, ia, W148caa_sigma_ooov_cooo.cptr(), W149caaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   76] -- 
  // |-- [    0] --| S2(i,j,k,a) += (   -1.00000000) V2(a,a3,c0,a2) W150(c0,j,i,a3,a2,k) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  FC_FUNC(g_if_sigma_ooov_cooo_no76_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO76_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W150caaaaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ia, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [   77] -- 
  // |-- [    0] --| S2(i,j,k,a) += (   -1.00000000) V2(a,a3,c0,a2) W151(c0,j,i,a3,a2,k) 
  int sa(s_eri);
  int ia(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ia]);
  FC_FUNC(g_if_sigma_ooov_cooo_no77_x0_type1_eri_v,G_IF_SIGMA_OOOV_COOO_NO77_X0_TYPE1_ERI_V)
    (sa, ia, V2_sym.cptr(), W151caaaaa_sigma_ooov_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  orz::mr::Femto::my_timer.push_back(boost::make_tuple("sigma_ooov_cooo", time_sigma.elapsed_cputime(), time_sigma.elapsed_wallclocktime()));
  orz::mr::Femto::file_timing << "* " << boost::format("%20s : %10.7f %10.7f ") % orz::mr::Femto::my_timer.back().get<0>() % orz::mr::Femto::my_timer.back().get<1>() % orz::mr::Femto::my_timer.back().get<2>() << endl;
  flush(orz::mr::Femto::file_timing);

  return retval; 
} 
