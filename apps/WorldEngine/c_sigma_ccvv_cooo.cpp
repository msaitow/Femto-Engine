                                                                                
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
#include <sci/icmr/Femto/elems/c_sigma_ccvv_cooo.h>                                  
                                                                                
using std::cout;                                                                
using std::endl;                                                                
                                                                                
#define FLOPCOUNT                                                               
                                                                                
using orz::mr::Femto::Fc0;                                                      
using orz::mr::Femto::h1_int;                                                   
using orz::mr::Femto::h6_int;                                                   
                                                                                
//  `7MM"""YMM                         mm               
//    MM    `7                         MM                  
//    MM   d  .gP"Ya `7MMpMMMb.pMMMb.mmMMmm ,pW"Wq.      
//    MM""MM ,M'   Yb  MM    MM    MM  MM  6W'   `Wb   
//    MM   Y 8M""""""  MM    MM    MM  MM  8M     M8 
//    MM     YM.    ,  MM    MM    MM  MM  YA.   ,A9       
//  .JMML.    `Mbmmd'.JMML  JMML  JMML.`Mbmo`Ybmd9'        

//                                   Generated date : Sun Apr 20 10:26:25 2014

                                                                                
// ***************************************************************************  
// orz::mr::mrci                                                                
// ***************************************************************************  
									    /*!        
   @brief CT input                                                              
									     */        
                                                                                
orz::mr::BareAmpPack orz::mr::Femto::sigma_ccvv_cooo(const orz::mr::Input &ctinp,                                    
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
  // -- Title : sigma_ccvv_cooo

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W4(x,a2) += (    1.00000000) T2(a1,x,a0,a2) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (    4.00000000) V2(a2,b,w,a) W4(x,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia2);
  orz::DTensor W4c_sigma_ccvv_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, sa2));
  FC_FUNC(g_if_sigma_ccvv_cooo_no0_x0_type1_eri_o,G_IF_SIGMA_CCVV_COOO_NO0_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), W4c_sigma_ccvv_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_cooo_no0_x1_type1_eri_o,G_IF_SIGMA_CCVV_COOO_NO0_X1_TYPE1_ERI_O)
      (sa2, ia2, sb, ib, V2_sym.cptr(), W4c_sigma_ccvv_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W5(w,a2) += (    1.00000000) T2(a1,w,a0,a2) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (   -2.00000000) V2(a2,b,x,a) W5(w,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia2);
  orz::DTensor W5c_sigma_ccvv_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, sa2));
  FC_FUNC(g_if_sigma_ccvv_cooo_no1_x0_type1_eri_o,G_IF_SIGMA_CCVV_COOO_NO1_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), W5c_sigma_ccvv_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_cooo_no1_x1_type1_eri_o,G_IF_SIGMA_CCVV_COOO_NO1_X1_TYPE1_ERI_O)
      (sa2, ia2, sb, ib, V2_sym.cptr(), W5c_sigma_ccvv_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W6(w,a2) += (    1.00000000) T2(a1,w,a0,a2) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (    4.00000000) V2(a2,a,x,b) W6(w,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia2);
  orz::DTensor W6c_sigma_ccvv_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, sa2));
  FC_FUNC(g_if_sigma_ccvv_cooo_no2_x0_type1_eri_o,G_IF_SIGMA_CCVV_COOO_NO2_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), W6c_sigma_ccvv_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_cooo_no2_x1_type1_eri_o,G_IF_SIGMA_CCVV_COOO_NO2_X1_TYPE1_ERI_O)
      (sa2, ia2, sb, ib, V2_sym.cptr(), W6c_sigma_ccvv_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W7(x,a2) += (    1.00000000) T2(a1,x,a0,a2) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (   -2.00000000) V2(a2,a,w,b) W7(x,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia2);
  orz::DTensor W7c_sigma_ccvv_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, sa2));
  FC_FUNC(g_if_sigma_ccvv_cooo_no3_x0_type1_eri_o,G_IF_SIGMA_CCVV_COOO_NO3_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), W7c_sigma_ccvv_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_cooo_no3_x1_type1_eri_o,G_IF_SIGMA_CCVV_COOO_NO3_X1_TYPE1_ERI_O)
      (sa2, ia2, sb, ib, V2_sym.cptr(), W7c_sigma_ccvv_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    4] -- 
  // |-- [    0] --| W8(w,a2) += (    1.00000000) T2(w,a1,a0,a2) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (    1.00000000) V2(a2,b,x,a) W8(w,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia2);
  orz::DTensor W8c_sigma_ccvv_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, sa2));
  FC_FUNC(g_if_sigma_ccvv_cooo_no4_x0_type1_eri_o,G_IF_SIGMA_CCVV_COOO_NO4_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), W8c_sigma_ccvv_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_cooo_no4_x1_type1_eri_o,G_IF_SIGMA_CCVV_COOO_NO4_X1_TYPE1_ERI_O)
      (sa2, ia2, sb, ib, V2_sym.cptr(), W8c_sigma_ccvv_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    5] -- 
  // |-- [    0] --| W9(x,a2) += (    1.00000000) T2(x,a1,a0,a2) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (   -2.00000000) V2(a2,b,w,a) W9(x,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia2);
  orz::DTensor W9c_sigma_ccvv_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, sa2));
  FC_FUNC(g_if_sigma_ccvv_cooo_no5_x0_type1_eri_o,G_IF_SIGMA_CCVV_COOO_NO5_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), W9c_sigma_ccvv_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_cooo_no5_x1_type1_eri_o,G_IF_SIGMA_CCVV_COOO_NO5_X1_TYPE1_ERI_O)
      (sa2, ia2, sb, ib, V2_sym.cptr(), W9c_sigma_ccvv_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    6] -- 
  // |-- [    0] --| W10(w,a2) += (    1.00000000) T2(w,a1,a0,a2) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (   -2.00000000) V2(a2,a,x,b) W10(w,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia2);
  orz::DTensor W10c_sigma_ccvv_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, sa2));
  FC_FUNC(g_if_sigma_ccvv_cooo_no6_x0_type1_eri_o,G_IF_SIGMA_CCVV_COOO_NO6_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), W10c_sigma_ccvv_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_cooo_no6_x1_type1_eri_o,G_IF_SIGMA_CCVV_COOO_NO6_X1_TYPE1_ERI_O)
      (sa2, ia2, sb, ib, V2_sym.cptr(), W10c_sigma_ccvv_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  // --> @[b, "virtual"] [notNeeded]
    retval.acc_amp2(ib, S2b);
  } // End ib
  } // End sb
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    7] -- 
  // |-- [    0] --| W11(x,a2) += (    1.00000000) T2(x,a1,a0,a2) D1(a1,a0) 
  // |-- [    1] --| S2(w,x,a,b) += (    1.00000000) V2(a2,a,w,b) W11(x,a2) 
  int sa2(s_eri);
  int ia2(i_eri);
  double flops(0); // Flop count  
// |------> 1 (alloc)
// allocSigma.first = 1
  // Pref: 0
  T2b = T2.get_amp2(ia2);
  orz::DTensor W11c_sigma_ccvv_cooo(orz::mr::sizeof_sympack_Xc(symblockinfo, sa2));
  FC_FUNC(g_if_sigma_ccvv_cooo_no7_x0_type1_eri_o,G_IF_SIGMA_CCVV_COOO_NO7_X0_TYPE1_ERI_O)
    (sa2, ia2, T2b.cptr(), W11c_sigma_ccvv_cooo.cptr(), nir, nsym, psym, &flops);
  for(int sb = 0;sb < nir;++sb){ 
  for(int ib = symblockinfo.psym()(sb,I_V,I_BEGIN);ib <= symblockinfo.psym()(sb,I_V,I_END);++ib){ 
    S2b = orz::DTensor(retval.namps_iamp()[ib]);
    FC_FUNC(g_if_sigma_ccvv_cooo_no7_x1_type1_eri_o,G_IF_SIGMA_CCVV_COOO_NO7_X1_TYPE1_ERI_O)
      (sa2, ia2, sb, ib, V2_sym.cptr(), W11c_sigma_ccvv_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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

//-@type(2).declaration(begin)
  // --  Title : sigma_ccvv_cooo
  //  >> Intermediates for the external contractions are defined here << 
  orz::DTensor W0ca_sigma_ccvv_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W1ca_sigma_ccvv_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W2ca_sigma_ccvv_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
  orz::DTensor W3ca_sigma_ccvv_cooo(orz::mr::sizeof_sympack_Xca(symblockinfo, 0));
//-@type(2).declaration(end)

  //*-- Entering to take the type 0 contractions --*//
//-@type(0).contraction(begin)

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| W0(w,a3) += (    1.00000000) T2(w,a2,a1,a0) D2(a3,a1,a2,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccvv_cooo_no0_x0_type0_eri_v,G_IF_SIGMA_CCVV_COOO_NO0_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W0ca_sigma_ccvv_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| W1(x,a3) += (    1.00000000) T2(x,a2,a1,a0) D2(a3,a1,a2,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccvv_cooo_no1_x0_type0_eri_v,G_IF_SIGMA_CCVV_COOO_NO1_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W1ca_sigma_ccvv_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| W2(w,a3) += (    1.00000000) T2(w,a2,a1,a0) D2(a3,a1,a2,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccvv_cooo_no2_x0_type0_eri_v,G_IF_SIGMA_CCVV_COOO_NO2_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W2ca_sigma_ccvv_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| W3(x,a3) += (    1.00000000) T2(x,a2,a1,a0) D2(a3,a1,a2,a0) 
  double flops(0); // Flop count  
  // Pref: 2
  for(int sa0 = 0;sa0 < nir;++sa0){ 
  for(int ia0 = symblockinfo.psym()(sa0,I_O,I_BEGIN);ia0 <= symblockinfo.psym()(sa0,I_O,I_END);++ia0){ 
    T2b = T2.get_amp2(ia0);
    FC_FUNC(g_if_sigma_ccvv_cooo_no3_x0_type0_eri_v,G_IF_SIGMA_CCVV_COOO_NO3_X0_TYPE0_ERI_V)
      (sa0, ia0, T2b.cptr(), W3ca_sigma_ccvv_cooo.cptr(), nir, nsym, psym, &flops);
  // --> @[a0, "active"] [notNeeded]
  } // End ia0
  } // End sa0
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
  // -- Title : sigma_ccvv_cooo

  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    0] -- 
  // |-- [    0] --| S2(w,x,a,b) += (   -2.00000000) V2(b,x,a3,a) W0(w,a3) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  FC_FUNC(g_if_sigma_ccvv_cooo_no0_x0_type1_eri_v,G_IF_SIGMA_CCVV_COOO_NO0_X0_TYPE1_ERI_V)
    (sb, ib, V2_sym.cptr(), W0ca_sigma_ccvv_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    1] -- 
  // |-- [    0] --| S2(w,x,a,b) += (    1.00000000) V2(b,w,a3,a) W1(x,a3) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  FC_FUNC(g_if_sigma_ccvv_cooo_no1_x0_type1_eri_v,G_IF_SIGMA_CCVV_COOO_NO1_X0_TYPE1_ERI_V)
    (sb, ib, V2_sym.cptr(), W1ca_sigma_ccvv_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    2] -- 
  // |-- [    0] --| S2(w,x,a,b) += (    1.00000000) V2(b,a3,x,a) W2(w,a3) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  FC_FUNC(g_if_sigma_ccvv_cooo_no2_x0_type1_eri_v,G_IF_SIGMA_CCVV_COOO_NO2_X0_TYPE1_ERI_V)
    (sb, ib, V2_sym.cptr(), W2ca_sigma_ccvv_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ib, S2b);
  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  // -- No. [    3] -- 
  // |-- [    0] --| S2(w,x,a,b) += (   -2.00000000) V2(b,a3,w,a) W3(x,a3) 
  int sb(s_eri);
  int ib(i_eri);
  double flops(0); // Flop count  
// |------> 0 (alloc)
// allocSigma.first = 0
  // Pref: 2
  S2b = orz::DTensor(retval.namps_iamp()[ib]);
  FC_FUNC(g_if_sigma_ccvv_cooo_no3_x0_type1_eri_v,G_IF_SIGMA_CCVV_COOO_NO3_X0_TYPE1_ERI_V)
    (sb, ib, V2_sym.cptr(), W3ca_sigma_ccvv_cooo.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
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
  orz::mr::Femto::my_timer.push_back(boost::make_tuple("sigma_ccvv_cooo", time_sigma.elapsed_cputime(), time_sigma.elapsed_wallclocktime()));
  orz::mr::Femto::file_timing << "* " << boost::format("%20s : %10.7f %10.7f ") % orz::mr::Femto::my_timer.back().get<0>() % orz::mr::Femto::my_timer.back().get<1>() % orz::mr::Femto::my_timer.back().get<2>() << endl;
  flush(orz::mr::Femto::file_timing);

  return retval; 
} 
