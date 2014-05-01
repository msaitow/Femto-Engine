                                                                                
#include <orz/orz.h>                                                            
#include <orz/openmp.h>                                                         
#include <orz/cblas.h>                                                          
#include <orz/clapack.h>                                                        
#include <tensor/tensor.h>                                                      
#include <sci/hint/para_disttools.h>                                            
#include <sci/ctnew2/ct.h>                                                      
#include <sci/ctnew2/ct_f.h>                                                    
#include <sci/ctnew2/ctclass_input.h>                                           
#include <sci/ctnew2/ctclass_symblock.h>                                        
#include <sci/ctnew2/ctclass_hintmo.h>                                          
#include <sci/ctnew2/ctclass_rdmpack.h>                                         
#include <sci/ctnew2/ctclass_bareamppack.h>                                     
#include <sci/ctnew2/ctclass_orthamppack.h>                                     
#include <sci/ctnew2/diaghessian.h>                                             
#include <sci/ctnew2/symamp2.h>                                                 
#include <sci/ctnew2/mrci.h>                                                    
#include <sci/ctnew2/c_sigma_oovv_ooov.h>                                            
                                                                                
using std::cout;                                                                
using std::endl;                                                                
                                                                                
#define FLOPCOUNT                                                               
                                                                                
//  __/\\\\\\\\\\\\\\\____________________________________________________________________                                   
//   _\/\\\///////////_____________________________________________________________________                                             
//    _\/\\\_______________________________________________________/\\\_____________________                                         
//     _\/\\\\\\\\\\\__________/\\\\\\\\______/\\\\\__/\\\\\_____/\\\\\\\\\\\______/\\\\\____ 
//      _\/\\\///////_________/\\\/////\\\___/\\\///\\\\\///\\\__\////\\\////_____/\\\///\\\__               
//       _\/\\\_______________/\\\\\\\\\\\___\/\\\_\//\\\__\/\\\_____\/\\\________/\\\__\//\\\_       
//        _\/\\\______________\//\\///////____\/\\\__\/\\\__\/\\\_____\/\\\_/\\___\//\\\__/\\\__            
//         _\/\\\_______________\//\\\\\\\\\\__\/\\\__\/\\\__\/\\\_____\//\\\\\_____\///\\\\\/___    
//          _\///_________________\//////////___\///___\///___\///_______\/////________\/////_____                                   

                                                                                
// ***************************************************************************  
// orz::ct::mrci                                                                
// ***************************************************************************  
									    /*!        
   @brief CT input                                                              
									     */        
                                                                                
orz::ct::BareAmpPack orz::ct::sigma_oovv_ooov(const orz::ct::Input &ctinp,                                    
                                  const orz::ct::SymBlockInfo &symblockinfo,                                
                                  const orz::ct::HintMO &hintmo,                                            
                                  const orz::ct::RdmPack &rdmPack_sym,                                      
                                  const orz::DTensor &rdm4,                                                 
                                  const orz::ct::BareAmpPack &T2,                             
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
                                                                                                                
  orz::DTensor moint1 = hintmo.int1(); // Setting up one-body integrals                                         
  const orz::DTensor moint1_sym = orz::ct::sympack_int1(symblockinfo, moint1); // moint1=(IR-COV index)         
  orz::DTensor V2(nmo,nmo,nmo);                                                                    
  double * const V2_ptr = V2.cptr();                                                  
                                                                                                                
  std::ostringstream stm;                                                                                       
  stm << num_sigma;                                                                                             
  std::string name_of_sigma = "S2" + stm.str() + "]"; // Name of the Sigma vector  
  orz::ct::BareAmpPack retval                                                                                   
    = orz::ct::BareAmpPack(ctinp, symblockinfo, name_of_sigma); // Sigma(a, a', e, e') tensor                   
                                                                                                                
  orz::DTensor S2b; // Container of S2_aae,[b] tensor                                   
                                                                                                                
  orz::DTensor T2b; // Container of T2_aae,[b] tensor                                             
  orz::DTensor rdm4_sym;                                                                                        


  {
  // No.0
  //* X(i,k,o2,a)  <--  (    1.00000000)  D3(i,o3,k,o2,o4,o1) T2(o1,o3,o4,a) 
  //* S2(i,k,a,c)  <--  (    1.00000000) X(i,k,o2,a) h(o2,c) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sa, X);
    FC_FUNC(g_if_sigma_oovv_ooov_no0_x0, G_IF_SIGMA_OOVV_OOOV_NO0_X0)
      (sa, ia, T2b.cptr(), Xaaa.cptr(), nir, nsym, psym);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      FC_FUNC(g_if_sigma_oovv_ooov_no1_x0, G_IF_SIGMA_OOVV_OOOV_NO1_X0)
        (sa, ia, sc, ic, Xaaa.cptr(), moint1_sym.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.1
  //* X(i,o2,k,c)  <--  (    1.00000000)  D3(i,o2,k,o3,o1,o4) T2(o4,o3,o1,c) 
  //* S2(i,k,a,c)  <--  (    1.00000000) X(i,o2,k,c) h(o2,a) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sc, X);
    FC_FUNC(g_if_sigma_oovv_ooov_no0_x1, G_IF_SIGMA_OOVV_OOOV_NO0_X1)
      (sc, ic, T2b.cptr(), Xaaa.cptr(), nir, nsym, psym);
    FC_FUNC(g_if_sigma_oovv_ooov_no1_x1, G_IF_SIGMA_OOVV_OOOV_NO1_X1)
      (sc, ic, Xaaa.cptr(), moint1_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.2
  //* X(i,k,o1,a)  <--  (    1.00000000)  D2(i,o3,k,o2) T2(o2,o3,o1,a) 
  //* S2(i,k,a,c)  <--  (    1.00000000) X(i,k,o1,a) h(o1,c) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sa, X);
    FC_FUNC(g_if_sigma_oovv_ooov_no0_x2, G_IF_SIGMA_OOVV_OOOV_NO0_X2)
      (sa, ia, T2b.cptr(), Xaaa.cptr(), nir, nsym, psym);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      FC_FUNC(g_if_sigma_oovv_ooov_no1_x2, G_IF_SIGMA_OOVV_OOOV_NO1_X2)
        (sa, ia, sc, ic, Xaaa.cptr(), moint1_sym.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.3
  //* X(i,k,o3,c)  <--  (    1.00000000)  D2(i,o1,k,o2) T2(o1,o2,o3,c) 
  //* S2(i,k,a,c)  <--  (    1.00000000) X(i,k,o3,c) h(o3,a) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sc, X);
    FC_FUNC(g_if_sigma_oovv_ooov_no0_x3, G_IF_SIGMA_OOVV_OOOV_NO0_X3)
      (sc, ic, T2b.cptr(), Xaaa.cptr(), nir, nsym, psym);
    FC_FUNC(g_if_sigma_oovv_ooov_no1_x3, G_IF_SIGMA_OOVV_OOOV_NO1_X3)
      (sc, ic, Xaaa.cptr(), moint1_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.4
  //* X(i,k,o2,a)  <--  (    1.00000000)  D3(i,o3,k,o2,o4,o1) T2(o1,o3,o4,a) 
  //* S2(i,k,a,c)  <--  (    2.00000000) X(i,k,o2,a) Y0(o2,c) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nocc, nvir);
  orz::DTensor Y0 = orz::ct::sympack_Xav(symblockinfo, 0, Y);
  for(int sc1 = 0;sc1 < nir;++sc1){ 
  for(int ic1 = symblockinfo.psym()(sc1,I_C,I_BEGIN);ic1 <= symblockinfo.psym()(sc1,I_C,I_END);++ic1){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(ic1);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, ic1, sc1, V2); // V2=(IR-COV index) 
    FC_FUNC(g_if_sigma_oovv_ooov_y4, G_IF_SIGMA_OOVV_OOOV_Y4)
      (sc1, ic1, V2_sym.cptr(), Y0.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sa, X);
    FC_FUNC(g_if_sigma_oovv_ooov_no0_x4, G_IF_SIGMA_OOVV_OOOV_NO0_X4)
      (sa, ia, T2b.cptr(), Xaaa.cptr(), nir, nsym, psym);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      FC_FUNC(g_if_sigma_oovv_ooov_no1_x4, G_IF_SIGMA_OOVV_OOOV_NO1_X4)
        (sa, ia, sc, ic, Xaaa.cptr(), Y0.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.5
  //* X(i,k,o2,a)  <--  (    1.00000000)  D3(i,o3,k,o2,o4,o1) T2(o1,o3,o4,a) 
  //* S2(i,k,a,c)  <--  (   -1.00000000) X(i,k,o2,a) Y1(o2,c) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nocc, nvir);
  orz::DTensor Y1 = orz::ct::sympack_Xav(symblockinfo, 0, Y);
  for(int sc1 = 0;sc1 < nir;++sc1){ 
  for(int ic1 = symblockinfo.psym()(sc1,I_C,I_BEGIN);ic1 <= symblockinfo.psym()(sc1,I_C,I_END);++ic1){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(ic1);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, ic1, sc1, V2); // V2=(IR-COV index) 
    FC_FUNC(g_if_sigma_oovv_ooov_y5, G_IF_SIGMA_OOVV_OOOV_Y5)
      (sc1, ic1, V2_sym.cptr(), Y1.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sa, X);
    FC_FUNC(g_if_sigma_oovv_ooov_no0_x5, G_IF_SIGMA_OOVV_OOOV_NO0_X5)
      (sa, ia, T2b.cptr(), Xaaa.cptr(), nir, nsym, psym);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      FC_FUNC(g_if_sigma_oovv_ooov_no1_x5, G_IF_SIGMA_OOVV_OOOV_NO1_X5)
        (sa, ia, sc, ic, Xaaa.cptr(), Y1.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.6
  //* X(i,o3,k,c)  <--  (    1.00000000)  D3(i,o3,k,o2,o4,o1) T2(o1,o2,o4,c) 
  //* S2(i,k,a,c)  <--  (    2.00000000) X(i,o3,k,c) Y2(o3,a) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nocc, nvir);
  orz::DTensor Y2 = orz::ct::sympack_Xav(symblockinfo, 0, Y);
  for(int sc1 = 0;sc1 < nir;++sc1){ 
  for(int ic1 = symblockinfo.psym()(sc1,I_C,I_BEGIN);ic1 <= symblockinfo.psym()(sc1,I_C,I_END);++ic1){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(ic1);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, ic1, sc1, V2); // V2=(IR-COV index) 
    FC_FUNC(g_if_sigma_oovv_ooov_y6, G_IF_SIGMA_OOVV_OOOV_Y6)
      (sc1, ic1, V2_sym.cptr(), Y2.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sc, X);
    FC_FUNC(g_if_sigma_oovv_ooov_no0_x6, G_IF_SIGMA_OOVV_OOOV_NO0_X6)
      (sc, ic, T2b.cptr(), Xaaa.cptr(), nir, nsym, psym);
    FC_FUNC(g_if_sigma_oovv_ooov_no1_x6, G_IF_SIGMA_OOVV_OOOV_NO1_X6)
      (sc, ic, Xaaa.cptr(), Y2.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.7
  //* X(i,o3,k,c)  <--  (    1.00000000)  D3(i,o3,k,o2,o4,o1) T2(o1,o2,o4,c) 
  //* S2(i,k,a,c)  <--  (   -1.00000000) X(i,o3,k,c) Y3(o3,a) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nocc, nvir);
  orz::DTensor Y3 = orz::ct::sympack_Xav(symblockinfo, 0, Y);
  for(int sc1 = 0;sc1 < nir;++sc1){ 
  for(int ic1 = symblockinfo.psym()(sc1,I_C,I_BEGIN);ic1 <= symblockinfo.psym()(sc1,I_C,I_END);++ic1){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(ic1);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, ic1, sc1, V2); // V2=(IR-COV index) 
    FC_FUNC(g_if_sigma_oovv_ooov_y7, G_IF_SIGMA_OOVV_OOOV_Y7)
      (sc1, ic1, V2_sym.cptr(), Y3.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sc, X);
    FC_FUNC(g_if_sigma_oovv_ooov_no0_x7, G_IF_SIGMA_OOVV_OOOV_NO0_X7)
      (sc, ic, T2b.cptr(), Xaaa.cptr(), nir, nsym, psym);
    FC_FUNC(g_if_sigma_oovv_ooov_no1_x7, G_IF_SIGMA_OOVV_OOOV_NO1_X7)
      (sc, ic, Xaaa.cptr(), Y3.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.8
  //* X(i,k,o1,a)  <--  (    1.00000000)  D2(i,o3,k,o2) T2(o2,o3,o1,a) 
  //* S2(i,k,a,c)  <--  (    2.00000000) X(i,k,o1,a) Y4(o1,c) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nocc, nvir);
  orz::DTensor Y4 = orz::ct::sympack_Xav(symblockinfo, 0, Y);
  for(int sc1 = 0;sc1 < nir;++sc1){ 
  for(int ic1 = symblockinfo.psym()(sc1,I_C,I_BEGIN);ic1 <= symblockinfo.psym()(sc1,I_C,I_END);++ic1){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(ic1);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, ic1, sc1, V2); // V2=(IR-COV index) 
    FC_FUNC(g_if_sigma_oovv_ooov_y8, G_IF_SIGMA_OOVV_OOOV_Y8)
      (sc1, ic1, V2_sym.cptr(), Y4.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sa, X);
    FC_FUNC(g_if_sigma_oovv_ooov_no0_x8, G_IF_SIGMA_OOVV_OOOV_NO0_X8)
      (sa, ia, T2b.cptr(), Xaaa.cptr(), nir, nsym, psym);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      FC_FUNC(g_if_sigma_oovv_ooov_no1_x8, G_IF_SIGMA_OOVV_OOOV_NO1_X8)
        (sa, ia, sc, ic, Xaaa.cptr(), Y4.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.9
  //* X(i,k,o1,a)  <--  (    1.00000000)  D2(i,o3,k,o2) T2(o2,o3,o1,a) 
  //* S2(i,k,a,c)  <--  (   -1.00000000) X(i,k,o1,a) Y5(o1,c) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nocc, nvir);
  orz::DTensor Y5 = orz::ct::sympack_Xav(symblockinfo, 0, Y);
  for(int sc1 = 0;sc1 < nir;++sc1){ 
  for(int ic1 = symblockinfo.psym()(sc1,I_C,I_BEGIN);ic1 <= symblockinfo.psym()(sc1,I_C,I_END);++ic1){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(ic1);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, ic1, sc1, V2); // V2=(IR-COV index) 
    FC_FUNC(g_if_sigma_oovv_ooov_y9, G_IF_SIGMA_OOVV_OOOV_Y9)
      (sc1, ic1, V2_sym.cptr(), Y5.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sa, X);
    FC_FUNC(g_if_sigma_oovv_ooov_no0_x9, G_IF_SIGMA_OOVV_OOOV_NO0_X9)
      (sa, ia, T2b.cptr(), Xaaa.cptr(), nir, nsym, psym);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      FC_FUNC(g_if_sigma_oovv_ooov_no1_x9, G_IF_SIGMA_OOVV_OOOV_NO1_X9)
        (sa, ia, sc, ic, Xaaa.cptr(), Y5.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.10
  //* X(i,k,o1,c)  <--  (    1.00000000)  D2(i,o3,k,o2) T2(o3,o2,o1,c) 
  //* S2(i,k,a,c)  <--  (    2.00000000) X(i,k,o1,c) Y6(o1,a) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nocc, nvir);
  orz::DTensor Y6 = orz::ct::sympack_Xav(symblockinfo, 0, Y);
  for(int sc1 = 0;sc1 < nir;++sc1){ 
  for(int ic1 = symblockinfo.psym()(sc1,I_C,I_BEGIN);ic1 <= symblockinfo.psym()(sc1,I_C,I_END);++ic1){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(ic1);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, ic1, sc1, V2); // V2=(IR-COV index) 
    FC_FUNC(g_if_sigma_oovv_ooov_y10, G_IF_SIGMA_OOVV_OOOV_Y10)
      (sc1, ic1, V2_sym.cptr(), Y6.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sc, X);
    FC_FUNC(g_if_sigma_oovv_ooov_no0_x10, G_IF_SIGMA_OOVV_OOOV_NO0_X10)
      (sc, ic, T2b.cptr(), Xaaa.cptr(), nir, nsym, psym);
    FC_FUNC(g_if_sigma_oovv_ooov_no1_x10, G_IF_SIGMA_OOVV_OOOV_NO1_X10)
      (sc, ic, Xaaa.cptr(), Y6.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.11
  //* X(i,k,o1,c)  <--  (    1.00000000)  D2(i,o3,k,o2) T2(o3,o2,o1,c) 
  //* S2(i,k,a,c)  <--  (   -1.00000000) X(i,k,o1,c) Y7(o1,a) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nocc, nvir);
  orz::DTensor Y7 = orz::ct::sympack_Xav(symblockinfo, 0, Y);
  for(int sc1 = 0;sc1 < nir;++sc1){ 
  for(int ic1 = symblockinfo.psym()(sc1,I_C,I_BEGIN);ic1 <= symblockinfo.psym()(sc1,I_C,I_END);++ic1){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(ic1);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, ic1, sc1, V2); // V2=(IR-COV index) 
    FC_FUNC(g_if_sigma_oovv_ooov_y11, G_IF_SIGMA_OOVV_OOOV_Y11)
      (sc1, ic1, V2_sym.cptr(), Y7.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sc, X);
    FC_FUNC(g_if_sigma_oovv_ooov_no0_x11, G_IF_SIGMA_OOVV_OOOV_NO0_X11)
      (sc, ic, T2b.cptr(), Xaaa.cptr(), nir, nsym, psym);
    FC_FUNC(g_if_sigma_oovv_ooov_no1_x11, G_IF_SIGMA_OOVV_OOOV_NO1_X11)
      (sc, ic, Xaaa.cptr(), Y7.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.12
  //* X(i,o4,k,o6,o2,c)  <--  (    1.00000000)  D4(i,o4,k,o3,o5,o1,o6,o2) V2(c,o3,o1,o5) 
  //* S2(i,k,a,c)  <--  (    1.00000000) T2(o2,o4,o6,a) X(i,o4,k,o6,o2,c) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(ic);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, ic, sc, V2); // V2=(IR-COV index) 
    for(int si = 0;si < nir;++si){ 
    for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
      for(int so4 = 0;so4 < nir;++so4){ 
      for(int io4 = symblockinfo.psym()(so4,I_O,I_BEGIN);io4 <= symblockinfo.psym()(so4,I_O,I_END);++io4){ 
        // Load D4 from disk, or GA ....                                                     
        int imoi = amo2imo[ii] - nclosed;                              
        int imoj = amo2imo[io4] - nclosed;                              
                                                                                             
        orz::DTensor rdm4_ij_sliced = rdm4(orz::Slice(imoi,imoi+1), orz::Slice(imoj,imoj+1), 
                                           orz::Slice(),            orz::Slice(),            
                                           orz::Slice(),            orz::Slice(),            
                                           orz::Slice(),            orz::Slice()).copy();    
        rdm4_sym = orz::ct::sympack_rdm4_2(symblockinfo, ii, si, io4, so4, rdm4_ij_sliced);    
        FC_FUNC(g_if_set_d4,G_IF_SET_D4)(si, so4, ii, io4, rdm4_sym.cptr(), nir, nsym, psym);  
        orz::DTensor X(nocc, nocc, nocc);
        orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, si^so4^sc, X);
        FC_FUNC(g_if_sigma_oovv_ooov_no0_x12, G_IF_SIGMA_OOVV_OOOV_NO0_X12)
          (sc, ic, si, ii, so4, io4, V2_sym.cptr(), Xaaa.cptr(), nir, nsym, psym);
        for(int sa = 0;sa < nir;++sa){ 
        for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
          T2b = T2.get_amp2(ia);
          FC_FUNC(g_if_sigma_oovv_ooov_no1_x12, G_IF_SIGMA_OOVV_OOOV_NO1_X12)
            (sa, ia, sc, ic, si, ii, so4, io4, T2b.cptr(), Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
        }
        }
        FC_FUNC(g_if_unset_d4,G_IF_UNSET_D4)();
      }
      }
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.13
  //* X(o5,o4,o1,o2,o6,o3,c,a)  <--  (    1.00000000)  T2(o5,o4,o1,c) V2(o2,o6,o3,a) 
  //* S2(i,k,a,c)  <--  (    1.00000000) D4(o2,o6,i,o3,k,o4,o1,o5) X(o5,o4,o1,o2,o6,o3,c,a) 
  for(int so2 = 0;so2 < nir;++so2){ 
  for(int io2 = symblockinfo.psym()(so2,I_O,I_BEGIN);io2 <= symblockinfo.psym()(so2,I_O,I_END);++io2){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(io2);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, io2, so2, V2); // V2=(IR-COV index) 
    for(int so6 = 0;so6 < nir;++so6){ 
    for(int io6 = symblockinfo.psym()(so6,I_O,I_BEGIN);io6 <= symblockinfo.psym()(so6,I_O,I_END);++io6){ 
      for(int sc = 0;sc < nir;++sc){ 
      for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
        S2b = orz::DTensor(retval.namps_iamp()[ic]);
        T2b = T2.get_amp2(ic);
        orz::DTensor X(nocc, nocc, nocc, nocc, nvir);
        orz::DTensor Xaaaav = orz::ct::sympack_Xaaaav(symblockinfo, so2^so6^sc, X);
        FC_FUNC(g_if_sigma_oovv_ooov_no0_x13, G_IF_SIGMA_OOVV_OOOV_NO0_X13)
          (sc, ic, so2, io2, so6, io6, T2b.cptr(), V2_sym.cptr(), Xaaaav.cptr(), nir, nsym, psym);
        // Load D4 from disk, or GA ....                                                     
        int imoi = amo2imo[io2] - nclosed;                              
        int imoj = amo2imo[io6] - nclosed;                              
                                                                                             
        orz::DTensor rdm4_ij_sliced = rdm4(orz::Slice(imoi,imoi+1), orz::Slice(imoj,imoj+1), 
                                           orz::Slice(),            orz::Slice(),            
                                           orz::Slice(),            orz::Slice(),            
                                           orz::Slice(),            orz::Slice()).copy();    
        rdm4_sym = orz::ct::sympack_rdm4_2(symblockinfo, io2, so2, io6, so6, rdm4_ij_sliced);    
        FC_FUNC(g_if_set_d4,G_IF_SET_D4)(so2, so6, io2, io6, rdm4_sym.cptr(), nir, nsym, psym);  
        FC_FUNC(g_if_sigma_oovv_ooov_no1_x13, G_IF_SIGMA_OOVV_OOOV_NO1_X13)
          (sc, ic, so2, io2, so6, io6, Xaaaav.cptr(), S2b.cptr(), nir, nsym, psym);
        retval.acc_amp2(ic, S2b);
        FC_FUNC(g_if_unset_d4,G_IF_UNSET_D4)();
      }
      }
    }
    }
  }
  }
  }


  {
  // No.14
  //* X(i,o4,k,o2,o1,c)  <--  (    1.00000000)  D3(i,o4,k,o3,o5,o2) V2(c,o3,o1,o5) 
  //* S2(i,k,a,c)  <--  (    1.00000000) T2(o2,o4,o1,a) X(i,o4,k,o2,o1,c) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(ic);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, ic, sc, V2); // V2=(IR-COV index) 
    orz::DTensor X(nocc, nocc, nocc, nocc, nocc);
    orz::DTensor Xaaaaa = orz::ct::sympack_Xaaaaa(symblockinfo, sc, X);
    FC_FUNC(g_if_sigma_oovv_ooov_no0_x14, G_IF_SIGMA_OOVV_OOOV_NO0_X14)
      (sc, ic, V2_sym.cptr(), Xaaaaa.cptr(), nir, nsym, psym);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      T2b = T2.get_amp2(ia);
      FC_FUNC(g_if_sigma_oovv_ooov_no1_x14, G_IF_SIGMA_OOVV_OOOV_NO1_X14)
        (sa, ia, sc, ic, T2b.cptr(), Xaaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.15
  //* X(i,k,o3,o2,o1,a)  <--  (    1.00000000)  D3(i,o4,k,o3,o5,o2) V2(a,o4,o1,o5) 
  //* S2(i,k,a,c)  <--  (    1.00000000) T2(o2,o3,o1,c) X(i,k,o3,o2,o1,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(ia);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, ia, sa, V2); // V2=(IR-COV index) 
    orz::DTensor X(nocc, nocc, nocc, nocc, nocc);
    orz::DTensor Xaaaaa = orz::ct::sympack_Xaaaaa(symblockinfo, sa, X);
    FC_FUNC(g_if_sigma_oovv_ooov_no0_x15, G_IF_SIGMA_OOVV_OOOV_NO0_X15)
      (sa, ia, V2_sym.cptr(), Xaaaaa.cptr(), nir, nsym, psym);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      T2b = T2.get_amp2(ic);
      FC_FUNC(g_if_sigma_oovv_ooov_no1_x15, G_IF_SIGMA_OOVV_OOOV_NO1_X15)
        (sa, ia, sc, ic, T2b.cptr(), Xaaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.16
  //* X(i,o4,k,o3,o1,c)  <--  (    1.00000000)  D3(i,o4,k,o3,o5,o2) V2(c,o1,o2,o5) 
  //* S2(i,k,a,c)  <--  (    1.00000000) T2(o3,o4,o1,a) X(i,o4,k,o3,o1,c) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(ic);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, ic, sc, V2); // V2=(IR-COV index) 
    orz::DTensor X(nocc, nocc, nocc, nocc, nocc);
    orz::DTensor Xaaaaa = orz::ct::sympack_Xaaaaa(symblockinfo, sc, X);
    FC_FUNC(g_if_sigma_oovv_ooov_no0_x16, G_IF_SIGMA_OOVV_OOOV_NO0_X16)
      (sc, ic, V2_sym.cptr(), Xaaaaa.cptr(), nir, nsym, psym);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      T2b = T2.get_amp2(ia);
      FC_FUNC(g_if_sigma_oovv_ooov_no1_x16, G_IF_SIGMA_OOVV_OOOV_NO1_X16)
        (sa, ia, sc, ic, T2b.cptr(), Xaaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.17
  //* X(i,o4,k,o3,o1,a)  <--  (    1.00000000)  D3(i,o4,k,o3,o5,o2) V2(a,o1,o2,o5) 
  //* S2(i,k,a,c)  <--  (    1.00000000) T2(o4,o3,o1,c) X(i,o4,k,o3,o1,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(ia);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, ia, sa, V2); // V2=(IR-COV index) 
    orz::DTensor X(nocc, nocc, nocc, nocc, nocc);
    orz::DTensor Xaaaaa = orz::ct::sympack_Xaaaaa(symblockinfo, sa, X);
    FC_FUNC(g_if_sigma_oovv_ooov_no0_x17, G_IF_SIGMA_OOVV_OOOV_NO0_X17)
      (sa, ia, V2_sym.cptr(), Xaaaaa.cptr(), nir, nsym, psym);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      T2b = T2.get_amp2(ic);
      FC_FUNC(g_if_sigma_oovv_ooov_no1_x17, G_IF_SIGMA_OOVV_OOOV_NO1_X17)
        (sa, ia, sc, ic, T2b.cptr(), Xaaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.18
  //* X(o1,o2,o4,o3,c,a)  <--  (    1.00000000)  T2(o1,o2,o4,v1) V2(c,v1,o3,a) 
  //* S2(i,k,a,c)  <--  (    1.00000000) D3(i,o3,k,o2,o4,o1) X(o1,o2,o4,o3,c,a) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(ic);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, ic, sc, V2); // V2=(IR-COV index) 
    orz::DTensor X(nocc, nocc, nocc, nocc, nvir);
    orz::DTensor Xaaaav = orz::ct::sympack_Xaaaav(symblockinfo, sc, X);
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_oovv_ooov_no0_x18, G_IF_SIGMA_OOVV_OOOV_NO0_X18)
        (sc, ic, sv1, iv1, T2b.cptr(), V2_sym.cptr(), Xaaaav.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_oovv_ooov_no1_x18, G_IF_SIGMA_OOVV_OOOV_NO1_X18)
      (sc, ic, Xaaaav.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.19
  //* X(o4,o2,o1,o3,c,a)  <--  (    1.00000000)  T2(o4,o2,o1,v1) V2(c,o3,a,v1) 
  //* S2(i,k,a,c)  <--  (    1.00000000) D3(i,o2,k,o3,o1,o4) X(o4,o2,o1,o3,c,a) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(ic);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, ic, sc, V2); // V2=(IR-COV index) 
    orz::DTensor X(nocc, nocc, nocc, nocc, nvir);
    orz::DTensor Xaaaav = orz::ct::sympack_Xaaaav(symblockinfo, sc, X);
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_oovv_ooov_no0_x19, G_IF_SIGMA_OOVV_OOOV_NO0_X19)
        (sc, ic, sv1, iv1, T2b.cptr(), V2_sym.cptr(), Xaaaav.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_oovv_ooov_no1_x19, G_IF_SIGMA_OOVV_OOOV_NO1_X19)
      (sc, ic, Xaaaav.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.20
  //* X(o3,o2,c,a)  <--  (    1.00000000)  T2(o3,o2,o1,v1) V2(c,v1,o1,a) 
  //* S2(i,k,a,c)  <--  (    1.00000000) D2(i,o3,k,o2) X(o3,o2,c,a) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(ic);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, ic, sc, V2); // V2=(IR-COV index) 
    orz::DTensor X(nocc, nocc, nvir);
    orz::DTensor Xaav = orz::ct::sympack_Xaav(symblockinfo, sc, X);
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_oovv_ooov_no0_x20, G_IF_SIGMA_OOVV_OOOV_NO0_X20)
        (sc, ic, sv1, iv1, T2b.cptr(), V2_sym.cptr(), Xaav.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_oovv_ooov_no1_x20, G_IF_SIGMA_OOVV_OOOV_NO1_X20)
      (sc, ic, Xaav.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.21
  //* X(o2,o1,c,a)  <--  (    1.00000000)  T2(o2,o1,o3,v1) V2(c,o3,a,v1) 
  //* S2(i,k,a,c)  <--  (    1.00000000) D2(i,o1,k,o2) X(o2,o1,c,a) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(ic);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, ic, sc, V2); // V2=(IR-COV index) 
    orz::DTensor X(nocc, nocc, nvir);
    orz::DTensor Xaav = orz::ct::sympack_Xaav(symblockinfo, sc, X);
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_oovv_ooov_no0_x21, G_IF_SIGMA_OOVV_OOOV_NO0_X21)
        (sc, ic, sv1, iv1, T2b.cptr(), V2_sym.cptr(), Xaav.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_oovv_ooov_no1_x21, G_IF_SIGMA_OOVV_OOOV_NO1_X21)
      (sc, ic, Xaav.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }

  for(int ssig = 0;ssig < nir;++ssig){                                                                                      
  for(int isig = symblockinfo.psym()(ssig,I_V,I_BEGIN);isig <= symblockinfo.psym()(ssig,I_V,I_END);++isig){                 
    S2b = retval.get_amp2(isig);                                                                    
    FC_FUNC(g_if_sigma_oovv_scale,G_IF_SIGMA_OOVV_SCALE) // to make hamiltonian symmetric (scale by some factor of overlap) 
      (ssig, isig, S2b.cptr(), nir, nsym, psym);                                                    
    retval.put_amp2(isig, S2b); // S2ija, [b] <<-- Sb                                               
  } // End isig                                                                                                             
  } // End ssig                                                                                                             

  return retval; 
} 
