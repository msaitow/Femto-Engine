                                                                                
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
#include <sci/ctnew2/c_sigma_ooov_ooov.h>                                            
                                                                                
using std::cout;                                                                
using std::endl;                                                                
                                                                                
#define FLOPCOUNT                                                               
                                                                                
//      ______                  __           
//     / ____/___   ____ ___   / /_ ____     
//    / /_   / _ \ / __ `__ \ / __// __ \ 
//   / __/  /  __// / / / / // /_ / /_/ /    
//  /_/     \___//_/ /_/ /_/ \__/ \____/  

                                                                                
// ***************************************************************************  
// orz::ct::mrci                                                                
// ***************************************************************************  
									    /*!        
   @brief CT input                                                              
									     */        
                                                                                
orz::ct::BareAmpPack orz::ct::sigma_ooov_ooov(const orz::ct::Input &ctinp,                                    
                                  const orz::ct::SymBlockInfo &symblockinfo,                                
                                  const orz::ct::HintMO &hintmo,                                            
                                  const orz::ct::RdmPack &rdmPack_sym,                                      
                                  const orz::DTensor &rdm4,                                                 
                                  const orz::ct::BareAmpPack &T2,                             
                                  const int num_sigma,
                                  const double Ecas) {
                                                                                                                
                                                                                                                
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
  //* X(o2,o3,o4,a)  <--  (    1.00000000)  T2(o2,o3,o1,a) h(o4,o1) 
  //* S2(i,k,m,a)  <--  (    1.00000000) D3(i,m,k,o3,o4,o2) X(o2,o3,o4,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    T2b = T2.get_amp2(ia);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sa, X);
    FC_FUNC(g_if_sigma_ooov_ooov_no0_x0, G_IF_SIGMA_OOOV_OOOV_NO0_X0)
      (sa, ia, T2b.cptr(), moint1_sym.cptr(), Xaaa.cptr(), nir, nsym, psym);
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x0, G_IF_SIGMA_OOOV_OOOV_NO1_X0)
      (sa, ia, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.1
  //* X(o3,o2,m,a)  <--  (    1.00000000)  T2(o3,o2,o1,a) h(m,o1) 
  //* S2(i,k,m,a)  <--  (    1.00000000) D2(i,o3,k,o2) X(o3,o2,m,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    T2b = T2.get_amp2(ia);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sa, X);
    FC_FUNC(g_if_sigma_ooov_ooov_no0_x1, G_IF_SIGMA_OOOV_OOOV_NO0_X1)
      (sa, ia, T2b.cptr(), moint1_sym.cptr(), Xaaa.cptr(), nir, nsym, psym);
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x1, G_IF_SIGMA_OOOV_OOOV_NO1_X1)
      (sa, ia, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.2
  //* X(o3,o4,o2,a)  <--  (    1.00000000)  T2(o1,o3,o4,a) h(o2,o1) 
  //* S2(i,k,m,a)  <--  (   -1.00000000) D3(i,m,k,o3,o4,o2) X(o3,o4,o2,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    T2b = T2.get_amp2(ia);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sa, X);
    FC_FUNC(g_if_sigma_ooov_ooov_no0_x2, G_IF_SIGMA_OOOV_OOOV_NO0_X2)
      (sa, ia, T2b.cptr(), moint1_sym.cptr(), Xaaa.cptr(), nir, nsym, psym);
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x2, G_IF_SIGMA_OOOV_OOOV_NO1_X2)
      (sa, ia, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.3
  //* X(o2,m,o3,a)  <--  (    1.00000000)  T2(o1,o2,m,a) h(o3,o1) 
  //* S2(i,k,m,a)  <--  (   -1.00000000) D2(i,o3,k,o2) X(o2,m,o3,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    T2b = T2.get_amp2(ia);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sa, X);
    FC_FUNC(g_if_sigma_ooov_ooov_no0_x3, G_IF_SIGMA_OOOV_OOOV_NO0_X3)
      (sa, ia, T2b.cptr(), moint1_sym.cptr(), Xaaa.cptr(), nir, nsym, psym);
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x3, G_IF_SIGMA_OOOV_OOOV_NO1_X3)
      (sa, ia, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.4
  //* X(o2,o4,o3,a)  <--  (    1.00000000)  T2(o2,o1,o4,a) h(o3,o1) 
  //* S2(i,k,m,a)  <--  (   -1.00000000) D3(i,m,k,o3,o4,o2) X(o2,o4,o3,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    T2b = T2.get_amp2(ia);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sa, X);
    FC_FUNC(g_if_sigma_ooov_ooov_no0_x4, G_IF_SIGMA_OOOV_OOOV_NO0_X4)
      (sa, ia, T2b.cptr(), moint1_sym.cptr(), Xaaa.cptr(), nir, nsym, psym);
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x4, G_IF_SIGMA_OOOV_OOOV_NO1_X4)
      (sa, ia, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.5
  //* X(o3,m,o2,a)  <--  (    1.00000000)  T2(o3,o1,m,a) h(o2,o1) 
  //* S2(i,k,m,a)  <--  (   -1.00000000) D2(i,o3,k,o2) X(o3,m,o2,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    T2b = T2.get_amp2(ia);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sa, X);
    FC_FUNC(g_if_sigma_ooov_ooov_no0_x5, G_IF_SIGMA_OOOV_OOOV_NO0_X5)
      (sa, ia, T2b.cptr(), moint1_sym.cptr(), Xaaa.cptr(), nir, nsym, psym);
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x5, G_IF_SIGMA_OOOV_OOOV_NO1_X5)
      (sa, ia, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.6
  //* X(o3,o2,o1,a)  <--  (    1.00000000)  T2(o3,o2,o1,v1) h(a,v1) 
  //* S2(i,k,m,a)  <--  (    1.00000000) D3(i,m,k,o2,o1,o3) X(o3,o2,o1,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sa, X);
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ooov_ooov_no0_x6, G_IF_SIGMA_OOOV_OOOV_NO0_X6)
        (sa, ia, sv1, iv1, T2b.cptr(), moint1_sym.cptr(), Xaaa.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x6, G_IF_SIGMA_OOOV_OOOV_NO1_X6)
      (sa, ia, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.7
  //* X(o1,o2,m,a)  <--  (    1.00000000)  T2(o1,o2,m,v1) h(a,v1) 
  //* S2(i,k,m,a)  <--  (    1.00000000) D2(i,o1,k,o2) X(o1,o2,m,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sa, X);
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ooov_ooov_no0_x7, G_IF_SIGMA_OOOV_OOOV_NO0_X7)
        (sa, ia, sv1, iv1, T2b.cptr(), moint1_sym.cptr(), Xaaa.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x7, G_IF_SIGMA_OOOV_OOOV_NO1_X7)
      (sa, ia, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.8
  //* X(o2,o3,o4,a)  <--  (    1.00000000)  T2(o2,o3,o1,a) Y0(o1,o4) 
  //* S2(i,k,m,a)  <--  (    2.00000000) D3(i,m,k,o3,o4,o2) X(o2,o3,o4,a) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nocc, nocc);
  orz::DTensor Y0 = orz::ct::sympack_Xaa(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_ooov_ooov_y8, G_IF_SIGMA_OOOV_OOOV_Y8)
      (sc1, ic1, V2_sym.cptr(), Y0.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    T2b = T2.get_amp2(ia);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sa, X);
    FC_FUNC(g_if_sigma_ooov_ooov_no0_x8, G_IF_SIGMA_OOOV_OOOV_NO0_X8)
      (sa, ia, T2b.cptr(), Y0.cptr(), Xaaa.cptr(), nir, nsym, psym);
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x8, G_IF_SIGMA_OOOV_OOOV_NO1_X8)
      (sa, ia, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.9
  //* X(o2,o3,o4,a)  <--  (    1.00000000)  T2(o2,o3,o1,a) Y1(o1,o4) 
  //* S2(i,k,m,a)  <--  (   -1.00000000) D3(i,m,k,o3,o4,o2) X(o2,o3,o4,a) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nocc, nocc);
  orz::DTensor Y1 = orz::ct::sympack_Xaa(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_ooov_ooov_y9, G_IF_SIGMA_OOOV_OOOV_Y9)
      (sc1, ic1, V2_sym.cptr(), Y1.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    T2b = T2.get_amp2(ia);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sa, X);
    FC_FUNC(g_if_sigma_ooov_ooov_no0_x9, G_IF_SIGMA_OOOV_OOOV_NO0_X9)
      (sa, ia, T2b.cptr(), Y1.cptr(), Xaaa.cptr(), nir, nsym, psym);
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x9, G_IF_SIGMA_OOOV_OOOV_NO1_X9)
      (sa, ia, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.10
  //* X(o3,o2,m,a)  <--  (    1.00000000)  T2(o3,o2,o1,a) Y2(m,o1) 
  //* S2(i,k,m,a)  <--  (    2.00000000) D2(i,o3,k,o2) X(o3,o2,m,a) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nocc, nocc);
  orz::DTensor Y2 = orz::ct::sympack_Xaa(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_ooov_ooov_y10, G_IF_SIGMA_OOOV_OOOV_Y10)
      (sc1, ic1, V2_sym.cptr(), Y2.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    T2b = T2.get_amp2(ia);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sa, X);
    FC_FUNC(g_if_sigma_ooov_ooov_no0_x10, G_IF_SIGMA_OOOV_OOOV_NO0_X10)
      (sa, ia, T2b.cptr(), Y2.cptr(), Xaaa.cptr(), nir, nsym, psym);
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x10, G_IF_SIGMA_OOOV_OOOV_NO1_X10)
      (sa, ia, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.11
  //* X(o3,o2,m,a)  <--  (    1.00000000)  T2(o3,o2,o1,a) Y3(m,o1) 
  //* S2(i,k,m,a)  <--  (   -1.00000000) D2(i,o3,k,o2) X(o3,o2,m,a) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nocc, nocc);
  orz::DTensor Y3 = orz::ct::sympack_Xaa(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_ooov_ooov_y11, G_IF_SIGMA_OOOV_OOOV_Y11)
      (sc1, ic1, V2_sym.cptr(), Y3.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    T2b = T2.get_amp2(ia);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sa, X);
    FC_FUNC(g_if_sigma_ooov_ooov_no0_x11, G_IF_SIGMA_OOOV_OOOV_NO0_X11)
      (sa, ia, T2b.cptr(), Y3.cptr(), Xaaa.cptr(), nir, nsym, psym);
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x11, G_IF_SIGMA_OOOV_OOOV_NO1_X11)
      (sa, ia, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.12
  //* X(o3,o4,o2,a)  <--  (    1.00000000)  T2(o1,o3,o4,a) Y4(o1,o2) 
  //* S2(i,k,m,a)  <--  (   -2.00000000) D3(i,m,k,o3,o4,o2) X(o3,o4,o2,a) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nocc, nocc);
  orz::DTensor Y4 = orz::ct::sympack_Xaa(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_ooov_ooov_y12, G_IF_SIGMA_OOOV_OOOV_Y12)
      (sc1, ic1, V2_sym.cptr(), Y4.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    T2b = T2.get_amp2(ia);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sa, X);
    FC_FUNC(g_if_sigma_ooov_ooov_no0_x12, G_IF_SIGMA_OOOV_OOOV_NO0_X12)
      (sa, ia, T2b.cptr(), Y4.cptr(), Xaaa.cptr(), nir, nsym, psym);
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x12, G_IF_SIGMA_OOOV_OOOV_NO1_X12)
      (sa, ia, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.13
  //* X(o3,o4,o2,a)  <--  (    1.00000000)  T2(o1,o3,o4,a) Y5(o1,o2) 
  //* S2(i,k,m,a)  <--  (    1.00000000) D3(i,m,k,o3,o4,o2) X(o3,o4,o2,a) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nocc, nocc);
  orz::DTensor Y5 = orz::ct::sympack_Xaa(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_ooov_ooov_y13, G_IF_SIGMA_OOOV_OOOV_Y13)
      (sc1, ic1, V2_sym.cptr(), Y5.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    T2b = T2.get_amp2(ia);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sa, X);
    FC_FUNC(g_if_sigma_ooov_ooov_no0_x13, G_IF_SIGMA_OOOV_OOOV_NO0_X13)
      (sa, ia, T2b.cptr(), Y5.cptr(), Xaaa.cptr(), nir, nsym, psym);
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x13, G_IF_SIGMA_OOOV_OOOV_NO1_X13)
      (sa, ia, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.14
  //* X(o2,m,o3,a)  <--  (    1.00000000)  T2(o1,o2,m,a) Y6(o1,o3) 
  //* S2(i,k,m,a)  <--  (   -2.00000000) D2(i,o3,k,o2) X(o2,m,o3,a) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nocc, nocc);
  orz::DTensor Y6 = orz::ct::sympack_Xaa(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_ooov_ooov_y14, G_IF_SIGMA_OOOV_OOOV_Y14)
      (sc1, ic1, V2_sym.cptr(), Y6.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    T2b = T2.get_amp2(ia);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sa, X);
    FC_FUNC(g_if_sigma_ooov_ooov_no0_x14, G_IF_SIGMA_OOOV_OOOV_NO0_X14)
      (sa, ia, T2b.cptr(), Y6.cptr(), Xaaa.cptr(), nir, nsym, psym);
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x14, G_IF_SIGMA_OOOV_OOOV_NO1_X14)
      (sa, ia, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.15
  //* X(o2,m,o3,a)  <--  (    1.00000000)  T2(o1,o2,m,a) Y7(o1,o3) 
  //* S2(i,k,m,a)  <--  (    1.00000000) D2(i,o3,k,o2) X(o2,m,o3,a) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nocc, nocc);
  orz::DTensor Y7 = orz::ct::sympack_Xaa(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_ooov_ooov_y15, G_IF_SIGMA_OOOV_OOOV_Y15)
      (sc1, ic1, V2_sym.cptr(), Y7.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    T2b = T2.get_amp2(ia);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sa, X);
    FC_FUNC(g_if_sigma_ooov_ooov_no0_x15, G_IF_SIGMA_OOOV_OOOV_NO0_X15)
      (sa, ia, T2b.cptr(), Y7.cptr(), Xaaa.cptr(), nir, nsym, psym);
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x15, G_IF_SIGMA_OOOV_OOOV_NO1_X15)
      (sa, ia, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.16
  //* X(o2,o4,o3,a)  <--  (    1.00000000)  T2(o2,o1,o4,a) Y8(o1,o3) 
  //* S2(i,k,m,a)  <--  (   -2.00000000) D3(i,m,k,o3,o4,o2) X(o2,o4,o3,a) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nocc, nocc);
  orz::DTensor Y8 = orz::ct::sympack_Xaa(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_ooov_ooov_y16, G_IF_SIGMA_OOOV_OOOV_Y16)
      (sc1, ic1, V2_sym.cptr(), Y8.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    T2b = T2.get_amp2(ia);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sa, X);
    FC_FUNC(g_if_sigma_ooov_ooov_no0_x16, G_IF_SIGMA_OOOV_OOOV_NO0_X16)
      (sa, ia, T2b.cptr(), Y8.cptr(), Xaaa.cptr(), nir, nsym, psym);
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x16, G_IF_SIGMA_OOOV_OOOV_NO1_X16)
      (sa, ia, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.17
  //* X(o2,o4,o3,a)  <--  (    1.00000000)  T2(o2,o1,o4,a) Y9(o1,o3) 
  //* S2(i,k,m,a)  <--  (    1.00000000) D3(i,m,k,o3,o4,o2) X(o2,o4,o3,a) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nocc, nocc);
  orz::DTensor Y9 = orz::ct::sympack_Xaa(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_ooov_ooov_y17, G_IF_SIGMA_OOOV_OOOV_Y17)
      (sc1, ic1, V2_sym.cptr(), Y9.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    T2b = T2.get_amp2(ia);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sa, X);
    FC_FUNC(g_if_sigma_ooov_ooov_no0_x17, G_IF_SIGMA_OOOV_OOOV_NO0_X17)
      (sa, ia, T2b.cptr(), Y9.cptr(), Xaaa.cptr(), nir, nsym, psym);
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x17, G_IF_SIGMA_OOOV_OOOV_NO1_X17)
      (sa, ia, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.18
  //* X(o3,m,o2,a)  <--  (    1.00000000)  T2(o3,o1,m,a) Y10(o1,o2) 
  //* S2(i,k,m,a)  <--  (   -2.00000000) D2(i,o3,k,o2) X(o3,m,o2,a) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nocc, nocc);
  orz::DTensor Y10 = orz::ct::sympack_Xaa(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_ooov_ooov_y18, G_IF_SIGMA_OOOV_OOOV_Y18)
      (sc1, ic1, V2_sym.cptr(), Y10.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    T2b = T2.get_amp2(ia);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sa, X);
    FC_FUNC(g_if_sigma_ooov_ooov_no0_x18, G_IF_SIGMA_OOOV_OOOV_NO0_X18)
      (sa, ia, T2b.cptr(), Y10.cptr(), Xaaa.cptr(), nir, nsym, psym);
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x18, G_IF_SIGMA_OOOV_OOOV_NO1_X18)
      (sa, ia, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.19
  //* X(o3,m,o2,a)  <--  (    1.00000000)  T2(o3,o1,m,a) Y11(o1,o2) 
  //* S2(i,k,m,a)  <--  (    1.00000000) D2(i,o3,k,o2) X(o3,m,o2,a) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nocc, nocc);
  orz::DTensor Y11 = orz::ct::sympack_Xaa(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_ooov_ooov_y19, G_IF_SIGMA_OOOV_OOOV_Y19)
      (sc1, ic1, V2_sym.cptr(), Y11.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    T2b = T2.get_amp2(ia);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sa, X);
    FC_FUNC(g_if_sigma_ooov_ooov_no0_x19, G_IF_SIGMA_OOOV_OOOV_NO0_X19)
      (sa, ia, T2b.cptr(), Y11.cptr(), Xaaa.cptr(), nir, nsym, psym);
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x19, G_IF_SIGMA_OOOV_OOOV_NO1_X19)
      (sa, ia, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.20
  //* X(o1,o2,o3,a)  <--  (    1.00000000)  T2(o1,o2,o3,v1) Y12(a,v1) 
  //* S2(i,k,m,a)  <--  (    2.00000000) D3(i,m,k,o2,o3,o1) X(o1,o2,o3,a) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nvir, nvir);
  orz::DTensor Y12 = orz::ct::sympack_Xvv(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_ooov_ooov_y20, G_IF_SIGMA_OOOV_OOOV_Y20)
      (sc1, ic1, V2_sym.cptr(), Y12.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sa, X);
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ooov_ooov_no0_x20, G_IF_SIGMA_OOOV_OOOV_NO0_X20)
        (sa, ia, sv1, iv1, T2b.cptr(), Y12.cptr(), Xaaa.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x20, G_IF_SIGMA_OOOV_OOOV_NO1_X20)
      (sa, ia, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.21
  //* X(o1,o2,o3,a)  <--  (    1.00000000)  T2(o1,o2,o3,v1) Y13(a,v1) 
  //* S2(i,k,m,a)  <--  (   -1.00000000) D3(i,m,k,o2,o3,o1) X(o1,o2,o3,a) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nvir, nvir);
  orz::DTensor Y13 = orz::ct::sympack_Xvv(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_ooov_ooov_y21, G_IF_SIGMA_OOOV_OOOV_Y21)
      (sc1, ic1, V2_sym.cptr(), Y13.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sa, X);
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ooov_ooov_no0_x21, G_IF_SIGMA_OOOV_OOOV_NO0_X21)
        (sa, ia, sv1, iv1, T2b.cptr(), Y13.cptr(), Xaaa.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x21, G_IF_SIGMA_OOOV_OOOV_NO1_X21)
      (sa, ia, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.22
  //* X(o2,o1,m,a)  <--  (    1.00000000)  T2(o2,o1,m,v1) Y14(a,v1) 
  //* S2(i,k,m,a)  <--  (    2.00000000) D2(i,o2,k,o1) X(o2,o1,m,a) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nvir, nvir);
  orz::DTensor Y14 = orz::ct::sympack_Xvv(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_ooov_ooov_y22, G_IF_SIGMA_OOOV_OOOV_Y22)
      (sc1, ic1, V2_sym.cptr(), Y14.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sa, X);
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ooov_ooov_no0_x22, G_IF_SIGMA_OOOV_OOOV_NO0_X22)
        (sa, ia, sv1, iv1, T2b.cptr(), Y14.cptr(), Xaaa.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x22, G_IF_SIGMA_OOOV_OOOV_NO1_X22)
      (sa, ia, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.23
  //* X(o2,o1,m,a)  <--  (    1.00000000)  T2(o2,o1,m,v1) Y15(a,v1) 
  //* S2(i,k,m,a)  <--  (   -1.00000000) D2(i,o2,k,o1) X(o2,o1,m,a) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nvir, nvir);
  orz::DTensor Y15 = orz::ct::sympack_Xvv(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_ooov_ooov_y23, G_IF_SIGMA_OOOV_OOOV_Y23)
      (sc1, ic1, V2_sym.cptr(), Y15.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sa, X);
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ooov_ooov_no0_x23, G_IF_SIGMA_OOOV_OOOV_NO0_X23)
        (sa, ia, sv1, iv1, T2b.cptr(), Y15.cptr(), Xaaa.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x23, G_IF_SIGMA_OOOV_OOOV_NO1_X23)
      (sa, ia, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.24
  //* X(o3,i,m,k,o4,o1)  <--  (    1.00000000)  D4(o6,o3,i,m,k,o4,o5,o2) V2(o6,o1,o2,o5) 
  //* S2(i,k,m,a)  <--  (    1.00000000) T2(o3,o4,o1,a) X(o3,i,m,k,o4,o1) 
  for(int so3 = 0;so3 < nir;++so3){ 
  for(int io3 = symblockinfo.psym()(so3,I_O,I_BEGIN);io3 <= symblockinfo.psym()(so3,I_O,I_END);++io3){ 
    orz::DTensor X(nocc, nocc, nocc, nocc, nocc);
    orz::DTensor Xaaaaa = orz::ct::sympack_Xaaaaa(symblockinfo, so3, X);
    for(int so6 = 0;so6 < nir;++so6){ 
    for(int io6 = symblockinfo.psym()(so6,I_O,I_BEGIN);io6 <= symblockinfo.psym()(so6,I_O,I_END);++io6){ 
      // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
      V2 <<= 0.0;                                                                          
      shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(io6);
      for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
        // Load a signle record of integals                                                             
        const int &imo2 = loadbuf_ptr->i0;                                                              
        const int &imo3 = loadbuf_ptr->i1;                                                              
        const int &imo4 = loadbuf_ptr->i2;                                                              
        const double &v = loadbuf_ptr->v;                                                               
        V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
      }                                                                                                 
      const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, io6, so6, V2); // V2=(IR-COV index) 
      // Load D4 from disk, or GA ....                                                     
      int imoi = amo2imo[io6] - nclosed;                              
      int imoj = amo2imo[io3] - nclosed;                              
                                                                                           
      orz::DTensor rdm4_ij_sliced = rdm4(orz::Slice(imoi,imoi+1), orz::Slice(imoj,imoj+1), 
                                         orz::Slice(),            orz::Slice(),            
                                         orz::Slice(),            orz::Slice(),            
                                         orz::Slice(),            orz::Slice()).copy();    
      rdm4_sym = orz::ct::sympack_rdm4_2(symblockinfo, io6, so6, io3, so3, rdm4_ij_sliced);    
      FC_FUNC(g_if_set_d4,G_IF_SET_D4)(so6, so3, io6, io3, rdm4_sym.cptr(), nir, nsym, psym);  
      FC_FUNC(g_if_sigma_ooov_ooov_no0_x24, G_IF_SIGMA_OOOV_OOOV_NO0_X24)
        (so3, io3, so6, io6, V2_sym.cptr(), Xaaaaa.cptr(), nir, nsym, psym);
      FC_FUNC(g_if_unset_d4,G_IF_UNSET_D4)();
    }
    }
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      S2b = orz::DTensor(retval.namps_iamp()[ia]);
      T2b = T2.get_amp2(ia);
      FC_FUNC(g_if_sigma_ooov_ooov_no1_x24, G_IF_SIGMA_OOOV_OOOV_NO1_X24)
        (sa, ia, so3, io3, T2b.cptr(), Xaaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ia, S2b);
    }
    }
  }
  }
  }


  {
  // No.25
  //* X(i,o4,k,o3,m,o1)  <--  (    1.00000000)  D3(i,o4,k,o3,o5,o2) V2(m,o1,o2,o5) 
  //* S2(i,k,m,a)  <--  (    1.00000000) T2(o4,o3,o1,a) X(i,o4,k,o3,m,o1) 
  for(int sm = 0;sm < nir;++sm){ 
  for(int im = symblockinfo.psym()(sm,I_O,I_BEGIN);im <= symblockinfo.psym()(sm,I_O,I_END);++im){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(im);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, im, sm, V2); // V2=(IR-COV index) 
    orz::DTensor X(nocc, nocc, nocc, nocc, nocc);
    orz::DTensor Xaaaaa = orz::ct::sympack_Xaaaaa(symblockinfo, sm, X);
    FC_FUNC(g_if_sigma_ooov_ooov_no0_x25, G_IF_SIGMA_OOOV_OOOV_NO0_X25)
      (sm, im, V2_sym.cptr(), Xaaaaa.cptr(), nir, nsym, psym);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      S2b = orz::DTensor(retval.namps_iamp()[ia]);
      T2b = T2.get_amp2(ia);
      FC_FUNC(g_if_sigma_ooov_ooov_no1_x25, G_IF_SIGMA_OOOV_OOOV_NO1_X25)
        (sa, ia, sm, im, T2b.cptr(), Xaaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ia, S2b);
    }
    }
  }
  }
  }


  {
  // No.26
  //* X(i,k,o3,o2,m,o1)  <--  (    1.00000000)  D3(i,o4,k,o3,o5,o2) V2(m,o4,o1,o5) 
  //* S2(i,k,m,a)  <--  (    1.00000000) T2(o2,o3,o1,a) X(i,k,o3,o2,m,o1) 
  for(int sm = 0;sm < nir;++sm){ 
  for(int im = symblockinfo.psym()(sm,I_O,I_BEGIN);im <= symblockinfo.psym()(sm,I_O,I_END);++im){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(im);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, im, sm, V2); // V2=(IR-COV index) 
    orz::DTensor X(nocc, nocc, nocc, nocc, nocc);
    orz::DTensor Xaaaaa = orz::ct::sympack_Xaaaaa(symblockinfo, sm, X);
    FC_FUNC(g_if_sigma_ooov_ooov_no0_x26, G_IF_SIGMA_OOOV_OOOV_NO0_X26)
      (sm, im, V2_sym.cptr(), Xaaaaa.cptr(), nir, nsym, psym);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      S2b = orz::DTensor(retval.namps_iamp()[ia]);
      T2b = T2.get_amp2(ia);
      FC_FUNC(g_if_sigma_ooov_ooov_no1_x26, G_IF_SIGMA_OOOV_OOOV_NO1_X26)
        (sa, ia, sm, im, T2b.cptr(), Xaaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ia, S2b);
    }
    }
  }
  }
  }


  {
  // No.27
  //* X(o6,m,i,o4,k,o1)  <--  (    1.00000000)  D4(o3,o6,m,i,o4,k,o2,o5) V2(o3,o1,o2,o5) 
  //* S2(i,k,m,a)  <--  (   -1.00000000) T2(o1,o4,o6,a) X(o6,m,i,o4,k,o1) 
  for(int so6 = 0;so6 < nir;++so6){ 
  for(int io6 = symblockinfo.psym()(so6,I_O,I_BEGIN);io6 <= symblockinfo.psym()(so6,I_O,I_END);++io6){ 
    orz::DTensor X(nocc, nocc, nocc, nocc, nocc);
    orz::DTensor Xaaaaa = orz::ct::sympack_Xaaaaa(symblockinfo, so6, X);
    for(int so3 = 0;so3 < nir;++so3){ 
    for(int io3 = symblockinfo.psym()(so3,I_O,I_BEGIN);io3 <= symblockinfo.psym()(so3,I_O,I_END);++io3){ 
      // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
      V2 <<= 0.0;                                                                          
      shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(io3);
      for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
        // Load a signle record of integals                                                             
        const int &imo2 = loadbuf_ptr->i0;                                                              
        const int &imo3 = loadbuf_ptr->i1;                                                              
        const int &imo4 = loadbuf_ptr->i2;                                                              
        const double &v = loadbuf_ptr->v;                                                               
        V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
      }                                                                                                 
      const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, io3, so3, V2); // V2=(IR-COV index) 
      // Load D4 from disk, or GA ....                                                     
      int imoi = amo2imo[io3] - nclosed;                              
      int imoj = amo2imo[io6] - nclosed;                              
                                                                                           
      orz::DTensor rdm4_ij_sliced = rdm4(orz::Slice(imoi,imoi+1), orz::Slice(imoj,imoj+1), 
                                         orz::Slice(),            orz::Slice(),            
                                         orz::Slice(),            orz::Slice(),            
                                         orz::Slice(),            orz::Slice()).copy();    
      rdm4_sym = orz::ct::sympack_rdm4_2(symblockinfo, io3, so3, io6, so6, rdm4_ij_sliced);    
      FC_FUNC(g_if_set_d4,G_IF_SET_D4)(so3, so6, io3, io6, rdm4_sym.cptr(), nir, nsym, psym);  
      FC_FUNC(g_if_sigma_ooov_ooov_no0_x27, G_IF_SIGMA_OOOV_OOOV_NO0_X27)
        (so3, io3, so6, io6, V2_sym.cptr(), Xaaaaa.cptr(), nir, nsym, psym);
      FC_FUNC(g_if_unset_d4,G_IF_UNSET_D4)();
    }
    }
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      S2b = orz::DTensor(retval.namps_iamp()[ia]);
      T2b = T2.get_amp2(ia);
      FC_FUNC(g_if_sigma_ooov_ooov_no1_x27, G_IF_SIGMA_OOOV_OOOV_NO1_X27)
        (sa, ia, so6, io6, T2b.cptr(), Xaaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ia, S2b);
    }
    }
  }
  }
  }


  {
  // No.28
  //* X(i,k,o3,o1)  <--  (    1.00000000)  D3(i,o4,k,o3,o5,o2) V2(o1,o4,o2,o5) 
  //* S2(i,k,m,a)  <--  (   -1.00000000) T2(o1,o3,m,a) X(i,k,o3,o1) 
  for(int so1 = 0;so1 < nir;++so1){ 
  for(int io1 = symblockinfo.psym()(so1,I_O,I_BEGIN);io1 <= symblockinfo.psym()(so1,I_O,I_END);++io1){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(io1);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, io1, so1, V2); // V2=(IR-COV index) 
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, so1, X);
    FC_FUNC(g_if_sigma_ooov_ooov_no0_x28, G_IF_SIGMA_OOOV_OOOV_NO0_X28)
      (so1, io1, V2_sym.cptr(), Xaaa.cptr(), nir, nsym, psym);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      S2b = orz::DTensor(retval.namps_iamp()[ia]);
      T2b = T2.get_amp2(ia);
      FC_FUNC(g_if_sigma_ooov_ooov_no1_x28, G_IF_SIGMA_OOOV_OOOV_NO1_X28)
        (sa, ia, so1, io1, T2b.cptr(), Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ia, S2b);
    }
    }
  }
  }
  }


  {
  // No.29
  //* X(i,k,o3,o5,m,o1)  <--  (    1.00000000)  D3(i,o4,k,o3,o5,o2) V2(m,o4,o1,o2) 
  //* S2(i,k,m,a)  <--  (   -1.00000000) T2(o1,o3,o5,a) X(i,k,o3,o5,m,o1) 
  for(int sm = 0;sm < nir;++sm){ 
  for(int im = symblockinfo.psym()(sm,I_O,I_BEGIN);im <= symblockinfo.psym()(sm,I_O,I_END);++im){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(im);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, im, sm, V2); // V2=(IR-COV index) 
    orz::DTensor X(nocc, nocc, nocc, nocc, nocc);
    orz::DTensor Xaaaaa = orz::ct::sympack_Xaaaaa(symblockinfo, sm, X);
    FC_FUNC(g_if_sigma_ooov_ooov_no0_x29, G_IF_SIGMA_OOOV_OOOV_NO0_X29)
      (sm, im, V2_sym.cptr(), Xaaaaa.cptr(), nir, nsym, psym);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      S2b = orz::DTensor(retval.namps_iamp()[ia]);
      T2b = T2.get_amp2(ia);
      FC_FUNC(g_if_sigma_ooov_ooov_no1_x29, G_IF_SIGMA_OOOV_OOOV_NO1_X29)
        (sa, ia, sm, im, T2b.cptr(), Xaaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ia, S2b);
    }
    }
  }
  }
  }


  {
  // No.30
  //* X(k,m,i,o3,o6,o1)  <--  (    1.00000000)  D4(o4,k,m,i,o2,o5,o3,o6) V2(o4,o1,o2,o5) 
  //* S2(i,k,m,a)  <--  (   -1.00000000) T2(o3,o1,o6,a) X(k,m,i,o3,o6,o1) 
  for(int sk = 0;sk < nir;++sk){ 
  for(int ik = symblockinfo.psym()(sk,I_O,I_BEGIN);ik <= symblockinfo.psym()(sk,I_O,I_END);++ik){ 
    orz::DTensor X(nocc, nocc, nocc, nocc, nocc);
    orz::DTensor Xaaaaa = orz::ct::sympack_Xaaaaa(symblockinfo, sk, X);
    for(int so4 = 0;so4 < nir;++so4){ 
    for(int io4 = symblockinfo.psym()(so4,I_O,I_BEGIN);io4 <= symblockinfo.psym()(so4,I_O,I_END);++io4){ 
      // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
      V2 <<= 0.0;                                                                          
      shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(io4);
      for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
        // Load a signle record of integals                                                             
        const int &imo2 = loadbuf_ptr->i0;                                                              
        const int &imo3 = loadbuf_ptr->i1;                                                              
        const int &imo4 = loadbuf_ptr->i2;                                                              
        const double &v = loadbuf_ptr->v;                                                               
        V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
      }                                                                                                 
      const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, io4, so4, V2); // V2=(IR-COV index) 
      // Load D4 from disk, or GA ....                                                     
      int imoi = amo2imo[io4] - nclosed;                              
      int imoj = amo2imo[ik] - nclosed;                              
                                                                                           
      orz::DTensor rdm4_ij_sliced = rdm4(orz::Slice(imoi,imoi+1), orz::Slice(imoj,imoj+1), 
                                         orz::Slice(),            orz::Slice(),            
                                         orz::Slice(),            orz::Slice(),            
                                         orz::Slice(),            orz::Slice()).copy();    
      rdm4_sym = orz::ct::sympack_rdm4_2(symblockinfo, io4, so4, ik, sk, rdm4_ij_sliced);    
      FC_FUNC(g_if_set_d4,G_IF_SET_D4)(so4, sk, io4, ik, rdm4_sym.cptr(), nir, nsym, psym);  
      FC_FUNC(g_if_sigma_ooov_ooov_no0_x30, G_IF_SIGMA_OOOV_OOOV_NO0_X30)
        (sk, ik, so4, io4, V2_sym.cptr(), Xaaaaa.cptr(), nir, nsym, psym);
      FC_FUNC(g_if_unset_d4,G_IF_UNSET_D4)();
    }
    }
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      S2b = orz::DTensor(retval.namps_iamp()[ia]);
      T2b = T2.get_amp2(ia);
      FC_FUNC(g_if_sigma_ooov_ooov_no1_x30, G_IF_SIGMA_OOOV_OOOV_NO1_X30)
        (sa, ia, sk, ik, T2b.cptr(), Xaaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ia, S2b);
    }
    }
  }
  }
  }


  {
  // No.31
  //* X(i,o4,k,o1)  <--  (    1.00000000)  D3(i,o4,k,o3,o5,o2) V2(o1,o3,o2,o5) 
  //* S2(i,k,m,a)  <--  (   -1.00000000) T2(o4,o1,m,a) X(i,o4,k,o1) 
  for(int so1 = 0;so1 < nir;++so1){ 
  for(int io1 = symblockinfo.psym()(so1,I_O,I_BEGIN);io1 <= symblockinfo.psym()(so1,I_O,I_END);++io1){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(io1);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, io1, so1, V2); // V2=(IR-COV index) 
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, so1, X);
    FC_FUNC(g_if_sigma_ooov_ooov_no0_x31, G_IF_SIGMA_OOOV_OOOV_NO0_X31)
      (so1, io1, V2_sym.cptr(), Xaaa.cptr(), nir, nsym, psym);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      S2b = orz::DTensor(retval.namps_iamp()[ia]);
      T2b = T2.get_amp2(ia);
      FC_FUNC(g_if_sigma_ooov_ooov_no1_x31, G_IF_SIGMA_OOOV_OOOV_NO1_X31)
        (sa, ia, so1, io1, T2b.cptr(), Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ia, S2b);
    }
    }
  }
  }
  }


  {
  // No.32
  //* X(i,k,o5,o2,m,o1)  <--  (    1.00000000)  D3(i,o4,k,o3,o5,o2) V2(m,o4,o1,o3) 
  //* S2(i,k,m,a)  <--  (   -1.00000000) T2(o2,o1,o5,a) X(i,k,o5,o2,m,o1) 
  for(int sm = 0;sm < nir;++sm){ 
  for(int im = symblockinfo.psym()(sm,I_O,I_BEGIN);im <= symblockinfo.psym()(sm,I_O,I_END);++im){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(im);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, im, sm, V2); // V2=(IR-COV index) 
    orz::DTensor X(nocc, nocc, nocc, nocc, nocc);
    orz::DTensor Xaaaaa = orz::ct::sympack_Xaaaaa(symblockinfo, sm, X);
    FC_FUNC(g_if_sigma_ooov_ooov_no0_x32, G_IF_SIGMA_OOOV_OOOV_NO0_X32)
      (sm, im, V2_sym.cptr(), Xaaaaa.cptr(), nir, nsym, psym);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      S2b = orz::DTensor(retval.namps_iamp()[ia]);
      T2b = T2.get_amp2(ia);
      FC_FUNC(g_if_sigma_ooov_ooov_no1_x32, G_IF_SIGMA_OOOV_OOOV_NO1_X32)
        (sa, ia, sm, im, T2b.cptr(), Xaaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ia, S2b);
    }
    }
  }
  }
  }


  {
  // No.33
  //* X(o5,o4,o3,a)  <--  (    1.00000000)  T2(o2,o1,o5,a) V2(o1,o4,o2,o3) 
  //* S2(i,k,m,a)  <--  (   -1.00000000) D3(i,m,k,o4,o5,o3) X(o5,o4,o3,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    T2b = T2.get_amp2(ia);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sa, X);
    for(int so1 = 0;so1 < nir;++so1){ 
    for(int io1 = symblockinfo.psym()(so1,I_O,I_BEGIN);io1 <= symblockinfo.psym()(so1,I_O,I_END);++io1){ 
      // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
      V2 <<= 0.0;                                                                          
      shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(io1);
      for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
        // Load a signle record of integals                                                             
        const int &imo2 = loadbuf_ptr->i0;                                                              
        const int &imo3 = loadbuf_ptr->i1;                                                              
        const int &imo4 = loadbuf_ptr->i2;                                                              
        const double &v = loadbuf_ptr->v;                                                               
        V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
      }                                                                                                 
      const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, io1, so1, V2); // V2=(IR-COV index) 
      FC_FUNC(g_if_sigma_ooov_ooov_no0_x33, G_IF_SIGMA_OOOV_OOOV_NO0_X33)
        (sa, ia, so1, io1, T2b.cptr(), V2_sym.cptr(), Xaaa.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x33, G_IF_SIGMA_OOOV_OOOV_NO1_X33)
      (sa, ia, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.34
  //* X(m,o3,o4,a)  <--  (    1.00000000)  T2(o2,o1,m,a) V2(o1,o3,o2,o4) 
  //* S2(i,k,m,a)  <--  (   -1.00000000) D2(i,o4,k,o3) X(m,o3,o4,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    T2b = T2.get_amp2(ia);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sa, X);
    for(int so1 = 0;so1 < nir;++so1){ 
    for(int io1 = symblockinfo.psym()(so1,I_O,I_BEGIN);io1 <= symblockinfo.psym()(so1,I_O,I_END);++io1){ 
      // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
      V2 <<= 0.0;                                                                          
      shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(io1);
      for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
        // Load a signle record of integals                                                             
        const int &imo2 = loadbuf_ptr->i0;                                                              
        const int &imo3 = loadbuf_ptr->i1;                                                              
        const int &imo4 = loadbuf_ptr->i2;                                                              
        const double &v = loadbuf_ptr->v;                                                               
        V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
      }                                                                                                 
      const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, io1, so1, V2); // V2=(IR-COV index) 
      FC_FUNC(g_if_sigma_ooov_ooov_no0_x34, G_IF_SIGMA_OOOV_OOOV_NO0_X34)
        (sa, ia, so1, io1, T2b.cptr(), V2_sym.cptr(), Xaaa.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x34, G_IF_SIGMA_OOOV_OOOV_NO1_X34)
      (sa, ia, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.35
  //* X(o2,o3,o5,o1,o4,a)  <--  (    1.00000000)  T2(o2,o3,o5,v1) V2(a,v1,o1,o4) 
  //* S2(i,k,m,a)  <--  (    1.00000000) D4(i,m,k,o3,o4,o1,o5,o2) X(o2,o3,o5,o1,o4,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
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
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ooov_ooov_no0_x35, G_IF_SIGMA_OOOV_OOOV_NO0_X35)
        (sa, ia, sv1, iv1, T2b.cptr(), V2_sym.cptr(), Xaaaaa.cptr(), nir, nsym, psym);
    }
    }
    for(int si = 0;si < nir;++si){ 
    for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
      for(int sm = 0;sm < nir;++sm){ 
      for(int im = symblockinfo.psym()(sm,I_O,I_BEGIN);im <= symblockinfo.psym()(sm,I_O,I_END);++im){ 
        // Load D4 from disk, or GA ....                                                     
        int imoi = amo2imo[ii] - nclosed;                              
        int imoj = amo2imo[im] - nclosed;                              
                                                                                             
        orz::DTensor rdm4_ij_sliced = rdm4(orz::Slice(imoi,imoi+1), orz::Slice(imoj,imoj+1), 
                                           orz::Slice(),            orz::Slice(),            
                                           orz::Slice(),            orz::Slice(),            
                                           orz::Slice(),            orz::Slice()).copy();    
        rdm4_sym = orz::ct::sympack_rdm4_2(symblockinfo, ii, si, im, sm, rdm4_ij_sliced);    
        FC_FUNC(g_if_set_d4,G_IF_SET_D4)(si, sm, ii, im, rdm4_sym.cptr(), nir, nsym, psym);  
        FC_FUNC(g_if_sigma_ooov_ooov_no1_x35, G_IF_SIGMA_OOOV_OOOV_NO1_X35)
          (sa, ia, si, ii, sm, im, Xaaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
        FC_FUNC(g_if_unset_d4,G_IF_UNSET_D4)();
      }
      }
    }
    }
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.36
  //* X(o3,o2,m,o1,o4,a)  <--  (    1.00000000)  T2(o3,o2,m,v1) V2(a,v1,o1,o4) 
  //* S2(i,k,m,a)  <--  (    1.00000000) D3(i,o3,k,o2,o4,o1) X(o3,o2,m,o1,o4,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
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
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ooov_ooov_no0_x36, G_IF_SIGMA_OOOV_OOOV_NO0_X36)
        (sa, ia, sv1, iv1, T2b.cptr(), V2_sym.cptr(), Xaaaaa.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x36, G_IF_SIGMA_OOOV_OOOV_NO1_X36)
      (sa, ia, Xaaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.37
  //* X(o1,o2,o4,m,o3,a)  <--  (    1.00000000)  T2(o1,o2,o4,v1) V2(a,v1,m,o3) 
  //* S2(i,k,m,a)  <--  (    1.00000000) D3(i,o3,k,o2,o4,o1) X(o1,o2,o4,m,o3,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
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
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ooov_ooov_no0_x37, G_IF_SIGMA_OOOV_OOOV_NO0_X37)
        (sa, ia, sv1, iv1, T2b.cptr(), V2_sym.cptr(), Xaaaaa.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x37, G_IF_SIGMA_OOOV_OOOV_NO1_X37)
      (sa, ia, Xaaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.38
  //* X(o2,o3,o4,a)  <--  (    1.00000000)  T2(o2,o3,o1,v1) V2(a,v1,o1,o4) 
  //* S2(i,k,m,a)  <--  (    1.00000000) D3(i,m,k,o3,o4,o2) X(o2,o3,o4,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
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
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sa, X);
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ooov_ooov_no0_x38, G_IF_SIGMA_OOOV_OOOV_NO0_X38)
        (sa, ia, sv1, iv1, T2b.cptr(), V2_sym.cptr(), Xaaa.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x38, G_IF_SIGMA_OOOV_OOOV_NO1_X38)
      (sa, ia, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.39
  //* X(o3,o2,m,a)  <--  (    1.00000000)  T2(o3,o2,o1,v1) V2(a,v1,m,o1) 
  //* S2(i,k,m,a)  <--  (    1.00000000) D2(i,o3,k,o2) X(o3,o2,m,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
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
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sa, X);
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ooov_ooov_no0_x39, G_IF_SIGMA_OOOV_OOOV_NO0_X39)
        (sa, ia, sv1, iv1, T2b.cptr(), V2_sym.cptr(), Xaaa.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x39, G_IF_SIGMA_OOOV_OOOV_NO1_X39)
      (sa, ia, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.40
  //* X(o4,o5,o1,o3,o2,a)  <--  (    1.00000000)  T2(o4,o5,o1,v1) V2(a,o3,o2,v1) 
  //* S2(i,k,m,a)  <--  (    1.00000000) D4(i,m,k,o3,o1,o4,o2,o5) X(o4,o5,o1,o3,o2,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
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
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ooov_ooov_no0_x40, G_IF_SIGMA_OOOV_OOOV_NO0_X40)
        (sa, ia, sv1, iv1, T2b.cptr(), V2_sym.cptr(), Xaaaaa.cptr(), nir, nsym, psym);
    }
    }
    for(int si = 0;si < nir;++si){ 
    for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
      for(int sm = 0;sm < nir;++sm){ 
      for(int im = symblockinfo.psym()(sm,I_O,I_BEGIN);im <= symblockinfo.psym()(sm,I_O,I_END);++im){ 
        // Load D4 from disk, or GA ....                                                     
        int imoi = amo2imo[ii] - nclosed;                              
        int imoj = amo2imo[im] - nclosed;                              
                                                                                             
        orz::DTensor rdm4_ij_sliced = rdm4(orz::Slice(imoi,imoi+1), orz::Slice(imoj,imoj+1), 
                                           orz::Slice(),            orz::Slice(),            
                                           orz::Slice(),            orz::Slice(),            
                                           orz::Slice(),            orz::Slice()).copy();    
        rdm4_sym = orz::ct::sympack_rdm4_2(symblockinfo, ii, si, im, sm, rdm4_ij_sliced);    
        FC_FUNC(g_if_set_d4,G_IF_SET_D4)(si, sm, ii, im, rdm4_sym.cptr(), nir, nsym, psym);  
        FC_FUNC(g_if_sigma_ooov_ooov_no1_x40, G_IF_SIGMA_OOOV_OOOV_NO1_X40)
          (sa, ia, si, ii, sm, im, Xaaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
        FC_FUNC(g_if_unset_d4,G_IF_UNSET_D4)();
      }
      }
    }
    }
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.41
  //* X(o3,o1,m,o2,o4,a)  <--  (    1.00000000)  T2(o3,o1,m,v1) V2(a,o2,o4,v1) 
  //* S2(i,k,m,a)  <--  (    1.00000000) D3(i,o3,k,o2,o4,o1) X(o3,o1,m,o2,o4,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
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
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ooov_ooov_no0_x41, G_IF_SIGMA_OOOV_OOOV_NO0_X41)
        (sa, ia, sv1, iv1, T2b.cptr(), V2_sym.cptr(), Xaaaaa.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x41, G_IF_SIGMA_OOOV_OOOV_NO1_X41)
      (sa, ia, Xaaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.42
  //* X(o1,o3,o4,o2,m,a)  <--  (    1.00000000)  T2(o1,o3,o4,v1) V2(a,o2,m,v1) 
  //* S2(i,k,m,a)  <--  (    1.00000000) D3(i,o3,k,o2,o4,o1) X(o1,o3,o4,o2,m,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
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
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ooov_ooov_no0_x42, G_IF_SIGMA_OOOV_OOOV_NO0_X42)
        (sa, ia, sv1, iv1, T2b.cptr(), V2_sym.cptr(), Xaaaaa.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x42, G_IF_SIGMA_OOOV_OOOV_NO1_X42)
      (sa, ia, Xaaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.43
  //* X(o2,o3,o1,a)  <--  (    1.00000000)  T2(o2,o3,o4,v1) V2(a,o4,o1,v1) 
  //* S2(i,k,m,a)  <--  (    1.00000000) D3(i,m,k,o2,o1,o3) X(o2,o3,o1,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
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
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sa, X);
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ooov_ooov_no0_x43, G_IF_SIGMA_OOOV_OOOV_NO0_X43)
        (sa, ia, sv1, iv1, T2b.cptr(), V2_sym.cptr(), Xaaa.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x43, G_IF_SIGMA_OOOV_OOOV_NO1_X43)
      (sa, ia, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.44
  //* X(o2,o1,m,a)  <--  (    1.00000000)  T2(o2,o1,o3,v1) V2(a,o3,m,v1) 
  //* S2(i,k,m,a)  <--  (    1.00000000) D2(i,o1,k,o2) X(o2,o1,m,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
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
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sa, X);
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ooov_ooov_no0_x44, G_IF_SIGMA_OOOV_OOOV_NO0_X44)
        (sa, ia, sv1, iv1, T2b.cptr(), V2_sym.cptr(), Xaaa.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x44, G_IF_SIGMA_OOOV_OOOV_NO1_X44)
      (sa, ia, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.45
  //* S2(i,k,m,a)  <--  (    1.00000000) Ecas D3(i,m,k,o2,o1,o3) T2(o3,o2,o1,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    T2b = T2.get_amp2(ia);
    FC_FUNC(g_if_sigma_ooov_ooov_no0_x45, G_IF_SIGMA_OOOV_OOOV_NO0_X45)
      (sa, ia, &Ecas, T2b.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.46
  //* S2(i,k,m,a)  <--  (    1.00000000) Ecas D2(i,o1,k,o2) T2(o1,o2,m,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    T2b = T2.get_amp2(ia);
    FC_FUNC(g_if_sigma_ooov_ooov_no0_x46, G_IF_SIGMA_OOOV_OOOV_NO0_X46)
      (sa, ia, &Ecas, T2b.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }

  return retval; 
} 
