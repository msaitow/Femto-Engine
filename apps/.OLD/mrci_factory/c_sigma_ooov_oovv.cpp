                                                                                
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
#include <sci/ctnew2/c_sigma_ooov_oovv.h>                                            
                                                                                
using std::cout;                                                                
using std::endl;                                                                
                                                                                
#define FLOPCOUNT                                                               
                                                                                
//       #                #########      #     #   # 
//  ########## ##########         #   #######  #   # 
//      #    #         #          #    # #     #   # 
//      #    #        #   ########     # #     #   # 
//     #     #     # #           #  ##########    #  
//    #   # #       #            #       #       #   
//   #     #         #    ########       #     ##    

                                                                                
// ***************************************************************************  
// orz::ct::mrci                                                                
// ***************************************************************************  
									    /*!        
   @brief CT input                                                              
									     */        
                                                                                
orz::ct::BareAmpPack orz::ct::sigma_ooov_oovv(const orz::ct::Input &ctinp,                                    
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
  //* X(o1,o2,o3,a)  <--  (    1.00000000)  T2(o1,o2,v1,a) h(o3,v1) 
  //* S2(i,k,m,a)  <--  (    1.00000000) D3(i,m,k,o2,o3,o1) X(o1,o2,o3,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    T2b = T2.get_amp2(ia);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sa, X);
    FC_FUNC(g_if_sigma_ooov_oovv_no0_x0, G_IF_SIGMA_OOOV_OOVV_NO0_X0)
      (sa, ia, T2b.cptr(), moint1_sym.cptr(), Xaaa.cptr(), nir, nsym, psym);
    FC_FUNC(g_if_sigma_ooov_oovv_no1_x0, G_IF_SIGMA_OOOV_OOVV_NO1_X0)
      (sa, ia, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.1
  //* X(o2,o1,m,a)  <--  (    1.00000000)  T2(o2,o1,v1,a) h(m,v1) 
  //* S2(i,k,m,a)  <--  (    1.00000000) D2(i,o2,k,o1) X(o2,o1,m,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    T2b = T2.get_amp2(ia);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sa, X);
    FC_FUNC(g_if_sigma_ooov_oovv_no0_x1, G_IF_SIGMA_OOOV_OOVV_NO0_X1)
      (sa, ia, T2b.cptr(), moint1_sym.cptr(), Xaaa.cptr(), nir, nsym, psym);
    FC_FUNC(g_if_sigma_ooov_oovv_no1_x1, G_IF_SIGMA_OOOV_OOVV_NO1_X1)
      (sa, ia, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.2
  //* X(o2,o3,o1,a)  <--  (    1.00000000)  T2(o2,o3,a,v1) h(o1,v1) 
  //* S2(i,k,m,a)  <--  (    1.00000000) D3(i,m,k,o2,o1,o3) X(o2,o3,o1,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sa, X);
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ooov_oovv_no0_x2, G_IF_SIGMA_OOOV_OOVV_NO0_X2)
        (sa, ia, sv1, iv1, T2b.cptr(), moint1_sym.cptr(), Xaaa.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_ooov_oovv_no1_x2, G_IF_SIGMA_OOOV_OOVV_NO1_X2)
      (sa, ia, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.3
  //* X(o2,o1,m,a)  <--  (    1.00000000)  T2(o2,o1,a,v1) h(m,v1) 
  //* S2(i,k,m,a)  <--  (    1.00000000) D2(i,o1,k,o2) X(o2,o1,m,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sa, X);
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ooov_oovv_no0_x3, G_IF_SIGMA_OOOV_OOVV_NO0_X3)
        (sa, ia, sv1, iv1, T2b.cptr(), moint1_sym.cptr(), Xaaa.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_ooov_oovv_no1_x3, G_IF_SIGMA_OOOV_OOVV_NO1_X3)
      (sa, ia, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.4
  //* X(o1,o2,o3,a)  <--  (    1.00000000)  T2(o1,o2,v1,a) Y0(o3,v1) 
  //* S2(i,k,m,a)  <--  (    2.00000000) D3(i,m,k,o2,o3,o1) X(o1,o2,o3,a) 
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
    FC_FUNC(g_if_sigma_ooov_oovv_y4, G_IF_SIGMA_OOOV_OOVV_Y4)
      (sc1, ic1, V2_sym.cptr(), Y0.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    T2b = T2.get_amp2(ia);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sa, X);
    FC_FUNC(g_if_sigma_ooov_oovv_no0_x4, G_IF_SIGMA_OOOV_OOVV_NO0_X4)
      (sa, ia, T2b.cptr(), Y0.cptr(), Xaaa.cptr(), nir, nsym, psym);
    FC_FUNC(g_if_sigma_ooov_oovv_no1_x4, G_IF_SIGMA_OOOV_OOVV_NO1_X4)
      (sa, ia, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.5
  //* X(o1,o2,o3,a)  <--  (    1.00000000)  T2(o1,o2,v1,a) Y1(o3,v1) 
  //* S2(i,k,m,a)  <--  (   -1.00000000) D3(i,m,k,o2,o3,o1) X(o1,o2,o3,a) 
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
    FC_FUNC(g_if_sigma_ooov_oovv_y5, G_IF_SIGMA_OOOV_OOVV_Y5)
      (sc1, ic1, V2_sym.cptr(), Y1.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    T2b = T2.get_amp2(ia);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sa, X);
    FC_FUNC(g_if_sigma_ooov_oovv_no0_x5, G_IF_SIGMA_OOOV_OOVV_NO0_X5)
      (sa, ia, T2b.cptr(), Y1.cptr(), Xaaa.cptr(), nir, nsym, psym);
    FC_FUNC(g_if_sigma_ooov_oovv_no1_x5, G_IF_SIGMA_OOOV_OOVV_NO1_X5)
      (sa, ia, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.6
  //* X(o2,o1,m,a)  <--  (    1.00000000)  T2(o2,o1,v1,a) Y2(m,v1) 
  //* S2(i,k,m,a)  <--  (    2.00000000) D2(i,o2,k,o1) X(o2,o1,m,a) 
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
    FC_FUNC(g_if_sigma_ooov_oovv_y6, G_IF_SIGMA_OOOV_OOVV_Y6)
      (sc1, ic1, V2_sym.cptr(), Y2.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    T2b = T2.get_amp2(ia);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sa, X);
    FC_FUNC(g_if_sigma_ooov_oovv_no0_x6, G_IF_SIGMA_OOOV_OOVV_NO0_X6)
      (sa, ia, T2b.cptr(), Y2.cptr(), Xaaa.cptr(), nir, nsym, psym);
    FC_FUNC(g_if_sigma_ooov_oovv_no1_x6, G_IF_SIGMA_OOOV_OOVV_NO1_X6)
      (sa, ia, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.7
  //* X(o2,o1,m,a)  <--  (    1.00000000)  T2(o2,o1,v1,a) Y3(m,v1) 
  //* S2(i,k,m,a)  <--  (   -1.00000000) D2(i,o2,k,o1) X(o2,o1,m,a) 
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
    FC_FUNC(g_if_sigma_ooov_oovv_y7, G_IF_SIGMA_OOOV_OOVV_Y7)
      (sc1, ic1, V2_sym.cptr(), Y3.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    T2b = T2.get_amp2(ia);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sa, X);
    FC_FUNC(g_if_sigma_ooov_oovv_no0_x7, G_IF_SIGMA_OOOV_OOVV_NO0_X7)
      (sa, ia, T2b.cptr(), Y3.cptr(), Xaaa.cptr(), nir, nsym, psym);
    FC_FUNC(g_if_sigma_ooov_oovv_no1_x7, G_IF_SIGMA_OOOV_OOVV_NO1_X7)
      (sa, ia, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.8
  //* X(o2,o1,o3,a)  <--  (    1.00000000)  T2(o2,o1,a,v1) Y4(o3,v1) 
  //* S2(i,k,m,a)  <--  (    2.00000000) D3(i,m,k,o2,o3,o1) X(o2,o1,o3,a) 
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
    FC_FUNC(g_if_sigma_ooov_oovv_y8, G_IF_SIGMA_OOOV_OOVV_Y8)
      (sc1, ic1, V2_sym.cptr(), Y4.cptr(), nir, nsym, psym);
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
      FC_FUNC(g_if_sigma_ooov_oovv_no0_x8, G_IF_SIGMA_OOOV_OOVV_NO0_X8)
        (sa, ia, sv1, iv1, T2b.cptr(), Y4.cptr(), Xaaa.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_ooov_oovv_no1_x8, G_IF_SIGMA_OOOV_OOVV_NO1_X8)
      (sa, ia, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.9
  //* X(o2,o1,o3,a)  <--  (    1.00000000)  T2(o2,o1,a,v1) Y5(o3,v1) 
  //* S2(i,k,m,a)  <--  (   -1.00000000) D3(i,m,k,o2,o3,o1) X(o2,o1,o3,a) 
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
    FC_FUNC(g_if_sigma_ooov_oovv_y9, G_IF_SIGMA_OOOV_OOVV_Y9)
      (sc1, ic1, V2_sym.cptr(), Y5.cptr(), nir, nsym, psym);
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
      FC_FUNC(g_if_sigma_ooov_oovv_no0_x9, G_IF_SIGMA_OOOV_OOVV_NO0_X9)
        (sa, ia, sv1, iv1, T2b.cptr(), Y5.cptr(), Xaaa.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_ooov_oovv_no1_x9, G_IF_SIGMA_OOOV_OOVV_NO1_X9)
      (sa, ia, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.10
  //* X(o1,o2,m,a)  <--  (    1.00000000)  T2(o1,o2,a,v1) Y6(m,v1) 
  //* S2(i,k,m,a)  <--  (    2.00000000) D2(i,o2,k,o1) X(o1,o2,m,a) 
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
    FC_FUNC(g_if_sigma_ooov_oovv_y10, G_IF_SIGMA_OOOV_OOVV_Y10)
      (sc1, ic1, V2_sym.cptr(), Y6.cptr(), nir, nsym, psym);
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
      FC_FUNC(g_if_sigma_ooov_oovv_no0_x10, G_IF_SIGMA_OOOV_OOVV_NO0_X10)
        (sa, ia, sv1, iv1, T2b.cptr(), Y6.cptr(), Xaaa.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_ooov_oovv_no1_x10, G_IF_SIGMA_OOOV_OOVV_NO1_X10)
      (sa, ia, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.11
  //* X(o1,o2,m,a)  <--  (    1.00000000)  T2(o1,o2,a,v1) Y7(m,v1) 
  //* S2(i,k,m,a)  <--  (   -1.00000000) D2(i,o2,k,o1) X(o1,o2,m,a) 
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
    FC_FUNC(g_if_sigma_ooov_oovv_y11, G_IF_SIGMA_OOOV_OOVV_Y11)
      (sc1, ic1, V2_sym.cptr(), Y7.cptr(), nir, nsym, psym);
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
      FC_FUNC(g_if_sigma_ooov_oovv_no0_x11, G_IF_SIGMA_OOOV_OOVV_NO0_X11)
        (sa, ia, sv1, iv1, T2b.cptr(), Y7.cptr(), Xaaa.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_ooov_oovv_no1_x11, G_IF_SIGMA_OOOV_OOVV_NO1_X11)
      (sa, ia, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.12
  //* X(o1,o3,o2,o5,o4,a)  <--  (    1.00000000)  T2(o1,o3,v1,a) V2(o2,o5,o4,v1) 
  //* S2(i,k,m,a)  <--  (    1.00000000) D4(o2,o5,m,i,o3,k,o1,o4) X(o1,o3,o2,o5,o4,a) 
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
    for(int so5 = 0;so5 < nir;++so5){ 
    for(int io5 = symblockinfo.psym()(so5,I_O,I_BEGIN);io5 <= symblockinfo.psym()(so5,I_O,I_END);++io5){ 
      for(int sa = 0;sa < nir;++sa){ 
      for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
        S2b = orz::DTensor(retval.namps_iamp()[ia]);
        T2b = T2.get_amp2(ia);
        orz::DTensor X(nocc, nocc, nocc);
        orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, so2^so5^sa, X);
        FC_FUNC(g_if_sigma_ooov_oovv_no0_x12, G_IF_SIGMA_OOOV_OOVV_NO0_X12)
          (sa, ia, so2, io2, so5, io5, T2b.cptr(), V2_sym.cptr(), Xaaa.cptr(), nir, nsym, psym);
        // Load D4 from disk, or GA ....                                                     
        int imoi = amo2imo[io2] - nclosed;                              
        int imoj = amo2imo[io5] - nclosed;                              
                                                                                             
        orz::DTensor rdm4_ij_sliced = rdm4(orz::Slice(imoi,imoi+1), orz::Slice(imoj,imoj+1), 
                                           orz::Slice(),            orz::Slice(),            
                                           orz::Slice(),            orz::Slice(),            
                                           orz::Slice(),            orz::Slice()).copy();    
        rdm4_sym = orz::ct::sympack_rdm4_2(symblockinfo, io2, so2, io5, so5, rdm4_ij_sliced);    
        FC_FUNC(g_if_set_d4,G_IF_SET_D4)(so2, so5, io2, io5, rdm4_sym.cptr(), nir, nsym, psym);  
        FC_FUNC(g_if_sigma_ooov_oovv_no1_x12, G_IF_SIGMA_OOOV_OOVV_NO1_X12)
          (sa, ia, so2, io2, so5, io5, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
        retval.acc_amp2(ia, S2b);
        FC_FUNC(g_if_unset_d4,G_IF_UNSET_D4)();
      }
      }
    }
    }
  }
  }
  }


  {
  // No.13
  //* X(o1,o2,o4,m,o3,a)  <--  (    1.00000000)  T2(o1,o2,v1,a) V2(v1,o4,m,o3) 
  //* S2(i,k,m,a)  <--  (    1.00000000) D3(i,o3,k,o2,o4,o1) X(o1,o2,o4,m,o3,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    T2b = T2.get_amp2(ia);
    orz::DTensor X(nocc, nocc, nocc, nocc, nocc);
    orz::DTensor Xaaaaa = orz::ct::sympack_Xaaaaa(symblockinfo, sa, X);
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
      V2 <<= 0.0;                                                                          
      shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iv1);
      for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
        // Load a signle record of integals                                                             
        const int &imo2 = loadbuf_ptr->i0;                                                              
        const int &imo3 = loadbuf_ptr->i1;                                                              
        const int &imo4 = loadbuf_ptr->i2;                                                              
        const double &v = loadbuf_ptr->v;                                                               
        V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
      }                                                                                                 
      const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iv1, sv1, V2); // V2=(IR-COV index) 
      FC_FUNC(g_if_sigma_ooov_oovv_no0_x13, G_IF_SIGMA_OOOV_OOVV_NO0_X13)
        (sa, ia, sv1, iv1, T2b.cptr(), V2_sym.cptr(), Xaaaaa.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_ooov_oovv_no1_x13, G_IF_SIGMA_OOOV_OOVV_NO1_X13)
      (sa, ia, Xaaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.14
  //* X(o3,o2,m,o1,o4,a)  <--  (    1.00000000)  T2(o3,o2,v1,a) V2(v1,m,o1,o4) 
  //* S2(i,k,m,a)  <--  (    1.00000000) D3(i,o3,k,o2,o4,o1) X(o3,o2,m,o1,o4,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    T2b = T2.get_amp2(ia);
    orz::DTensor X(nocc, nocc, nocc, nocc, nocc);
    orz::DTensor Xaaaaa = orz::ct::sympack_Xaaaaa(symblockinfo, sa, X);
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
      V2 <<= 0.0;                                                                          
      shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iv1);
      for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
        // Load a signle record of integals                                                             
        const int &imo2 = loadbuf_ptr->i0;                                                              
        const int &imo3 = loadbuf_ptr->i1;                                                              
        const int &imo4 = loadbuf_ptr->i2;                                                              
        const double &v = loadbuf_ptr->v;                                                               
        V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
      }                                                                                                 
      const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iv1, sv1, V2); // V2=(IR-COV index) 
      FC_FUNC(g_if_sigma_ooov_oovv_no0_x14, G_IF_SIGMA_OOOV_OOVV_NO0_X14)
        (sa, ia, sv1, iv1, T2b.cptr(), V2_sym.cptr(), Xaaaaa.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_ooov_oovv_no1_x14, G_IF_SIGMA_OOOV_OOVV_NO1_X14)
      (sa, ia, Xaaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.15
  //* X(o3,o5,o2,o1,o4,a)  <--  (    1.00000000)  T2(o3,o5,a,v1) V2(v1,o2,o1,o4) 
  //* S2(i,k,m,a)  <--  (    1.00000000) D4(i,m,k,o3,o1,o4,o2,o5) X(o3,o5,o2,o1,o4,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor X(nocc, nocc, nocc, nocc, nocc);
    orz::DTensor Xaaaaa = orz::ct::sympack_Xaaaaa(symblockinfo, sa, X);
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
      V2 <<= 0.0;                                                                          
      shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iv1);
      for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
        // Load a signle record of integals                                                             
        const int &imo2 = loadbuf_ptr->i0;                                                              
        const int &imo3 = loadbuf_ptr->i1;                                                              
        const int &imo4 = loadbuf_ptr->i2;                                                              
        const double &v = loadbuf_ptr->v;                                                               
        V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
      }                                                                                                 
      const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iv1, sv1, V2); // V2=(IR-COV index) 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ooov_oovv_no0_x15, G_IF_SIGMA_OOOV_OOVV_NO0_X15)
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
        FC_FUNC(g_if_sigma_ooov_oovv_no1_x15, G_IF_SIGMA_OOOV_OOVV_NO1_X15)
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
  // No.16
  //* X(o2,o1,o4,m,o3,a)  <--  (    1.00000000)  T2(o2,o1,a,v1) V2(v1,o4,m,o3) 
  //* S2(i,k,m,a)  <--  (    1.00000000) D3(i,o3,k,o2,o4,o1) X(o2,o1,o4,m,o3,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor X(nocc, nocc, nocc, nocc, nocc);
    orz::DTensor Xaaaaa = orz::ct::sympack_Xaaaaa(symblockinfo, sa, X);
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
      V2 <<= 0.0;                                                                          
      shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iv1);
      for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
        // Load a signle record of integals                                                             
        const int &imo2 = loadbuf_ptr->i0;                                                              
        const int &imo3 = loadbuf_ptr->i1;                                                              
        const int &imo4 = loadbuf_ptr->i2;                                                              
        const double &v = loadbuf_ptr->v;                                                               
        V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
      }                                                                                                 
      const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iv1, sv1, V2); // V2=(IR-COV index) 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ooov_oovv_no0_x16, G_IF_SIGMA_OOOV_OOVV_NO0_X16)
        (sa, ia, sv1, iv1, T2b.cptr(), V2_sym.cptr(), Xaaaaa.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_ooov_oovv_no1_x16, G_IF_SIGMA_OOOV_OOVV_NO1_X16)
      (sa, ia, Xaaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.17
  //* X(o2,o3,m,o1,o4,a)  <--  (    1.00000000)  T2(o2,o3,a,v1) V2(v1,m,o1,o4) 
  //* S2(i,k,m,a)  <--  (    1.00000000) D3(i,o3,k,o2,o4,o1) X(o2,o3,m,o1,o4,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor X(nocc, nocc, nocc, nocc, nocc);
    orz::DTensor Xaaaaa = orz::ct::sympack_Xaaaaa(symblockinfo, sa, X);
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
      V2 <<= 0.0;                                                                          
      shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iv1);
      for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
        // Load a signle record of integals                                                             
        const int &imo2 = loadbuf_ptr->i0;                                                              
        const int &imo3 = loadbuf_ptr->i1;                                                              
        const int &imo4 = loadbuf_ptr->i2;                                                              
        const double &v = loadbuf_ptr->v;                                                               
        V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
      }                                                                                                 
      const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iv1, sv1, V2); // V2=(IR-COV index) 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ooov_oovv_no0_x17, G_IF_SIGMA_OOOV_OOVV_NO0_X17)
        (sa, ia, sv1, iv1, T2b.cptr(), V2_sym.cptr(), Xaaaaa.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_ooov_oovv_no1_x17, G_IF_SIGMA_OOOV_OOVV_NO1_X17)
      (sa, ia, Xaaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.18
  //* X(o1,o2,o3,a)  <--  (    1.00000000)  T2(o1,o2,v2,v1) V2(a,v1,o3,v2) 
  //* S2(i,k,m,a)  <--  (    1.00000000) D3(i,m,k,o2,o3,o1) X(o1,o2,o3,a) 
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
      FC_FUNC(g_if_sigma_ooov_oovv_no0_x18, G_IF_SIGMA_OOOV_OOVV_NO0_X18)
        (sa, ia, sv1, iv1, T2b.cptr(), V2_sym.cptr(), Xaaa.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_ooov_oovv_no1_x18, G_IF_SIGMA_OOOV_OOVV_NO1_X18)
      (sa, ia, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.19
  //* X(o2,o1,m,a)  <--  (    1.00000000)  T2(o2,o1,v2,v1) V2(a,v1,m,v2) 
  //* S2(i,k,m,a)  <--  (    1.00000000) D2(i,o2,k,o1) X(o2,o1,m,a) 
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
      FC_FUNC(g_if_sigma_ooov_oovv_no0_x19, G_IF_SIGMA_OOOV_OOVV_NO0_X19)
        (sa, ia, sv1, iv1, T2b.cptr(), V2_sym.cptr(), Xaaa.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_ooov_oovv_no1_x19, G_IF_SIGMA_OOOV_OOVV_NO1_X19)
      (sa, ia, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.20
  //* X(o2,o3,o1,a)  <--  (    1.00000000)  T2(o2,o3,v1,v2) V2(a,v1,o1,v2) 
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
    for(int sv2 = 0;sv2 < nir;++sv2){ 
    for(int iv2 = symblockinfo.psym()(sv2,I_V,I_BEGIN);iv2 <= symblockinfo.psym()(sv2,I_V,I_END);++iv2){ 
      T2b = T2.get_amp2(iv2);
      FC_FUNC(g_if_sigma_ooov_oovv_no0_x20, G_IF_SIGMA_OOOV_OOVV_NO0_X20)
        (sa, ia, sv2, iv2, T2b.cptr(), V2_sym.cptr(), Xaaa.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_ooov_oovv_no1_x20, G_IF_SIGMA_OOOV_OOVV_NO1_X20)
      (sa, ia, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.21
  //* X(o2,o1,m,a)  <--  (    1.00000000)  T2(o2,o1,v1,v2) V2(a,v1,m,v2) 
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
    for(int sv2 = 0;sv2 < nir;++sv2){ 
    for(int iv2 = symblockinfo.psym()(sv2,I_V,I_BEGIN);iv2 <= symblockinfo.psym()(sv2,I_V,I_END);++iv2){ 
      T2b = T2.get_amp2(iv2);
      FC_FUNC(g_if_sigma_ooov_oovv_no0_x21, G_IF_SIGMA_OOOV_OOVV_NO0_X21)
        (sa, ia, sv2, iv2, T2b.cptr(), V2_sym.cptr(), Xaaa.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_ooov_oovv_no1_x21, G_IF_SIGMA_OOOV_OOVV_NO1_X21)
      (sa, ia, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }

  return retval; 
} 
