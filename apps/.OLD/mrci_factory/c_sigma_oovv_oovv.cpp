                                                                                
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
#include <sci/ctnew2/c_sigma_oovv_oovv.h>                                            
                                                                                
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
                                                                                
orz::ct::BareAmpPack orz::ct::sigma_oovv_oovv(const orz::ct::Input &ctinp,                                    
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
  //* S2(i,k,a,c)  <--  (    2.00000000) Y0 D2(i,o2,k,o1) T2(o2,o1,a,c) 
  // The effective tensor is detected .... 
  double Y0 = 0;
  FC_FUNC(g_if_sigma_oovv_oovv_y0, G_IF_SIGMA_OOVV_OOVV_Y0)
    (moint1_sym.cptr(), &Y0, nir, nsym, psym);
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    FC_FUNC(g_if_sigma_oovv_oovv_no0_x0, G_IF_SIGMA_OOVV_OOVV_NO0_X0)
      (sc, ic, &Y0, T2b.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.1
  //* S2(i,k,a,c)  <--  (    2.00000000) Y1 D2(i,o2,k,o1) T2(o2,o1,a,c) 
  // The effective tensor is detected .... 
  double Y1 = 0;
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
    FC_FUNC(g_if_sigma_oovv_oovv_y1, G_IF_SIGMA_OOVV_OOVV_Y1)
      (sc1, ic1, V2_sym.cptr(), &Y1, nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    FC_FUNC(g_if_sigma_oovv_oovv_no0_x1, G_IF_SIGMA_OOVV_OOVV_NO0_X1)
      (sc, ic, &Y1, T2b.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.2
  //* S2(i,k,a,c)  <--  (   -1.00000000) Y2 D2(i,o2,k,o1) T2(o2,o1,a,c) 
  // The effective tensor is detected .... 
  double Y2 = 0;
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
    FC_FUNC(g_if_sigma_oovv_oovv_y2, G_IF_SIGMA_OOVV_OOVV_Y2)
      (sc1, ic1, V2_sym.cptr(), &Y2, nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    FC_FUNC(g_if_sigma_oovv_oovv_no0_x2, G_IF_SIGMA_OOVV_OOVV_NO0_X2)
      (sc, ic, &Y2, T2b.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.3
  //* S2(i,k,a,c)  <--  (    2.00000000) Y3 D2(i,o2,k,o1) T2(o1,o2,c,a) 
  // The effective tensor is detected .... 
  double Y3 = 0;
  FC_FUNC(g_if_sigma_oovv_oovv_y3, G_IF_SIGMA_OOVV_OOVV_Y3)
    (moint1_sym.cptr(), &Y3, nir, nsym, psym);
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      T2b = T2.get_amp2(ia);
      FC_FUNC(g_if_sigma_oovv_oovv_no0_x3, G_IF_SIGMA_OOVV_OOVV_NO0_X3)
        (sa, ia, sc, ic, &Y3, T2b.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.4
  //* S2(i,k,a,c)  <--  (    2.00000000) Y4 D2(i,o2,k,o1) T2(o1,o2,c,a) 
  // The effective tensor is detected .... 
  double Y4 = 0;
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
    FC_FUNC(g_if_sigma_oovv_oovv_y4, G_IF_SIGMA_OOVV_OOVV_Y4)
      (sc1, ic1, V2_sym.cptr(), &Y4, nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      T2b = T2.get_amp2(ia);
      FC_FUNC(g_if_sigma_oovv_oovv_no0_x4, G_IF_SIGMA_OOVV_OOVV_NO0_X4)
        (sa, ia, sc, ic, &Y4, T2b.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.5
  //* S2(i,k,a,c)  <--  (   -1.00000000) Y5 D2(i,o2,k,o1) T2(o1,o2,c,a) 
  // The effective tensor is detected .... 
  double Y5 = 0;
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
    FC_FUNC(g_if_sigma_oovv_oovv_y5, G_IF_SIGMA_OOVV_OOVV_Y5)
      (sc1, ic1, V2_sym.cptr(), &Y5, nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      T2b = T2.get_amp2(ia);
      FC_FUNC(g_if_sigma_oovv_oovv_no0_x5, G_IF_SIGMA_OOVV_OOVV_NO0_X5)
        (sa, ia, sc, ic, &Y5, T2b.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.6
  //* X(i,o3,k,o2)  <--  (    1.00000000)  D3(i,o3,k,o2,o4,o1) h(o4,o1) 
  //* S2(i,k,a,c)  <--  (    1.00000000) T2(o3,o2,a,c) X(i,o3,k,o2) 
  orz::DTensor X(nocc, nocc, nocc, nocc);
  orz::DTensor Xaaaa = orz::ct::sympack_Xaaaa(symblockinfo, 0, X);
  FC_FUNC(g_if_sigma_oovv_oovv_no0_x6, G_IF_SIGMA_OOVV_OOVV_NO0_X6)
    (moint1_sym.cptr(), Xaaaa.cptr(), nir, nsym, psym);
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    FC_FUNC(g_if_sigma_oovv_oovv_no1_x6, G_IF_SIGMA_OOVV_OOVV_NO1_X6)
      (sc, ic, T2b.cptr(), Xaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.7
  //* X(i,o2,k,o3)  <--  (    1.00000000)  D3(i,o2,k,o3,o1,o4) h(o1,o4) 
  //* S2(i,k,a,c)  <--  (    1.00000000) T2(o3,o2,c,a) X(i,o2,k,o3) 
  orz::DTensor X(nocc, nocc, nocc, nocc);
  orz::DTensor Xaaaa = orz::ct::sympack_Xaaaa(symblockinfo, 0, X);
  FC_FUNC(g_if_sigma_oovv_oovv_no0_x7, G_IF_SIGMA_OOVV_OOVV_NO0_X7)
    (moint1_sym.cptr(), Xaaaa.cptr(), nir, nsym, psym);
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      T2b = T2.get_amp2(ia);
      FC_FUNC(g_if_sigma_oovv_oovv_no1_x7, G_IF_SIGMA_OOVV_OOVV_NO1_X7)
        (sa, ia, sc, ic, T2b.cptr(), Xaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.8
  //* X(o1,o2,a,c)  <--  (    1.00000000)  T2(o1,o2,v1,a) h(c,v1) 
  //* S2(i,k,a,c)  <--  (    1.00000000) D2(i,o2,k,o1) X(o1,o2,a,c) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      orz::DTensor X(nocc, nocc);
      orz::DTensor Xaa = orz::ct::sympack_Xaa(symblockinfo, sa^sc, X);
      FC_FUNC(g_if_sigma_oovv_oovv_no0_x8, G_IF_SIGMA_OOVV_OOVV_NO0_X8)
        (sa, ia, sc, ic, T2b.cptr(), moint1_sym.cptr(), Xaa.cptr(), nir, nsym, psym);
      FC_FUNC(g_if_sigma_oovv_oovv_no1_x8, G_IF_SIGMA_OOVV_OOVV_NO1_X8)
        (sa, ia, sc, ic, Xaa.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.9
  //* X(o2,o1,c,a)  <--  (    1.00000000)  T2(o2,o1,v1,c) h(a,v1) 
  //* S2(i,k,a,c)  <--  (    1.00000000) D2(i,o2,k,o1) X(o2,o1,c,a) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    orz::DTensor X(nocc, nocc, nvir);
    orz::DTensor Xaav = orz::ct::sympack_Xaav(symblockinfo, sc, X);
    FC_FUNC(g_if_sigma_oovv_oovv_no0_x9, G_IF_SIGMA_OOVV_OOVV_NO0_X9)
      (sc, ic, T2b.cptr(), moint1_sym.cptr(), Xaav.cptr(), nir, nsym, psym);
    FC_FUNC(g_if_sigma_oovv_oovv_no1_x9, G_IF_SIGMA_OOVV_OOVV_NO1_X9)
      (sc, ic, Xaav.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.10
  //* X(o2,o1,a,c)  <--  (    1.00000000)  T2(o2,o1,a,v1) h(c,v1) 
  //* S2(i,k,a,c)  <--  (    1.00000000) D2(i,o2,k,o1) X(o2,o1,a,c) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    orz::DTensor X(nocc, nocc, nvir);
    orz::DTensor Xaav = orz::ct::sympack_Xaav(symblockinfo, sc, X);
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_oovv_oovv_no0_x10, G_IF_SIGMA_OOVV_OOVV_NO0_X10)
        (sc, ic, sv1, iv1, T2b.cptr(), moint1_sym.cptr(), Xaav.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_oovv_oovv_no1_x10, G_IF_SIGMA_OOVV_OOVV_NO1_X10)
      (sc, ic, Xaav.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.11
  //* X(o2,o1,c,a)  <--  (    1.00000000)  T2(o2,o1,c,v1) h(a,v1) 
  //* S2(i,k,a,c)  <--  (    1.00000000) D2(i,o1,k,o2) X(o2,o1,c,a) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    orz::DTensor X(nocc, nocc, nvir);
    orz::DTensor Xaav = orz::ct::sympack_Xaav(symblockinfo, sc, X);
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_oovv_oovv_no0_x11, G_IF_SIGMA_OOVV_OOVV_NO0_X11)
        (sc, ic, sv1, iv1, T2b.cptr(), moint1_sym.cptr(), Xaav.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_oovv_oovv_no1_x11, G_IF_SIGMA_OOVV_OOVV_NO1_X11)
      (sc, ic, Xaav.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.12
  //* X(i,o3,k,o2)  <--  (    1.00000000)  D3(i,o3,k,o2,o4,o1) Y6(o1,o4) 
  //* S2(i,k,a,c)  <--  (    2.00000000) T2(o3,o2,a,c) X(i,o3,k,o2) 
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
    FC_FUNC(g_if_sigma_oovv_oovv_y12, G_IF_SIGMA_OOVV_OOVV_Y12)
      (sc1, ic1, V2_sym.cptr(), Y6.cptr(), nir, nsym, psym);
  }
  }
  orz::DTensor X(nocc, nocc, nocc, nocc);
  orz::DTensor Xaaaa = orz::ct::sympack_Xaaaa(symblockinfo, 0, X);
  FC_FUNC(g_if_sigma_oovv_oovv_no0_x12, G_IF_SIGMA_OOVV_OOVV_NO0_X12)
    (Y6.cptr(), Xaaaa.cptr(), nir, nsym, psym);
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    FC_FUNC(g_if_sigma_oovv_oovv_no1_x12, G_IF_SIGMA_OOVV_OOVV_NO1_X12)
      (sc, ic, T2b.cptr(), Xaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.13
  //* X(i,o3,k,o2)  <--  (    1.00000000)  D3(i,o3,k,o2,o4,o1) Y7(o1,o4) 
  //* S2(i,k,a,c)  <--  (   -1.00000000) T2(o3,o2,a,c) X(i,o3,k,o2) 
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
    FC_FUNC(g_if_sigma_oovv_oovv_y13, G_IF_SIGMA_OOVV_OOVV_Y13)
      (sc1, ic1, V2_sym.cptr(), Y7.cptr(), nir, nsym, psym);
  }
  }
  orz::DTensor X(nocc, nocc, nocc, nocc);
  orz::DTensor Xaaaa = orz::ct::sympack_Xaaaa(symblockinfo, 0, X);
  FC_FUNC(g_if_sigma_oovv_oovv_no0_x13, G_IF_SIGMA_OOVV_OOVV_NO0_X13)
    (Y7.cptr(), Xaaaa.cptr(), nir, nsym, psym);
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    FC_FUNC(g_if_sigma_oovv_oovv_no1_x13, G_IF_SIGMA_OOVV_OOVV_NO1_X13)
      (sc, ic, T2b.cptr(), Xaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.14
  //* X(i,o3,k,o2)  <--  (    1.00000000)  D3(i,o3,k,o2,o4,o1) Y8(o1,o4) 
  //* S2(i,k,a,c)  <--  (    2.00000000) T2(o2,o3,c,a) X(i,o3,k,o2) 
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
    FC_FUNC(g_if_sigma_oovv_oovv_y14, G_IF_SIGMA_OOVV_OOVV_Y14)
      (sc1, ic1, V2_sym.cptr(), Y8.cptr(), nir, nsym, psym);
  }
  }
  orz::DTensor X(nocc, nocc, nocc, nocc);
  orz::DTensor Xaaaa = orz::ct::sympack_Xaaaa(symblockinfo, 0, X);
  FC_FUNC(g_if_sigma_oovv_oovv_no0_x14, G_IF_SIGMA_OOVV_OOVV_NO0_X14)
    (Y8.cptr(), Xaaaa.cptr(), nir, nsym, psym);
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      T2b = T2.get_amp2(ia);
      FC_FUNC(g_if_sigma_oovv_oovv_no1_x14, G_IF_SIGMA_OOVV_OOVV_NO1_X14)
        (sa, ia, sc, ic, T2b.cptr(), Xaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.15
  //* X(i,o3,k,o2)  <--  (    1.00000000)  D3(i,o3,k,o2,o4,o1) Y9(o1,o4) 
  //* S2(i,k,a,c)  <--  (   -1.00000000) T2(o2,o3,c,a) X(i,o3,k,o2) 
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
    FC_FUNC(g_if_sigma_oovv_oovv_y15, G_IF_SIGMA_OOVV_OOVV_Y15)
      (sc1, ic1, V2_sym.cptr(), Y9.cptr(), nir, nsym, psym);
  }
  }
  orz::DTensor X(nocc, nocc, nocc, nocc);
  orz::DTensor Xaaaa = orz::ct::sympack_Xaaaa(symblockinfo, 0, X);
  FC_FUNC(g_if_sigma_oovv_oovv_no0_x15, G_IF_SIGMA_OOVV_OOVV_NO0_X15)
    (Y9.cptr(), Xaaaa.cptr(), nir, nsym, psym);
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      T2b = T2.get_amp2(ia);
      FC_FUNC(g_if_sigma_oovv_oovv_no1_x15, G_IF_SIGMA_OOVV_OOVV_NO1_X15)
        (sa, ia, sc, ic, T2b.cptr(), Xaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.16
  //* X(o1,o2,a,c)  <--  (    1.00000000)  T2(o1,o2,v1,a) Y10(c,v1) 
  //* S2(i,k,a,c)  <--  (    2.00000000) D2(i,o2,k,o1) X(o1,o2,a,c) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nvir, nvir);
  orz::DTensor Y10 = orz::ct::sympack_Xvv(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_oovv_oovv_y16, G_IF_SIGMA_OOVV_OOVV_Y16)
      (sc1, ic1, V2_sym.cptr(), Y10.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      orz::DTensor X(nocc, nocc);
      orz::DTensor Xaa = orz::ct::sympack_Xaa(symblockinfo, sa^sc, X);
      FC_FUNC(g_if_sigma_oovv_oovv_no0_x16, G_IF_SIGMA_OOVV_OOVV_NO0_X16)
        (sa, ia, sc, ic, T2b.cptr(), Y10.cptr(), Xaa.cptr(), nir, nsym, psym);
      FC_FUNC(g_if_sigma_oovv_oovv_no1_x16, G_IF_SIGMA_OOVV_OOVV_NO1_X16)
        (sa, ia, sc, ic, Xaa.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.17
  //* X(o1,o2,a,c)  <--  (    1.00000000)  T2(o1,o2,v1,a) Y11(c,v1) 
  //* S2(i,k,a,c)  <--  (   -1.00000000) D2(i,o2,k,o1) X(o1,o2,a,c) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nvir, nvir);
  orz::DTensor Y11 = orz::ct::sympack_Xvv(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_oovv_oovv_y17, G_IF_SIGMA_OOVV_OOVV_Y17)
      (sc1, ic1, V2_sym.cptr(), Y11.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      orz::DTensor X(nocc, nocc);
      orz::DTensor Xaa = orz::ct::sympack_Xaa(symblockinfo, sa^sc, X);
      FC_FUNC(g_if_sigma_oovv_oovv_no0_x17, G_IF_SIGMA_OOVV_OOVV_NO0_X17)
        (sa, ia, sc, ic, T2b.cptr(), Y11.cptr(), Xaa.cptr(), nir, nsym, psym);
      FC_FUNC(g_if_sigma_oovv_oovv_no1_x17, G_IF_SIGMA_OOVV_OOVV_NO1_X17)
        (sa, ia, sc, ic, Xaa.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.18
  //* X(o2,o1,c,a)  <--  (    1.00000000)  T2(o2,o1,v1,c) Y12(a,v1) 
  //* S2(i,k,a,c)  <--  (    2.00000000) D2(i,o2,k,o1) X(o2,o1,c,a) 
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
    FC_FUNC(g_if_sigma_oovv_oovv_y18, G_IF_SIGMA_OOVV_OOVV_Y18)
      (sc1, ic1, V2_sym.cptr(), Y12.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    orz::DTensor X(nocc, nocc, nvir);
    orz::DTensor Xaav = orz::ct::sympack_Xaav(symblockinfo, sc, X);
    FC_FUNC(g_if_sigma_oovv_oovv_no0_x18, G_IF_SIGMA_OOVV_OOVV_NO0_X18)
      (sc, ic, T2b.cptr(), Y12.cptr(), Xaav.cptr(), nir, nsym, psym);
    FC_FUNC(g_if_sigma_oovv_oovv_no1_x18, G_IF_SIGMA_OOVV_OOVV_NO1_X18)
      (sc, ic, Xaav.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.19
  //* X(o2,o1,c,a)  <--  (    1.00000000)  T2(o2,o1,v1,c) Y13(a,v1) 
  //* S2(i,k,a,c)  <--  (   -1.00000000) D2(i,o2,k,o1) X(o2,o1,c,a) 
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
    FC_FUNC(g_if_sigma_oovv_oovv_y19, G_IF_SIGMA_OOVV_OOVV_Y19)
      (sc1, ic1, V2_sym.cptr(), Y13.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    orz::DTensor X(nocc, nocc, nvir);
    orz::DTensor Xaav = orz::ct::sympack_Xaav(symblockinfo, sc, X);
    FC_FUNC(g_if_sigma_oovv_oovv_no0_x19, G_IF_SIGMA_OOVV_OOVV_NO0_X19)
      (sc, ic, T2b.cptr(), Y13.cptr(), Xaav.cptr(), nir, nsym, psym);
    FC_FUNC(g_if_sigma_oovv_oovv_no1_x19, G_IF_SIGMA_OOVV_OOVV_NO1_X19)
      (sc, ic, Xaav.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.20
  //* X(o2,o1,a,c)  <--  (    1.00000000)  T2(o2,o1,a,v1) Y14(c,v1) 
  //* S2(i,k,a,c)  <--  (    2.00000000) D2(i,o2,k,o1) X(o2,o1,a,c) 
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
    FC_FUNC(g_if_sigma_oovv_oovv_y20, G_IF_SIGMA_OOVV_OOVV_Y20)
      (sc1, ic1, V2_sym.cptr(), Y14.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    orz::DTensor X(nocc, nocc, nvir);
    orz::DTensor Xaav = orz::ct::sympack_Xaav(symblockinfo, sc, X);
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_oovv_oovv_no0_x20, G_IF_SIGMA_OOVV_OOVV_NO0_X20)
        (sc, ic, sv1, iv1, T2b.cptr(), Y14.cptr(), Xaav.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_oovv_oovv_no1_x20, G_IF_SIGMA_OOVV_OOVV_NO1_X20)
      (sc, ic, Xaav.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.21
  //* X(o2,o1,a,c)  <--  (    1.00000000)  T2(o2,o1,a,v1) Y15(c,v1) 
  //* S2(i,k,a,c)  <--  (   -1.00000000) D2(i,o2,k,o1) X(o2,o1,a,c) 
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
    FC_FUNC(g_if_sigma_oovv_oovv_y21, G_IF_SIGMA_OOVV_OOVV_Y21)
      (sc1, ic1, V2_sym.cptr(), Y15.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    orz::DTensor X(nocc, nocc, nvir);
    orz::DTensor Xaav = orz::ct::sympack_Xaav(symblockinfo, sc, X);
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_oovv_oovv_no0_x21, G_IF_SIGMA_OOVV_OOVV_NO0_X21)
        (sc, ic, sv1, iv1, T2b.cptr(), Y15.cptr(), Xaav.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_oovv_oovv_no1_x21, G_IF_SIGMA_OOVV_OOVV_NO1_X21)
      (sc, ic, Xaav.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.22
  //* X(o1,o2,c,a)  <--  (    1.00000000)  T2(o1,o2,c,v1) Y16(a,v1) 
  //* S2(i,k,a,c)  <--  (    2.00000000) D2(i,o2,k,o1) X(o1,o2,c,a) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nvir, nvir);
  orz::DTensor Y16 = orz::ct::sympack_Xvv(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_oovv_oovv_y22, G_IF_SIGMA_OOVV_OOVV_Y22)
      (sc1, ic1, V2_sym.cptr(), Y16.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    orz::DTensor X(nocc, nocc, nvir);
    orz::DTensor Xaav = orz::ct::sympack_Xaav(symblockinfo, sc, X);
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_oovv_oovv_no0_x22, G_IF_SIGMA_OOVV_OOVV_NO0_X22)
        (sc, ic, sv1, iv1, T2b.cptr(), Y16.cptr(), Xaav.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_oovv_oovv_no1_x22, G_IF_SIGMA_OOVV_OOVV_NO1_X22)
      (sc, ic, Xaav.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.23
  //* X(o1,o2,c,a)  <--  (    1.00000000)  T2(o1,o2,c,v1) Y17(a,v1) 
  //* S2(i,k,a,c)  <--  (   -1.00000000) D2(i,o2,k,o1) X(o1,o2,c,a) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nvir, nvir);
  orz::DTensor Y17 = orz::ct::sympack_Xvv(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_oovv_oovv_y23, G_IF_SIGMA_OOVV_OOVV_Y23)
      (sc1, ic1, V2_sym.cptr(), Y17.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    orz::DTensor X(nocc, nocc, nvir);
    orz::DTensor Xaav = orz::ct::sympack_Xaav(symblockinfo, sc, X);
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_oovv_oovv_no0_x23, G_IF_SIGMA_OOVV_OOVV_NO0_X23)
        (sc, ic, sv1, iv1, T2b.cptr(), Y17.cptr(), Xaav.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_oovv_oovv_no1_x23, G_IF_SIGMA_OOVV_OOVV_NO1_X23)
      (sc, ic, Xaav.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.24
  //* X(o4,i,o3,k)  <--  (    1.00000000)  D4(o1,o5,o4,i,o3,k,o2,o6) V2(o1,o5,o2,o6) 
  //* S2(i,k,a,c)  <--  (    0.50000000) T2(o4,o3,a,c) X(o4,i,o3,k) 
  orz::DTensor X(nocc, nocc, nocc, nocc);
  orz::DTensor Xaaaa = orz::ct::sympack_Xaaaa(symblockinfo, 0, X);
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
    for(int so5 = 0;so5 < nir;++so5){ 
    for(int io5 = symblockinfo.psym()(so5,I_O,I_BEGIN);io5 <= symblockinfo.psym()(so5,I_O,I_END);++io5){ 
      // Load D4 from disk, or GA ....                                                     
      int imoi = amo2imo[io1] - nclosed;                              
      int imoj = amo2imo[io5] - nclosed;                              
                                                                                           
      orz::DTensor rdm4_ij_sliced = rdm4(orz::Slice(imoi,imoi+1), orz::Slice(imoj,imoj+1), 
                                         orz::Slice(),            orz::Slice(),            
                                         orz::Slice(),            orz::Slice(),            
                                         orz::Slice(),            orz::Slice()).copy();    
      rdm4_sym = orz::ct::sympack_rdm4_2(symblockinfo, io1, so1, io5, so5, rdm4_ij_sliced);    
      FC_FUNC(g_if_set_d4,G_IF_SET_D4)(so1, so5, io1, io5, rdm4_sym.cptr(), nir, nsym, psym);  
      FC_FUNC(g_if_sigma_oovv_oovv_no0_x24, G_IF_SIGMA_OOVV_OOVV_NO0_X24)
        (so1, io1, so5, io5, V2_sym.cptr(), Xaaaa.cptr(), nir, nsym, psym);
      FC_FUNC(g_if_unset_d4,G_IF_UNSET_D4)();
    }
    }
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    FC_FUNC(g_if_sigma_oovv_oovv_no1_x24, G_IF_SIGMA_OOVV_OOVV_NO1_X24)
      (sc, ic, T2b.cptr(), Xaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.25
  //* X(i,o3,k,o4)  <--  (    1.00000000)  D4(o1,o5,i,o3,k,o4,o2,o6) V2(o1,o5,o2,o6) 
  //* S2(i,k,a,c)  <--  (    0.50000000) T2(o4,o3,c,a) X(i,o3,k,o4) 
  orz::DTensor X(nocc, nocc, nocc, nocc);
  orz::DTensor Xaaaa = orz::ct::sympack_Xaaaa(symblockinfo, 0, X);
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
    for(int so5 = 0;so5 < nir;++so5){ 
    for(int io5 = symblockinfo.psym()(so5,I_O,I_BEGIN);io5 <= symblockinfo.psym()(so5,I_O,I_END);++io5){ 
      // Load D4 from disk, or GA ....                                                     
      int imoi = amo2imo[io1] - nclosed;                              
      int imoj = amo2imo[io5] - nclosed;                              
                                                                                           
      orz::DTensor rdm4_ij_sliced = rdm4(orz::Slice(imoi,imoi+1), orz::Slice(imoj,imoj+1), 
                                         orz::Slice(),            orz::Slice(),            
                                         orz::Slice(),            orz::Slice(),            
                                         orz::Slice(),            orz::Slice()).copy();    
      rdm4_sym = orz::ct::sympack_rdm4_2(symblockinfo, io1, so1, io5, so5, rdm4_ij_sliced);    
      FC_FUNC(g_if_set_d4,G_IF_SET_D4)(so1, so5, io1, io5, rdm4_sym.cptr(), nir, nsym, psym);  
      FC_FUNC(g_if_sigma_oovv_oovv_no0_x25, G_IF_SIGMA_OOVV_OOVV_NO0_X25)
        (so1, io1, so5, io5, V2_sym.cptr(), Xaaaa.cptr(), nir, nsym, psym);
      FC_FUNC(g_if_unset_d4,G_IF_UNSET_D4)();
    }
    }
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      T2b = T2.get_amp2(ia);
      FC_FUNC(g_if_sigma_oovv_oovv_no1_x25, G_IF_SIGMA_OOVV_OOVV_NO1_X25)
        (sa, ia, sc, ic, T2b.cptr(), Xaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.26
  //* X(o2,o3,o1,o4,a,c)  <--  (    1.00000000)  T2(o2,o3,v1,a) V2(c,v1,o1,o4) 
  //* S2(i,k,a,c)  <--  (    1.00000000) D3(i,o3,k,o2,o4,o1) X(o2,o3,o1,o4,a,c) 
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
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      T2b = T2.get_amp2(ia);
      orz::DTensor X(nocc, nocc, nocc, nocc);
      orz::DTensor Xaaaa = orz::ct::sympack_Xaaaa(symblockinfo, sa^sc, X);
      FC_FUNC(g_if_sigma_oovv_oovv_no0_x26, G_IF_SIGMA_OOVV_OOVV_NO0_X26)
        (sa, ia, sc, ic, T2b.cptr(), V2_sym.cptr(), Xaaaa.cptr(), nir, nsym, psym);
      FC_FUNC(g_if_sigma_oovv_oovv_no1_x26, G_IF_SIGMA_OOVV_OOVV_NO1_X26)
        (sa, ia, sc, ic, Xaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.27
  //* X(o3,o2,o1,o4,c,a)  <--  (    1.00000000)  T2(o3,o2,v1,c) V2(a,v1,o1,o4) 
  //* S2(i,k,a,c)  <--  (    1.00000000) D3(i,o3,k,o2,o4,o1) X(o3,o2,o1,o4,c,a) 
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
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      T2b = T2.get_amp2(ic);
      orz::DTensor X(nocc, nocc, nocc, nocc);
      orz::DTensor Xaaaa = orz::ct::sympack_Xaaaa(symblockinfo, sc^sa, X);
      FC_FUNC(g_if_sigma_oovv_oovv_no0_x27, G_IF_SIGMA_OOVV_OOVV_NO0_X27)
        (sa, ia, sc, ic, T2b.cptr(), V2_sym.cptr(), Xaaaa.cptr(), nir, nsym, psym);
      FC_FUNC(g_if_sigma_oovv_oovv_no1_x27, G_IF_SIGMA_OOVV_OOVV_NO1_X27)
        (sa, ia, sc, ic, Xaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.28
  //* X(o3,o2,o1,o4,a,c)  <--  (    1.00000000)  T2(o3,o2,a,v1) V2(c,v1,o1,o4) 
  //* S2(i,k,a,c)  <--  (    1.00000000) D3(i,o3,k,o2,o4,o1) X(o3,o2,o1,o4,a,c) 
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
      FC_FUNC(g_if_sigma_oovv_oovv_no0_x28, G_IF_SIGMA_OOVV_OOVV_NO0_X28)
        (sc, ic, sv1, iv1, T2b.cptr(), V2_sym.cptr(), Xaaaav.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_oovv_oovv_no1_x28, G_IF_SIGMA_OOVV_OOVV_NO1_X28)
      (sc, ic, Xaaaav.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.29
  //* X(o2,o3,o1,o4,c,a)  <--  (    1.00000000)  T2(o2,o3,c,v1) V2(v1,a,o1,o4) 
  //* S2(i,k,a,c)  <--  (    1.00000000) D3(i,o3,k,o2,o4,o1) X(o2,o3,o1,o4,c,a) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    orz::DTensor X(nocc, nocc, nocc, nocc, nvir);
    orz::DTensor Xaaaav = orz::ct::sympack_Xaaaav(symblockinfo, sc, X);
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
      FC_FUNC(g_if_sigma_oovv_oovv_no0_x29, G_IF_SIGMA_OOVV_OOVV_NO0_X29)
        (sc, ic, sv1, iv1, T2b.cptr(), V2_sym.cptr(), Xaaaav.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_oovv_oovv_no1_x29, G_IF_SIGMA_OOVV_OOVV_NO1_X29)
      (sc, ic, Xaaaav.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.30
  //* X(o1,o3,o2,o4,a,c)  <--  (    1.00000000)  T2(o1,o3,v1,a) V2(c,o2,o4,v1) 
  //* S2(i,k,a,c)  <--  (    1.00000000) D3(i,o3,k,o2,o4,o1) X(o1,o3,o2,o4,a,c) 
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
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      T2b = T2.get_amp2(ia);
      orz::DTensor X(nocc, nocc, nocc, nocc);
      orz::DTensor Xaaaa = orz::ct::sympack_Xaaaa(symblockinfo, sa^sc, X);
      FC_FUNC(g_if_sigma_oovv_oovv_no0_x30, G_IF_SIGMA_OOVV_OOVV_NO0_X30)
        (sa, ia, sc, ic, T2b.cptr(), V2_sym.cptr(), Xaaaa.cptr(), nir, nsym, psym);
      FC_FUNC(g_if_sigma_oovv_oovv_no1_x30, G_IF_SIGMA_OOVV_OOVV_NO1_X30)
        (sa, ia, sc, ic, Xaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.31
  //* X(o1,o2,o3,o4,c,a)  <--  (    1.00000000)  T2(o1,o2,v1,c) V2(a,o3,o4,v1) 
  //* S2(i,k,a,c)  <--  (    1.00000000) D3(i,o3,k,o2,o4,o1) X(o1,o2,o3,o4,c,a) 
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
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      T2b = T2.get_amp2(ic);
      orz::DTensor X(nocc, nocc, nocc, nocc);
      orz::DTensor Xaaaa = orz::ct::sympack_Xaaaa(symblockinfo, sc^sa, X);
      FC_FUNC(g_if_sigma_oovv_oovv_no0_x31, G_IF_SIGMA_OOVV_OOVV_NO0_X31)
        (sa, ia, sc, ic, T2b.cptr(), V2_sym.cptr(), Xaaaa.cptr(), nir, nsym, psym);
      FC_FUNC(g_if_sigma_oovv_oovv_no1_x31, G_IF_SIGMA_OOVV_OOVV_NO1_X31)
        (sa, ia, sc, ic, Xaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.32
  //* X(o3,o1,o2,o4,a,c)  <--  (    1.00000000)  T2(o3,o1,a,v1) V2(c,o2,o4,v1) 
  //* S2(i,k,a,c)  <--  (    1.00000000) D3(i,o3,k,o2,o4,o1) X(o3,o1,o2,o4,a,c) 
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
      FC_FUNC(g_if_sigma_oovv_oovv_no0_x32, G_IF_SIGMA_OOVV_OOVV_NO0_X32)
        (sc, ic, sv1, iv1, T2b.cptr(), V2_sym.cptr(), Xaaaav.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_oovv_oovv_no1_x32, G_IF_SIGMA_OOVV_OOVV_NO1_X32)
      (sc, ic, Xaaaav.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.33
  //* X(o3,o4,o1,o2,c,a)  <--  (    1.00000000)  T2(o3,o4,c,v1) V2(v1,o1,o2,a) 
  //* S2(i,k,a,c)  <--  (    1.00000000) D3(i,o2,k,o3,o1,o4) X(o3,o4,o1,o2,c,a) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    orz::DTensor X(nocc, nocc, nocc, nocc, nvir);
    orz::DTensor Xaaaav = orz::ct::sympack_Xaaaav(symblockinfo, sc, X);
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
      FC_FUNC(g_if_sigma_oovv_oovv_no0_x33, G_IF_SIGMA_OOVV_OOVV_NO0_X33)
        (sc, ic, sv1, iv1, T2b.cptr(), V2_sym.cptr(), Xaaaav.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_oovv_oovv_no1_x33, G_IF_SIGMA_OOVV_OOVV_NO1_X33)
      (sc, ic, Xaaaav.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.34
  //* X(o2,o1,c,a)  <--  (    1.00000000)  T2(o2,o1,v2,v1) V2(c,v1,a,v2) 
  //* S2(i,k,a,c)  <--  (    1.00000000) D2(i,o2,k,o1) X(o2,o1,c,a) 
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
      FC_FUNC(g_if_sigma_oovv_oovv_no0_x34, G_IF_SIGMA_OOVV_OOVV_NO0_X34)
        (sc, ic, sv1, iv1, T2b.cptr(), V2_sym.cptr(), Xaav.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_oovv_oovv_no1_x34, G_IF_SIGMA_OOVV_OOVV_NO1_X34)
      (sc, ic, Xaav.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.35
  //* X(o2,o1,c,a)  <--  (    1.00000000)  T2(o2,o1,v1,v2) V2(c,v1,a,v2) 
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
    for(int sv2 = 0;sv2 < nir;++sv2){ 
    for(int iv2 = symblockinfo.psym()(sv2,I_V,I_BEGIN);iv2 <= symblockinfo.psym()(sv2,I_V,I_END);++iv2){ 
      T2b = T2.get_amp2(iv2);
      FC_FUNC(g_if_sigma_oovv_oovv_no0_x35, G_IF_SIGMA_OOVV_OOVV_NO0_X35)
        (sc, ic, sv2, iv2, T2b.cptr(), V2_sym.cptr(), Xaav.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_oovv_oovv_no1_x35, G_IF_SIGMA_OOVV_OOVV_NO1_X35)
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
