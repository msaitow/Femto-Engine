                                                                                
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
#include <sci/ctnew2/c_sigma_ccvv_ccvv.h>                                            
                                                                                
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
                                                                                
orz::ct::BareAmpPack orz::ct::sigma_ccvv_ccvv(const orz::ct::Input &ctinp,                                    
                                  const orz::ct::SymBlockInfo &symblockinfo,                                
                                  const orz::ct::HintMO &hintmo,                                            
                                  const orz::ct::RdmPack &rdmPack_sym,                                      
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


  {
  // No.0
  //* S2(w,y,a,c)  <--  (    8.00000000) Y0 T2(w,y,a,c) 
  // The effective tensor is detected .... 
  double Y0 = 0;
  FC_FUNC(g_if_sigma_ccvv_ccvv_y0, G_IF_SIGMA_CCVV_CCVV_Y0)
    (moint1_sym.cptr(), &Y0, nir, nsym, psym);
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x0, G_IF_SIGMA_CCVV_CCVV_NO0_X0)
      (sc, ic, &Y0, T2b.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.1
  //* S2(w,y,a,c)  <--  (   -4.00000000) Y1 T2(y,w,a,c) 
  // The effective tensor is detected .... 
  double Y1 = 0;
  FC_FUNC(g_if_sigma_ccvv_ccvv_y1, G_IF_SIGMA_CCVV_CCVV_Y1)
    (moint1_sym.cptr(), &Y1, nir, nsym, psym);
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x1, G_IF_SIGMA_CCVV_CCVV_NO0_X1)
      (sc, ic, &Y1, T2b.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.2
  //* S2(w,y,a,c)  <--  (   -4.00000000) T2(c1,y,a,c) h(c1,w) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x2, G_IF_SIGMA_CCVV_CCVV_NO0_X2)
      (sc, ic, T2b.cptr(), moint1_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.3
  //* S2(w,y,a,c)  <--  (    2.00000000) T2(y,c1,a,c) h(c1,w) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x3, G_IF_SIGMA_CCVV_CCVV_NO0_X3)
      (sc, ic, T2b.cptr(), moint1_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.4
  //* S2(w,y,a,c)  <--  (    2.00000000) T2(c1,w,a,c) h(c1,y) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x4, G_IF_SIGMA_CCVV_CCVV_NO0_X4)
      (sc, ic, T2b.cptr(), moint1_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.5
  //* S2(w,y,a,c)  <--  (   -4.00000000) T2(w,c1,a,c) h(c1,y) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x5, G_IF_SIGMA_CCVV_CCVV_NO0_X5)
      (sc, ic, T2b.cptr(), moint1_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.6
  //* S2(w,y,a,c)  <--  (    8.00000000) Y2 T2(y,w,c,a) 
  // The effective tensor is detected .... 
  double Y2 = 0;
  FC_FUNC(g_if_sigma_ccvv_ccvv_y6, G_IF_SIGMA_CCVV_CCVV_Y6)
    (moint1_sym.cptr(), &Y2, nir, nsym, psym);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x6, G_IF_SIGMA_CCVV_CCVV_NO0_X6)
        (sa, ia, sc, ic, &Y2, T2b.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.7
  //* S2(w,y,a,c)  <--  (   -4.00000000) Y3 T2(w,y,c,a) 
  // The effective tensor is detected .... 
  double Y3 = 0;
  FC_FUNC(g_if_sigma_ccvv_ccvv_y7, G_IF_SIGMA_CCVV_CCVV_Y7)
    (moint1_sym.cptr(), &Y3, nir, nsym, psym);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x7, G_IF_SIGMA_CCVV_CCVV_NO0_X7)
        (sa, ia, sc, ic, &Y3, T2b.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.8
  //* S2(w,y,a,c)  <--  (   -4.00000000) T2(y,c1,c,a) h(c1,w) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      T2b = T2.get_amp2(ia);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x8, G_IF_SIGMA_CCVV_CCVV_NO0_X8)
        (sa, ia, sc, ic, T2b.cptr(), moint1_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.9
  //* S2(w,y,a,c)  <--  (    2.00000000) T2(c1,y,c,a) h(c1,w) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      T2b = T2.get_amp2(ia);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x9, G_IF_SIGMA_CCVV_CCVV_NO0_X9)
        (sa, ia, sc, ic, T2b.cptr(), moint1_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.10
  //* S2(w,y,a,c)  <--  (    2.00000000) T2(w,c1,c,a) h(c1,y) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      T2b = T2.get_amp2(ia);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x10, G_IF_SIGMA_CCVV_CCVV_NO0_X10)
        (sa, ia, sc, ic, T2b.cptr(), moint1_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.11
  //* S2(w,y,a,c)  <--  (   -4.00000000) T2(c1,w,c,a) h(c1,y) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      T2b = T2.get_amp2(ia);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x11, G_IF_SIGMA_CCVV_CCVV_NO0_X11)
        (sa, ia, sc, ic, T2b.cptr(), moint1_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.12
  //* X()  <--  (    1.00000000)  D1(o1,o2) h(o2,o1) 
  //* S2(w,y,a,c)  <--  (    4.00000000) X T2(w,y,a,c) 
  double X = 0;
  FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x12, G_IF_SIGMA_CCVV_CCVV_NO0_X12)
    (moint1_sym.cptr(), &X, nir, nsym, psym);
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x12, G_IF_SIGMA_CCVV_CCVV_NO1_X12)
      (sc, ic, &X, T2b.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.13
  //* X()  <--  (    1.00000000)  D1(o1,o2) h(o2,o1) 
  //* S2(w,y,a,c)  <--  (   -2.00000000) X T2(y,w,a,c) 
  double X = 0;
  FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x13, G_IF_SIGMA_CCVV_CCVV_NO0_X13)
    (moint1_sym.cptr(), &X, nir, nsym, psym);
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x13, G_IF_SIGMA_CCVV_CCVV_NO1_X13)
      (sc, ic, &X, T2b.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.14
  //* X()  <--  (    1.00000000)  D1(o1,o2) Y4(o1,o2) 
  //* S2(w,y,a,c)  <--  (   -4.00000000) X T2(y,w,a,c) 
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y14, G_IF_SIGMA_CCVV_CCVV_Y14)
      (sc1, ic1, V2_sym.cptr(), Y4.cptr(), nir, nsym, psym);
  }
  }
  double X = 0;
  FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x14, G_IF_SIGMA_CCVV_CCVV_NO0_X14)
    (Y4.cptr(), &X, nir, nsym, psym);
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x14, G_IF_SIGMA_CCVV_CCVV_NO1_X14)
      (sc, ic, &X, T2b.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.15
  //* X()  <--  (    1.00000000)  D1(o1,o2) Y5(o1,o2) 
  //* S2(w,y,a,c)  <--  (    2.00000000) X T2(y,w,a,c) 
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y15, G_IF_SIGMA_CCVV_CCVV_Y15)
      (sc1, ic1, V2_sym.cptr(), Y5.cptr(), nir, nsym, psym);
  }
  }
  double X = 0;
  FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x15, G_IF_SIGMA_CCVV_CCVV_NO0_X15)
    (Y5.cptr(), &X, nir, nsym, psym);
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x15, G_IF_SIGMA_CCVV_CCVV_NO1_X15)
      (sc, ic, &X, T2b.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.16
  //* X()  <--  (    1.00000000)  D1(o1,o2) h(o2,o1) 
  //* S2(w,y,a,c)  <--  (    4.00000000) X T2(y,w,c,a) 
  double X = 0;
  FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x16, G_IF_SIGMA_CCVV_CCVV_NO0_X16)
    (moint1_sym.cptr(), &X, nir, nsym, psym);
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      T2b = T2.get_amp2(ia);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x16, G_IF_SIGMA_CCVV_CCVV_NO1_X16)
        (sa, ia, sc, ic, &X, T2b.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.17
  //* X()  <--  (    1.00000000)  D1(o1,o2) h(o1,o2) 
  //* S2(w,y,a,c)  <--  (   -2.00000000) X T2(w,y,c,a) 
  double X = 0;
  FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x17, G_IF_SIGMA_CCVV_CCVV_NO0_X17)
    (moint1_sym.cptr(), &X, nir, nsym, psym);
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      T2b = T2.get_amp2(ia);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x17, G_IF_SIGMA_CCVV_CCVV_NO1_X17)
        (sa, ia, sc, ic, &X, T2b.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.18
  //* X()  <--  (    1.00000000)  D1(o1,o2) Y6(o1,o2) 
  //* S2(w,y,a,c)  <--  (   -4.00000000) X T2(w,y,c,a) 
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y18, G_IF_SIGMA_CCVV_CCVV_Y18)
      (sc1, ic1, V2_sym.cptr(), Y6.cptr(), nir, nsym, psym);
  }
  }
  double X = 0;
  FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x18, G_IF_SIGMA_CCVV_CCVV_NO0_X18)
    (Y6.cptr(), &X, nir, nsym, psym);
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      T2b = T2.get_amp2(ia);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x18, G_IF_SIGMA_CCVV_CCVV_NO1_X18)
        (sa, ia, sc, ic, &X, T2b.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.19
  //* X()  <--  (    1.00000000)  D1(o1,o2) Y7(o1,o2) 
  //* S2(w,y,a,c)  <--  (    2.00000000) X T2(w,y,c,a) 
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y19, G_IF_SIGMA_CCVV_CCVV_Y19)
      (sc1, ic1, V2_sym.cptr(), Y7.cptr(), nir, nsym, psym);
  }
  }
  double X = 0;
  FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x19, G_IF_SIGMA_CCVV_CCVV_NO0_X19)
    (Y7.cptr(), &X, nir, nsym, psym);
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      T2b = T2.get_amp2(ia);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x19, G_IF_SIGMA_CCVV_CCVV_NO1_X19)
        (sa, ia, sc, ic, &X, T2b.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.20
  //* S2(w,y,a,c)  <--  (    4.00000000) T2(y,w,v1,a) h(c,v1) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x20, G_IF_SIGMA_CCVV_CCVV_NO0_X20)
        (sa, ia, sc, ic, T2b.cptr(), moint1_sym.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.21
  //* S2(w,y,a,c)  <--  (   -2.00000000) T2(w,y,v1,a) h(c,v1) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x21, G_IF_SIGMA_CCVV_CCVV_NO0_X21)
        (sa, ia, sc, ic, T2b.cptr(), moint1_sym.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.22
  //* S2(w,y,a,c)  <--  (   -4.00000000) T2(w,y,v1,a) Y8(c,v1) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nvir, nvir);
  orz::DTensor Y8 = orz::ct::sympack_Xvv(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y22, G_IF_SIGMA_CCVV_CCVV_Y22)
      (sc1, ic1, V2_sym.cptr(), Y8.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x22, G_IF_SIGMA_CCVV_CCVV_NO0_X22)
        (sa, ia, sc, ic, T2b.cptr(), Y8.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.23
  //* S2(w,y,a,c)  <--  (    2.00000000) T2(w,y,v1,a) Y9(c,v1) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nvir, nvir);
  orz::DTensor Y9 = orz::ct::sympack_Xvv(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y23, G_IF_SIGMA_CCVV_CCVV_Y23)
      (sc1, ic1, V2_sym.cptr(), Y9.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x23, G_IF_SIGMA_CCVV_CCVV_NO0_X23)
        (sa, ia, sc, ic, T2b.cptr(), Y9.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.24
  //* S2(w,y,a,c)  <--  (    4.00000000) T2(w,y,v1,c) h(a,v1) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x24, G_IF_SIGMA_CCVV_CCVV_NO0_X24)
      (sc, ic, T2b.cptr(), moint1_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.25
  //* S2(w,y,a,c)  <--  (   -2.00000000) T2(y,w,v1,c) h(a,v1) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x25, G_IF_SIGMA_CCVV_CCVV_NO0_X25)
      (sc, ic, T2b.cptr(), moint1_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.26
  //* S2(w,y,a,c)  <--  (   -4.00000000) T2(y,w,v1,c) Y10(a,v1) 
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y26, G_IF_SIGMA_CCVV_CCVV_Y26)
      (sc1, ic1, V2_sym.cptr(), Y10.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x26, G_IF_SIGMA_CCVV_CCVV_NO0_X26)
      (sc, ic, T2b.cptr(), Y10.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.27
  //* S2(w,y,a,c)  <--  (    2.00000000) T2(y,w,v1,c) Y11(a,v1) 
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y27, G_IF_SIGMA_CCVV_CCVV_Y27)
      (sc1, ic1, V2_sym.cptr(), Y11.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x27, G_IF_SIGMA_CCVV_CCVV_NO0_X27)
      (sc, ic, T2b.cptr(), Y11.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.28
  //* S2(w,y,a,c)  <--  (    4.00000000) T2(w,y,a,v1) h(c,v1) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x28, G_IF_SIGMA_CCVV_CCVV_NO0_X28)
        (sc, ic, sv1, iv1, T2b.cptr(), moint1_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.29
  //* S2(w,y,a,c)  <--  (   -2.00000000) T2(y,w,a,v1) h(c,v1) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x29, G_IF_SIGMA_CCVV_CCVV_NO0_X29)
        (sc, ic, sv1, iv1, T2b.cptr(), moint1_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.30
  //* S2(w,y,a,c)  <--  (   -4.00000000) T2(y,w,a,v1) Y12(c,v1) 
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y30, G_IF_SIGMA_CCVV_CCVV_Y30)
      (sc1, ic1, V2_sym.cptr(), Y12.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x30, G_IF_SIGMA_CCVV_CCVV_NO0_X30)
        (sc, ic, sv1, iv1, T2b.cptr(), Y12.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.31
  //* S2(w,y,a,c)  <--  (    2.00000000) T2(y,w,a,v1) Y13(c,v1) 
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y31, G_IF_SIGMA_CCVV_CCVV_Y31)
      (sc1, ic1, V2_sym.cptr(), Y13.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x31, G_IF_SIGMA_CCVV_CCVV_NO0_X31)
        (sc, ic, sv1, iv1, T2b.cptr(), Y13.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.32
  //* S2(w,y,a,c)  <--  (    4.00000000) T2(y,w,c,v1) h(a,v1) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x32, G_IF_SIGMA_CCVV_CCVV_NO0_X32)
        (sc, ic, sv1, iv1, T2b.cptr(), moint1_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.33
  //* S2(w,y,a,c)  <--  (   -2.00000000) T2(w,y,c,v1) h(a,v1) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x33, G_IF_SIGMA_CCVV_CCVV_NO0_X33)
        (sc, ic, sv1, iv1, T2b.cptr(), moint1_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.34
  //* S2(w,y,a,c)  <--  (   -4.00000000) T2(w,y,c,v1) Y14(a,v1) 
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y34, G_IF_SIGMA_CCVV_CCVV_Y34)
      (sc1, ic1, V2_sym.cptr(), Y14.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x34, G_IF_SIGMA_CCVV_CCVV_NO0_X34)
        (sc, ic, sv1, iv1, T2b.cptr(), Y14.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.35
  //* S2(w,y,a,c)  <--  (    2.00000000) T2(w,y,c,v1) Y15(a,v1) 
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y35, G_IF_SIGMA_CCVV_CCVV_Y35)
      (sc1, ic1, V2_sym.cptr(), Y15.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x35, G_IF_SIGMA_CCVV_CCVV_NO0_X35)
        (sc, ic, sv1, iv1, T2b.cptr(), Y15.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.36
  //* S2(w,y,a,c)  <--  (    8.00000000) Y16 T2(w,y,a,c) 
  // The effective tensor is detected .... 
  double Y16 = 0;
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y36, G_IF_SIGMA_CCVV_CCVV_Y36)
      (sc1, ic1, V2_sym.cptr(), &Y16, nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x36, G_IF_SIGMA_CCVV_CCVV_NO0_X36)
      (sc, ic, &Y16, T2b.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.37
  //* S2(w,y,a,c)  <--  (   -4.00000000) Y17 T2(w,y,a,c) 
  // The effective tensor is detected .... 
  double Y17 = 0;
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y37, G_IF_SIGMA_CCVV_CCVV_Y37)
      (sc1, ic1, V2_sym.cptr(), &Y17, nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x37, G_IF_SIGMA_CCVV_CCVV_NO0_X37)
      (sc, ic, &Y17, T2b.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.38
  //* S2(w,y,a,c)  <--  (   -8.00000000) T2(w,c2,a,c) Y18(c2,y) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nclosed, nclosed);
  orz::DTensor Y18 = orz::ct::sympack_Xcc(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y38, G_IF_SIGMA_CCVV_CCVV_Y38)
      (sc1, ic1, V2_sym.cptr(), Y18.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x38, G_IF_SIGMA_CCVV_CCVV_NO0_X38)
      (sc, ic, T2b.cptr(), Y18.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.39
  //* S2(w,y,a,c)  <--  (    4.00000000) T2(w,c2,a,c) Y19(c2,y) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nclosed, nclosed);
  orz::DTensor Y19 = orz::ct::sympack_Xcc(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y39, G_IF_SIGMA_CCVV_CCVV_Y39)
      (sc1, ic1, V2_sym.cptr(), Y19.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x39, G_IF_SIGMA_CCVV_CCVV_NO0_X39)
      (sc, ic, T2b.cptr(), Y19.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.40
  //* S2(w,y,a,c)  <--  (   -4.00000000) Y20 T2(y,w,a,c) 
  // The effective tensor is detected .... 
  double Y20 = 0;
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y40, G_IF_SIGMA_CCVV_CCVV_Y40)
      (sc1, ic1, V2_sym.cptr(), &Y20, nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x40, G_IF_SIGMA_CCVV_CCVV_NO0_X40)
      (sc, ic, &Y20, T2b.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.41
  //* S2(w,y,a,c)  <--  (    2.00000000) Y21 T2(y,w,a,c) 
  // The effective tensor is detected .... 
  double Y21 = 0;
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y41, G_IF_SIGMA_CCVV_CCVV_Y41)
      (sc1, ic1, V2_sym.cptr(), &Y21, nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x41, G_IF_SIGMA_CCVV_CCVV_NO0_X41)
      (sc, ic, &Y21, T2b.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.42
  //* S2(w,y,a,c)  <--  (    4.00000000) T2(c2,w,a,c) Y22(c2,y) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nclosed, nclosed);
  orz::DTensor Y22 = orz::ct::sympack_Xcc(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y42, G_IF_SIGMA_CCVV_CCVV_Y42)
      (sc1, ic1, V2_sym.cptr(), Y22.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x42, G_IF_SIGMA_CCVV_CCVV_NO0_X42)
      (sc, ic, T2b.cptr(), Y22.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.43
  //* S2(w,y,a,c)  <--  (   -2.00000000) T2(c2,w,a,c) Y23(c2,y) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nclosed, nclosed);
  orz::DTensor Y23 = orz::ct::sympack_Xcc(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y43, G_IF_SIGMA_CCVV_CCVV_Y43)
      (sc1, ic1, V2_sym.cptr(), Y23.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x43, G_IF_SIGMA_CCVV_CCVV_NO0_X43)
      (sc, ic, T2b.cptr(), Y23.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.44
  //* S2(w,y,a,c)  <--  (    4.00000000) T2(y,c2,a,c) Y24(c2,w) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nclosed, nclosed);
  orz::DTensor Y24 = orz::ct::sympack_Xcc(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y44, G_IF_SIGMA_CCVV_CCVV_Y44)
      (sc1, ic1, V2_sym.cptr(), Y24.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x44, G_IF_SIGMA_CCVV_CCVV_NO0_X44)
      (sc, ic, T2b.cptr(), Y24.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.45
  //* S2(w,y,a,c)  <--  (   -2.00000000) T2(y,c2,a,c) Y25(c2,w) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nclosed, nclosed);
  orz::DTensor Y25 = orz::ct::sympack_Xcc(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y45, G_IF_SIGMA_CCVV_CCVV_Y45)
      (sc1, ic1, V2_sym.cptr(), Y25.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x45, G_IF_SIGMA_CCVV_CCVV_NO0_X45)
      (sc, ic, T2b.cptr(), Y25.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.46
  //* S2(w,y,a,c)  <--  (   -8.00000000) T2(c2,y,a,c) Y26(c2,w) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nclosed, nclosed);
  orz::DTensor Y26 = orz::ct::sympack_Xcc(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y46, G_IF_SIGMA_CCVV_CCVV_Y46)
      (sc1, ic1, V2_sym.cptr(), Y26.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x46, G_IF_SIGMA_CCVV_CCVV_NO0_X46)
      (sc, ic, T2b.cptr(), Y26.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.47
  //* S2(w,y,a,c)  <--  (    4.00000000) T2(c2,y,a,c) Y27(c2,w) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nclosed, nclosed);
  orz::DTensor Y27 = orz::ct::sympack_Xcc(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y47, G_IF_SIGMA_CCVV_CCVV_Y47)
      (sc1, ic1, V2_sym.cptr(), Y27.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x47, G_IF_SIGMA_CCVV_CCVV_NO0_X47)
      (sc, ic, T2b.cptr(), Y27.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.48
  //* S2(w,y,a,c)  <--  (    4.00000000) T2(c2,c1,a,c) V2(c1,y,c2,w) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
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
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x48, G_IF_SIGMA_CCVV_CCVV_NO0_X48)
        (sc, ic, sc1, ic1, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.49
  //* S2(w,y,a,c)  <--  (   -2.00000000) T2(c2,c1,a,c) V2(c1,w,c2,y) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
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
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x49, G_IF_SIGMA_CCVV_CCVV_NO0_X49)
        (sc, ic, sc1, ic1, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.50
  //* S2(w,y,a,c)  <--  (    8.00000000) Y28 T2(y,w,c,a) 
  // The effective tensor is detected .... 
  double Y28 = 0;
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y50, G_IF_SIGMA_CCVV_CCVV_Y50)
      (sc1, ic1, V2_sym.cptr(), &Y28, nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x50, G_IF_SIGMA_CCVV_CCVV_NO0_X50)
        (sa, ia, sc, ic, &Y28, T2b.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.51
  //* S2(w,y,a,c)  <--  (   -4.00000000) Y29 T2(y,w,c,a) 
  // The effective tensor is detected .... 
  double Y29 = 0;
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y51, G_IF_SIGMA_CCVV_CCVV_Y51)
      (sc1, ic1, V2_sym.cptr(), &Y29, nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x51, G_IF_SIGMA_CCVV_CCVV_NO0_X51)
        (sa, ia, sc, ic, &Y29, T2b.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.52
  //* S2(w,y,a,c)  <--  (   -8.00000000) T2(y,c2,c,a) Y30(c2,w) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nclosed, nclosed);
  orz::DTensor Y30 = orz::ct::sympack_Xcc(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y52, G_IF_SIGMA_CCVV_CCVV_Y52)
      (sc1, ic1, V2_sym.cptr(), Y30.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      T2b = T2.get_amp2(ia);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x52, G_IF_SIGMA_CCVV_CCVV_NO0_X52)
        (sa, ia, sc, ic, T2b.cptr(), Y30.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.53
  //* S2(w,y,a,c)  <--  (    4.00000000) T2(y,c2,c,a) Y31(c2,w) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nclosed, nclosed);
  orz::DTensor Y31 = orz::ct::sympack_Xcc(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y53, G_IF_SIGMA_CCVV_CCVV_Y53)
      (sc1, ic1, V2_sym.cptr(), Y31.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      T2b = T2.get_amp2(ia);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x53, G_IF_SIGMA_CCVV_CCVV_NO0_X53)
        (sa, ia, sc, ic, T2b.cptr(), Y31.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.54
  //* S2(w,y,a,c)  <--  (   -4.00000000) Y32 T2(w,y,c,a) 
  // The effective tensor is detected .... 
  double Y32 = 0;
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y54, G_IF_SIGMA_CCVV_CCVV_Y54)
      (sc1, ic1, V2_sym.cptr(), &Y32, nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x54, G_IF_SIGMA_CCVV_CCVV_NO0_X54)
        (sa, ia, sc, ic, &Y32, T2b.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.55
  //* S2(w,y,a,c)  <--  (    2.00000000) Y33 T2(w,y,c,a) 
  // The effective tensor is detected .... 
  double Y33 = 0;
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y55, G_IF_SIGMA_CCVV_CCVV_Y55)
      (sc1, ic1, V2_sym.cptr(), &Y33, nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x55, G_IF_SIGMA_CCVV_CCVV_NO0_X55)
        (sa, ia, sc, ic, &Y33, T2b.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.56
  //* S2(w,y,a,c)  <--  (    4.00000000) T2(c2,y,c,a) Y34(c2,w) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nclosed, nclosed);
  orz::DTensor Y34 = orz::ct::sympack_Xcc(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y56, G_IF_SIGMA_CCVV_CCVV_Y56)
      (sc1, ic1, V2_sym.cptr(), Y34.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      T2b = T2.get_amp2(ia);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x56, G_IF_SIGMA_CCVV_CCVV_NO0_X56)
        (sa, ia, sc, ic, T2b.cptr(), Y34.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.57
  //* S2(w,y,a,c)  <--  (   -2.00000000) T2(c2,y,c,a) Y35(c2,w) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nclosed, nclosed);
  orz::DTensor Y35 = orz::ct::sympack_Xcc(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y57, G_IF_SIGMA_CCVV_CCVV_Y57)
      (sc1, ic1, V2_sym.cptr(), Y35.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      T2b = T2.get_amp2(ia);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x57, G_IF_SIGMA_CCVV_CCVV_NO0_X57)
        (sa, ia, sc, ic, T2b.cptr(), Y35.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.58
  //* S2(w,y,a,c)  <--  (    4.00000000) T2(w,c2,c,a) Y36(c2,y) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nclosed, nclosed);
  orz::DTensor Y36 = orz::ct::sympack_Xcc(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y58, G_IF_SIGMA_CCVV_CCVV_Y58)
      (sc1, ic1, V2_sym.cptr(), Y36.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      T2b = T2.get_amp2(ia);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x58, G_IF_SIGMA_CCVV_CCVV_NO0_X58)
        (sa, ia, sc, ic, T2b.cptr(), Y36.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.59
  //* S2(w,y,a,c)  <--  (   -2.00000000) T2(w,c2,c,a) Y37(c2,y) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nclosed, nclosed);
  orz::DTensor Y37 = orz::ct::sympack_Xcc(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y59, G_IF_SIGMA_CCVV_CCVV_Y59)
      (sc1, ic1, V2_sym.cptr(), Y37.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      T2b = T2.get_amp2(ia);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x59, G_IF_SIGMA_CCVV_CCVV_NO0_X59)
        (sa, ia, sc, ic, T2b.cptr(), Y37.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.60
  //* S2(w,y,a,c)  <--  (   -8.00000000) T2(c2,w,c,a) Y38(c2,y) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nclosed, nclosed);
  orz::DTensor Y38 = orz::ct::sympack_Xcc(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y60, G_IF_SIGMA_CCVV_CCVV_Y60)
      (sc1, ic1, V2_sym.cptr(), Y38.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      T2b = T2.get_amp2(ia);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x60, G_IF_SIGMA_CCVV_CCVV_NO0_X60)
        (sa, ia, sc, ic, T2b.cptr(), Y38.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.61
  //* S2(w,y,a,c)  <--  (    4.00000000) T2(c2,w,c,a) Y39(c2,y) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nclosed, nclosed);
  orz::DTensor Y39 = orz::ct::sympack_Xcc(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y61, G_IF_SIGMA_CCVV_CCVV_Y61)
      (sc1, ic1, V2_sym.cptr(), Y39.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      T2b = T2.get_amp2(ia);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x61, G_IF_SIGMA_CCVV_CCVV_NO0_X61)
        (sa, ia, sc, ic, T2b.cptr(), Y39.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.62
  //* S2(w,y,a,c)  <--  (    4.00000000) T2(c2,c1,c,a) V2(c1,w,c2,y) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      T2b = T2.get_amp2(ia);
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
        FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x62, G_IF_SIGMA_CCVV_CCVV_NO0_X62)
          (sa, ia, sc, ic, sc1, ic1, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym);
      }
      }
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.63
  //* S2(w,y,a,c)  <--  (   -2.00000000) T2(c2,c1,c,a) V2(c1,y,c2,w) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      T2b = T2.get_amp2(ia);
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
        FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x63, G_IF_SIGMA_CCVV_CCVV_NO0_X63)
          (sa, ia, sc, ic, sc1, ic1, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym);
      }
      }
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.64
  //* X()  <--  (    1.00000000)  D1(o1,o2) Y40(o1,o2) 
  //* S2(w,y,a,c)  <--  (    8.00000000) X T2(w,y,a,c) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nocc, nocc);
  orz::DTensor Y40 = orz::ct::sympack_Xaa(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y64, G_IF_SIGMA_CCVV_CCVV_Y64)
      (sc1, ic1, V2_sym.cptr(), Y40.cptr(), nir, nsym, psym);
  }
  }
  double X = 0;
  FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x64, G_IF_SIGMA_CCVV_CCVV_NO0_X64)
    (Y40.cptr(), &X, nir, nsym, psym);
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x64, G_IF_SIGMA_CCVV_CCVV_NO1_X64)
      (sc, ic, &X, T2b.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.65
  //* X(c1,y)  <--  (    1.00000000)  D1(o1,o2) V2(c1,y,o1,o2) 
  //* S2(w,y,a,c)  <--  (   -4.00000000) T2(w,c1,a,c) X(c1,y) 
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
    orz::DTensor X(nclosed);
    orz::DTensor Xc = orz::ct::sympack_Xc(symblockinfo, sc1, X);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x65, G_IF_SIGMA_CCVV_CCVV_NO0_X65)
      (sc1, ic1, V2_sym.cptr(), Xc.cptr(), nir, nsym, psym);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      T2b = T2.get_amp2(ic);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x65, G_IF_SIGMA_CCVV_CCVV_NO1_X65)
        (sc, ic, sc1, ic1, T2b.cptr(), Xc.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.66
  //* X(c1,y)  <--  (    1.00000000)  D1(o1,o2) V2(c1,y,o1,o2) 
  //* S2(w,y,a,c)  <--  (    2.00000000) T2(c1,w,a,c) X(c1,y) 
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
    orz::DTensor X(nclosed);
    orz::DTensor Xc = orz::ct::sympack_Xc(symblockinfo, sc1, X);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x66, G_IF_SIGMA_CCVV_CCVV_NO0_X66)
      (sc1, ic1, V2_sym.cptr(), Xc.cptr(), nir, nsym, psym);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      T2b = T2.get_amp2(ic);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x66, G_IF_SIGMA_CCVV_CCVV_NO1_X66)
        (sc, ic, sc1, ic1, T2b.cptr(), Xc.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.67
  //* X(c1,y)  <--  (    1.00000000)  D1(o1,o2) V2(c1,o1,y,o2) 
  //* S2(w,y,a,c)  <--  (   -1.00000000) T2(c1,w,a,c) X(c1,y) 
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
    orz::DTensor X(nclosed);
    orz::DTensor Xc = orz::ct::sympack_Xc(symblockinfo, sc1, X);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x67, G_IF_SIGMA_CCVV_CCVV_NO0_X67)
      (sc1, ic1, V2_sym.cptr(), Xc.cptr(), nir, nsym, psym);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      T2b = T2.get_amp2(ic);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x67, G_IF_SIGMA_CCVV_CCVV_NO1_X67)
        (sc, ic, sc1, ic1, T2b.cptr(), Xc.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.68
  //* X(c1,w)  <--  (    1.00000000)  D1(o1,o2) V2(c1,w,o1,o2) 
  //* S2(w,y,a,c)  <--  (    2.00000000) T2(y,c1,a,c) X(c1,w) 
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
    orz::DTensor X(nclosed);
    orz::DTensor Xc = orz::ct::sympack_Xc(symblockinfo, sc1, X);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x68, G_IF_SIGMA_CCVV_CCVV_NO0_X68)
      (sc1, ic1, V2_sym.cptr(), Xc.cptr(), nir, nsym, psym);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      T2b = T2.get_amp2(ic);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x68, G_IF_SIGMA_CCVV_CCVV_NO1_X68)
        (sc, ic, sc1, ic1, T2b.cptr(), Xc.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.69
  //* X(c1,w)  <--  (    1.00000000)  D1(o1,o2) V2(c1,o1,w,o2) 
  //* S2(w,y,a,c)  <--  (   -1.00000000) T2(y,c1,a,c) X(c1,w) 
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
    orz::DTensor X(nclosed);
    orz::DTensor Xc = orz::ct::sympack_Xc(symblockinfo, sc1, X);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x69, G_IF_SIGMA_CCVV_CCVV_NO0_X69)
      (sc1, ic1, V2_sym.cptr(), Xc.cptr(), nir, nsym, psym);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      T2b = T2.get_amp2(ic);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x69, G_IF_SIGMA_CCVV_CCVV_NO1_X69)
        (sc, ic, sc1, ic1, T2b.cptr(), Xc.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.70
  //* X(c1,w)  <--  (    1.00000000)  D1(o1,o2) V2(c1,w,o1,o2) 
  //* S2(w,y,a,c)  <--  (   -4.00000000) T2(c1,y,a,c) X(c1,w) 
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
    orz::DTensor X(nclosed);
    orz::DTensor Xc = orz::ct::sympack_Xc(symblockinfo, sc1, X);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x70, G_IF_SIGMA_CCVV_CCVV_NO0_X70)
      (sc1, ic1, V2_sym.cptr(), Xc.cptr(), nir, nsym, psym);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      T2b = T2.get_amp2(ic);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x70, G_IF_SIGMA_CCVV_CCVV_NO1_X70)
        (sc, ic, sc1, ic1, T2b.cptr(), Xc.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.71
  //* X()  <--  (    1.00000000)  D1(o1,o2) Y41(o1,o2) 
  //* S2(w,y,a,c)  <--  (    8.00000000) X T2(y,w,c,a) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nocc, nocc);
  orz::DTensor Y41 = orz::ct::sympack_Xaa(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y71, G_IF_SIGMA_CCVV_CCVV_Y71)
      (sc1, ic1, V2_sym.cptr(), Y41.cptr(), nir, nsym, psym);
  }
  }
  double X = 0;
  FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x71, G_IF_SIGMA_CCVV_CCVV_NO0_X71)
    (Y41.cptr(), &X, nir, nsym, psym);
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      T2b = T2.get_amp2(ia);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x71, G_IF_SIGMA_CCVV_CCVV_NO1_X71)
        (sa, ia, sc, ic, &X, T2b.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.72
  //* X(c1,w)  <--  (    1.00000000)  D1(o1,o2) V2(c1,w,o1,o2) 
  //* S2(w,y,a,c)  <--  (   -4.00000000) T2(y,c1,c,a) X(c1,w) 
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
    orz::DTensor X(nclosed);
    orz::DTensor Xc = orz::ct::sympack_Xc(symblockinfo, sc1, X);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x72, G_IF_SIGMA_CCVV_CCVV_NO0_X72)
      (sc1, ic1, V2_sym.cptr(), Xc.cptr(), nir, nsym, psym);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      for(int sa = 0;sa < nir;++sa){ 
      for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
        T2b = T2.get_amp2(ia);
        FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x72, G_IF_SIGMA_CCVV_CCVV_NO1_X72)
          (sa, ia, sc, ic, sc1, ic1, T2b.cptr(), Xc.cptr(), S2b.cptr(), nir, nsym, psym);
      }
      }
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.73
  //* X(c1,w)  <--  (    1.00000000)  D1(o1,o2) V2(c1,w,o1,o2) 
  //* S2(w,y,a,c)  <--  (    2.00000000) T2(c1,y,c,a) X(c1,w) 
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
    orz::DTensor X(nclosed);
    orz::DTensor Xc = orz::ct::sympack_Xc(symblockinfo, sc1, X);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x73, G_IF_SIGMA_CCVV_CCVV_NO0_X73)
      (sc1, ic1, V2_sym.cptr(), Xc.cptr(), nir, nsym, psym);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      for(int sa = 0;sa < nir;++sa){ 
      for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
        T2b = T2.get_amp2(ia);
        FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x73, G_IF_SIGMA_CCVV_CCVV_NO1_X73)
          (sa, ia, sc, ic, sc1, ic1, T2b.cptr(), Xc.cptr(), S2b.cptr(), nir, nsym, psym);
      }
      }
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.74
  //* X(c1,w)  <--  (    1.00000000)  D1(o1,o2) V2(c1,o1,w,o2) 
  //* S2(w,y,a,c)  <--  (   -1.00000000) T2(c1,y,c,a) X(c1,w) 
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
    orz::DTensor X(nclosed);
    orz::DTensor Xc = orz::ct::sympack_Xc(symblockinfo, sc1, X);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x74, G_IF_SIGMA_CCVV_CCVV_NO0_X74)
      (sc1, ic1, V2_sym.cptr(), Xc.cptr(), nir, nsym, psym);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      for(int sa = 0;sa < nir;++sa){ 
      for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
        T2b = T2.get_amp2(ia);
        FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x74, G_IF_SIGMA_CCVV_CCVV_NO1_X74)
          (sa, ia, sc, ic, sc1, ic1, T2b.cptr(), Xc.cptr(), S2b.cptr(), nir, nsym, psym);
      }
      }
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.75
  //* X(c1,y)  <--  (    1.00000000)  D1(o1,o2) V2(c1,y,o1,o2) 
  //* S2(w,y,a,c)  <--  (    2.00000000) T2(w,c1,c,a) X(c1,y) 
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
    orz::DTensor X(nclosed);
    orz::DTensor Xc = orz::ct::sympack_Xc(symblockinfo, sc1, X);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x75, G_IF_SIGMA_CCVV_CCVV_NO0_X75)
      (sc1, ic1, V2_sym.cptr(), Xc.cptr(), nir, nsym, psym);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      for(int sa = 0;sa < nir;++sa){ 
      for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
        T2b = T2.get_amp2(ia);
        FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x75, G_IF_SIGMA_CCVV_CCVV_NO1_X75)
          (sa, ia, sc, ic, sc1, ic1, T2b.cptr(), Xc.cptr(), S2b.cptr(), nir, nsym, psym);
      }
      }
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.76
  //* X(c1,y)  <--  (    1.00000000)  D1(o1,o2) V2(c1,o1,y,o2) 
  //* S2(w,y,a,c)  <--  (   -1.00000000) T2(w,c1,c,a) X(c1,y) 
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
    orz::DTensor X(nclosed);
    orz::DTensor Xc = orz::ct::sympack_Xc(symblockinfo, sc1, X);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x76, G_IF_SIGMA_CCVV_CCVV_NO0_X76)
      (sc1, ic1, V2_sym.cptr(), Xc.cptr(), nir, nsym, psym);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      for(int sa = 0;sa < nir;++sa){ 
      for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
        T2b = T2.get_amp2(ia);
        FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x76, G_IF_SIGMA_CCVV_CCVV_NO1_X76)
          (sa, ia, sc, ic, sc1, ic1, T2b.cptr(), Xc.cptr(), S2b.cptr(), nir, nsym, psym);
      }
      }
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.77
  //* X(c1,y)  <--  (    1.00000000)  D1(o1,o2) V2(c1,y,o1,o2) 
  //* S2(w,y,a,c)  <--  (   -4.00000000) T2(c1,w,c,a) X(c1,y) 
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
    orz::DTensor X(nclosed);
    orz::DTensor Xc = orz::ct::sympack_Xc(symblockinfo, sc1, X);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x77, G_IF_SIGMA_CCVV_CCVV_NO0_X77)
      (sc1, ic1, V2_sym.cptr(), Xc.cptr(), nir, nsym, psym);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      for(int sa = 0;sa < nir;++sa){ 
      for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
        T2b = T2.get_amp2(ia);
        FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x77, G_IF_SIGMA_CCVV_CCVV_NO1_X77)
          (sa, ia, sc, ic, sc1, ic1, T2b.cptr(), Xc.cptr(), S2b.cptr(), nir, nsym, psym);
      }
      }
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.78
  //* X()  <--  (    1.00000000)  D1(o1,o2) Y42(o1,o2) 
  //* S2(w,y,a,c)  <--  (   -4.00000000) X T2(w,y,a,c) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nocc, nocc);
  orz::DTensor Y42 = orz::ct::sympack_Xaa(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y78, G_IF_SIGMA_CCVV_CCVV_Y78)
      (sc1, ic1, V2_sym.cptr(), Y42.cptr(), nir, nsym, psym);
  }
  }
  double X = 0;
  FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x78, G_IF_SIGMA_CCVV_CCVV_NO0_X78)
    (Y42.cptr(), &X, nir, nsym, psym);
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x78, G_IF_SIGMA_CCVV_CCVV_NO1_X78)
      (sc, ic, &X, T2b.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.79
  //* X(c1,y)  <--  (    1.00000000)  D1(o1,o2) V2(c1,o1,y,o2) 
  //* S2(w,y,a,c)  <--  (    2.00000000) T2(w,c1,a,c) X(c1,y) 
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
    orz::DTensor X(nclosed);
    orz::DTensor Xc = orz::ct::sympack_Xc(symblockinfo, sc1, X);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x79, G_IF_SIGMA_CCVV_CCVV_NO0_X79)
      (sc1, ic1, V2_sym.cptr(), Xc.cptr(), nir, nsym, psym);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      T2b = T2.get_amp2(ic);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x79, G_IF_SIGMA_CCVV_CCVV_NO1_X79)
        (sc, ic, sc1, ic1, T2b.cptr(), Xc.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.80
  //* X(c1,w)  <--  (    1.00000000)  D1(o1,o2) V2(c1,o1,w,o2) 
  //* S2(w,y,a,c)  <--  (    2.00000000) T2(c1,y,a,c) X(c1,w) 
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
    orz::DTensor X(nclosed);
    orz::DTensor Xc = orz::ct::sympack_Xc(symblockinfo, sc1, X);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x80, G_IF_SIGMA_CCVV_CCVV_NO0_X80)
      (sc1, ic1, V2_sym.cptr(), Xc.cptr(), nir, nsym, psym);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      T2b = T2.get_amp2(ic);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x80, G_IF_SIGMA_CCVV_CCVV_NO1_X80)
        (sc, ic, sc1, ic1, T2b.cptr(), Xc.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.81
  //* X()  <--  (    1.00000000)  D1(o1,o2) Y43(o1,o2) 
  //* S2(w,y,a,c)  <--  (   -4.00000000) X T2(y,w,c,a) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nocc, nocc);
  orz::DTensor Y43 = orz::ct::sympack_Xaa(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y81, G_IF_SIGMA_CCVV_CCVV_Y81)
      (sc1, ic1, V2_sym.cptr(), Y43.cptr(), nir, nsym, psym);
  }
  }
  double X = 0;
  FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x81, G_IF_SIGMA_CCVV_CCVV_NO0_X81)
    (Y43.cptr(), &X, nir, nsym, psym);
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      T2b = T2.get_amp2(ia);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x81, G_IF_SIGMA_CCVV_CCVV_NO1_X81)
        (sa, ia, sc, ic, &X, T2b.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.82
  //* X(c1,w)  <--  (    1.00000000)  D1(o1,o2) V2(c1,o1,w,o2) 
  //* S2(w,y,a,c)  <--  (    2.00000000) T2(y,c1,c,a) X(c1,w) 
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
    orz::DTensor X(nclosed);
    orz::DTensor Xc = orz::ct::sympack_Xc(symblockinfo, sc1, X);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x82, G_IF_SIGMA_CCVV_CCVV_NO0_X82)
      (sc1, ic1, V2_sym.cptr(), Xc.cptr(), nir, nsym, psym);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      for(int sa = 0;sa < nir;++sa){ 
      for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
        T2b = T2.get_amp2(ia);
        FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x82, G_IF_SIGMA_CCVV_CCVV_NO1_X82)
          (sa, ia, sc, ic, sc1, ic1, T2b.cptr(), Xc.cptr(), S2b.cptr(), nir, nsym, psym);
      }
      }
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.83
  //* X(c1,y)  <--  (    1.00000000)  D1(o1,o2) V2(c1,o1,y,o2) 
  //* S2(w,y,a,c)  <--  (    2.00000000) T2(c1,w,c,a) X(c1,y) 
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
    orz::DTensor X(nclosed);
    orz::DTensor Xc = orz::ct::sympack_Xc(symblockinfo, sc1, X);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x83, G_IF_SIGMA_CCVV_CCVV_NO0_X83)
      (sc1, ic1, V2_sym.cptr(), Xc.cptr(), nir, nsym, psym);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      for(int sa = 0;sa < nir;++sa){ 
      for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
        T2b = T2.get_amp2(ia);
        FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x83, G_IF_SIGMA_CCVV_CCVV_NO1_X83)
          (sa, ia, sc, ic, sc1, ic1, T2b.cptr(), Xc.cptr(), S2b.cptr(), nir, nsym, psym);
      }
      }
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.84
  //* S2(w,y,a,c)  <--  (    8.00000000) T2(y,w,v1,a) Y44(c,v1) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nvir, nvir);
  orz::DTensor Y44 = orz::ct::sympack_Xvv(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y84, G_IF_SIGMA_CCVV_CCVV_Y84)
      (sc1, ic1, V2_sym.cptr(), Y44.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x84, G_IF_SIGMA_CCVV_CCVV_NO0_X84)
        (sa, ia, sc, ic, T2b.cptr(), Y44.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.85
  //* S2(w,y,a,c)  <--  (   -4.00000000) T2(y,c1,v1,a) V2(c,v1,c1,w) 
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
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x85, G_IF_SIGMA_CCVV_CCVV_NO0_X85)
        (sa, ia, sc, ic, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.86
  //* S2(w,y,a,c)  <--  (    2.00000000) T2(c1,y,v1,a) V2(c,v1,c1,w) 
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
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x86, G_IF_SIGMA_CCVV_CCVV_NO0_X86)
        (sa, ia, sc, ic, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.87
  //* S2(w,y,a,c)  <--  (   -4.00000000) T2(c1,y,v1,a) V2(c,w,c1,v1) 
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
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x87, G_IF_SIGMA_CCVV_CCVV_NO0_X87)
        (sa, ia, sc, ic, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.88
  //* S2(w,y,a,c)  <--  (    2.00000000) T2(w,c1,v1,a) V2(c,v1,c1,y) 
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
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x88, G_IF_SIGMA_CCVV_CCVV_NO0_X88)
        (sa, ia, sc, ic, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.89
  //* S2(w,y,a,c)  <--  (   -4.00000000) T2(w,c1,v1,a) V2(c,y,c1,v1) 
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
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x89, G_IF_SIGMA_CCVV_CCVV_NO0_X89)
        (sa, ia, sc, ic, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.90
  //* S2(w,y,a,c)  <--  (   -4.00000000) T2(c1,w,v1,a) V2(c,v1,c1,y) 
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
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x90, G_IF_SIGMA_CCVV_CCVV_NO0_X90)
        (sa, ia, sc, ic, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.91
  //* S2(w,y,a,c)  <--  (    8.00000000) T2(w,y,v1,c) Y45(a,v1) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nvir, nvir);
  orz::DTensor Y45 = orz::ct::sympack_Xvv(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y91, G_IF_SIGMA_CCVV_CCVV_Y91)
      (sc1, ic1, V2_sym.cptr(), Y45.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x91, G_IF_SIGMA_CCVV_CCVV_NO0_X91)
      (sc, ic, T2b.cptr(), Y45.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.92
  //* S2(w,y,a,c)  <--  (   -4.00000000) T2(w,c1,v1,c) V2(a,v1,c1,y) 
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
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x92, G_IF_SIGMA_CCVV_CCVV_NO0_X92)
        (sa, ia, sc, ic, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.93
  //* S2(w,y,a,c)  <--  (    2.00000000) T2(c1,w,v1,c) V2(a,v1,c1,y) 
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
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x93, G_IF_SIGMA_CCVV_CCVV_NO0_X93)
        (sa, ia, sc, ic, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.94
  //* S2(w,y,a,c)  <--  (   -4.00000000) T2(c1,w,v1,c) V2(v1,c1,y,a) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
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
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x94, G_IF_SIGMA_CCVV_CCVV_NO0_X94)
        (sc, ic, sv1, iv1, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.95
  //* S2(w,y,a,c)  <--  (    2.00000000) T2(y,c1,v1,c) V2(a,v1,c1,w) 
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
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x95, G_IF_SIGMA_CCVV_CCVV_NO0_X95)
        (sa, ia, sc, ic, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.96
  //* S2(w,y,a,c)  <--  (   -4.00000000) T2(y,c1,v1,c) V2(v1,c1,w,a) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
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
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x96, G_IF_SIGMA_CCVV_CCVV_NO0_X96)
        (sc, ic, sv1, iv1, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.97
  //* S2(w,y,a,c)  <--  (   -4.00000000) T2(c1,y,v1,c) V2(a,v1,c1,w) 
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
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x97, G_IF_SIGMA_CCVV_CCVV_NO0_X97)
        (sa, ia, sc, ic, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.98
  //* S2(w,y,a,c)  <--  (    8.00000000) T2(w,y,a,v1) Y46(c,v1) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nvir, nvir);
  orz::DTensor Y46 = orz::ct::sympack_Xvv(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y98, G_IF_SIGMA_CCVV_CCVV_Y98)
      (sc1, ic1, V2_sym.cptr(), Y46.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x98, G_IF_SIGMA_CCVV_CCVV_NO0_X98)
        (sc, ic, sv1, iv1, T2b.cptr(), Y46.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.99
  //* S2(w,y,a,c)  <--  (   -4.00000000) T2(w,c1,a,v1) V2(c,v1,c1,y) 
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
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x99, G_IF_SIGMA_CCVV_CCVV_NO0_X99)
        (sc, ic, sv1, iv1, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.100
  //* S2(w,y,a,c)  <--  (    2.00000000) T2(c1,w,a,v1) V2(c,v1,c1,y) 
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
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x100, G_IF_SIGMA_CCVV_CCVV_NO0_X100)
        (sc, ic, sv1, iv1, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.101
  //* S2(w,y,a,c)  <--  (   -4.00000000) T2(c1,w,a,v1) V2(c,y,c1,v1) 
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
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x101, G_IF_SIGMA_CCVV_CCVV_NO0_X101)
        (sc, ic, sv1, iv1, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.102
  //* S2(w,y,a,c)  <--  (    2.00000000) T2(y,c1,a,v1) V2(c,v1,c1,w) 
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
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x102, G_IF_SIGMA_CCVV_CCVV_NO0_X102)
        (sc, ic, sv1, iv1, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.103
  //* S2(w,y,a,c)  <--  (   -4.00000000) T2(y,c1,a,v1) V2(c,w,c1,v1) 
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
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x103, G_IF_SIGMA_CCVV_CCVV_NO0_X103)
        (sc, ic, sv1, iv1, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.104
  //* S2(w,y,a,c)  <--  (   -4.00000000) T2(c1,y,a,v1) V2(c,v1,c1,w) 
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
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x104, G_IF_SIGMA_CCVV_CCVV_NO0_X104)
        (sc, ic, sv1, iv1, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.105
  //* S2(w,y,a,c)  <--  (    8.00000000) T2(y,w,c,v1) Y47(a,v1) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nvir, nvir);
  orz::DTensor Y47 = orz::ct::sympack_Xvv(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y105, G_IF_SIGMA_CCVV_CCVV_Y105)
      (sc1, ic1, V2_sym.cptr(), Y47.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x105, G_IF_SIGMA_CCVV_CCVV_NO0_X105)
        (sc, ic, sv1, iv1, T2b.cptr(), Y47.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.106
  //* S2(w,y,a,c)  <--  (   -4.00000000) T2(y,c1,c,v1) V2(v1,a,c1,w) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
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
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x106, G_IF_SIGMA_CCVV_CCVV_NO0_X106)
        (sc, ic, sv1, iv1, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.107
  //* S2(w,y,a,c)  <--  (    2.00000000) T2(c1,y,c,v1) V2(v1,a,c1,w) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
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
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x107, G_IF_SIGMA_CCVV_CCVV_NO0_X107)
        (sc, ic, sv1, iv1, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.108
  //* S2(w,y,a,c)  <--  (   -4.00000000) T2(c1,y,c,v1) V2(v1,c1,w,a) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
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
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x108, G_IF_SIGMA_CCVV_CCVV_NO0_X108)
        (sc, ic, sv1, iv1, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.109
  //* S2(w,y,a,c)  <--  (    2.00000000) T2(w,c1,c,v1) V2(v1,a,c1,y) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
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
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x109, G_IF_SIGMA_CCVV_CCVV_NO0_X109)
        (sc, ic, sv1, iv1, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.110
  //* S2(w,y,a,c)  <--  (   -4.00000000) T2(w,c1,c,v1) V2(v1,c1,y,a) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
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
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x110, G_IF_SIGMA_CCVV_CCVV_NO0_X110)
        (sc, ic, sv1, iv1, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.111
  //* S2(w,y,a,c)  <--  (   -4.00000000) T2(c1,w,c,v1) V2(v1,a,c1,y) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
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
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x111, G_IF_SIGMA_CCVV_CCVV_NO0_X111)
        (sc, ic, sv1, iv1, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.112
  //* S2(w,y,a,c)  <--  (    8.00000000) T2(c1,w,v1,a) V2(c,y,c1,v1) 
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
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x112, G_IF_SIGMA_CCVV_CCVV_NO0_X112)
        (sa, ia, sc, ic, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.113
  //* S2(w,y,a,c)  <--  (    2.00000000) T2(y,c1,v1,a) V2(c,w,c1,v1) 
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
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x113, G_IF_SIGMA_CCVV_CCVV_NO0_X113)
        (sa, ia, sc, ic, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.114
  //* S2(w,y,a,c)  <--  (   -4.00000000) T2(y,w,v1,a) Y48(c,v1) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nvir, nvir);
  orz::DTensor Y48 = orz::ct::sympack_Xvv(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y114, G_IF_SIGMA_CCVV_CCVV_Y114)
      (sc1, ic1, V2_sym.cptr(), Y48.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    T2b = T2.get_amp2(ia);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x114, G_IF_SIGMA_CCVV_CCVV_NO0_X114)
        (sa, ia, sc, ic, T2b.cptr(), Y48.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.115
  //* S2(w,y,a,c)  <--  (    8.00000000) T2(c1,y,v1,c) V2(v1,c1,w,a) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
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
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x115, G_IF_SIGMA_CCVV_CCVV_NO0_X115)
        (sc, ic, sv1, iv1, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.116
  //* S2(w,y,a,c)  <--  (    2.00000000) T2(w,c1,v1,c) V2(v1,c1,y,a) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
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
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x116, G_IF_SIGMA_CCVV_CCVV_NO0_X116)
        (sc, ic, sv1, iv1, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.117
  //* S2(w,y,a,c)  <--  (   -4.00000000) T2(w,y,v1,c) Y49(a,v1) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nvir, nvir);
  orz::DTensor Y49 = orz::ct::sympack_Xvv(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y117, G_IF_SIGMA_CCVV_CCVV_Y117)
      (sc1, ic1, V2_sym.cptr(), Y49.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x117, G_IF_SIGMA_CCVV_CCVV_NO0_X117)
      (sc, ic, T2b.cptr(), Y49.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.118
  //* S2(w,y,a,c)  <--  (    8.00000000) T2(w,c1,a,v1) V2(c,y,c1,v1) 
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
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x118, G_IF_SIGMA_CCVV_CCVV_NO0_X118)
        (sc, ic, sv1, iv1, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.119
  //* S2(w,y,a,c)  <--  (   -4.00000000) T2(w,y,a,v1) Y50(c,v1) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nvir, nvir);
  orz::DTensor Y50 = orz::ct::sympack_Xvv(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y119, G_IF_SIGMA_CCVV_CCVV_Y119)
      (sc1, ic1, V2_sym.cptr(), Y50.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x119, G_IF_SIGMA_CCVV_CCVV_NO0_X119)
        (sc, ic, sv1, iv1, T2b.cptr(), Y50.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.120
  //* S2(w,y,a,c)  <--  (    2.00000000) T2(c1,y,a,v1) V2(c,w,c1,v1) 
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
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x120, G_IF_SIGMA_CCVV_CCVV_NO0_X120)
        (sc, ic, sv1, iv1, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.121
  //* S2(w,y,a,c)  <--  (    8.00000000) T2(y,c1,c,v1) V2(v1,c1,w,a) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
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
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x121, G_IF_SIGMA_CCVV_CCVV_NO0_X121)
        (sc, ic, sv1, iv1, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.122
  //* S2(w,y,a,c)  <--  (   -4.00000000) T2(y,w,c,v1) Y51(a,v1) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nvir, nvir);
  orz::DTensor Y51 = orz::ct::sympack_Xvv(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_y122, G_IF_SIGMA_CCVV_CCVV_Y122)
      (sc1, ic1, V2_sym.cptr(), Y51.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x122, G_IF_SIGMA_CCVV_CCVV_NO0_X122)
        (sc, ic, sv1, iv1, T2b.cptr(), Y51.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.123
  //* S2(w,y,a,c)  <--  (    2.00000000) T2(c1,w,c,v1) V2(v1,c1,y,a) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
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
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x123, G_IF_SIGMA_CCVV_CCVV_NO0_X123)
        (sc, ic, sv1, iv1, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.124
  //* X()  <--  (    1.00000000)  D2(o1,o3,o2,o4) V2(o1,o3,o2,o4) 
  //* S2(w,y,a,c)  <--  (    2.00000000) X T2(w,y,a,c) 
  double X = 0;
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x124, G_IF_SIGMA_CCVV_CCVV_NO0_X124)
      (so1, io1, V2_sym.cptr(), &X, nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x124, G_IF_SIGMA_CCVV_CCVV_NO1_X124)
      (sc, ic, &X, T2b.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.125
  //* X()  <--  (    1.00000000)  D2(o1,o3,o2,o4) V2(o1,o3,o2,o4) 
  //* S2(w,y,a,c)  <--  (   -1.00000000) X T2(y,w,a,c) 
  double X = 0;
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x125, G_IF_SIGMA_CCVV_CCVV_NO0_X125)
      (so1, io1, V2_sym.cptr(), &X, nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x125, G_IF_SIGMA_CCVV_CCVV_NO1_X125)
      (sc, ic, &X, T2b.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.126
  //* X()  <--  (    1.00000000)  D2(o1,o3,o2,o4) V2(o1,o3,o2,o4) 
  //* S2(w,y,a,c)  <--  (    2.00000000) X T2(y,w,c,a) 
  double X = 0;
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x126, G_IF_SIGMA_CCVV_CCVV_NO0_X126)
      (so1, io1, V2_sym.cptr(), &X, nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      T2b = T2.get_amp2(ia);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x126, G_IF_SIGMA_CCVV_CCVV_NO1_X126)
        (sa, ia, sc, ic, &X, T2b.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.127
  //* X()  <--  (    1.00000000)  D2(o1,o3,o2,o4) V2(o1,o3,o2,o4) 
  //* S2(w,y,a,c)  <--  (   -1.00000000) X T2(w,y,c,a) 
  double X = 0;
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
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x127, G_IF_SIGMA_CCVV_CCVV_NO0_X127)
      (so1, io1, V2_sym.cptr(), &X, nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      T2b = T2.get_amp2(ia);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x127, G_IF_SIGMA_CCVV_CCVV_NO1_X127)
        (sa, ia, sc, ic, &X, T2b.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.128
  //* X(c,v1)  <--  (    1.00000000)  D1(o1,o2) V2(c,v1,o1,o2) 
  //* S2(w,y,a,c)  <--  (    4.00000000) T2(y,w,v1,a) X(c,v1) 
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
    orz::DTensor X(nvir);
    orz::DTensor Xv = orz::ct::sympack_Xv(symblockinfo, sc, X);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x128, G_IF_SIGMA_CCVV_CCVV_NO0_X128)
      (sc, ic, V2_sym.cptr(), Xv.cptr(), nir, nsym, psym);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      T2b = T2.get_amp2(ia);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x128, G_IF_SIGMA_CCVV_CCVV_NO1_X128)
        (sa, ia, sc, ic, T2b.cptr(), Xv.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.129
  //* X(c,v1)  <--  (    1.00000000)  D1(o1,o2) V2(c,v1,o1,o2) 
  //* S2(w,y,a,c)  <--  (   -2.00000000) T2(w,y,v1,a) X(c,v1) 
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
    orz::DTensor X(nvir);
    orz::DTensor Xv = orz::ct::sympack_Xv(symblockinfo, sc, X);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x129, G_IF_SIGMA_CCVV_CCVV_NO0_X129)
      (sc, ic, V2_sym.cptr(), Xv.cptr(), nir, nsym, psym);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      T2b = T2.get_amp2(ia);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x129, G_IF_SIGMA_CCVV_CCVV_NO1_X129)
        (sa, ia, sc, ic, T2b.cptr(), Xv.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.130
  //* X(c,v1)  <--  (    1.00000000)  D1(o1,o2) V2(c,o1,o2,v1) 
  //* S2(w,y,a,c)  <--  (    1.00000000) T2(w,y,v1,a) X(c,v1) 
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
    orz::DTensor X(nvir);
    orz::DTensor Xv = orz::ct::sympack_Xv(symblockinfo, sc, X);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x130, G_IF_SIGMA_CCVV_CCVV_NO0_X130)
      (sc, ic, V2_sym.cptr(), Xv.cptr(), nir, nsym, psym);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      T2b = T2.get_amp2(ia);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x130, G_IF_SIGMA_CCVV_CCVV_NO1_X130)
        (sa, ia, sc, ic, T2b.cptr(), Xv.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.131
  //* X(a,v1)  <--  (    1.00000000)  D1(o1,o2) V2(a,v1,o1,o2) 
  //* S2(w,y,a,c)  <--  (    4.00000000) T2(w,y,v1,c) X(a,v1) 
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
    orz::DTensor X(nvir);
    orz::DTensor Xv = orz::ct::sympack_Xv(symblockinfo, sa, X);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x131, G_IF_SIGMA_CCVV_CCVV_NO0_X131)
      (sa, ia, V2_sym.cptr(), Xv.cptr(), nir, nsym, psym);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      T2b = T2.get_amp2(ic);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x131, G_IF_SIGMA_CCVV_CCVV_NO1_X131)
        (sa, ia, sc, ic, T2b.cptr(), Xv.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.132
  //* X(a,v1)  <--  (    1.00000000)  D1(o1,o2) V2(a,v1,o1,o2) 
  //* S2(w,y,a,c)  <--  (   -2.00000000) T2(y,w,v1,c) X(a,v1) 
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
    orz::DTensor X(nvir);
    orz::DTensor Xv = orz::ct::sympack_Xv(symblockinfo, sa, X);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x132, G_IF_SIGMA_CCVV_CCVV_NO0_X132)
      (sa, ia, V2_sym.cptr(), Xv.cptr(), nir, nsym, psym);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      T2b = T2.get_amp2(ic);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x132, G_IF_SIGMA_CCVV_CCVV_NO1_X132)
        (sa, ia, sc, ic, T2b.cptr(), Xv.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.133
  //* X(a,v1)  <--  (    1.00000000)  D1(o1,o2) V2(a,o1,o2,v1) 
  //* S2(w,y,a,c)  <--  (    1.00000000) T2(y,w,v1,c) X(a,v1) 
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
    orz::DTensor X(nvir);
    orz::DTensor Xv = orz::ct::sympack_Xv(symblockinfo, sa, X);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x133, G_IF_SIGMA_CCVV_CCVV_NO0_X133)
      (sa, ia, V2_sym.cptr(), Xv.cptr(), nir, nsym, psym);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      T2b = T2.get_amp2(ic);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x133, G_IF_SIGMA_CCVV_CCVV_NO1_X133)
        (sa, ia, sc, ic, T2b.cptr(), Xv.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.134
  //* X(c,v1)  <--  (    1.00000000)  D1(o1,o2) V2(c,v1,o1,o2) 
  //* S2(w,y,a,c)  <--  (    4.00000000) T2(w,y,a,v1) X(c,v1) 
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
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      double X = 0;
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x134, G_IF_SIGMA_CCVV_CCVV_NO0_X134)
        (sc, ic, sv1, iv1, V2_sym.cptr(), &X, nir, nsym, psym);
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x134, G_IF_SIGMA_CCVV_CCVV_NO1_X134)
        (sc, ic, sv1, iv1, T2b.cptr(), &X, S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.135
  //* X(c,v1)  <--  (    1.00000000)  D1(o1,o2) V2(c,v1,o1,o2) 
  //* S2(w,y,a,c)  <--  (   -2.00000000) T2(y,w,a,v1) X(c,v1) 
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
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      double X = 0;
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x135, G_IF_SIGMA_CCVV_CCVV_NO0_X135)
        (sc, ic, sv1, iv1, V2_sym.cptr(), &X, nir, nsym, psym);
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x135, G_IF_SIGMA_CCVV_CCVV_NO1_X135)
        (sc, ic, sv1, iv1, T2b.cptr(), &X, S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.136
  //* X(c,v1)  <--  (    1.00000000)  D1(o1,o2) V2(c,o1,o2,v1) 
  //* S2(w,y,a,c)  <--  (    1.00000000) T2(y,w,a,v1) X(c,v1) 
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
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      double X = 0;
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x136, G_IF_SIGMA_CCVV_CCVV_NO0_X136)
        (sc, ic, sv1, iv1, V2_sym.cptr(), &X, nir, nsym, psym);
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x136, G_IF_SIGMA_CCVV_CCVV_NO1_X136)
        (sc, ic, sv1, iv1, T2b.cptr(), &X, S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.137
  //* X(v1,a)  <--  (    1.00000000)  D1(o1,o2) V2(v1,a,o1,o2) 
  //* S2(w,y,a,c)  <--  (    4.00000000) T2(y,w,c,v1) X(v1,a) 
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
    orz::DTensor X(nvir);
    orz::DTensor Xv = orz::ct::sympack_Xv(symblockinfo, sv1, X);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x137, G_IF_SIGMA_CCVV_CCVV_NO0_X137)
      (sv1, iv1, V2_sym.cptr(), Xv.cptr(), nir, nsym, psym);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x137, G_IF_SIGMA_CCVV_CCVV_NO1_X137)
        (sc, ic, sv1, iv1, T2b.cptr(), Xv.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.138
  //* X(v1,a)  <--  (    1.00000000)  D1(o1,o2) V2(v1,a,o1,o2) 
  //* S2(w,y,a,c)  <--  (   -2.00000000) T2(w,y,c,v1) X(v1,a) 
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
    orz::DTensor X(nvir);
    orz::DTensor Xv = orz::ct::sympack_Xv(symblockinfo, sv1, X);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x138, G_IF_SIGMA_CCVV_CCVV_NO0_X138)
      (sv1, iv1, V2_sym.cptr(), Xv.cptr(), nir, nsym, psym);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x138, G_IF_SIGMA_CCVV_CCVV_NO1_X138)
        (sc, ic, sv1, iv1, T2b.cptr(), Xv.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.139
  //* X(v1,a)  <--  (    1.00000000)  D1(o1,o2) V2(v1,o2,o1,a) 
  //* S2(w,y,a,c)  <--  (    1.00000000) T2(w,y,c,v1) X(v1,a) 
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
    orz::DTensor X(nvir);
    orz::DTensor Xv = orz::ct::sympack_Xv(symblockinfo, sv1, X);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x139, G_IF_SIGMA_CCVV_CCVV_NO0_X139)
      (sv1, iv1, V2_sym.cptr(), Xv.cptr(), nir, nsym, psym);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x139, G_IF_SIGMA_CCVV_CCVV_NO1_X139)
        (sc, ic, sv1, iv1, T2b.cptr(), Xv.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.140
  //* X(c,v1)  <--  (    1.00000000)  D1(o1,o2) V2(c,o1,o2,v1) 
  //* S2(w,y,a,c)  <--  (   -2.00000000) T2(y,w,v1,a) X(c,v1) 
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
    orz::DTensor X(nvir);
    orz::DTensor Xv = orz::ct::sympack_Xv(symblockinfo, sc, X);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x140, G_IF_SIGMA_CCVV_CCVV_NO0_X140)
      (sc, ic, V2_sym.cptr(), Xv.cptr(), nir, nsym, psym);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      T2b = T2.get_amp2(ia);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x140, G_IF_SIGMA_CCVV_CCVV_NO1_X140)
        (sa, ia, sc, ic, T2b.cptr(), Xv.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.141
  //* X(a,v1)  <--  (    1.00000000)  D1(o1,o2) V2(a,o1,o2,v1) 
  //* S2(w,y,a,c)  <--  (   -2.00000000) T2(w,y,v1,c) X(a,v1) 
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
    orz::DTensor X(nvir);
    orz::DTensor Xv = orz::ct::sympack_Xv(symblockinfo, sa, X);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x141, G_IF_SIGMA_CCVV_CCVV_NO0_X141)
      (sa, ia, V2_sym.cptr(), Xv.cptr(), nir, nsym, psym);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      T2b = T2.get_amp2(ic);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x141, G_IF_SIGMA_CCVV_CCVV_NO1_X141)
        (sa, ia, sc, ic, T2b.cptr(), Xv.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.142
  //* X(c,v1)  <--  (    1.00000000)  D1(o1,o2) V2(c,o1,o2,v1) 
  //* S2(w,y,a,c)  <--  (   -2.00000000) T2(w,y,a,v1) X(c,v1) 
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
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      double X = 0;
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x142, G_IF_SIGMA_CCVV_CCVV_NO0_X142)
        (sc, ic, sv1, iv1, V2_sym.cptr(), &X, nir, nsym, psym);
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x142, G_IF_SIGMA_CCVV_CCVV_NO1_X142)
        (sc, ic, sv1, iv1, T2b.cptr(), &X, S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.143
  //* X(v1,a)  <--  (    1.00000000)  D1(o1,o2) V2(v1,o1,o2,a) 
  //* S2(w,y,a,c)  <--  (   -2.00000000) T2(y,w,c,v1) X(v1,a) 
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
    orz::DTensor X(nvir);
    orz::DTensor Xv = orz::ct::sympack_Xv(symblockinfo, sv1, X);
    FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x143, G_IF_SIGMA_CCVV_CCVV_NO0_X143)
      (sv1, iv1, V2_sym.cptr(), Xv.cptr(), nir, nsym, psym);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no1_x143, G_IF_SIGMA_CCVV_CCVV_NO1_X143)
        (sc, ic, sv1, iv1, T2b.cptr(), Xv.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.144
  //* S2(w,y,a,c)  <--  (    4.00000000) T2(w,y,v2,v1) V2(c,v1,a,v2) 
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
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x144, G_IF_SIGMA_CCVV_CCVV_NO0_X144)
        (sc, ic, sv1, iv1, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.145
  //* S2(w,y,a,c)  <--  (   -2.00000000) T2(y,w,v2,v1) V2(c,v1,a,v2) 
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
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x145, G_IF_SIGMA_CCVV_CCVV_NO0_X145)
        (sc, ic, sv1, iv1, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.146
  //* S2(w,y,a,c)  <--  (    4.00000000) T2(y,w,v2,v1) V2(c,v2,a,v1) 
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
    for(int sv1 = 0;sv1 < nir;++sv1){ 
    for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
      T2b = T2.get_amp2(iv1);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x146, G_IF_SIGMA_CCVV_CCVV_NO0_X146)
        (sc, ic, sv1, iv1, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }


  {
  // No.147
  //* S2(w,y,a,c)  <--  (   -2.00000000) T2(w,y,v1,v2) V2(c,v1,a,v2) 
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
    for(int sv2 = 0;sv2 < nir;++sv2){ 
    for(int iv2 = symblockinfo.psym()(sv2,I_V,I_BEGIN);iv2 <= symblockinfo.psym()(sv2,I_V,I_END);++iv2){ 
      T2b = T2.get_amp2(iv2);
      FC_FUNC(g_if_sigma_ccvv_ccvv_no0_x147, G_IF_SIGMA_CCVV_CCVV_NO0_X147)
        (sc, ic, sv2, iv2, T2b.cptr(), V2_sym.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(ic, S2b);
  }
  }
  }

  for(int ssig = 0;ssig < nir;++ssig){                                                                                      
  for(int isig = symblockinfo.psym()(ssig,I_V,I_BEGIN);isig <= symblockinfo.psym()(ssig,I_V,I_END);++isig){                 
    S2b = retval.get_amp2(isig);                                                                    
    FC_FUNC(g_if_sigma_ccvv_scale,G_IF_SIGMA_CCVV_SCALE) // to make hamiltonian symmetric (scale by some factor of overlap) 
      (ssig, isig, S2b.cptr(), nir, nsym, psym);                                                    
    retval.put_amp2(isig, S2b); // S2ija, [b] <<-- Sb                                               
  } // End isig                                                                                                             
  } // End ssig                                                                                                             

  return retval; 
} 
