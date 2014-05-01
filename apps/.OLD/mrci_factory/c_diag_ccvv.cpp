                                                                                
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
#include <sci/ctnew2/c_diag_ccvv.h>                                            
                                                                                
using std::cout;                                                                
using std::endl;                                                                
                                                                                
#define FLOPCOUNT                                                               
                                                                                
//  8888888888                     888                  
//  888                            888                  
//  888                            888                  
//  8888888  .d88b.  88888b.d88b.  888888  .d88b.       
//  888     d8P  Y8b 888 "888 "88b 888    d88""88b  
//  888     88888888 888  888  888 888    888  888      
//  888     Y8b.     888  888  888 Y88b.  Y88..88P      
//  888      "Y8888  888  888  888  "Y888  "Y88P"   

                                                                                
// ***************************************************************************  
// orz::ct::mrci                                                                
// ***************************************************************************  
									    /*!        
   @brief CT input                                                              
									     */        
                                                                                
orz::ct::BareAmpPack orz::ct::diag_ccvv(const orz::ct::Input &ctinp,                                    
					 const orz::ct::SymBlockInfo &symblockinfo,                                
					 const orz::ct::HintMO &hintmo,                                            
					 const orz::ct::RdmPack &rdmPack_sym,                                      
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
  std::string name_of_sigma = "Hdiag" + stm.str() + "]"; // Name of the Sigma vector  
  orz::ct::BareAmpPack retval                                                                                   
    = orz::ct::BareAmpPack(ctinp, symblockinfo, name_of_sigma); // Sigma(a, a', e, e') tensor                   
                                                                                                                
  orz::DTensor Hdiagb; // Container of S2_aae,[b] tensor                                   
                                                                                                                


  {
  // No.0
  //* Hdiag(w,y,a,c)  <--  (    8.00000000) Y0 
  // The effective tensor is detected .... 
  double Y0 = 0;
  FC_FUNC(g_if_diag_ccvv_y0, G_IF_DIAG_CCVV_Y0)
    (moint1_sym.cptr(), &Y0, nir, nsym, psym);
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
    FC_FUNC(g_if_diag_ccvv_no0_x0, G_IF_DIAG_CCVV_NO0_X0)
      (sc, ic, &Y0, Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.1
  //* Hdiag(w,w,a,c)  <--  (   -4.00000000) Y1 
  // The effective tensor is detected .... 
  double Y1 = 0;
  FC_FUNC(g_if_diag_ccvv_y1, G_IF_DIAG_CCVV_Y1)
    (moint1_sym.cptr(), &Y1, nir, nsym, psym);
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
    FC_FUNC(g_if_diag_ccvv_no0_x1, G_IF_DIAG_CCVV_NO0_X1)
      (sc, ic, &Y1, Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.2
  //* Hdiag(w,y,a,c)  <--  (   -4.00000000) h(w,w) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
    FC_FUNC(g_if_diag_ccvv_no0_x2, G_IF_DIAG_CCVV_NO0_X2)
      (sc, ic, moint1_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.3
  //* Hdiag(w,w,a,c)  <--  (    4.00000000) h(w,w) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
    FC_FUNC(g_if_diag_ccvv_no0_x3, G_IF_DIAG_CCVV_NO0_X3)
      (sc, ic, moint1_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.4
  //* Hdiag(w,y,a,c)  <--  (   -4.00000000) h(y,y) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
    FC_FUNC(g_if_diag_ccvv_no0_x4, G_IF_DIAG_CCVV_NO0_X4)
      (sc, ic, moint1_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.5
  //* Hdiag(w,w,a,a)  <--  (    8.00000000) Y2 
  // The effective tensor is detected .... 
  double Y2 = 0;
  FC_FUNC(g_if_diag_ccvv_y5, G_IF_DIAG_CCVV_Y5)
    (moint1_sym.cptr(), &Y2, nir, nsym, psym);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ccvv_no0_x5, G_IF_DIAG_CCVV_NO0_X5)
      (sa, ia, &Y2, Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.6
  //* Hdiag(w,y,a,a)  <--  (   -4.00000000) Y3 
  // The effective tensor is detected .... 
  double Y3 = 0;
  FC_FUNC(g_if_diag_ccvv_y6, G_IF_DIAG_CCVV_Y6)
    (moint1_sym.cptr(), &Y3, nir, nsym, psym);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ccvv_no0_x6, G_IF_DIAG_CCVV_NO0_X6)
      (sa, ia, &Y3, Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.7
  //* Hdiag(w,w,a,a)  <--  (   -8.00000000) h(w,w) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ccvv_no0_x7, G_IF_DIAG_CCVV_NO0_X7)
      (sa, ia, moint1_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.8
  //* Hdiag(w,y,a,a)  <--  (    2.00000000) h(w,w) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ccvv_no0_x8, G_IF_DIAG_CCVV_NO0_X8)
      (sa, ia, moint1_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.9
  //* Hdiag(w,y,a,a)  <--  (    2.00000000) h(y,y) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ccvv_no0_x9, G_IF_DIAG_CCVV_NO0_X9)
      (sa, ia, moint1_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.10
  //* Hdiag(w,y,a,c)  <--  (    4.00000000) D1(o1,o2) h(o1,o2) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
    FC_FUNC(g_if_diag_ccvv_no0_x10, G_IF_DIAG_CCVV_NO0_X10)
      (sc, ic, moint1_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.11
  //* Hdiag(w,w,a,c)  <--  (   -2.00000000) D1(o1,o2) h(o2,o1) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
    FC_FUNC(g_if_diag_ccvv_no0_x11, G_IF_DIAG_CCVV_NO0_X11)
      (sc, ic, moint1_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.12
  //* Hdiag(w,w,a,a)  <--  (    4.00000000) D1(o1,o2) h(o1,o2) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ccvv_no0_x12, G_IF_DIAG_CCVV_NO0_X12)
      (sa, ia, moint1_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.13
  //* Hdiag(w,y,a,a)  <--  (   -2.00000000) D1(o1,o2) h(o1,o2) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ccvv_no0_x13, G_IF_DIAG_CCVV_NO0_X13)
      (sa, ia, moint1_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.14
  //* Hdiag(w,w,a,a)  <--  (    8.00000000) h(a,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ccvv_no0_x14, G_IF_DIAG_CCVV_NO0_X14)
      (sa, ia, moint1_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.15
  //* Hdiag(w,y,a,a)  <--  (   -4.00000000) h(a,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ccvv_no0_x15, G_IF_DIAG_CCVV_NO0_X15)
      (sa, ia, moint1_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.16
  //* Hdiag(w,y,a,c)  <--  (    4.00000000) h(a,a) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
    FC_FUNC(g_if_diag_ccvv_no0_x16, G_IF_DIAG_CCVV_NO0_X16)
      (sc, ic, moint1_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.17
  //* Hdiag(w,w,a,c)  <--  (   -2.00000000) h(a,a) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
    FC_FUNC(g_if_diag_ccvv_no0_x17, G_IF_DIAG_CCVV_NO0_X17)
      (sc, ic, moint1_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.18
  //* Hdiag(w,w,a,c)  <--  (   -2.00000000) h(c,c) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
    FC_FUNC(g_if_diag_ccvv_no0_x18, G_IF_DIAG_CCVV_NO0_X18)
      (sc, ic, moint1_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.19
  //* Hdiag(w,w,a,c)  <--  (   -4.00000000) Y4 
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
    FC_FUNC(g_if_diag_ccvv_y19, G_IF_DIAG_CCVV_Y19)
      (sc1, ic1, V2_sym.cptr(), &Y4, nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
    FC_FUNC(g_if_diag_ccvv_no0_x19, G_IF_DIAG_CCVV_NO0_X19)
      (sc, ic, &Y4, Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.20
  //* Hdiag(w,w,a,c)  <--  (    2.00000000) Y5 
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
    FC_FUNC(g_if_diag_ccvv_y20, G_IF_DIAG_CCVV_Y20)
      (sc1, ic1, V2_sym.cptr(), &Y5, nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
    FC_FUNC(g_if_diag_ccvv_no0_x20, G_IF_DIAG_CCVV_NO0_X20)
      (sc, ic, &Y5, Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.21
  //* Hdiag(w,w,a,c)  <--  (   -4.00000000) Y6(a,a) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nvir, nvir);
  orz::DTensor Y6 = orz::ct::sympack_Xvv(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_diag_ccvv_y21, G_IF_DIAG_CCVV_Y21)
      (sc1, ic1, V2_sym.cptr(), Y6.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
    FC_FUNC(g_if_diag_ccvv_no0_x21, G_IF_DIAG_CCVV_NO0_X21)
      (sc, ic, Y6.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.22
  //* Hdiag(w,w,a,c)  <--  (   -4.00000000) Y7(c,c) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nvir, nvir);
  orz::DTensor Y7 = orz::ct::sympack_Xvv(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_diag_ccvv_y22, G_IF_DIAG_CCVV_Y22)
      (sc1, ic1, V2_sym.cptr(), Y7.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
    FC_FUNC(g_if_diag_ccvv_no0_x22, G_IF_DIAG_CCVV_NO0_X22)
      (sc, ic, Y7.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.23
  //* Hdiag(w,w,a,c)  <--  (    2.00000000) Y8(a,a) 
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
    FC_FUNC(g_if_diag_ccvv_y23, G_IF_DIAG_CCVV_Y23)
      (sc1, ic1, V2_sym.cptr(), Y8.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
    FC_FUNC(g_if_diag_ccvv_no0_x23, G_IF_DIAG_CCVV_NO0_X23)
      (sc, ic, Y8.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.24
  //* Hdiag(w,w,a,c)  <--  (    2.00000000) Y9(c,c) 
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
    FC_FUNC(g_if_diag_ccvv_y24, G_IF_DIAG_CCVV_Y24)
      (sc1, ic1, V2_sym.cptr(), Y9.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
    FC_FUNC(g_if_diag_ccvv_no0_x24, G_IF_DIAG_CCVV_NO0_X24)
      (sc, ic, Y9.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.25
  //* Hdiag(w,w,a,c)  <--  (   -2.00000000) V2(c,c,a,a) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
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
    FC_FUNC(g_if_diag_ccvv_no0_x25, G_IF_DIAG_CCVV_NO0_X25)
      (sc, ic, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.26
  //* Hdiag(w,w,a,c)  <--  (    4.00000000) V2(c,a,a,c) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
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
    FC_FUNC(g_if_diag_ccvv_no0_x26, G_IF_DIAG_CCVV_NO0_X26)
      (sc, ic, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.27
  //* Hdiag(w,y,a,c)  <--  (    4.00000000) h(c,c) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
    FC_FUNC(g_if_diag_ccvv_no0_x27, G_IF_DIAG_CCVV_NO0_X27)
      (sc, ic, moint1_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.28
  //* Hdiag(w,y,a,c)  <--  (    8.00000000) Y10 
  // The effective tensor is detected .... 
  double Y10 = 0;
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
    FC_FUNC(g_if_diag_ccvv_y28, G_IF_DIAG_CCVV_Y28)
      (sc1, ic1, V2_sym.cptr(), &Y10, nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
    FC_FUNC(g_if_diag_ccvv_no0_x28, G_IF_DIAG_CCVV_NO0_X28)
      (sc, ic, &Y10, Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.29
  //* Hdiag(w,y,a,c)  <--  (   -8.00000000) Y11(w,w) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nclosed, nclosed);
  orz::DTensor Y11 = orz::ct::sympack_Xcc(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_diag_ccvv_y29, G_IF_DIAG_CCVV_Y29)
      (sc1, ic1, V2_sym.cptr(), Y11.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
    FC_FUNC(g_if_diag_ccvv_no0_x29, G_IF_DIAG_CCVV_NO0_X29)
      (sc, ic, Y11.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.30
  //* Hdiag(w,w,a,c)  <--  (    8.00000000) Y12(w,w) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nclosed, nclosed);
  orz::DTensor Y12 = orz::ct::sympack_Xcc(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_diag_ccvv_y30, G_IF_DIAG_CCVV_Y30)
      (sc1, ic1, V2_sym.cptr(), Y12.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
    FC_FUNC(g_if_diag_ccvv_no0_x30, G_IF_DIAG_CCVV_NO0_X30)
      (sc, ic, Y12.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.31
  //* Hdiag(w,y,a,c)  <--  (   -8.00000000) Y13(y,y) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nclosed, nclosed);
  orz::DTensor Y13 = orz::ct::sympack_Xcc(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_diag_ccvv_y31, G_IF_DIAG_CCVV_Y31)
      (sc1, ic1, V2_sym.cptr(), Y13.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
    FC_FUNC(g_if_diag_ccvv_no0_x31, G_IF_DIAG_CCVV_NO0_X31)
      (sc, ic, Y13.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.32
  //* Hdiag(w,y,a,c)  <--  (   -4.00000000) Y14 
  // The effective tensor is detected .... 
  double Y14 = 0;
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
    FC_FUNC(g_if_diag_ccvv_y32, G_IF_DIAG_CCVV_Y32)
      (sc1, ic1, V2_sym.cptr(), &Y14, nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
    FC_FUNC(g_if_diag_ccvv_no0_x32, G_IF_DIAG_CCVV_NO0_X32)
      (sc, ic, &Y14, Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.33
  //* Hdiag(w,y,a,c)  <--  (    4.00000000) Y15(w,w) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nclosed, nclosed);
  orz::DTensor Y15 = orz::ct::sympack_Xcc(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_diag_ccvv_y33, G_IF_DIAG_CCVV_Y33)
      (sc1, ic1, V2_sym.cptr(), Y15.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
    FC_FUNC(g_if_diag_ccvv_no0_x33, G_IF_DIAG_CCVV_NO0_X33)
      (sc, ic, Y15.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.34
  //* Hdiag(w,w,a,c)  <--  (   -4.00000000) Y16(w,w) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nclosed, nclosed);
  orz::DTensor Y16 = orz::ct::sympack_Xcc(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_diag_ccvv_y34, G_IF_DIAG_CCVV_Y34)
      (sc1, ic1, V2_sym.cptr(), Y16.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
    FC_FUNC(g_if_diag_ccvv_no0_x34, G_IF_DIAG_CCVV_NO0_X34)
      (sc, ic, Y16.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.35
  //* Hdiag(w,y,a,c)  <--  (    4.00000000) Y17(y,y) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nclosed, nclosed);
  orz::DTensor Y17 = orz::ct::sympack_Xcc(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_diag_ccvv_y35, G_IF_DIAG_CCVV_Y35)
      (sc1, ic1, V2_sym.cptr(), Y17.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
    FC_FUNC(g_if_diag_ccvv_no0_x35, G_IF_DIAG_CCVV_NO0_X35)
      (sc, ic, Y17.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.36
  //* Hdiag(w,y,a,c)  <--  (    4.00000000) V2(w,w,y,y) 
  for(int sw = 0;sw < nir;++sw){ 
  for(int iw = symblockinfo.psym()(sw,I_C,I_BEGIN);iw <= symblockinfo.psym()(sw,I_C,I_END);++iw){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iw);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iw, sw, V2); // V2=(IR-COV index) 
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
      FC_FUNC(g_if_diag_ccvv_no0_x36, G_IF_DIAG_CCVV_NO0_X36)
        (sc, ic, sw, iw, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.37
  //* Hdiag(w,y,a,c)  <--  (   -2.00000000) V2(w,y,w,y) 
  for(int sw = 0;sw < nir;++sw){ 
  for(int iw = symblockinfo.psym()(sw,I_C,I_BEGIN);iw <= symblockinfo.psym()(sw,I_C,I_END);++iw){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iw);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iw, sw, V2); // V2=(IR-COV index) 
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
      FC_FUNC(g_if_diag_ccvv_no0_x37, G_IF_DIAG_CCVV_NO0_X37)
        (sc, ic, sw, iw, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.38
  //* Hdiag(w,y,a,c)  <--  (    8.00000000) Y18(a,a) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nvir, nvir);
  orz::DTensor Y18 = orz::ct::sympack_Xvv(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_diag_ccvv_y38, G_IF_DIAG_CCVV_Y38)
      (sc1, ic1, V2_sym.cptr(), Y18.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
    FC_FUNC(g_if_diag_ccvv_no0_x38, G_IF_DIAG_CCVV_NO0_X38)
      (sc, ic, Y18.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.39
  //* Hdiag(w,y,a,c)  <--  (   -4.00000000) V2(a,a,w,w) 
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
      Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
      FC_FUNC(g_if_diag_ccvv_no0_x39, G_IF_DIAG_CCVV_NO0_X39)
        (sa, ia, sc, ic, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.40
  //* Hdiag(w,y,a,c)  <--  (   -4.00000000) V2(a,a,y,y) 
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
      Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
      FC_FUNC(g_if_diag_ccvv_no0_x40, G_IF_DIAG_CCVV_NO0_X40)
        (sa, ia, sc, ic, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.41
  //* Hdiag(w,y,a,c)  <--  (    8.00000000) Y19(c,c) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nvir, nvir);
  orz::DTensor Y19 = orz::ct::sympack_Xvv(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_diag_ccvv_y41, G_IF_DIAG_CCVV_Y41)
      (sc1, ic1, V2_sym.cptr(), Y19.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
    FC_FUNC(g_if_diag_ccvv_no0_x41, G_IF_DIAG_CCVV_NO0_X41)
      (sc, ic, Y19.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.42
  //* Hdiag(w,y,a,c)  <--  (   -4.00000000) V2(c,c,w,w) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
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
    FC_FUNC(g_if_diag_ccvv_no0_x42, G_IF_DIAG_CCVV_NO0_X42)
      (sc, ic, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.43
  //* Hdiag(w,y,a,c)  <--  (   -4.00000000) V2(c,c,y,y) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
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
    FC_FUNC(g_if_diag_ccvv_no0_x43, G_IF_DIAG_CCVV_NO0_X43)
      (sc, ic, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.44
  //* Hdiag(w,y,a,c)  <--  (    8.00000000) V2(a,w,w,a) 
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
      Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
      FC_FUNC(g_if_diag_ccvv_no0_x44, G_IF_DIAG_CCVV_NO0_X44)
        (sa, ia, sc, ic, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.45
  //* Hdiag(w,y,a,c)  <--  (   -4.00000000) Y20(a,a) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nvir, nvir);
  orz::DTensor Y20 = orz::ct::sympack_Xvv(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_diag_ccvv_y45, G_IF_DIAG_CCVV_Y45)
      (sc1, ic1, V2_sym.cptr(), Y20.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
    FC_FUNC(g_if_diag_ccvv_no0_x45, G_IF_DIAG_CCVV_NO0_X45)
      (sc, ic, Y20.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.46
  //* Hdiag(w,y,a,c)  <--  (    2.00000000) V2(a,y,y,a) 
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
      Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
      FC_FUNC(g_if_diag_ccvv_no0_x46, G_IF_DIAG_CCVV_NO0_X46)
        (sa, ia, sc, ic, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.47
  //* Hdiag(w,y,a,c)  <--  (    8.00000000) V2(c,y,y,c) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
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
    FC_FUNC(g_if_diag_ccvv_no0_x47, G_IF_DIAG_CCVV_NO0_X47)
      (sc, ic, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.48
  //* Hdiag(w,y,a,c)  <--  (    2.00000000) V2(c,w,w,c) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
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
    FC_FUNC(g_if_diag_ccvv_no0_x48, G_IF_DIAG_CCVV_NO0_X48)
      (sc, ic, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.49
  //* Hdiag(w,y,a,c)  <--  (   -4.00000000) Y21(c,c) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nvir, nvir);
  orz::DTensor Y21 = orz::ct::sympack_Xvv(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_diag_ccvv_y49, G_IF_DIAG_CCVV_Y49)
      (sc1, ic1, V2_sym.cptr(), Y21.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
    FC_FUNC(g_if_diag_ccvv_no0_x49, G_IF_DIAG_CCVV_NO0_X49)
      (sc, ic, Y21.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.50
  //* Hdiag(w,y,a,c)  <--  (    4.00000000) V2(c,c,a,a) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
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
    FC_FUNC(g_if_diag_ccvv_no0_x50, G_IF_DIAG_CCVV_NO0_X50)
      (sc, ic, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.51
  //* Hdiag(w,y,a,c)  <--  (   -2.00000000) V2(c,a,a,c) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
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
    FC_FUNC(g_if_diag_ccvv_no0_x51, G_IF_DIAG_CCVV_NO0_X51)
      (sc, ic, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.52
  //* Hdiag(w,w,a,a)  <--  (    8.00000000) Y22 
  // The effective tensor is detected .... 
  double Y22 = 0;
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
    FC_FUNC(g_if_diag_ccvv_y52, G_IF_DIAG_CCVV_Y52)
      (sc1, ic1, V2_sym.cptr(), &Y22, nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ccvv_no0_x52, G_IF_DIAG_CCVV_NO0_X52)
      (sa, ia, &Y22, Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.53
  //* Hdiag(w,y,a,a)  <--  (   -4.00000000) Y23 
  // The effective tensor is detected .... 
  double Y23 = 0;
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
    FC_FUNC(g_if_diag_ccvv_y53, G_IF_DIAG_CCVV_Y53)
      (sc1, ic1, V2_sym.cptr(), &Y23, nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ccvv_no0_x53, G_IF_DIAG_CCVV_NO0_X53)
      (sa, ia, &Y23, Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.54
  //* Hdiag(w,w,a,a)  <--  (  -16.00000000) Y24(w,w) 
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
    FC_FUNC(g_if_diag_ccvv_y54, G_IF_DIAG_CCVV_Y54)
      (sc1, ic1, V2_sym.cptr(), Y24.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ccvv_no0_x54, G_IF_DIAG_CCVV_NO0_X54)
      (sa, ia, Y24.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.55
  //* Hdiag(w,y,a,a)  <--  (    4.00000000) Y25(w,w) 
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
    FC_FUNC(g_if_diag_ccvv_y55, G_IF_DIAG_CCVV_Y55)
      (sc1, ic1, V2_sym.cptr(), Y25.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ccvv_no0_x55, G_IF_DIAG_CCVV_NO0_X55)
      (sa, ia, Y25.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.56
  //* Hdiag(w,y,a,a)  <--  (    4.00000000) Y26(y,y) 
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
    FC_FUNC(g_if_diag_ccvv_y56, G_IF_DIAG_CCVV_Y56)
      (sc1, ic1, V2_sym.cptr(), Y26.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ccvv_no0_x56, G_IF_DIAG_CCVV_NO0_X56)
      (sa, ia, Y26.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.57
  //* Hdiag(w,w,a,a)  <--  (   -4.00000000) Y27 
  // The effective tensor is detected .... 
  double Y27 = 0;
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
    FC_FUNC(g_if_diag_ccvv_y57, G_IF_DIAG_CCVV_Y57)
      (sc1, ic1, V2_sym.cptr(), &Y27, nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ccvv_no0_x57, G_IF_DIAG_CCVV_NO0_X57)
      (sa, ia, &Y27, Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.58
  //* Hdiag(w,y,a,a)  <--  (    2.00000000) Y28 
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
    FC_FUNC(g_if_diag_ccvv_y58, G_IF_DIAG_CCVV_Y58)
      (sc1, ic1, V2_sym.cptr(), &Y28, nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ccvv_no0_x58, G_IF_DIAG_CCVV_NO0_X58)
      (sa, ia, &Y28, Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.59
  //* Hdiag(w,y,a,a)  <--  (   -2.00000000) Y29(w,w) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nclosed, nclosed);
  orz::DTensor Y29 = orz::ct::sympack_Xcc(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_diag_ccvv_y59, G_IF_DIAG_CCVV_Y59)
      (sc1, ic1, V2_sym.cptr(), Y29.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ccvv_no0_x59, G_IF_DIAG_CCVV_NO0_X59)
      (sa, ia, Y29.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.60
  //* Hdiag(w,y,a,a)  <--  (   -2.00000000) Y30(y,y) 
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
    FC_FUNC(g_if_diag_ccvv_y60, G_IF_DIAG_CCVV_Y60)
      (sc1, ic1, V2_sym.cptr(), Y30.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ccvv_no0_x60, G_IF_DIAG_CCVV_NO0_X60)
      (sa, ia, Y30.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.61
  //* Hdiag(w,y,a,a)  <--  (    4.00000000) V2(w,y,w,y) 
  for(int sw = 0;sw < nir;++sw){ 
  for(int iw = symblockinfo.psym()(sw,I_C,I_BEGIN);iw <= symblockinfo.psym()(sw,I_C,I_END);++iw){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iw);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iw, sw, V2); // V2=(IR-COV index) 
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
      FC_FUNC(g_if_diag_ccvv_no0_x61, G_IF_DIAG_CCVV_NO0_X61)
        (sa, ia, sw, iw, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(ia, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.62
  //* Hdiag(w,y,a,a)  <--  (   -2.00000000) V2(w,w,y,y) 
  for(int sw = 0;sw < nir;++sw){ 
  for(int iw = symblockinfo.psym()(sw,I_C,I_BEGIN);iw <= symblockinfo.psym()(sw,I_C,I_END);++iw){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iw);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iw, sw, V2); // V2=(IR-COV index) 
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
      FC_FUNC(g_if_diag_ccvv_no0_x62, G_IF_DIAG_CCVV_NO0_X62)
        (sa, ia, sw, iw, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(ia, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.63
  //* Hdiag(w,w,a,a)  <--  (    8.00000000) Y31(w,w) 
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
    FC_FUNC(g_if_diag_ccvv_y63, G_IF_DIAG_CCVV_Y63)
      (sc1, ic1, V2_sym.cptr(), Y31.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ccvv_no0_x63, G_IF_DIAG_CCVV_NO0_X63)
      (sa, ia, Y31.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.64
  //* Hdiag(w,y,a,c)  <--  (    8.00000000) D1(o1,o2) Y32(o1,o2) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nocc, nocc);
  orz::DTensor Y32 = orz::ct::sympack_Xaa(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_diag_ccvv_y64, G_IF_DIAG_CCVV_Y64)
      (sc1, ic1, V2_sym.cptr(), Y32.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
    FC_FUNC(g_if_diag_ccvv_no0_x64, G_IF_DIAG_CCVV_NO0_X64)
      (sc, ic, Y32.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.65
  //* Hdiag(w,w,a,c)  <--  (   -4.00000000) D1(o1,o2) Y33(o1,o2) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nocc, nocc);
  orz::DTensor Y33 = orz::ct::sympack_Xaa(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_diag_ccvv_y65, G_IF_DIAG_CCVV_Y65)
      (sc1, ic1, V2_sym.cptr(), Y33.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
    FC_FUNC(g_if_diag_ccvv_no0_x65, G_IF_DIAG_CCVV_NO0_X65)
      (sc, ic, Y33.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.66
  //* Hdiag(w,y,a,c)  <--  (   -4.00000000) D1(o1,o2) V2(w,w,o1,o2) 
  for(int sw = 0;sw < nir;++sw){ 
  for(int iw = symblockinfo.psym()(sw,I_C,I_BEGIN);iw <= symblockinfo.psym()(sw,I_C,I_END);++iw){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iw);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iw, sw, V2); // V2=(IR-COV index) 
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
      FC_FUNC(g_if_diag_ccvv_no0_x66, G_IF_DIAG_CCVV_NO0_X66)
        (sc, ic, sw, iw, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.67
  //* Hdiag(w,w,a,c)  <--  (    4.00000000) D1(o1,o2) V2(w,w,o1,o2) 
  for(int sw = 0;sw < nir;++sw){ 
  for(int iw = symblockinfo.psym()(sw,I_C,I_BEGIN);iw <= symblockinfo.psym()(sw,I_C,I_END);++iw){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iw);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iw, sw, V2); // V2=(IR-COV index) 
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
      FC_FUNC(g_if_diag_ccvv_no0_x67, G_IF_DIAG_CCVV_NO0_X67)
        (sc, ic, sw, iw, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.68
  //* Hdiag(w,y,a,c)  <--  (   -4.00000000) D1(o1,o2) V2(y,y,o1,o2) 
  for(int sy = 0;sy < nir;++sy){ 
  for(int iy = symblockinfo.psym()(sy,I_C,I_BEGIN);iy <= symblockinfo.psym()(sy,I_C,I_END);++iy){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iy);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iy, sy, V2); // V2=(IR-COV index) 
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
      FC_FUNC(g_if_diag_ccvv_no0_x68, G_IF_DIAG_CCVV_NO0_X68)
        (sc, ic, sy, iy, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.69
  //* Hdiag(w,w,a,a)  <--  (    4.00000000) D1(o1,o2) Y34(o1,o2) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nocc, nocc);
  orz::DTensor Y34 = orz::ct::sympack_Xaa(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_diag_ccvv_y69, G_IF_DIAG_CCVV_Y69)
      (sc1, ic1, V2_sym.cptr(), Y34.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ccvv_no0_x69, G_IF_DIAG_CCVV_NO0_X69)
      (sa, ia, Y34.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.70
  //* Hdiag(w,y,a,a)  <--  (   -2.00000000) D1(o1,o2) Y35(o1,o2) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nocc, nocc);
  orz::DTensor Y35 = orz::ct::sympack_Xaa(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_diag_ccvv_y70, G_IF_DIAG_CCVV_Y70)
      (sc1, ic1, V2_sym.cptr(), Y35.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ccvv_no0_x70, G_IF_DIAG_CCVV_NO0_X70)
      (sa, ia, Y35.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.71
  //* Hdiag(w,y,a,a)  <--  (    1.00000000) D1(o1,o2) V2(w,w,o1,o2) 
  for(int sw = 0;sw < nir;++sw){ 
  for(int iw = symblockinfo.psym()(sw,I_C,I_BEGIN);iw <= symblockinfo.psym()(sw,I_C,I_END);++iw){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iw);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iw, sw, V2); // V2=(IR-COV index) 
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
      FC_FUNC(g_if_diag_ccvv_no0_x71, G_IF_DIAG_CCVV_NO0_X71)
        (sa, ia, sw, iw, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(ia, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.72
  //* Hdiag(w,y,a,a)  <--  (    1.00000000) D1(o1,o2) V2(y,y,o1,o2) 
  for(int sy = 0;sy < nir;++sy){ 
  for(int iy = symblockinfo.psym()(sy,I_C,I_BEGIN);iy <= symblockinfo.psym()(sy,I_C,I_END);++iy){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iy);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iy, sy, V2); // V2=(IR-COV index) 
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
      FC_FUNC(g_if_diag_ccvv_no0_x72, G_IF_DIAG_CCVV_NO0_X72)
        (sa, ia, sy, iy, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(ia, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.73
  //* Hdiag(w,w,a,a)  <--  (   -4.00000000) D1(o1,o2) V2(w,w,o1,o2) 
  for(int sw = 0;sw < nir;++sw){ 
  for(int iw = symblockinfo.psym()(sw,I_C,I_BEGIN);iw <= symblockinfo.psym()(sw,I_C,I_END);++iw){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iw);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iw, sw, V2); // V2=(IR-COV index) 
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
      FC_FUNC(g_if_diag_ccvv_no0_x73, G_IF_DIAG_CCVV_NO0_X73)
        (sa, ia, sw, iw, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(ia, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.74
  //* Hdiag(w,w,a,c)  <--  (   -2.00000000) D1(o1,o2) V2(w,o2,w,o1) 
  for(int sw = 0;sw < nir;++sw){ 
  for(int iw = symblockinfo.psym()(sw,I_C,I_BEGIN);iw <= symblockinfo.psym()(sw,I_C,I_END);++iw){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iw);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iw, sw, V2); // V2=(IR-COV index) 
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
      FC_FUNC(g_if_diag_ccvv_no0_x74, G_IF_DIAG_CCVV_NO0_X74)
        (sc, ic, sw, iw, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.75
  //* Hdiag(w,y,a,c)  <--  (    2.00000000) D1(o1,o2) V2(w,o1,w,o2) 
  for(int sw = 0;sw < nir;++sw){ 
  for(int iw = symblockinfo.psym()(sw,I_C,I_BEGIN);iw <= symblockinfo.psym()(sw,I_C,I_END);++iw){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iw);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iw, sw, V2); // V2=(IR-COV index) 
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
      FC_FUNC(g_if_diag_ccvv_no0_x75, G_IF_DIAG_CCVV_NO0_X75)
        (sc, ic, sw, iw, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.76
  //* Hdiag(w,y,a,c)  <--  (    2.00000000) D1(o1,o2) V2(y,o1,y,o2) 
  for(int sy = 0;sy < nir;++sy){ 
  for(int iy = symblockinfo.psym()(sy,I_C,I_BEGIN);iy <= symblockinfo.psym()(sy,I_C,I_END);++iy){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iy);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iy, sy, V2); // V2=(IR-COV index) 
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
      FC_FUNC(g_if_diag_ccvv_no0_x76, G_IF_DIAG_CCVV_NO0_X76)
        (sc, ic, sy, iy, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.77
  //* Hdiag(w,y,a,c)  <--  (   -4.00000000) D1(o1,o2) Y36(o1,o2) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nocc, nocc);
  orz::DTensor Y36 = orz::ct::sympack_Xaa(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_diag_ccvv_y77, G_IF_DIAG_CCVV_Y77)
      (sc1, ic1, V2_sym.cptr(), Y36.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
    FC_FUNC(g_if_diag_ccvv_no0_x77, G_IF_DIAG_CCVV_NO0_X77)
      (sc, ic, Y36.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.78
  //* Hdiag(w,y,a,c)  <--  (    4.00000000) D1(o1,o2) V2(a,a,o1,o2) 
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
      Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
      FC_FUNC(g_if_diag_ccvv_no0_x78, G_IF_DIAG_CCVV_NO0_X78)
        (sa, ia, sc, ic, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.79
  //* Hdiag(w,y,a,c)  <--  (    4.00000000) D1(o1,o2) V2(c,c,o1,o2) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
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
    FC_FUNC(g_if_diag_ccvv_no0_x79, G_IF_DIAG_CCVV_NO0_X79)
      (sc, ic, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.80
  //* Hdiag(w,y,a,c)  <--  (   -2.00000000) D1(o1,o2) V2(a,o1,o2,a) 
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
      Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
      FC_FUNC(g_if_diag_ccvv_no0_x80, G_IF_DIAG_CCVV_NO0_X80)
        (sa, ia, sc, ic, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.81
  //* Hdiag(w,y,a,c)  <--  (   -2.00000000) D1(o1,o2) V2(c,o1,o2,c) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
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
    FC_FUNC(g_if_diag_ccvv_no0_x81, G_IF_DIAG_CCVV_NO0_X81)
      (sc, ic, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.82
  //* Hdiag(w,w,a,c)  <--  (    2.00000000) D1(o1,o2) Y37(o1,o2) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nocc, nocc);
  orz::DTensor Y37 = orz::ct::sympack_Xaa(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_diag_ccvv_y82, G_IF_DIAG_CCVV_Y82)
      (sc1, ic1, V2_sym.cptr(), Y37.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
    FC_FUNC(g_if_diag_ccvv_no0_x82, G_IF_DIAG_CCVV_NO0_X82)
      (sc, ic, Y37.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.83
  //* Hdiag(w,w,a,c)  <--  (   -2.00000000) D1(o1,o2) V2(a,a,o1,o2) 
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
      Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
      FC_FUNC(g_if_diag_ccvv_no0_x83, G_IF_DIAG_CCVV_NO0_X83)
        (sa, ia, sc, ic, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.84
  //* Hdiag(w,w,a,c)  <--  (   -2.00000000) D1(o1,o2) V2(c,c,o1,o2) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
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
    FC_FUNC(g_if_diag_ccvv_no0_x84, G_IF_DIAG_CCVV_NO0_X84)
      (sc, ic, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.85
  //* Hdiag(w,w,a,c)  <--  (    1.00000000) D1(o1,o2) V2(a,o1,o2,a) 
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
      Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
      FC_FUNC(g_if_diag_ccvv_no0_x85, G_IF_DIAG_CCVV_NO0_X85)
        (sa, ia, sc, ic, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.86
  //* Hdiag(w,w,a,c)  <--  (    1.00000000) D1(o1,o2) V2(c,o1,o2,c) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
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
    FC_FUNC(g_if_diag_ccvv_no0_x86, G_IF_DIAG_CCVV_NO0_X86)
      (sc, ic, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.87
  //* Hdiag(w,y,a,a)  <--  (   -0.50000000) D1(o1,o2) V2(y,o1,y,o2) 
  for(int sy = 0;sy < nir;++sy){ 
  for(int iy = symblockinfo.psym()(sy,I_C,I_BEGIN);iy <= symblockinfo.psym()(sy,I_C,I_END);++iy){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iy);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iy, sy, V2); // V2=(IR-COV index) 
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
      FC_FUNC(g_if_diag_ccvv_no0_x87, G_IF_DIAG_CCVV_NO0_X87)
        (sa, ia, sy, iy, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(ia, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.88
  //* Hdiag(w,w,a,a)  <--  (    2.00000000) D1(o1,o2) V2(w,o1,w,o2) 
  for(int sw = 0;sw < nir;++sw){ 
  for(int iw = symblockinfo.psym()(sw,I_C,I_BEGIN);iw <= symblockinfo.psym()(sw,I_C,I_END);++iw){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iw);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iw, sw, V2); // V2=(IR-COV index) 
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
      FC_FUNC(g_if_diag_ccvv_no0_x88, G_IF_DIAG_CCVV_NO0_X88)
        (sa, ia, sw, iw, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(ia, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.89
  //* Hdiag(w,w,a,a)  <--  (   -2.00000000) D1(o1,o2) Y38(o1,o2) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nocc, nocc);
  orz::DTensor Y38 = orz::ct::sympack_Xaa(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_diag_ccvv_y89, G_IF_DIAG_CCVV_Y89)
      (sc1, ic1, V2_sym.cptr(), Y38.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ccvv_no0_x89, G_IF_DIAG_CCVV_NO0_X89)
      (sa, ia, Y38.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.90
  //* Hdiag(w,y,a,a)  <--  (   -0.50000000) D1(o1,o2) V2(w,o1,w,o2) 
  for(int sw = 0;sw < nir;++sw){ 
  for(int iw = symblockinfo.psym()(sw,I_C,I_BEGIN);iw <= symblockinfo.psym()(sw,I_C,I_END);++iw){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iw);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iw, sw, V2); // V2=(IR-COV index) 
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
      FC_FUNC(g_if_diag_ccvv_no0_x90, G_IF_DIAG_CCVV_NO0_X90)
        (sa, ia, sw, iw, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(ia, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.91
  //* Hdiag(w,y,a,a)  <--  (    1.00000000) D1(o1,o2) Y39(o1,o2) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nocc, nocc);
  orz::DTensor Y39 = orz::ct::sympack_Xaa(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_diag_ccvv_y91, G_IF_DIAG_CCVV_Y91)
      (sc1, ic1, V2_sym.cptr(), Y39.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ccvv_no0_x91, G_IF_DIAG_CCVV_NO0_X91)
      (sa, ia, Y39.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.92
  //* Hdiag(w,w,a,a)  <--  (    8.00000000) Y40(a,a) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nvir, nvir);
  orz::DTensor Y40 = orz::ct::sympack_Xvv(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_diag_ccvv_y92, G_IF_DIAG_CCVV_Y92)
      (sc1, ic1, V2_sym.cptr(), Y40.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ccvv_no0_x92, G_IF_DIAG_CCVV_NO0_X92)
      (sa, ia, Y40.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.93
  //* Hdiag(w,y,a,a)  <--  (   -4.00000000) Y41(a,a) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nvir, nvir);
  orz::DTensor Y41 = orz::ct::sympack_Xvv(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_diag_ccvv_y93, G_IF_DIAG_CCVV_Y93)
      (sc1, ic1, V2_sym.cptr(), Y41.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ccvv_no0_x93, G_IF_DIAG_CCVV_NO0_X93)
      (sa, ia, Y41.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.94
  //* Hdiag(w,w,a,a)  <--  (   -8.00000000) V2(a,a,w,w) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
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
    FC_FUNC(g_if_diag_ccvv_no0_x94, G_IF_DIAG_CCVV_NO0_X94)
      (sa, ia, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.95
  //* Hdiag(w,y,a,a)  <--  (    2.00000000) V2(a,a,w,w) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
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
    FC_FUNC(g_if_diag_ccvv_no0_x95, G_IF_DIAG_CCVV_NO0_X95)
      (sa, ia, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.96
  //* Hdiag(w,y,a,a)  <--  (    2.00000000) V2(a,a,y,y) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
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
    FC_FUNC(g_if_diag_ccvv_no0_x96, G_IF_DIAG_CCVV_NO0_X96)
      (sa, ia, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.97
  //* Hdiag(w,w,a,c)  <--  (    4.00000000) V2(a,a,w,w) 
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
      Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
      FC_FUNC(g_if_diag_ccvv_no0_x97, G_IF_DIAG_CCVV_NO0_X97)
        (sa, ia, sc, ic, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.98
  //* Hdiag(w,w,a,c)  <--  (    4.00000000) V2(c,c,w,w) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
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
    FC_FUNC(g_if_diag_ccvv_no0_x98, G_IF_DIAG_CCVV_NO0_X98)
      (sc, ic, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.99
  //* Hdiag(w,w,a,a)  <--  (    8.00000000) V2(a,w,w,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
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
    FC_FUNC(g_if_diag_ccvv_no0_x99, G_IF_DIAG_CCVV_NO0_X99)
      (sa, ia, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.100
  //* Hdiag(w,y,a,a)  <--  (   -4.00000000) V2(a,w,w,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
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
    FC_FUNC(g_if_diag_ccvv_no0_x100, G_IF_DIAG_CCVV_NO0_X100)
      (sa, ia, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.101
  //* Hdiag(w,y,a,a)  <--  (   -4.00000000) V2(a,y,y,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
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
    FC_FUNC(g_if_diag_ccvv_no0_x101, G_IF_DIAG_CCVV_NO0_X101)
      (sa, ia, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.102
  //* Hdiag(w,y,a,a)  <--  (    2.00000000) Y42(a,a) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nvir, nvir);
  orz::DTensor Y42 = orz::ct::sympack_Xvv(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_diag_ccvv_y102, G_IF_DIAG_CCVV_Y102)
      (sc1, ic1, V2_sym.cptr(), Y42.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ccvv_no0_x102, G_IF_DIAG_CCVV_NO0_X102)
      (sa, ia, Y42.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.103
  //* Hdiag(w,w,a,a)  <--  (    2.00000000) V2(a,w,w,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
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
    FC_FUNC(g_if_diag_ccvv_no0_x103, G_IF_DIAG_CCVV_NO0_X103)
      (sa, ia, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.104
  //* Hdiag(w,w,a,a)  <--  (   -4.00000000) Y43(a,a) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nvir, nvir);
  orz::DTensor Y43 = orz::ct::sympack_Xvv(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_diag_ccvv_y104, G_IF_DIAG_CCVV_Y104)
      (sc1, ic1, V2_sym.cptr(), Y43.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ccvv_no0_x104, G_IF_DIAG_CCVV_NO0_X104)
      (sa, ia, Y43.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.105
  //* Hdiag(w,w,a,c)  <--  (   -8.00000000) V2(a,w,w,a) 
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
      Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
      FC_FUNC(g_if_diag_ccvv_no0_x105, G_IF_DIAG_CCVV_NO0_X105)
        (sa, ia, sc, ic, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.106
  //* Hdiag(w,w,a,c)  <--  (   -8.00000000) V2(c,w,w,c) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
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
    FC_FUNC(g_if_diag_ccvv_no0_x106, G_IF_DIAG_CCVV_NO0_X106)
      (sc, ic, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.107
  //* Hdiag(w,y,a,a)  <--  (   -0.50000000) D1(o1,o2) V2(y,o1,y,o2) 
  for(int sy = 0;sy < nir;++sy){ 
  for(int iy = symblockinfo.psym()(sy,I_C,I_BEGIN);iy <= symblockinfo.psym()(sy,I_C,I_END);++iy){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iy);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iy, sy, V2); // V2=(IR-COV index) 
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
      FC_FUNC(g_if_diag_ccvv_no0_x107, G_IF_DIAG_CCVV_NO0_X107)
        (sa, ia, sy, iy, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(ia, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.108
  //* Hdiag(w,w,a,a)  <--  (    2.00000000) D1(o1,o2) V2(w,o1,w,o2) 
  for(int sw = 0;sw < nir;++sw){ 
  for(int iw = symblockinfo.psym()(sw,I_C,I_BEGIN);iw <= symblockinfo.psym()(sw,I_C,I_END);++iw){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iw);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iw, sw, V2); // V2=(IR-COV index) 
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
      FC_FUNC(g_if_diag_ccvv_no0_x108, G_IF_DIAG_CCVV_NO0_X108)
        (sa, ia, sw, iw, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(ia, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.109
  //* Hdiag(w,w,a,a)  <--  (   -2.00000000) D1(o1,o2) Y44(o1,o2) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nocc, nocc);
  orz::DTensor Y44 = orz::ct::sympack_Xaa(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_diag_ccvv_y109, G_IF_DIAG_CCVV_Y109)
      (sc1, ic1, V2_sym.cptr(), Y44.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ccvv_no0_x109, G_IF_DIAG_CCVV_NO0_X109)
      (sa, ia, Y44.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.110
  //* Hdiag(w,y,a,a)  <--  (   -0.50000000) D1(o1,o2) V2(w,o1,w,o2) 
  for(int sw = 0;sw < nir;++sw){ 
  for(int iw = symblockinfo.psym()(sw,I_C,I_BEGIN);iw <= symblockinfo.psym()(sw,I_C,I_END);++iw){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iw);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iw, sw, V2); // V2=(IR-COV index) 
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
      FC_FUNC(g_if_diag_ccvv_no0_x110, G_IF_DIAG_CCVV_NO0_X110)
        (sa, ia, sw, iw, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(ia, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.111
  //* Hdiag(w,y,a,a)  <--  (    1.00000000) D1(o1,o2) Y45(o1,o2) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nocc, nocc);
  orz::DTensor Y45 = orz::ct::sympack_Xaa(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_diag_ccvv_y111, G_IF_DIAG_CCVV_Y111)
      (sc1, ic1, V2_sym.cptr(), Y45.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ccvv_no0_x111, G_IF_DIAG_CCVV_NO0_X111)
      (sa, ia, Y45.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.112
  //* Hdiag(w,w,a,a)  <--  (    4.00000000) D1(o1,o2) Y46(o1,o2) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nocc, nocc);
  orz::DTensor Y46 = orz::ct::sympack_Xaa(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_diag_ccvv_y112, G_IF_DIAG_CCVV_Y112)
      (sc1, ic1, V2_sym.cptr(), Y46.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ccvv_no0_x112, G_IF_DIAG_CCVV_NO0_X112)
      (sa, ia, Y46.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.113
  //* Hdiag(w,y,a,a)  <--  (   -2.00000000) D1(o1,o2) Y47(o1,o2) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nocc, nocc);
  orz::DTensor Y47 = orz::ct::sympack_Xaa(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_diag_ccvv_y113, G_IF_DIAG_CCVV_Y113)
      (sc1, ic1, V2_sym.cptr(), Y47.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ccvv_no0_x113, G_IF_DIAG_CCVV_NO0_X113)
      (sa, ia, Y47.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.114
  //* Hdiag(w,y,a,a)  <--  (    1.00000000) D1(o1,o2) V2(w,w,o1,o2) 
  for(int sw = 0;sw < nir;++sw){ 
  for(int iw = symblockinfo.psym()(sw,I_C,I_BEGIN);iw <= symblockinfo.psym()(sw,I_C,I_END);++iw){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iw);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iw, sw, V2); // V2=(IR-COV index) 
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
      FC_FUNC(g_if_diag_ccvv_no0_x114, G_IF_DIAG_CCVV_NO0_X114)
        (sa, ia, sw, iw, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(ia, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.115
  //* Hdiag(w,y,a,a)  <--  (    1.00000000) D1(o1,o2) V2(y,y,o1,o2) 
  for(int sy = 0;sy < nir;++sy){ 
  for(int iy = symblockinfo.psym()(sy,I_C,I_BEGIN);iy <= symblockinfo.psym()(sy,I_C,I_END);++iy){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iy);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iy, sy, V2); // V2=(IR-COV index) 
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
      FC_FUNC(g_if_diag_ccvv_no0_x115, G_IF_DIAG_CCVV_NO0_X115)
        (sa, ia, sy, iy, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(ia, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.116
  //* Hdiag(w,w,a,a)  <--  (   -4.00000000) D1(o1,o2) V2(w,w,o1,o2) 
  for(int sw = 0;sw < nir;++sw){ 
  for(int iw = symblockinfo.psym()(sw,I_C,I_BEGIN);iw <= symblockinfo.psym()(sw,I_C,I_END);++iw){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iw);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iw, sw, V2); // V2=(IR-COV index) 
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
      FC_FUNC(g_if_diag_ccvv_no0_x116, G_IF_DIAG_CCVV_NO0_X116)
        (sa, ia, sw, iw, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(ia, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.117
  //* Hdiag(w,y,a,c)  <--  (    2.00000000) D2(o1,o3,o2,o4) V2(o1,o3,o2,o4) 
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
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
      FC_FUNC(g_if_diag_ccvv_no0_x117, G_IF_DIAG_CCVV_NO0_X117)
        (sc, ic, so1, io1, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.118
  //* Hdiag(w,w,a,c)  <--  (   -1.00000000) D2(o1,o3,o2,o4) V2(o1,o3,o2,o4) 
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
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
      FC_FUNC(g_if_diag_ccvv_no0_x118, G_IF_DIAG_CCVV_NO0_X118)
        (sc, ic, so1, io1, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.119
  //* Hdiag(w,w,a,a)  <--  (    2.00000000) D2(o1,o3,o2,o4) V2(o1,o3,o2,o4) 
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
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
      FC_FUNC(g_if_diag_ccvv_no0_x119, G_IF_DIAG_CCVV_NO0_X119)
        (sa, ia, so1, io1, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(ia, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.120
  //* Hdiag(w,y,a,a)  <--  (   -1.00000000) D2(o1,o3,o2,o4) V2(o1,o3,o2,o4) 
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
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
      FC_FUNC(g_if_diag_ccvv_no0_x120, G_IF_DIAG_CCVV_NO0_X120)
        (sa, ia, so1, io1, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(ia, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.121
  //* Hdiag(w,w,a,a)  <--  (    4.00000000) D1(o1,o2) V2(a,a,o1,o2) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
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
    FC_FUNC(g_if_diag_ccvv_no0_x121, G_IF_DIAG_CCVV_NO0_X121)
      (sa, ia, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.122
  //* Hdiag(w,y,a,a)  <--  (   -2.00000000) D1(o1,o2) V2(a,a,o1,o2) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
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
    FC_FUNC(g_if_diag_ccvv_no0_x122, G_IF_DIAG_CCVV_NO0_X122)
      (sa, ia, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.123
  //* Hdiag(w,w,a,a)  <--  (   -2.00000000) D1(o1,o2) V2(a,o1,o2,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
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
    FC_FUNC(g_if_diag_ccvv_no0_x123, G_IF_DIAG_CCVV_NO0_X123)
      (sa, ia, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.124
  //* Hdiag(w,y,a,a)  <--  (    1.00000000) D1(o1,o2) V2(a,o1,o2,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
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
    FC_FUNC(g_if_diag_ccvv_no0_x124, G_IF_DIAG_CCVV_NO0_X124)
      (sa, ia, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.125
  //* Hdiag(w,w,a,a)  <--  (    8.00000000) V2(a,w,w,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
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
    FC_FUNC(g_if_diag_ccvv_no0_x125, G_IF_DIAG_CCVV_NO0_X125)
      (sa, ia, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.126
  //* Hdiag(w,y,a,a)  <--  (   -4.00000000) V2(a,w,w,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
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
    FC_FUNC(g_if_diag_ccvv_no0_x126, G_IF_DIAG_CCVV_NO0_X126)
      (sa, ia, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.127
  //* Hdiag(w,y,a,a)  <--  (   -4.00000000) V2(a,y,y,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
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
    FC_FUNC(g_if_diag_ccvv_no0_x127, G_IF_DIAG_CCVV_NO0_X127)
      (sa, ia, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.128
  //* Hdiag(w,y,a,a)  <--  (    2.00000000) Y48(a,a) 
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
    FC_FUNC(g_if_diag_ccvv_y128, G_IF_DIAG_CCVV_Y128)
      (sc1, ic1, V2_sym.cptr(), Y48.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ccvv_no0_x128, G_IF_DIAG_CCVV_NO0_X128)
      (sa, ia, Y48.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.129
  //* Hdiag(w,w,a,a)  <--  (    2.00000000) V2(a,w,w,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
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
    FC_FUNC(g_if_diag_ccvv_no0_x129, G_IF_DIAG_CCVV_NO0_X129)
      (sa, ia, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.130
  //* Hdiag(w,w,a,a)  <--  (   -4.00000000) Y49(a,a) 
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
    FC_FUNC(g_if_diag_ccvv_y130, G_IF_DIAG_CCVV_Y130)
      (sc1, ic1, V2_sym.cptr(), Y49.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ccvv_no0_x130, G_IF_DIAG_CCVV_NO0_X130)
      (sa, ia, Y49.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.131
  //* Hdiag(w,w,a,a)  <--  (    8.00000000) Y50(a,a) 
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
    FC_FUNC(g_if_diag_ccvv_y131, G_IF_DIAG_CCVV_Y131)
      (sc1, ic1, V2_sym.cptr(), Y50.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ccvv_no0_x131, G_IF_DIAG_CCVV_NO0_X131)
      (sa, ia, Y50.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.132
  //* Hdiag(w,y,a,a)  <--  (   -4.00000000) Y51(a,a) 
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
    FC_FUNC(g_if_diag_ccvv_y132, G_IF_DIAG_CCVV_Y132)
      (sc1, ic1, V2_sym.cptr(), Y51.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ccvv_no0_x132, G_IF_DIAG_CCVV_NO0_X132)
      (sa, ia, Y51.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.133
  //* Hdiag(w,w,a,a)  <--  (   -8.00000000) V2(a,a,w,w) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
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
    FC_FUNC(g_if_diag_ccvv_no0_x133, G_IF_DIAG_CCVV_NO0_X133)
      (sa, ia, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.134
  //* Hdiag(w,y,a,a)  <--  (    2.00000000) V2(a,a,w,w) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
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
    FC_FUNC(g_if_diag_ccvv_no0_x134, G_IF_DIAG_CCVV_NO0_X134)
      (sa, ia, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.135
  //* Hdiag(w,y,a,a)  <--  (    2.00000000) V2(a,a,y,y) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
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
    FC_FUNC(g_if_diag_ccvv_no0_x135, G_IF_DIAG_CCVV_NO0_X135)
      (sa, ia, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.136
  //* Hdiag(w,w,a,a)  <--  (   -2.00000000) D1(o1,o2) V2(a,o1,o2,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
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
    FC_FUNC(g_if_diag_ccvv_no0_x136, G_IF_DIAG_CCVV_NO0_X136)
      (sa, ia, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.137
  //* Hdiag(w,y,a,a)  <--  (    1.00000000) D1(o1,o2) V2(a,o1,o2,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
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
    FC_FUNC(g_if_diag_ccvv_no0_x137, G_IF_DIAG_CCVV_NO0_X137)
      (sa, ia, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.138
  //* Hdiag(w,w,a,a)  <--  (    4.00000000) D1(o1,o2) V2(a,a,o1,o2) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
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
    FC_FUNC(g_if_diag_ccvv_no0_x138, G_IF_DIAG_CCVV_NO0_X138)
      (sa, ia, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.139
  //* Hdiag(w,y,a,a)  <--  (   -2.00000000) D1(o1,o2) V2(a,a,o1,o2) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
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
    FC_FUNC(g_if_diag_ccvv_no0_x139, G_IF_DIAG_CCVV_NO0_X139)
      (sa, ia, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }

  for(int ssig = 0;ssig < nir;++ssig){                                                                                      
  for(int isig = symblockinfo.psym()(ssig,I_V,I_BEGIN);isig <= symblockinfo.psym()(ssig,I_V,I_END);++isig){                 
    Hdiagb = retval.get_amp2(isig);                                                                    
    FC_FUNC(g_if_sigma_ccvv_scale,G_IF_SIGMA_CCVV_SCALE) // to make hamiltonian symmetric (scale by some factor of overlap) 
      (ssig, isig, Hdiagb.cptr(), nir, nsym, psym);                                                    
    retval.put_amp2(isig, Hdiagb); // S2ija, [b] <<-- Sb                                               
  } // End isig                                                                                                             
  } // End ssig                                                                                                             

  return retval; 
} 
