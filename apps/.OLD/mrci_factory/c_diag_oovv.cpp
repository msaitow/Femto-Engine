                                                                                
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
#include <sci/ctnew2/c_diag_oovv.h>                                            
                                                                                
using std::cout;                                                                
using std::endl;                                                                
                                                                                
#define FLOPCOUNT                                                               
                                                                                
//  8 8888888888   8 8888888888            ,8.       ,8.    8888888 8888888888 ,o888888o.     
//  8 8888         8 8888                 ,888.     ,888.         8 8888    . 8888     `88.   
//  8 8888         8 8888                .`8888.   .`8888.        8 8888   ,8 8888       `8b  
//  8 8888         8 8888               ,8.`8888. ,8.`8888.       8 8888   88 8888        `8b 
//  8 888888888888 8 888888888888      ,8'8.`8888,8^8.`8888.      8 8888   88 8888         88 
//  8 8888         8 8888             ,8' `8.`8888' `8.`8888.     8 8888   88 8888         88 
//  8 8888         8 8888            ,8'   `8.`88'   `8.`8888.    8 8888   88 8888        ,8P 
//  8 8888         8 8888           ,8'     `8.`'     `8.`8888.   8 8888   `8 8888       ,8P  
//  8 8888         8 8888          ,8'       `8        `8.`8888.  8 8888    ` 8888     ,88'   
//  8 8888         8 888888888888 ,8'         `         `8.`8888. 8 8888       `8888888P'     

                                                                                
// ***************************************************************************  
// orz::ct::mrci                                                                
// ***************************************************************************  
									    /*!        
   @brief CT input                                                              
									     */        
                                                                                
orz::ct::BareAmpPack orz::ct::diag_oovv(const orz::ct::Input &ctinp,                                    
					 const orz::ct::SymBlockInfo &symblockinfo,                                
					 const orz::ct::HintMO &hintmo,                                            
					 const orz::ct::RdmPack &rdmPack_sym,                                      
					 const orz::DTensor &rdm4,                                                 
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
                                                                                                                
  orz::DTensor rdm4_sym;                                                                                        


  {
  // No.0
  //* Hdiag(i,k,a,c)  <--  (    2.00000000) Y0 D2(i,i,k,k) 
  // The effective tensor is detected .... 
  double Y0 = 0;
  FC_FUNC(g_if_diag_oovv_y0, G_IF_DIAG_OOVV_Y0)
    (moint1_sym.cptr(), &Y0, nir, nsym, psym);
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
    FC_FUNC(g_if_diag_oovv_no0_x0, G_IF_DIAG_OOVV_NO0_X0)
      (sc, ic, &Y0, Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.1
  //* Hdiag(i,k,a,c)  <--  (    1.00000000) D2(i,i,k,k) h(a,a) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
    FC_FUNC(g_if_diag_oovv_no0_x1, G_IF_DIAG_OOVV_NO0_X1)
      (sc, ic, moint1_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.2
  //* Hdiag(i,k,a,c)  <--  (    1.00000000) D2(i,i,k,k) h(c,c) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
    FC_FUNC(g_if_diag_oovv_no0_x2, G_IF_DIAG_OOVV_NO0_X2)
      (sc, ic, moint1_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.3
  //* Hdiag(i,k,a,c)  <--  (    2.00000000) Y1 D2(i,i,k,k) 
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
    FC_FUNC(g_if_diag_oovv_y3, G_IF_DIAG_OOVV_Y3)
      (sc1, ic1, V2_sym.cptr(), &Y1, nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
    FC_FUNC(g_if_diag_oovv_no0_x3, G_IF_DIAG_OOVV_NO0_X3)
      (sc, ic, &Y1, Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.4
  //* Hdiag(i,k,a,c)  <--  (   -1.00000000) Y2 D2(i,i,k,k) 
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
    FC_FUNC(g_if_diag_oovv_y4, G_IF_DIAG_OOVV_Y4)
      (sc1, ic1, V2_sym.cptr(), &Y2, nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
    FC_FUNC(g_if_diag_oovv_no0_x4, G_IF_DIAG_OOVV_NO0_X4)
      (sc, ic, &Y2, Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.5
  //* Hdiag(i,k,a,c)  <--  (    2.00000000) D2(i,i,k,k) Y3(a,a) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nvir, nvir);
  orz::DTensor Y3 = orz::ct::sympack_Xvv(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_diag_oovv_y5, G_IF_DIAG_OOVV_Y5)
      (sc1, ic1, V2_sym.cptr(), Y3.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
    FC_FUNC(g_if_diag_oovv_no0_x5, G_IF_DIAG_OOVV_NO0_X5)
      (sc, ic, Y3.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.6
  //* Hdiag(i,k,a,c)  <--  (    2.00000000) D2(i,i,k,k) Y4(c,c) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nvir, nvir);
  orz::DTensor Y4 = orz::ct::sympack_Xvv(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_diag_oovv_y6, G_IF_DIAG_OOVV_Y6)
      (sc1, ic1, V2_sym.cptr(), Y4.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
    FC_FUNC(g_if_diag_oovv_no0_x6, G_IF_DIAG_OOVV_NO0_X6)
      (sc, ic, Y4.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.7
  //* Hdiag(i,k,a,c)  <--  (   -1.00000000) D2(i,i,k,k) Y5(a,a) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nvir, nvir);
  orz::DTensor Y5 = orz::ct::sympack_Xvv(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_diag_oovv_y7, G_IF_DIAG_OOVV_Y7)
      (sc1, ic1, V2_sym.cptr(), Y5.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
    FC_FUNC(g_if_diag_oovv_no0_x7, G_IF_DIAG_OOVV_NO0_X7)
      (sc, ic, Y5.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.8
  //* Hdiag(i,k,a,c)  <--  (   -1.00000000) D2(i,i,k,k) Y6(c,c) 
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
    FC_FUNC(g_if_diag_oovv_y8, G_IF_DIAG_OOVV_Y8)
      (sc1, ic1, V2_sym.cptr(), Y6.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
    FC_FUNC(g_if_diag_oovv_no0_x8, G_IF_DIAG_OOVV_NO0_X8)
      (sc, ic, Y6.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.9
  //* Hdiag(i,k,a,c)  <--  (    1.00000000) D2(i,i,k,k) V2(c,c,a,a) 
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
    FC_FUNC(g_if_diag_oovv_no0_x9, G_IF_DIAG_OOVV_NO0_X9)
      (sc, ic, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.10
  //* Hdiag(i,k,a,a)  <--  (    2.00000000) Y7 D2(i,k,k,i) 
  // The effective tensor is detected .... 
  double Y7 = 0;
  FC_FUNC(g_if_diag_oovv_y10, G_IF_DIAG_OOVV_Y10)
    (moint1_sym.cptr(), &Y7, nir, nsym, psym);
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_oovv_no0_x10, G_IF_DIAG_OOVV_NO0_X10)
      (sa, ia, &Y7, Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.11
  //* Hdiag(i,k,a,c)  <--  (    1.00000000) D3(i,i,k,k,o1,o2) h(o1,o2) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
    FC_FUNC(g_if_diag_oovv_no0_x11, G_IF_DIAG_OOVV_NO0_X11)
      (sc, ic, moint1_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.12
  //* Hdiag(i,k,a,a)  <--  (    1.00000000) D3(i,k,k,i,o1,o2) h(o1,o2) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_oovv_no0_x12, G_IF_DIAG_OOVV_NO0_X12)
      (sa, ia, moint1_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.13
  //* Hdiag(i,k,a,a)  <--  (    2.00000000) D2(i,k,k,i) h(a,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_oovv_no0_x13, G_IF_DIAG_OOVV_NO0_X13)
      (sa, ia, moint1_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.14
  //* Hdiag(i,k,a,a)  <--  (    2.00000000) Y8 D2(i,k,k,i) 
  // The effective tensor is detected .... 
  double Y8 = 0;
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
    FC_FUNC(g_if_diag_oovv_y14, G_IF_DIAG_OOVV_Y14)
      (sc1, ic1, V2_sym.cptr(), &Y8, nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_oovv_no0_x14, G_IF_DIAG_OOVV_NO0_X14)
      (sa, ia, &Y8, Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.15
  //* Hdiag(i,k,a,a)  <--  (   -1.00000000) Y9 D2(i,k,k,i) 
  // The effective tensor is detected .... 
  double Y9 = 0;
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
    FC_FUNC(g_if_diag_oovv_y15, G_IF_DIAG_OOVV_Y15)
      (sc1, ic1, V2_sym.cptr(), &Y9, nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_oovv_no0_x15, G_IF_DIAG_OOVV_NO0_X15)
      (sa, ia, &Y9, Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.16
  //* Hdiag(i,k,a,c)  <--  (    2.00000000) D3(i,i,k,k,o1,o2) Y10(o1,o2) 
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
    FC_FUNC(g_if_diag_oovv_y16, G_IF_DIAG_OOVV_Y16)
      (sc1, ic1, V2_sym.cptr(), Y10.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
    FC_FUNC(g_if_diag_oovv_no0_x16, G_IF_DIAG_OOVV_NO0_X16)
      (sc, ic, Y10.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.17
  //* Hdiag(i,k,a,c)  <--  (   -1.00000000) D3(i,i,k,k,o1,o2) Y11(o1,o2) 
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
    FC_FUNC(g_if_diag_oovv_y17, G_IF_DIAG_OOVV_Y17)
      (sc1, ic1, V2_sym.cptr(), Y11.cptr(), nir, nsym, psym);
  }
  }
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
    FC_FUNC(g_if_diag_oovv_no0_x17, G_IF_DIAG_OOVV_NO0_X17)
      (sc, ic, Y11.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.18
  //* Hdiag(i,k,a,c)  <--  (    1.00000000) D3(i,i,k,k,o1,o2) V2(a,a,o1,o2) 
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
      FC_FUNC(g_if_diag_oovv_no0_x18, G_IF_DIAG_OOVV_NO0_X18)
        (sa, ia, sc, ic, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.19
  //* Hdiag(i,k,a,c)  <--  (    1.00000000) D3(i,i,k,k,o1,o2) V2(c,c,o1,o2) 
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
    FC_FUNC(g_if_diag_oovv_no0_x19, G_IF_DIAG_OOVV_NO0_X19)
      (sc, ic, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.20
  //* Hdiag(i,k,a,a)  <--  (    1.00000000) D3(i,k,k,i,o1,o2) Y12(o1,o2) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nocc, nocc);
  orz::DTensor Y12 = orz::ct::sympack_Xaa(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_diag_oovv_y20, G_IF_DIAG_OOVV_Y20)
      (sc1, ic1, V2_sym.cptr(), Y12.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_oovv_no0_x20, G_IF_DIAG_OOVV_NO0_X20)
      (sa, ia, Y12.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.21
  //* Hdiag(i,k,a,a)  <--  (   -0.50000000) D3(i,k,k,i,o1,o2) Y13(o1,o2) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nocc, nocc);
  orz::DTensor Y13 = orz::ct::sympack_Xaa(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_diag_oovv_y21, G_IF_DIAG_OOVV_Y21)
      (sc1, ic1, V2_sym.cptr(), Y13.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_oovv_no0_x21, G_IF_DIAG_OOVV_NO0_X21)
      (sa, ia, Y13.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.22
  //* Hdiag(i,k,a,a)  <--  (    2.00000000) D2(i,k,k,i) Y14(a,a) 
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
    FC_FUNC(g_if_diag_oovv_y22, G_IF_DIAG_OOVV_Y22)
      (sc1, ic1, V2_sym.cptr(), Y14.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_oovv_no0_x22, G_IF_DIAG_OOVV_NO0_X22)
      (sa, ia, Y14.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.23
  //* Hdiag(i,k,a,a)  <--  (   -1.00000000) D2(i,k,k,i) Y15(a,a) 
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
    FC_FUNC(g_if_diag_oovv_y23, G_IF_DIAG_OOVV_Y23)
      (sc1, ic1, V2_sym.cptr(), Y15.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_oovv_no0_x23, G_IF_DIAG_OOVV_NO0_X23)
      (sa, ia, Y15.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.24
  //* Hdiag(i,k,a,a)  <--  (   -0.50000000) D3(i,k,k,i,o1,o2) Y16(o1,o2) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nocc, nocc);
  orz::DTensor Y16 = orz::ct::sympack_Xaa(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_diag_oovv_y24, G_IF_DIAG_OOVV_Y24)
      (sc1, ic1, V2_sym.cptr(), Y16.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_oovv_no0_x24, G_IF_DIAG_OOVV_NO0_X24)
      (sa, ia, Y16.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.25
  //* Hdiag(i,k,a,a)  <--  (    1.00000000) D3(i,k,k,i,o1,o2) Y17(o1,o2) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nocc, nocc);
  orz::DTensor Y17 = orz::ct::sympack_Xaa(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_diag_oovv_y25, G_IF_DIAG_OOVV_Y25)
      (sc1, ic1, V2_sym.cptr(), Y17.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_oovv_no0_x25, G_IF_DIAG_OOVV_NO0_X25)
      (sa, ia, Y17.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.26
  //* Hdiag(i,k,a,c)  <--  (    0.50000000) D4(o1,o3,i,i,k,k,o2,o4) V2(o1,o3,o2,o4) 
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
    for(int so3 = 0;so3 < nir;++so3){ 
    for(int io3 = symblockinfo.psym()(so3,I_O,I_BEGIN);io3 <= symblockinfo.psym()(so3,I_O,I_END);++io3){ 
      for(int sc = 0;sc < nir;++sc){ 
      for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
        Hdiagb = orz::DTensor(retval.namps_iamp()[ic]);
        // Load D4 from disk, or GA ....                                                     
        int imoi = amo2imo[io1] - nclosed;                              
        int imoj = amo2imo[io3] - nclosed;                              
                                                                                             
        orz::DTensor rdm4_ij_sliced = rdm4(orz::Slice(imoi,imoi+1), orz::Slice(imoj,imoj+1), 
                                           orz::Slice(),            orz::Slice(),            
                                           orz::Slice(),            orz::Slice(),            
                                           orz::Slice(),            orz::Slice()).copy();    
        rdm4_sym = orz::ct::sympack_rdm4_2(symblockinfo, io1, so1, io3, so3, rdm4_ij_sliced);    
        FC_FUNC(g_if_set_d4,G_IF_SET_D4)(so1, so3, io1, io3, rdm4_sym.cptr(), nir, nsym, psym);  
        FC_FUNC(g_if_diag_oovv_no0_x26, G_IF_DIAG_OOVV_NO0_X26)
          (sc, ic, so1, io1, so3, io3, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
        retval.acc_amp2(ic, Hdiagb);
        FC_FUNC(g_if_unset_d4,G_IF_UNSET_D4)();
      }
      }
    }
    }
  }
  }
  }


  {
  // No.27
  //* Hdiag(i,k,a,a)  <--  (    0.50000000) D4(o1,o3,i,k,k,i,o2,o4) V2(o1,o3,o2,o4) 
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
    for(int so3 = 0;so3 < nir;++so3){ 
    for(int io3 = symblockinfo.psym()(so3,I_O,I_BEGIN);io3 <= symblockinfo.psym()(so3,I_O,I_END);++io3){ 
      for(int sa = 0;sa < nir;++sa){ 
      for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
        Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
        // Load D4 from disk, or GA ....                                                     
        int imoi = amo2imo[io1] - nclosed;                              
        int imoj = amo2imo[io3] - nclosed;                              
                                                                                             
        orz::DTensor rdm4_ij_sliced = rdm4(orz::Slice(imoi,imoi+1), orz::Slice(imoj,imoj+1), 
                                           orz::Slice(),            orz::Slice(),            
                                           orz::Slice(),            orz::Slice(),            
                                           orz::Slice(),            orz::Slice()).copy();    
        rdm4_sym = orz::ct::sympack_rdm4_2(symblockinfo, io1, so1, io3, so3, rdm4_ij_sliced);    
        FC_FUNC(g_if_set_d4,G_IF_SET_D4)(so1, so3, io1, io3, rdm4_sym.cptr(), nir, nsym, psym);  
        FC_FUNC(g_if_diag_oovv_no0_x27, G_IF_DIAG_OOVV_NO0_X27)
          (sa, ia, so1, io1, so3, io3, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
        retval.acc_amp2(ia, Hdiagb);
        FC_FUNC(g_if_unset_d4,G_IF_UNSET_D4)();
      }
      }
    }
    }
  }
  }
  }


  {
  // No.28
  //* Hdiag(i,k,a,a)  <--  (    1.00000000) D3(i,k,k,i,o1,o2) V2(a,a,o1,o2) 
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
    FC_FUNC(g_if_diag_oovv_no0_x28, G_IF_DIAG_OOVV_NO0_X28)
      (sa, ia, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.29
  //* Hdiag(i,k,a,a)  <--  (    1.00000000) D3(i,k,k,o1,o2,i) V2(a,o1,o2,a) 
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
    FC_FUNC(g_if_diag_oovv_no0_x29, G_IF_DIAG_OOVV_NO0_X29)
      (sa, ia, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.30
  //* Hdiag(i,k,a,c)  <--  (    1.00000000) D3(i,o1,k,k,o2,i) V2(a,o1,o2,a) 
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
      FC_FUNC(g_if_diag_oovv_no0_x30, G_IF_DIAG_OOVV_NO0_X30)
        (sa, ia, sc, ic, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.31
  //* Hdiag(i,k,a,c)  <--  (    1.00000000) D3(i,i,k,o1,o2,k) V2(c,o1,o2,c) 
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
    FC_FUNC(g_if_diag_oovv_no0_x31, G_IF_DIAG_OOVV_NO0_X31)
      (sc, ic, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }


  {
  // No.32
  //* Hdiag(i,k,a,a)  <--  (   -1.00000000) D2(i,k,k,i) Y18(a,a) 
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
    FC_FUNC(g_if_diag_oovv_y32, G_IF_DIAG_OOVV_Y32)
      (sc1, ic1, V2_sym.cptr(), Y18.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_oovv_no0_x32, G_IF_DIAG_OOVV_NO0_X32)
      (sa, ia, Y18.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.33
  //* Hdiag(i,k,a,a)  <--  (    2.00000000) D2(i,k,k,i) Y19(a,a) 
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
    FC_FUNC(g_if_diag_oovv_y33, G_IF_DIAG_OOVV_Y33)
      (sc1, ic1, V2_sym.cptr(), Y19.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_oovv_no0_x33, G_IF_DIAG_OOVV_NO0_X33)
      (sa, ia, Y19.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.34
  //* Hdiag(i,k,a,a)  <--  (    1.00000000) D3(i,k,k,o1,o2,i) V2(a,o1,o2,a) 
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
    FC_FUNC(g_if_diag_oovv_no0_x34, G_IF_DIAG_OOVV_NO0_X34)
      (sa, ia, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.35
  //* Hdiag(i,k,a,a)  <--  (    1.00000000) D3(i,k,k,i,o1,o2) V2(a,a,o1,o2) 
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
    FC_FUNC(g_if_diag_oovv_no0_x35, G_IF_DIAG_OOVV_NO0_X35)
      (sa, ia, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.36
  //* Hdiag(i,k,a,c)  <--  (    1.00000000) D2(i,k,k,i) V2(c,a,a,c) 
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
    FC_FUNC(g_if_diag_oovv_no0_x36, G_IF_DIAG_OOVV_NO0_X36)
      (sc, ic, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, Hdiagb);
  }
  }
  }

  for(int ssig = 0;ssig < nir;++ssig){                                                                                      
  for(int isig = symblockinfo.psym()(ssig,I_V,I_BEGIN);isig <= symblockinfo.psym()(ssig,I_V,I_END);++isig){                 
    Hdiagb = retval.get_amp2(isig);                                                                    
    FC_FUNC(g_if_sigma_oovv_scale,G_IF_SIGMA_OOVV_SCALE) // to make hamiltonian symmetric (scale by some factor of overlap) 
      (ssig, isig, Hdiagb.cptr(), nir, nsym, psym);                                                    
    retval.put_amp2(isig, Hdiagb); // S2ija, [b] <<-- Sb                                               
  } // End isig                                                                                                             
  } // End ssig                                                                                                             

  return retval; 
} 
