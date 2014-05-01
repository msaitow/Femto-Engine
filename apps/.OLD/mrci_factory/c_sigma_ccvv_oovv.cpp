                                                                                
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
#include <sci/ctnew2/c_sigma_ccvv_oovv.h>                                            
                                                                                
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
                                                                                
orz::ct::BareAmpPack orz::ct::sigma_ccvv_oovv(const orz::ct::Input &ctinp,                                    
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
  //* X(w,y,o1,o2)  <--  (    1.00000000)  D2(o1,o3,o2,o4) V2(w,o3,y,o4) 
  //* S2(w,y,a,c)  <--  (    1.00000000) T2(o1,o2,a,c) X(w,y,o1,o2) 
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
    orz::DTensor X(nclosed, nocc, nocc);
    orz::DTensor Xcaa = orz::ct::sympack_Xcaa(symblockinfo, sw, X);
    FC_FUNC(g_if_sigma_ccvv_oovv_no0_x0, G_IF_SIGMA_CCVV_OOVV_NO0_X0)
      (sw, iw, V2_sym.cptr(), Xcaa.cptr(), nir, nsym, psym);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      T2b = T2.get_amp2(ic);
      FC_FUNC(g_if_sigma_ccvv_oovv_no1_x0, G_IF_SIGMA_CCVV_OOVV_NO1_X0)
        (sc, ic, sw, iw, T2b.cptr(), Xcaa.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  }


  {
  // No.1
  //* X(w,y,o3,o4)  <--  (    1.00000000)  D2(o1,o3,o2,o4) V2(w,o2,y,o1) 
  //* S2(w,y,a,c)  <--  (    1.00000000) T2(o3,o4,c,a) X(w,y,o3,o4) 
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
    orz::DTensor X(nclosed, nocc, nocc);
    orz::DTensor Xcaa = orz::ct::sympack_Xcaa(symblockinfo, sw, X);
    FC_FUNC(g_if_sigma_ccvv_oovv_no0_x1, G_IF_SIGMA_CCVV_OOVV_NO0_X1)
      (sw, iw, V2_sym.cptr(), Xcaa.cptr(), nir, nsym, psym);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      for(int sa = 0;sa < nir;++sa){ 
      for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
        T2b = T2.get_amp2(ia);
        FC_FUNC(g_if_sigma_ccvv_oovv_no1_x1, G_IF_SIGMA_CCVV_OOVV_NO1_X1)
          (sa, ia, sc, ic, sw, iw, T2b.cptr(), Xcaa.cptr(), S2b.cptr(), nir, nsym, psym);
      }
      }
      retval.acc_amp2(ic, S2b);
    }
    }
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
