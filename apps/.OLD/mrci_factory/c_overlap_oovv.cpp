                                                                                
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
#include <sci/ctnew2/c_overlap_oovv.h>                                            
                                                                                
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
                                                                                
orz::ct::BareAmpPack orz::ct::overlap_oovv(const orz::ct::Input &ctinp,                                    
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
  std::string name_of_sigma = "O2" + stm.str() + "]"; // Name of the Sigma vector  
  orz::ct::BareAmpPack retval                                                                                   
    = orz::ct::BareAmpPack(ctinp, symblockinfo, name_of_sigma); // Sigma(a, a', e, e') tensor                   
                                                                                                                
  orz::DTensor O2b; // Container of S2_aae,[b] tensor                                   
                                                                                                                
  orz::DTensor T2b; // Container of T2_aae,[b] tensor                                             


  {
  // No.0
  //* O2(i,k,a,c)  <--  (    1.00000000) D2(i,o2,k,o1) T2(o2,o1,a,c) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    O2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    FC_FUNC(g_if_overlap_oovv_no0_x0, G_IF_OVERLAP_OOVV_NO0_X0)
      (sc, ic, T2b.cptr(), O2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ic, O2b);
  }
  }
  }


  {
  // No.1
  //* O2(i,k,a,c)  <--  (    1.00000000) D2(i,o1,k,o2) T2(o2,o1,c,a) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    O2b = orz::DTensor(retval.namps_iamp()[ic]);
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      T2b = T2.get_amp2(ia);
      FC_FUNC(g_if_overlap_oovv_no0_x1, G_IF_OVERLAP_OOVV_NO0_X1)
        (sa, ia, sc, ic, T2b.cptr(), O2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(ic, O2b);
    }
    }
  }
  }
  }

  for(int ssig = 0;ssig < nir;++ssig){                                                                                      
  for(int isig = symblockinfo.psym()(ssig,I_V,I_BEGIN);isig <= symblockinfo.psym()(ssig,I_V,I_END);++isig){                 
    O2b = retval.get_amp2(isig);                                                                    
    FC_FUNC(g_if_sigma_oovv_scale,G_IF_SIGMA_OOVV_SCALE) // to make hamiltonian symmetric (scale by some factor of overlap) 
      (ssig, isig, O2b.cptr(), nir, nsym, psym);                                                    
    retval.put_amp2(isig, O2b); // S2ija, [b] <<-- Sb                                               
  } // End isig                                                                                                             
  } // End ssig                                                                                                             

  return retval; 
} 
