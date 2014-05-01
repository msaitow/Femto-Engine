                                                                                
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
#include <sci/ctnew2/c_sigma_g_oovv.h>                                            
                                                                                
using std::cout;                                                                
using std::endl;                                                                
                                                                                
#define FLOPCOUNT                                                               
                                                                                
//                                                              
//   _______________                                  ______    
//  |          |                 .'. .`. `````|`````.~      ~.  
//  |______    |______         .'   `   `.    |    |          | 
//  |          |             .'           `.  |    |          | 
//  |          |___________.'               `.|     `.______.'  
//                                                              

                                                                                
// ***************************************************************************  
// orz::ct::mrci                                                                
// ***************************************************************************  
									    /*!        
   @brief CT input                                                              
									     */        
                                                                                
double orz::ct::sigma_g_oovv(const orz::ct::Input &ctinp,                                                  
                                  const orz::ct::SymBlockInfo &symblockinfo,                                
                                  const orz::ct::HintMO &hintmo,                                            
                                  const orz::ct::RdmPack &rdmPack_sym,                                      
                                  const double init_value,                                                  
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
                                                                                                                
  double S0 = init_value;                                                               
  orz::DTensor T2b; // Container of T2_aae,[b] tensor                                             


  {
  // No.0
  //* X(o1,o2,o4,o3)  <--  (    1.00000000)  T2(o1,o2,v1,v2) V2(v2,o4,o3,v1) 
  //* S0()  <--  (    1.00000000) D2(o1,o3,o2,o4) X(o1,o2,o4,o3) 
  orz::DTensor X(nocc, nocc, nocc, nocc);
  orz::DTensor Xaaaa = orz::ct::sympack_Xaaaa(symblockinfo, 0, X);
  for(int sv2 = 0;sv2 < nir;++sv2){ 
  for(int iv2 = symblockinfo.psym()(sv2,I_V,I_BEGIN);iv2 <= symblockinfo.psym()(sv2,I_V,I_END);++iv2){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iv2);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iv2, sv2, V2); // V2=(IR-COV index) 
    T2b = T2.get_amp2(iv2);
    FC_FUNC(g_if_sigma_g_oovv_no0_x0, G_IF_SIGMA_G_OOVV_NO0_X0)
      (sv2, iv2, T2b.cptr(), V2_sym.cptr(), Xaaaa.cptr(), nir, nsym, psym);
  }
  }
  FC_FUNC(g_if_sigma_g_oovv_no1_x0, G_IF_SIGMA_G_OOVV_NO1_X0)
    (Xaaaa.cptr(), &S0, nir, nsym, psym);
  }

  return  S0;
} 
