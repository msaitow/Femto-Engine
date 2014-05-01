                                                                                
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
#include <sci/ctnew2/c_sigma_ooov_ccvv.h>                                            
                                                                                
using std::cout;                                                                
using std::endl;                                                                
                                                                                
#define FLOPCOUNT                                                               
                                                                                
//  `MMMMMMM                                         
//   MM    \                         /              
//   MM       ____  ___  __    __   /M      _____    
//   MM   ,  6MMMMb `MM 6MMb  6MMb /MMMMM  6MMMMMb   
//   MMMMMM 6M'  `Mb MM69 `MM69 `Mb MM    6M'   `Mb  
//   MM   ` MM    MM MM'   MM'   MM MM    MM     MM  
//   MM     MMMMMMMM MM    MM    MM MM    MM     MM  
//   MM     MM       MM    MM    MM MM    MM     MM  
//   MM     YM    d9 MM    MM    MM YM.  ,YM.   ,M9  
//  _MM_     YMMMM9 _MM_  _MM_  _MM_ YMMM9 YMMMMM9   

                                                                                
// ***************************************************************************  
// orz::ct::mrci                                                                
// ***************************************************************************  
									    /*!        
   @brief CT input                                                              
									     */        
                                                                                
orz::ct::BareAmpPack orz::ct::sigma_ooov_ccvv(const orz::ct::Input &ctinp,                                    
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
  //* X(o1,a)  <--  (    1.00000000)  T2(c2,c1,v1,a) V2(v1,c2,c1,o1) 
  //* S2(i,k,m,a)  <--  (   -2.00000000) D2(i,m,k,o1) X(o1,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    T2b = T2.get_amp2(ia);
    orz::DTensor X(nocc);
    orz::DTensor Xa = orz::ct::sympack_Xa(symblockinfo, sa, X);
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
      FC_FUNC(g_if_sigma_ooov_ccvv_no0_x0, G_IF_SIGMA_OOOV_CCVV_NO0_X0)
        (sa, ia, sv1, iv1, T2b.cptr(), V2_sym.cptr(), Xa.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_ooov_ccvv_no1_x0, G_IF_SIGMA_OOOV_CCVV_NO1_X0)
      (sa, ia, Xa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.1
  //* X(o1,a)  <--  (    1.00000000)  T2(c2,c1,v1,a) V2(v1,c1,c2,o1) 
  //* S2(i,k,m,a)  <--  (    1.00000000) D2(i,m,k,o1) X(o1,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    T2b = T2.get_amp2(ia);
    orz::DTensor X(nocc);
    orz::DTensor Xa = orz::ct::sympack_Xa(symblockinfo, sa, X);
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
      FC_FUNC(g_if_sigma_ooov_ccvv_no0_x1, G_IF_SIGMA_OOOV_CCVV_NO0_X1)
        (sa, ia, sv1, iv1, T2b.cptr(), V2_sym.cptr(), Xa.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_ooov_ccvv_no1_x1, G_IF_SIGMA_OOOV_CCVV_NO1_X1)
      (sa, ia, Xa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.2
  //* X(o1,a)  <--  (    1.00000000)  T2(c2,c1,a,v1) V2(v1,c2,c1,o1) 
  //* S2(i,k,m,a)  <--  (    1.00000000) D2(i,m,k,o1) X(o1,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor X(nocc);
    orz::DTensor Xa = orz::ct::sympack_Xa(symblockinfo, sa, X);
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
      FC_FUNC(g_if_sigma_ooov_ccvv_no0_x2, G_IF_SIGMA_OOOV_CCVV_NO0_X2)
        (sa, ia, sv1, iv1, T2b.cptr(), V2_sym.cptr(), Xa.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_ooov_ccvv_no1_x2, G_IF_SIGMA_OOOV_CCVV_NO1_X2)
      (sa, ia, Xa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }


  {
  // No.3
  //* X(o1,a)  <--  (    1.00000000)  T2(c1,c2,a,v1) V2(v1,c2,c1,o1) 
  //* S2(i,k,m,a)  <--  (   -2.00000000) D2(i,m,k,o1) X(o1,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    S2b = orz::DTensor(retval.namps_iamp()[ia]);
    orz::DTensor X(nocc);
    orz::DTensor Xa = orz::ct::sympack_Xa(symblockinfo, sa, X);
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
      FC_FUNC(g_if_sigma_ooov_ccvv_no0_x3, G_IF_SIGMA_OOOV_CCVV_NO0_X3)
        (sa, ia, sv1, iv1, T2b.cptr(), V2_sym.cptr(), Xa.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_ooov_ccvv_no1_x3, G_IF_SIGMA_OOOV_CCVV_NO1_X3)
      (sa, ia, Xa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, S2b);
  }
  }
  }

  return retval; 
} 
