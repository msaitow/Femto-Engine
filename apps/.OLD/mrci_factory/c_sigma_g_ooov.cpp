                                                                                
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
#include <sci/ctnew2/c_sigma_g_ooov.h>                                            
                                                                                
using std::cout;                                                                
using std::endl;                                                                
                                                                                
#define FLOPCOUNT                                                               
                                                                                
//    o__ __o__/_                            o                         
//   <|    v                                <|>                        
//   < >                                    < >                        
//    |         o__  __o   \o__ __o__ __o    |        o__ __o         
//    o__/_    /v      |>   |     |     |>   o__/_   /v     v\        
//    |       />      //   / \   / \   / \   |      />       <\    
//   <o>      \o    o/     \o/   \o/   \o/   |      \         /   
//    |        v\  /v __o   |     |     |    o       o       o        
//   / \        <\/> __/>  / \   / \   / \   <\__    <\__ __/>  

                                                                                
// ***************************************************************************  
// orz::ct::mrci                                                                
// ***************************************************************************  
									    /*!        
   @brief CT input                                                              
									     */        
                                                                                
double orz::ct::sigma_g_ooov(const orz::ct::Input &ctinp,                                                  
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
  //* X(o2,o1,o4,o3)  <--  (    1.00000000)  T2(o2,o1,o4,v1) h(o3,v1) 
  //* S0()  <--  (    1.00000000) D2(o1,o3,o2,o4) X(o2,o1,o4,o3) 
  orz::DTensor X(nocc, nocc, nocc, nocc);
  orz::DTensor Xaaaa = orz::ct::sympack_Xaaaa(symblockinfo, 0, X);
  for(int sv1 = 0;sv1 < nir;++sv1){ 
  for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
    T2b = T2.get_amp2(iv1);
    FC_FUNC(g_if_sigma_g_ooov_no0_x0, G_IF_SIGMA_G_OOOV_NO0_X0)
      (sv1, iv1, T2b.cptr(), moint1_sym.cptr(), Xaaaa.cptr(), nir, nsym, psym);
  }
  }
  FC_FUNC(g_if_sigma_g_ooov_no1_x0, G_IF_SIGMA_G_OOOV_NO1_X0)
    (Xaaaa.cptr(), &S0, nir, nsym, psym);
  }


  {
  // No.1
  //* X(o3,o4,o1,o2)  <--  (    1.00000000)  T2(o3,o4,o1,v1) Y0(o2,v1) 
  //* S0()  <--  (    2.00000000) D2(o1,o3,o2,o4) X(o3,o4,o1,o2) 
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
    FC_FUNC(g_if_sigma_g_ooov_y1, G_IF_SIGMA_G_OOOV_Y1)
      (sc1, ic1, V2_sym.cptr(), Y0.cptr(), nir, nsym, psym);
  }
  }
  orz::DTensor X(nocc, nocc, nocc, nocc);
  orz::DTensor Xaaaa = orz::ct::sympack_Xaaaa(symblockinfo, 0, X);
  for(int sv1 = 0;sv1 < nir;++sv1){ 
  for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
    T2b = T2.get_amp2(iv1);
    FC_FUNC(g_if_sigma_g_ooov_no0_x1, G_IF_SIGMA_G_OOOV_NO0_X1)
      (sv1, iv1, T2b.cptr(), Y0.cptr(), Xaaaa.cptr(), nir, nsym, psym);
  }
  }
  FC_FUNC(g_if_sigma_g_ooov_no1_x1, G_IF_SIGMA_G_OOOV_NO1_X1)
    (Xaaaa.cptr(), &S0, nir, nsym, psym);
  }


  {
  // No.2
  //* X(o4,o3,o2,o1)  <--  (    1.00000000)  T2(o4,o3,o2,v1) Y1(o1,v1) 
  //* S0()  <--  (   -1.00000000) D2(o1,o3,o2,o4) X(o4,o3,o2,o1) 
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
    FC_FUNC(g_if_sigma_g_ooov_y2, G_IF_SIGMA_G_OOOV_Y2)
      (sc1, ic1, V2_sym.cptr(), Y1.cptr(), nir, nsym, psym);
  }
  }
  orz::DTensor X(nocc, nocc, nocc, nocc);
  orz::DTensor Xaaaa = orz::ct::sympack_Xaaaa(symblockinfo, 0, X);
  for(int sv1 = 0;sv1 < nir;++sv1){ 
  for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
    T2b = T2.get_amp2(iv1);
    FC_FUNC(g_if_sigma_g_ooov_no0_x2, G_IF_SIGMA_G_OOOV_NO0_X2)
      (sv1, iv1, T2b.cptr(), Y1.cptr(), Xaaaa.cptr(), nir, nsym, psym);
  }
  }
  FC_FUNC(g_if_sigma_g_ooov_no1_x2, G_IF_SIGMA_G_OOOV_NO1_X2)
    (Xaaaa.cptr(), &S0, nir, nsym, psym);
  }


  {
  // No.3
  //* X(o1,o2,o4,o5,o3,o6)  <--  (    1.00000000)  T2(o1,o2,o4,v1) V2(v1,o5,o3,o6) 
  //* S0()  <--  (    1.00000000) D3(o1,o4,o2,o5,o3,o6) X(o1,o2,o4,o5,o3,o6) 
  orz::DTensor X(nocc, nocc, nocc, nocc, nocc, nocc);
  orz::DTensor Xaaaaaa = orz::ct::sympack_Xaaaaaa(symblockinfo, 0, X);
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
    FC_FUNC(g_if_sigma_g_ooov_no0_x3, G_IF_SIGMA_G_OOOV_NO0_X3)
      (sv1, iv1, T2b.cptr(), V2_sym.cptr(), Xaaaaaa.cptr(), nir, nsym, psym);
  }
  }
  FC_FUNC(g_if_sigma_g_ooov_no1_x3, G_IF_SIGMA_G_OOOV_NO1_X3)
    (Xaaaaaa.cptr(), &S0, nir, nsym, psym);
  }


  {
  // No.4
  //* X(o1,o2,o4,o3)  <--  (    1.00000000)  T2(o1,o2,o5,v1) V2(v1,o4,o3,o5) 
  //* S0()  <--  (    1.00000000) D2(o1,o3,o2,o4) X(o1,o2,o4,o3) 
  orz::DTensor X(nocc, nocc, nocc, nocc);
  orz::DTensor Xaaaa = orz::ct::sympack_Xaaaa(symblockinfo, 0, X);
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
    FC_FUNC(g_if_sigma_g_ooov_no0_x4, G_IF_SIGMA_G_OOOV_NO0_X4)
      (sv1, iv1, T2b.cptr(), V2_sym.cptr(), Xaaaa.cptr(), nir, nsym, psym);
  }
  }
  FC_FUNC(g_if_sigma_g_ooov_no1_x4, G_IF_SIGMA_G_OOOV_NO1_X4)
    (Xaaaa.cptr(), &S0, nir, nsym, psym);
  }

  return  S0;
} 
