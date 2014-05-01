                                                                                
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
#include <sci/ctnew2/c_diag_ooov.h>                                            
                                                                                
using std::cout;                                                                
using std::endl;                                                                
                                                                                
#define FLOPCOUNT                                                               
                                                                                
//  ___________                __               
//  \_   _____/____    _____ _/  |_  ____      
//   |    __)_/ __ \  /     \\   __\/  _ \ 
//   |     \ \  ___/ |  Y Y  \|  | (  <_> )  
//   \___  /  \___  >|__|_|  /|__|  \____/   
//       \/       \/       \/                

                                                                                
// ***************************************************************************  
// orz::ct::mrci                                                                
// ***************************************************************************  
									    /*!        
   @brief CT input                                                              
									     */        
                                                                                
orz::ct::BareAmpPack orz::ct::diag_ooov(const orz::ct::Input &ctinp,                                    
					 const orz::ct::SymBlockInfo &symblockinfo,                                
					 const orz::ct::HintMO &hintmo,                                            
					 const orz::ct::RdmPack &rdmPack_sym,                                      
					 const orz::DTensor &rdm4,                                                 
					 const int num_sigma,
					 const double Ecas) {
                                                                                                                
                                                                                                                
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
  //* Hdiag(i,k,m,a)  <--  (    1.00000000) D3(i,m,k,k,o1,i) h(m,o1) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ooov_no0_x0, G_IF_DIAG_OOOV_NO0_X0)
      (sa, ia, moint1_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.1
  //* Hdiag(i,k,m,a)  <--  (    1.00000000) D2(i,i,k,k) h(m,m) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ooov_no0_x1, G_IF_DIAG_OOOV_NO0_X1)
      (sa, ia, moint1_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.2
  //* Hdiag(i,k,m,a)  <--  (   -1.00000000) D3(i,m,k,k,m,o1) h(i,o1) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ooov_no0_x2, G_IF_DIAG_OOOV_NO0_X2)
      (sa, ia, moint1_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.3
  //* Hdiag(i,k,m,a)  <--  (   -1.00000000) D2(i,o1,k,k) h(i,o1) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ooov_no0_x3, G_IF_DIAG_OOOV_NO0_X3)
      (sa, ia, moint1_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.4
  //* Hdiag(i,k,m,a)  <--  (   -1.00000000) D3(i,m,k,o1,m,i) h(k,o1) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ooov_no0_x4, G_IF_DIAG_OOOV_NO0_X4)
      (sa, ia, moint1_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.5
  //* Hdiag(i,k,m,a)  <--  (   -1.00000000) D2(i,i,k,o1) h(k,o1) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ooov_no0_x5, G_IF_DIAG_OOOV_NO0_X5)
      (sa, ia, moint1_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.6
  //* Hdiag(i,k,m,a)  <--  (    1.00000000) D3(i,m,k,k,m,i) h(a,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ooov_no0_x6, G_IF_DIAG_OOOV_NO0_X6)
      (sa, ia, moint1_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.7
  //* Hdiag(i,k,m,a)  <--  (    1.00000000) D2(i,i,k,k) h(a,a) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ooov_no0_x7, G_IF_DIAG_OOOV_NO0_X7)
      (sa, ia, moint1_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.8
  //* Hdiag(i,k,m,a)  <--  (    2.00000000) D3(i,m,k,k,o1,i) Y0(m,o1) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nocc, nocc);
  orz::DTensor Y0 = orz::ct::sympack_Xaa(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_diag_ooov_y8, G_IF_DIAG_OOOV_Y8)
      (sc1, ic1, V2_sym.cptr(), Y0.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ooov_no0_x8, G_IF_DIAG_OOOV_NO0_X8)
      (sa, ia, Y0.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.9
  //* Hdiag(i,k,m,a)  <--  (   -1.00000000) D3(i,m,k,k,o1,i) Y1(m,o1) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nocc, nocc);
  orz::DTensor Y1 = orz::ct::sympack_Xaa(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_diag_ooov_y9, G_IF_DIAG_OOOV_Y9)
      (sc1, ic1, V2_sym.cptr(), Y1.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ooov_no0_x9, G_IF_DIAG_OOOV_NO0_X9)
      (sa, ia, Y1.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.10
  //* Hdiag(i,k,m,a)  <--  (    2.00000000) D3(i,m,k,k,o1,i) V2(a,a,m,o1) 
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
    FC_FUNC(g_if_diag_ooov_no0_x10, G_IF_DIAG_OOOV_NO0_X10)
      (sa, ia, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.11
  //* Hdiag(i,k,m,a)  <--  (    2.00000000) D2(i,i,k,k) Y2(m,m) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nocc, nocc);
  orz::DTensor Y2 = orz::ct::sympack_Xaa(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_diag_ooov_y11, G_IF_DIAG_OOOV_Y11)
      (sc1, ic1, V2_sym.cptr(), Y2.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ooov_no0_x11, G_IF_DIAG_OOOV_NO0_X11)
      (sa, ia, Y2.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.12
  //* Hdiag(i,k,m,a)  <--  (   -1.00000000) D2(i,i,k,k) Y3(m,m) 
  // The effective tensor is detected .... 
  orz::DTensor Y(nocc, nocc);
  orz::DTensor Y3 = orz::ct::sympack_Xaa(symblockinfo, 0, Y);
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
    FC_FUNC(g_if_diag_ooov_y12, G_IF_DIAG_OOOV_Y12)
      (sc1, ic1, V2_sym.cptr(), Y3.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ooov_no0_x12, G_IF_DIAG_OOOV_NO0_X12)
      (sa, ia, Y3.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.13
  //* Hdiag(i,k,m,a)  <--  (    2.00000000) D2(i,i,k,k) Y4(a,a) 
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
    FC_FUNC(g_if_diag_ooov_y13, G_IF_DIAG_OOOV_Y13)
      (sc1, ic1, V2_sym.cptr(), Y4.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ooov_no0_x13, G_IF_DIAG_OOOV_NO0_X13)
      (sa, ia, Y4.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.14
  //* Hdiag(i,k,m,a)  <--  (   -1.00000000) D2(i,i,k,k) Y5(a,a) 
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
    FC_FUNC(g_if_diag_ooov_y14, G_IF_DIAG_OOOV_Y14)
      (sc1, ic1, V2_sym.cptr(), Y5.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ooov_no0_x14, G_IF_DIAG_OOOV_NO0_X14)
      (sa, ia, Y5.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.15
  //* Hdiag(i,k,m,a)  <--  (    1.00000000) D2(i,i,k,k) V2(a,a,m,m) 
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
    FC_FUNC(g_if_diag_ooov_no0_x15, G_IF_DIAG_OOOV_NO0_X15)
      (sa, ia, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.16
  //* Hdiag(i,k,m,a)  <--  (   -2.00000000) D3(i,m,k,k,m,o1) Y6(i,o1) 
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
    FC_FUNC(g_if_diag_ooov_y16, G_IF_DIAG_OOOV_Y16)
      (sc1, ic1, V2_sym.cptr(), Y6.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ooov_no0_x16, G_IF_DIAG_OOOV_NO0_X16)
      (sa, ia, Y6.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.17
  //* Hdiag(i,k,m,a)  <--  (    1.00000000) D3(i,m,k,k,m,o1) Y7(i,o1) 
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
    FC_FUNC(g_if_diag_ooov_y17, G_IF_DIAG_OOOV_Y17)
      (sc1, ic1, V2_sym.cptr(), Y7.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ooov_no0_x17, G_IF_DIAG_OOOV_NO0_X17)
      (sa, ia, Y7.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.18
  //* Hdiag(i,k,m,a)  <--  (   -2.00000000) D2(i,o1,k,k) Y8(i,o1) 
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
    FC_FUNC(g_if_diag_ooov_y18, G_IF_DIAG_OOOV_Y18)
      (sc1, ic1, V2_sym.cptr(), Y8.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ooov_no0_x18, G_IF_DIAG_OOOV_NO0_X18)
      (sa, ia, Y8.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.19
  //* Hdiag(i,k,m,a)  <--  (    1.00000000) D2(i,o1,k,k) Y9(i,o1) 
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
    FC_FUNC(g_if_diag_ooov_y19, G_IF_DIAG_OOOV_Y19)
      (sc1, ic1, V2_sym.cptr(), Y9.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ooov_no0_x19, G_IF_DIAG_OOOV_NO0_X19)
      (sa, ia, Y9.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.20
  //* Hdiag(i,k,m,a)  <--  (   -2.00000000) D3(i,m,k,o1,m,i) Y10(k,o1) 
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
    FC_FUNC(g_if_diag_ooov_y20, G_IF_DIAG_OOOV_Y20)
      (sc1, ic1, V2_sym.cptr(), Y10.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ooov_no0_x20, G_IF_DIAG_OOOV_NO0_X20)
      (sa, ia, Y10.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.21
  //* Hdiag(i,k,m,a)  <--  (    1.00000000) D3(i,m,k,o1,m,i) Y11(k,o1) 
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
    FC_FUNC(g_if_diag_ooov_y21, G_IF_DIAG_OOOV_Y21)
      (sc1, ic1, V2_sym.cptr(), Y11.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ooov_no0_x21, G_IF_DIAG_OOOV_NO0_X21)
      (sa, ia, Y11.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.22
  //* Hdiag(i,k,m,a)  <--  (   -2.00000000) D2(i,i,k,o1) Y12(k,o1) 
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
    FC_FUNC(g_if_diag_ooov_y22, G_IF_DIAG_OOOV_Y22)
      (sc1, ic1, V2_sym.cptr(), Y12.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ooov_no0_x22, G_IF_DIAG_OOOV_NO0_X22)
      (sa, ia, Y12.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.23
  //* Hdiag(i,k,m,a)  <--  (    1.00000000) D2(i,i,k,o1) Y13(k,o1) 
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
    FC_FUNC(g_if_diag_ooov_y23, G_IF_DIAG_OOOV_Y23)
      (sc1, ic1, V2_sym.cptr(), Y13.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ooov_no0_x23, G_IF_DIAG_OOOV_NO0_X23)
      (sa, ia, Y13.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.24
  //* Hdiag(i,k,m,a)  <--  (    2.00000000) D3(i,m,k,k,m,i) Y14(a,a) 
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
    FC_FUNC(g_if_diag_ooov_y24, G_IF_DIAG_OOOV_Y24)
      (sc1, ic1, V2_sym.cptr(), Y14.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ooov_no0_x24, G_IF_DIAG_OOOV_NO0_X24)
      (sa, ia, Y14.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.25
  //* Hdiag(i,k,m,a)  <--  (   -1.00000000) D3(i,m,k,k,m,i) Y15(a,a) 
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
    FC_FUNC(g_if_diag_ooov_y25, G_IF_DIAG_OOOV_Y25)
      (sc1, ic1, V2_sym.cptr(), Y15.cptr(), nir, nsym, psym);
  }
  }
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ooov_no0_x25, G_IF_DIAG_OOOV_NO0_X25)
      (sa, ia, Y15.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.26
  //* Hdiag(i,k,m,a)  <--  (    1.00000000) D4(m,i,k,k,o3,o1,i,o2) V2(m,o2,o1,o3) 
  for(int sm = 0;sm < nir;++sm){ 
  for(int im = symblockinfo.psym()(sm,I_O,I_BEGIN);im <= symblockinfo.psym()(sm,I_O,I_END);++im){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(im);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, im, sm, V2); // V2=(IR-COV index) 
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
      for(int si = 0;si < nir;++si){ 
      for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
        // Load D4 from disk, or GA ....                                                     
        int imoi = amo2imo[im] - nclosed;                              
        int imoj = amo2imo[ii] - nclosed;                              
                                                                                             
        orz::DTensor rdm4_ij_sliced = rdm4(orz::Slice(imoi,imoi+1), orz::Slice(imoj,imoj+1), 
                                           orz::Slice(),            orz::Slice(),            
                                           orz::Slice(),            orz::Slice(),            
                                           orz::Slice(),            orz::Slice()).copy();    
        rdm4_sym = orz::ct::sympack_rdm4_2(symblockinfo, im, sm, ii, si, rdm4_ij_sliced);    
        FC_FUNC(g_if_set_d4,G_IF_SET_D4)(sm, si, im, ii, rdm4_sym.cptr(), nir, nsym, psym);  
        FC_FUNC(g_if_diag_ooov_no0_x26, G_IF_DIAG_OOOV_NO0_X26)
          (sa, ia, si, ii, sm, im, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
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
  // No.27
  //* Hdiag(i,k,m,a)  <--  (    1.00000000) D3(i,i,k,k,o1,o2) V2(m,m,o1,o2) 
  for(int sm = 0;sm < nir;++sm){ 
  for(int im = symblockinfo.psym()(sm,I_O,I_BEGIN);im <= symblockinfo.psym()(sm,I_O,I_END);++im){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(im);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, im, sm, V2); // V2=(IR-COV index) 
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
      FC_FUNC(g_if_diag_ooov_no0_x27, G_IF_DIAG_OOOV_NO0_X27)
        (sa, ia, sm, im, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(ia, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.28
  //* Hdiag(i,k,m,a)  <--  (    1.00000000) D3(i,o1,k,k,o2,i) V2(m,o1,m,o2) 
  for(int sm = 0;sm < nir;++sm){ 
  for(int im = symblockinfo.psym()(sm,I_O,I_BEGIN);im <= symblockinfo.psym()(sm,I_O,I_END);++im){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(im);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, im, sm, V2); // V2=(IR-COV index) 
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
      FC_FUNC(g_if_diag_ooov_no0_x28, G_IF_DIAG_OOOV_NO0_X28)
        (sa, ia, sm, im, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(ia, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.29
  //* Hdiag(i,k,m,a)  <--  (   -1.00000000) D4(i,m,k,k,m,o2,o3,o1) V2(i,o2,o1,o3) 
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(ii);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, ii, si, V2); // V2=(IR-COV index) 
    for(int sm = 0;sm < nir;++sm){ 
    for(int im = symblockinfo.psym()(sm,I_O,I_BEGIN);im <= symblockinfo.psym()(sm,I_O,I_END);++im){ 
      for(int sa = 0;sa < nir;++sa){ 
      for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
        Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
        // Load D4 from disk, or GA ....                                                     
        int imoi = amo2imo[ii] - nclosed;                              
        int imoj = amo2imo[im] - nclosed;                              
                                                                                             
        orz::DTensor rdm4_ij_sliced = rdm4(orz::Slice(imoi,imoi+1), orz::Slice(imoj,imoj+1), 
                                           orz::Slice(),            orz::Slice(),            
                                           orz::Slice(),            orz::Slice(),            
                                           orz::Slice(),            orz::Slice()).copy();    
        rdm4_sym = orz::ct::sympack_rdm4_2(symblockinfo, ii, si, im, sm, rdm4_ij_sliced);    
        FC_FUNC(g_if_set_d4,G_IF_SET_D4)(si, sm, ii, im, rdm4_sym.cptr(), nir, nsym, psym);  
        FC_FUNC(g_if_diag_ooov_no0_x29, G_IF_DIAG_OOOV_NO0_X29)
          (sa, ia, si, ii, sm, im, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
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
  // No.30
  //* Hdiag(i,k,m,a)  <--  (   -1.00000000) D3(i,o2,k,k,o3,o1) V2(i,o2,o1,o3) 
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(ii);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, ii, si, V2); // V2=(IR-COV index) 
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
      FC_FUNC(g_if_diag_ooov_no0_x30, G_IF_DIAG_OOOV_NO0_X30)
        (sa, ia, si, ii, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(ia, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.31
  //* Hdiag(i,k,m,a)  <--  (   -1.00000000) D3(i,o2,k,k,m,o1) V2(i,o1,m,o2) 
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(ii);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, ii, si, V2); // V2=(IR-COV index) 
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
      FC_FUNC(g_if_diag_ooov_no0_x31, G_IF_DIAG_OOOV_NO0_X31)
        (sa, ia, si, ii, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(ia, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.32
  //* Hdiag(i,k,m,a)  <--  (   -1.00000000) D4(k,o2,i,m,m,i,o3,o1) V2(k,o2,o1,o3) 
  for(int sk = 0;sk < nir;++sk){ 
  for(int ik = symblockinfo.psym()(sk,I_O,I_BEGIN);ik <= symblockinfo.psym()(sk,I_O,I_END);++ik){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(ik);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, ik, sk, V2); // V2=(IR-COV index) 
    for(int so2 = 0;so2 < nir;++so2){ 
    for(int io2 = symblockinfo.psym()(so2,I_O,I_BEGIN);io2 <= symblockinfo.psym()(so2,I_O,I_END);++io2){ 
      for(int sa = 0;sa < nir;++sa){ 
      for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
        Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
        // Load D4 from disk, or GA ....                                                     
        int imoi = amo2imo[ik] - nclosed;                              
        int imoj = amo2imo[io2] - nclosed;                              
                                                                                             
        orz::DTensor rdm4_ij_sliced = rdm4(orz::Slice(imoi,imoi+1), orz::Slice(imoj,imoj+1), 
                                           orz::Slice(),            orz::Slice(),            
                                           orz::Slice(),            orz::Slice(),            
                                           orz::Slice(),            orz::Slice()).copy();    
        rdm4_sym = orz::ct::sympack_rdm4_2(symblockinfo, ik, sk, io2, so2, rdm4_ij_sliced);    
        FC_FUNC(g_if_set_d4,G_IF_SET_D4)(sk, so2, ik, io2, rdm4_sym.cptr(), nir, nsym, psym);  
        FC_FUNC(g_if_diag_ooov_no0_x32, G_IF_DIAG_OOOV_NO0_X32)
          (sa, ia, sk, ik, so2, io2, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
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
  // No.33
  //* Hdiag(i,k,m,a)  <--  (   -1.00000000) D3(i,i,k,o2,o3,o1) V2(k,o2,o1,o3) 
  for(int sk = 0;sk < nir;++sk){ 
  for(int ik = symblockinfo.psym()(sk,I_O,I_BEGIN);ik <= symblockinfo.psym()(sk,I_O,I_END);++ik){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(ik);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, ik, sk, V2); // V2=(IR-COV index) 
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
      FC_FUNC(g_if_diag_ooov_no0_x33, G_IF_DIAG_OOOV_NO0_X33)
        (sa, ia, sk, ik, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(ia, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.34
  //* Hdiag(i,k,m,a)  <--  (   -1.00000000) D3(i,m,o1,k,o2,i) V2(k,o1,m,o2) 
  for(int sk = 0;sk < nir;++sk){ 
  for(int ik = symblockinfo.psym()(sk,I_O,I_BEGIN);ik <= symblockinfo.psym()(sk,I_O,I_END);++ik){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(ik);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, ik, sk, V2); // V2=(IR-COV index) 
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
      FC_FUNC(g_if_diag_ooov_no0_x34, G_IF_DIAG_OOOV_NO0_X34)
        (sa, ia, sk, ik, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(ia, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.35
  //* Hdiag(i,k,m,a)  <--  (   -1.00000000) D3(i,m,k,o2,m,o1) V2(i,o1,k,o2) 
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(ii);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, ii, si, V2); // V2=(IR-COV index) 
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
      FC_FUNC(g_if_diag_ooov_no0_x35, G_IF_DIAG_OOOV_NO0_X35)
        (sa, ia, si, ii, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(ia, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.36
  //* Hdiag(i,k,m,a)  <--  (   -1.00000000) D2(i,o2,k,o1) V2(i,o2,k,o1) 
  for(int si = 0;si < nir;++si){ 
  for(int ii = symblockinfo.psym()(si,I_O,I_BEGIN);ii <= symblockinfo.psym()(si,I_O,I_END);++ii){ 
    // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
    V2 <<= 0.0;                                                                          
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(ii);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
      // Load a signle record of integals                                                             
      const int &imo2 = loadbuf_ptr->i0;                                                              
      const int &imo3 = loadbuf_ptr->i1;                                                              
      const int &imo4 = loadbuf_ptr->i2;                                                              
      const double &v = loadbuf_ptr->v;                                                               
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
    }                                                                                                 
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, ii, si, V2); // V2=(IR-COV index) 
    for(int sa = 0;sa < nir;++sa){ 
    for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
      Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
      FC_FUNC(g_if_diag_ooov_no0_x36, G_IF_DIAG_OOOV_NO0_X36)
        (sa, ia, si, ii, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(ia, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.37
  //* Hdiag(i,k,m,a)  <--  (    1.00000000) D4(i,m,k,k,m,i,o1,o2) V2(a,a,o1,o2) 
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
        FC_FUNC(g_if_diag_ooov_no0_x37, G_IF_DIAG_OOOV_NO0_X37)
          (sa, ia, si, ii, sm, im, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
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
  // No.38
  //* Hdiag(i,k,m,a)  <--  (    1.00000000) D3(i,i,k,k,o1,o2) V2(a,a,o1,o2) 
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
    FC_FUNC(g_if_diag_ooov_no0_x38, G_IF_DIAG_OOOV_NO0_X38)
      (sa, ia, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.39
  //* Hdiag(i,k,m,a)  <--  (    1.00000000) D4(i,m,k,o1,m,i,o2,k) V2(a,o1,o2,a) 
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
        FC_FUNC(g_if_diag_ooov_no0_x39, G_IF_DIAG_OOOV_NO0_X39)
          (sa, ia, si, ii, sm, im, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
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
  // No.40
  //* Hdiag(i,k,m,a)  <--  (    1.00000000) D3(i,i,k,o1,o2,k) V2(a,o1,o2,a) 
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
    FC_FUNC(g_if_diag_ooov_no0_x40, G_IF_DIAG_OOOV_NO0_X40)
      (sa, ia, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.41
  //* Hdiag(i,k,m,a)  <--  (    2.00000000) D3(i,k,k,o1,m,i) V2(a,m,o1,a) 
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
    FC_FUNC(g_if_diag_ooov_no0_x41, G_IF_DIAG_OOOV_NO0_X41)
      (sa, ia, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.42
  //* Hdiag(i,k,m,a)  <--  (    1.00000000) D2(i,k,k,i) V2(a,m,m,a) 
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
    FC_FUNC(g_if_diag_ooov_no0_x42, G_IF_DIAG_OOOV_NO0_X42)
      (sa, ia, V2_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.43
  //* Hdiag(i,k,m,a)  <--  (    1.00000000) Ecas D3(i,m,k,k,m,i) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ooov_no0_x43, G_IF_DIAG_OOOV_NO0_X43)
      (sa, ia, &Ecas, Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }


  {
  // No.44
  //* Hdiag(i,k,m,a)  <--  (    1.00000000) Ecas D2(i,i,k,k) 
  for(int sa = 0;sa < nir;++sa){ 
  for(int ia = symblockinfo.psym()(sa,I_V,I_BEGIN);ia <= symblockinfo.psym()(sa,I_V,I_END);++ia){ 
    Hdiagb = orz::DTensor(retval.namps_iamp()[ia]);
    FC_FUNC(g_if_diag_ooov_no0_x44, G_IF_DIAG_OOOV_NO0_X44)
      (sa, ia, &Ecas, Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(ia, Hdiagb);
  }
  }
  }

  return retval; 
} 
