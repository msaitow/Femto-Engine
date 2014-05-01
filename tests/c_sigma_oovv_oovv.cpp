                                                                                
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
                                                                                
//Timing object                                                                 
extern std::vector<boost::tuple<std::string, double, double, double> > my_timer;
extern double Fc0;                                                              

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
                                                                                
orz::ct::BareAmpPack orz::ct::sigma_oovv_oovv(const orz::ct::Input &ctinp,                                    
                                  const orz::ct::SymBlockInfo &symblockinfo,                                
                                  const orz::ct::HintMO &hintmo,                                            
                                  const int alloc_type,                                            
                                  const orz::ct::RdmPack &rdmPack,                                      
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
                                                                                                                 
  std::ostringstream stm;                                                                                        
  stm << num_sigma;                                                                                              
  std::string name_of_sigma = "S2" + stm.str() + "]"; // Name of the Sigma vector  
  orz::ct::BareAmpPack retval                                                                                    
    = orz::ct::BareAmpPack(ctinp, symblockinfo, name_of_sigma, alloc_type); // Sigma(a, a', e, e') tensor        
                                                                                                                 
  orz::DTensor S2b; // Container of S2_aae,[b] tensor                                   
                                                                                                                 
  orz::DTensor T2b; // Container of T2_aae,[b] tensor                                             
  orz::DTensor rdm4_sym;                                                                                         
  orz::DTensor rdm4_ij_sliced(ctinp.use_d4cum_of() ? nocc*nocc*nocc*nocc*nocc*nocc : 0);                         
  // set nproc, myrank                      
  const int nproc = orz::world().size();    
  const int myrank = orz::world().rank();   

  orz::DTensor moint1 = hintmo.int1(); // Setting up one-body integrals                                         
  const orz::DTensor moint1_sym = (myrank == 0) ? orz::ct::sympack_int1(symblockinfo, moint1) : orz::DTensor(); // moint1=(IR-COV index)
  orz::DTensor V2(nmo,nmo,nmo);                                                                    
  double * const V2_ptr = V2.cptr();                                                  

  //*-- FEMTO begins --//*
  // Label : noeri
  {

  //*-- Entering to take the type 1 contractions --*//
  { 
  // No. 0, [1]
  // S2(i,k,a,c) += (    2.00000000) Fc0 D2(i,o1,k,o2) T2(o1,o2,a,c) 
  double flops = 0; // Flop count
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    T2b = T2.get_amp2(ic);
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    FC_FUNC(g_if_sigma_oovv_oovv_no0_x0_type1_noeri, G_IF_SIGMA_OOVV_OOVV_NO0_X0_TYPE1_NOERI)
      (sc, ic, &Fc0, T2b.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
    retval.acc_amp2(ic, S2b);
  }
  }
  } // End scope

  { 
  // No. 1, [2]
  // W0(i,o1,k,o2) += (    1.00000000) D3(i,o1,k,o2,o3,o4) Fc1(o3,o4) 
  // S2(i,k,a,c) += (    2.00000000) T2(o1,o2,a,c) W0(i,o1,k,o2) 
  double flops = 0; // Flop count
  orz::DTensor W0aaaa(orz::ct::sizeof_sympack_Xaaaa(symblockinfo, 0));
  FC_FUNC(g_if_sigma_oovv_oovv_no0_x1_type1_noeri, G_IF_SIGMA_OOVV_OOVV_NO0_X1_TYPE1_NOERI)
    (W0aaaa.cptr(), nir, nsym, psym, &flops);
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    FC_FUNC(g_if_sigma_oovv_oovv_no1_x1_type1_noeri, G_IF_SIGMA_OOVV_OOVV_NO1_X1_TYPE1_NOERI)
      (sc, ic, T2b.cptr(), W0aaaa.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
    retval.acc_amp2(ic, S2b);
  }
  }
  } // End scope

  { 
  // No. 2, [2]
  // W1(i,k,a,v1) += (    1.00000000) D2(i,o1,k,o2) T2(o1,o2,a,v1) 
  // S2(i,k,a,c) += (    2.00000000) Fc1(c,v1) W1(i,k,a,v1) 
  double flops = 0; // Flop count
  for(int sv1 = 0;sv1 < nir;++sv1){ 
  for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
    T2b = T2.get_amp2(iv1);
    orz::DTensor W1aav(orz::ct::sizeof_sympack_Xaav(symblockinfo, sv1));
    FC_FUNC(g_if_sigma_oovv_oovv_no0_x2_type1_noeri, G_IF_SIGMA_OOVV_OOVV_NO0_X2_TYPE1_NOERI)
      (sv1, iv1, T2b.cptr(), W1aav.cptr(), nir, nsym, psym, &flops);
    for(int sc = 0;sc < nir;++sc){ 
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
      S2b = orz::DTensor(retval.namps_iamp()[ic]);
      FC_FUNC(g_if_sigma_oovv_oovv_no1_x2_type1_noeri, G_IF_SIGMA_OOVV_OOVV_NO1_X2_TYPE1_NOERI)
        (sc, ic, sv1, iv1, W1aav.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
      retval.acc_amp2(ic, S2b);
    }
    }
  }
  }
  } // End scope

  { 
  // No. 3, [2]
  // W2(i,k,v1,c) += (    1.00000000) D2(i,o1,k,o2) T2(o1,o2,v1,c) 
  // S2(i,k,a,c) += (    2.00000000) Fc1(a,v1) W2(i,k,v1,c) 
  double flops = 0; // Flop count
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    T2b = T2.get_amp2(ic);
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    orz::DTensor W2aav(orz::ct::sizeof_sympack_Xaav(symblockinfo, sc));
    FC_FUNC(g_if_sigma_oovv_oovv_no0_x3_type1_noeri, G_IF_SIGMA_OOVV_OOVV_NO0_X3_TYPE1_NOERI)
      (sc, ic, T2b.cptr(), W2aav.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_sigma_oovv_oovv_no1_x3_type1_noeri, G_IF_SIGMA_OOVV_OOVV_NO1_X3_TYPE1_NOERI)
      (sc, ic, W2aav.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
    retval.acc_amp2(ic, S2b);
  }
  }
  } // End scope


  } // End femto
  //*-- FEMTO ends --//*

  //*-- FEMTO begins --//*
  // Label : eri_o
  {

  //  >> Intermediates for the type 2 contractions are defined here << 
  orz::DTensor W0aaaa(orz::ct::sizeof_sympack_Xaaaa(symblockinfo, 0));
  for(int s_eri = 0;s_eri < nir;++s_eri){ 
  for(int i_eri = symblockinfo.psym()(s_eri,I_O,I_BEGIN);i_eri <= symblockinfo.psym()(s_eri,I_O,I_END);++i_eri){ 
  if(hintmo.iproc_havingimo()[i_eri] == myrank) {           
  // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
  V2 <<= 0.0;                                                                          
  shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(i_eri);
  for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
    // Load a signle record of integals                                                             
    const int &imo2 = loadbuf_ptr->i0;                                                              
    const int &imo3 = loadbuf_ptr->i1;                                                              
    const int &imo4 = loadbuf_ptr->i2;                                                              
    const double &v = loadbuf_ptr->v;                                                               
    V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
  }                                                                                                 
  const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, i_eri, s_eri, V2); // V2=(IR-COV index) 

  //*-- Entering to take the type 1 contractions --*//
  { 
  // No. 0, [1]
  // W0(i,o3,k,o4) += (    1.00000000) D4(o1,o5,i,o3,k,o4,o2,o6) V2(o1,o5,o2,o6) 
  double flops = 0; // Flop count
  int so1 = s_eri;
  int io1 = i_eri;
  for(int so5 = 0;so5 < nir;++so5){ 
  for(int io5 = symblockinfo.psym()(so5,I_O,I_BEGIN);io5 <= symblockinfo.psym()(so5,I_O,I_END);++io5){ 

    int imoi = amo2imo[io1] - nclosed;                              
    int imoj = amo2imo[io5] - nclosed;                              
                                                                                         
    // Generate D4 by cumulant expansion ....                                            
    if(ctinp.use_d4cum_of()){
    FC_FUNC(f_mrci_rdm4_cumulant_partial_opt, F_MRCI_RDM4_CUMULANT_PARTIAL_OPT)                  
      (nocc, 0, rdmPack.rdm1().cptr(), rdmPack.rdm2().cptr(), rdmPack.cum2().cptr(), rdmPack.rdm3().cptr(),     
       rdm4_ij_sliced.cptr(), imoi, imoj);                                               
    rdm4_sym = orz::ct::sympack_rdm4_2(symblockinfo, io1, so1, io5, so5, rdm4_ij_sliced);    
    flops += nocc*nocc*nocc*nocc*nocc*nocc*70;
    } // End if
    // Slice the already existing 8-index 4-RDM ....                                            
    else{
    const double* rdm4_ij_sliced = rdm4.cptr() + (imoi*nocc+imoj)*nocc*nocc*nocc*nocc*nocc*nocc;
    rdm4_sym = orz::ct::sympack_rdm4_2x(symblockinfo, io1, so1, io5, so5, rdm4_ij_sliced);    
    }
    FC_FUNC(g_if_set_d4,G_IF_SET_D4)(so1, so5, io1, io5, rdm4_sym.cptr(), nir, nsym, psym);  
    FC_FUNC(g_if_sigma_oovv_oovv_no0_x0_type1_eri_o, G_IF_SIGMA_OOVV_OOVV_NO0_X0_TYPE1_ERI_O)
      (so1, io1, so5, io5, V2_sym.cptr(), W0aaaa.cptr(), nir, nsym, psym, &flops);
    FC_FUNC(g_if_unset_d4,G_IF_UNSET_D4)();
  }
  }
  } // End scope

  } // End my_rank
  }
  }

  //*-- Entering to take the type 2 contractions --*//
  { 
  // No. 0, [1]
  // S2(i,k,a,c) += (    1.00000000) T2(o3,o4,a,c) W0(i,o3,k,o4) 
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    FC_FUNC(g_if_sigma_oovv_oovv_no0_x0_type2_eri_o, G_IF_SIGMA_OOVV_OOVV_NO0_X0_TYPE2_ERI_O)
      (sc, ic, T2b.cptr(), W0aaaa.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
    retval.acc_amp2(ic, S2b);
  }
  }
  } // End scope


  } // End femto
  //*-- FEMTO ends --//*

  //*-- FEMTO begins --//*
  // Label : eri_v
  {

  for(int s_eri = 0;s_eri < nir;++s_eri){ 
  for(int i_eri = symblockinfo.psym()(s_eri,I_V,I_BEGIN);i_eri <= symblockinfo.psym()(s_eri,I_V,I_END);++i_eri){ 
  if(hintmo.iproc_havingimo()[i_eri] == myrank) {           
  // Load ERIs from somewhere, e.g. disk, GA, etc..                                                 
  V2 <<= 0.0;                                                                          
  shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(i_eri);
  for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           
    // Load a signle record of integals                                                             
    const int &imo2 = loadbuf_ptr->i0;                                                              
    const int &imo3 = loadbuf_ptr->i1;                                                              
    const int &imo4 = loadbuf_ptr->i2;                                                              
    const double &v = loadbuf_ptr->v;                                                               
    V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       
  }                                                                                                 
  const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, i_eri, s_eri, V2); // V2=(IR-COV index) 

  //*-- Entering to take the type 1 contractions --*//
  { 
  // No. 0, [2]
  // W0(i,o3,k,o2,c,v1) += (    1.00000000) D3(i,o3,k,o2,o4,o1) V2(c,v1,o1,o4) 
  // S2(i,k,a,c) += (    2.00000000) T2(o3,o2,a,v1) W0(i,o3,k,o2,c,v1) 
  double flops = 0; // Flop count
  int sc = s_eri;
  int ic = i_eri;
  S2b = orz::DTensor(retval.namps_iamp()[ic]);
  for(int sv1 = 0;sv1 < nir;++sv1){ 
  for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
    orz::DTensor W0aaaa(orz::ct::sizeof_sympack_Xaaaa(symblockinfo, sc^sv1));
    FC_FUNC(g_if_sigma_oovv_oovv_no0_x0_type1_eri_v, G_IF_SIGMA_OOVV_OOVV_NO0_X0_TYPE1_ERI_V)
      (sc, ic, sv1, iv1, V2_sym.cptr(), W0aaaa.cptr(), nir, nsym, psym, &flops);
    T2b = T2.get_amp2(iv1);
    FC_FUNC(g_if_sigma_oovv_oovv_no1_x0_type1_eri_v, G_IF_SIGMA_OOVV_OOVV_NO1_X0_TYPE1_ERI_V)
      (sc, ic, sv1, iv1, T2b.cptr(), W0aaaa.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  }
  }
  retval.acc_amp2(ic, S2b);
  } // End scope

  { 
  // No. 1, [2]
  // W1(i,o2,k,o3,v1,a) += (    1.00000000) D3(i,o2,k,o3,o1,o4) V2(v1,a,o1,o4) 
  // S2(i,k,a,c) += (    2.00000000) T2(o3,o2,c,v1) W1(i,o2,k,o3,v1,a) 
  double flops = 0; // Flop count
  int sv1 = s_eri;
  int iv1 = i_eri;
  orz::DTensor W1aaaav(orz::ct::sizeof_sympack_Xaaaav(symblockinfo, sv1));
  FC_FUNC(g_if_sigma_oovv_oovv_no0_x1_type1_eri_v, G_IF_SIGMA_OOVV_OOVV_NO0_X1_TYPE1_ERI_V)
    (sv1, iv1, V2_sym.cptr(), W1aaaav.cptr(), nir, nsym, psym, &flops);
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(iv1);
    FC_FUNC(g_if_sigma_oovv_oovv_no1_x1_type1_eri_v, G_IF_SIGMA_OOVV_OOVV_NO1_X1_TYPE1_ERI_V)
      (sc, ic, sv1, iv1, T2b.cptr(), W1aaaav.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
    retval.acc_amp2(ic, S2b);
  }
  }
  } // End scope

  { 
  // No. 2, [2]
  // W2(i,o3,k,o1,c,v1) += (    1.00000000) D3(i,o3,k,o2,o4,o1) V2(c,o2,o4,v1) 
  // S2(i,k,a,c) += (    2.00000000) T2(o3,o1,a,v1) W2(i,o3,k,o1,c,v1) 
  double flops = 0; // Flop count
  int sc = s_eri;
  int ic = i_eri;
  S2b = orz::DTensor(retval.namps_iamp()[ic]);
  for(int sv1 = 0;sv1 < nir;++sv1){ 
  for(int iv1 = symblockinfo.psym()(sv1,I_V,I_BEGIN);iv1 <= symblockinfo.psym()(sv1,I_V,I_END);++iv1){ 
    orz::DTensor W2aaaa(orz::ct::sizeof_sympack_Xaaaa(symblockinfo, sc^sv1));
    FC_FUNC(g_if_sigma_oovv_oovv_no0_x2_type1_eri_v, G_IF_SIGMA_OOVV_OOVV_NO0_X2_TYPE1_ERI_V)
      (sc, ic, sv1, iv1, V2_sym.cptr(), W2aaaa.cptr(), nir, nsym, psym, &flops);
    T2b = T2.get_amp2(iv1);
    FC_FUNC(g_if_sigma_oovv_oovv_no1_x2_type1_eri_v, G_IF_SIGMA_OOVV_OOVV_NO1_X2_TYPE1_ERI_V)
      (sc, ic, sv1, iv1, T2b.cptr(), W2aaaa.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  }
  }
  retval.acc_amp2(ic, S2b);
  } // End scope

  { 
  // No. 3, [2]
  // W3(i,k,o3,o4,v1,a) += (    1.00000000) D3(i,o2,k,o3,o1,o4) V2(v1,o1,o2,a) 
  // S2(i,k,a,c) += (    2.00000000) T2(o4,o3,v1,c) W3(i,k,o3,o4,v1,a) 
  double flops = 0; // Flop count
  int sv1 = s_eri;
  int iv1 = i_eri;
  orz::DTensor W3aaaav(orz::ct::sizeof_sympack_Xaaaav(symblockinfo, sv1));
  FC_FUNC(g_if_sigma_oovv_oovv_no0_x3_type1_eri_v, G_IF_SIGMA_OOVV_OOVV_NO0_X3_TYPE1_ERI_V)
    (sv1, iv1, V2_sym.cptr(), W3aaaav.cptr(), nir, nsym, psym, &flops);
  for(int sc = 0;sc < nir;++sc){ 
  for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){ 
    S2b = orz::DTensor(retval.namps_iamp()[ic]);
    T2b = T2.get_amp2(ic);
    FC_FUNC(g_if_sigma_oovv_oovv_no1_x3_type1_eri_v, G_IF_SIGMA_OOVV_OOVV_NO1_X3_TYPE1_ERI_V)
      (sc, ic, sv1, iv1, T2b.cptr(), W3aaaav.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
    retval.acc_amp2(ic, S2b);
  }
  }
  } // End scope

  { 
  // No. 4, [2]
  // W4(o1,o2,c,a) += (    1.00000000) T2(o1,o2,v1,v2) V2(c,v2,a,v1) 
  // S2(i,k,a,c) += (    2.00000000) D2(i,o1,k,o2) W4(o1,o2,c,a) 
  double flops = 0; // Flop count
  int sc = s_eri;
  int ic = i_eri;
  S2b = orz::DTensor(retval.namps_iamp()[ic]);
  orz::DTensor W4aav(orz::ct::sizeof_sympack_Xaav(symblockinfo, sc));
  for(int sv2 = 0;sv2 < nir;++sv2){ 
  for(int iv2 = symblockinfo.psym()(sv2,I_V,I_BEGIN);iv2 <= symblockinfo.psym()(sv2,I_V,I_END);++iv2){ 
    T2b = T2.get_amp2(iv2);
    FC_FUNC(g_if_sigma_oovv_oovv_no0_x4_type1_eri_v, G_IF_SIGMA_OOVV_OOVV_NO0_X4_TYPE1_ERI_V)
      (sc, ic, sv2, iv2, T2b.cptr(), V2_sym.cptr(), W4aav.cptr(), nir, nsym, psym, &flops);
  }
  }
  FC_FUNC(g_if_sigma_oovv_oovv_no1_x4_type1_eri_v, G_IF_SIGMA_OOVV_OOVV_NO1_X4_TYPE1_ERI_V)
    (sc, ic, W4aav.cptr(), S2b.cptr(), nir, nsym, psym, &flops);
  retval.acc_amp2(ic, S2b);
  } // End scope

  } // End my_rank
  }
  }

  } // End femto
  //*-- FEMTO ends --//*


  if(myrank == 0 || retval.alloc_type() != 2){                                                                                                          
  for(int ssig = 0;ssig < nir;++ssig){                                                                                      
  for(int isig = symblockinfo.psym()(ssig,I_V,I_BEGIN);isig <= symblockinfo.psym()(ssig,I_V,I_END);++isig){                 
    S2b = retval.get_amp2(isig);                                                                   
    FC_FUNC(g_if_sigma_oovv_scale,G_IF_SIGMA_OOVV_SCALE) // to make hamiltonian symmetric (scale by some factor of overlap) 
      (ssig, isig, S2b.cptr(), nir, nsym, psym);                                                   
    retval.put_amp2(isig, S2b); // S2ija, [b] <<-- Sb                                              
  } // End isig                                                                                                             
  } // End ssig                                                                                                             
  } // End if                                                                                                               
  orz::world().barrier();

  return retval; 
} 
