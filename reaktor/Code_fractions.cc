//
//  Code_fractions.cc
//  
//
//  Created by Masaaki Saitow on 12/08/02.
//  Copyright (c) 2012 Masaaki Saitow. All rights reserved.
//

#include <fstream>
#include <SQreaktor.hpp>

using namespace std;
using namespace Femto;
using namespace Femto::Core;

namespace Femto { namespace Reaktor {

  // *********************************************************
  // 
  // *********************************************************
    void SQreaktor::CLoop(const string s, const SQindex &i, ofstream &f)
  {
    string itype;
    if     (i.get_char() == (char_state)0) itype = "I_C"; 
    else if(i.get_char() == (char_state)1) itype = "I_O"; 
    else if(i.get_char() == (char_state)2) itype = "I_V";

    f << "  " + s + "for(int s" + i.get_index() + " = 0;s" + i.get_index() + " < nir;++s" + i.get_index() + "){ " << endl;
    f << "  " + s + "for(int i" + i.get_index() + " = symblockinfo.psym()(s" + i.get_index() + "," + itype + ",I_BEGIN);i" + i.get_index() + " <= symblockinfo.psym()(s" + i.get_index() + "," + itype + ",I_END);++i" + i.get_index() + "){ " << endl;
  } 


  // *********************************************************
  // 
  // *********************************************************
  void SQreaktor::MPIiproc(const string s, const SQindex &i, ofstream &f)
  {
    f << "  " + s + "if(hintmo.iproc_havingimo()[i" + i.get_index() + "] == myrank) {           " << endl;           
  }

  // *********************************************************
  // 
  // *********************************************************
  void SQreaktor::ReadAmp(const string s, const SQtensor &t, ofstream &f)
  {
    if(t.get_name() == LTensor_.get_name() && isBareLHS_) ReadRetval(s, t, f);
    else f << "  " + s + t.get_name() + "b = " + t.get_name() + ".get_amp2(i" + t.get_indices()[extamp_]->get_index() + ");" << endl;
  }


  // *********************************************************
  // 
  // *********************************************************
  void SQreaktor::ReadRetval(const string s, const SQtensor &t, ofstream &f)
  {
    f << "  " + s + t.get_name() + "b = orz::DTensor(retval.namps_iamp()[i" + t.get_indices()[extamp_]->get_index() + "]);" << endl;
  }


  // *********************************************************
  // 
  // *********************************************************
  void SQreaktor::ReadD4C(const string s, const SQtensor &t, ofstream &f)
  {
    f << "  " + s + "orz::DTensor " + t.get_name() + ";" << endl;        
    f << "  " + s + "orz::LoadBin(ctinp.dir()/(format(\"D4C_g[%d]\")%i" + t.get_indices()[extd4c_]->get_index() + ").str()) >> " + t.get_name() +";" << endl;
  }


#ifdef _VERBOSE_MODE
  // *********************************************************
  // ReadERI for verbose mode 
  // *********************************************************
  void SQreaktor::ReadERI(const string s, const SQtensor &t, ofstream &f)
  {
    
    f << "  " + s + "if(hintmo.iproc_havingimo()[i" + t.get_indices()[exth2_]->get_index() + "] == myrank) {           " << endl;
    f << "  " + s + "// Load ERIs from somewhere, e.g. disk, GA, etc..                                                 " << endl;  
    f << "  " + s + t.get_name() + " <<= 0.0;                                                                          " << endl;      
    f <<        s + "#ifdef _VERBOSE                                                                                   " << endl;
#ifndef _CHRONO
    f << "  " + s + "double start_read(ticktack());                                                                    " << endl;
#else
    f << "  " + s + "boost::chrono::high_resolution_clock::time_point t_start(boost::chrono::high_resolution_clock::now());" << endl;
#endif
    f <<        s + "#endif                                                                                            " << endl;
    f << "  " + s + "shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(i" + t.get_indices()[exth2_]->get_index() + ");" << endl;         
    f << "  " + s + "for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           " << endl;            
    f << "  " + s + "  // Load a signle record of integals                                                             " << endl;               
    f << "  " + s + "  const int &imo2 = loadbuf_ptr->i0;                                                              " << endl;                 
    f << "  " + s + "  const int &imo3 = loadbuf_ptr->i1;                                                              " << endl;            
    f << "  " + s + "  const int &imo4 = loadbuf_ptr->i2;                                                              " << endl;                 
    f << "  " + s + "  const double &v = loadbuf_ptr->v;                                                               " << endl;                   
    f << "  " + s + "  " + t.get_name() + "_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       " << endl;                      
    f << "  " + s + "}                                                                                                 " << endl;
    f <<        s + "#ifdef _VERBOSE                                                                                   " << endl;
#ifndef _CHRONO
    f << "  " + s + "t_readERI += (ticktack()-start_read);                                                             " << endl;
#else
    f << "  " + s + "boost::chrono::high_resolution_clock::time_point t_stop(boost::chrono::high_resolution_clock::now());" << endl;
    f << "  " + s + "boost::chrono::duration<double> elapsed = t_stop - t_start;                                       " << endl;
    f << "  " + s + "t_readERI += elapsed.count();                                                                     " << endl;
#endif
    f <<        s + "#endif                                                                                            " << endl;                        
    f << "  " + s + "const orz::DTensor " + t.get_name() + "_sym = orz::mr::sympack_int2(symblockinfo, i" + t.get_indices()[exth2_]->get_index() + ", s" + t.get_indices()[exth2_]->get_index() + ", " + t.get_name() + "); // V2=(IR-COV index) " << endl;
  }


#else
  // *********************************************************
  // 
  // *********************************************************
  void SQreaktor::ReadERI(const string s, const SQtensor &t, ofstream &f)
  {
    
    f << "  " + s + "if(hintmo.iproc_havingimo()[i" + t.get_indices()[exth2_]->get_index() + "] == myrank) {           " << endl;
    f << "  " + s + "// Load ERIs from somewhere, e.g. disk, GA, etc..                                                 " << endl;  
    f << "  " + s + t.get_name() + " <<= 0.0;                                                                          " << endl;      
    f << "  " + s + "shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(i" + t.get_indices()[exth2_]->get_index() + ");" << endl;         
    f << "  " + s + "for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           " << endl;            
    f << "  " + s + "  // Load a signle record of integals                                                             " << endl;               
    f << "  " + s + "  const int &imo2 = loadbuf_ptr->i0;                                                              " << endl;                 
    f << "  " + s + "  const int &imo3 = loadbuf_ptr->i1;                                                              " << endl;            
    f << "  " + s + "  const int &imo4 = loadbuf_ptr->i2;                                                              " << endl;                 
    f << "  " + s + "  const double &v = loadbuf_ptr->v;                                                               " << endl;                   
    f << "  " + s + "  " + t.get_name() + "_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       " << endl;                      
    f << "  " + s + "}                                                                                                 " << endl;                        
    f << "  " + s + "const orz::DTensor " + t.get_name() + "_sym = orz::mr::sympack_int2(symblockinfo, i" + t.get_indices()[exth2_]->get_index() + ", s" + t.get_indices()[exth2_]->get_index() + ", " + t.get_name() + "); // V2=(IR-COV index) " << endl;
  }
#endif

//*OLD_MPI*   // *********************************************************
//*OLD_MPI*   // 
//*OLD_MPI*   // *********************************************************
//*OLD_MPI*   void SQreaktor::ReadERI(const string s, const SQtensor &t, ofstream &f)
//*OLD_MPI*   {
//*OLD_MPI*      f << "  " + s + "if(my_imo <= i" + t.get_indices()[exth2_]->get_index() + " && i" + t.get_indices()[exth2_]->get_index() +" < (my_imo + my_nmo)) {                                                  " << endl;          
//*OLD_MPI*      f << "  " + s + "// Load ERIs from somewhere, e.g. disk, GA, etc..                                                 " << endl;  
//*OLD_MPI*      f << "  " + s + t.get_name() + " <<= 0.0;                                                                          " << endl;      
//*OLD_MPI*      f << "  " + s + "shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(i" + t.get_indices()[exth2_]->get_index() + ");" << endl;         
//*OLD_MPI*      f << "  " + s + "for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           " << endl;            
//*OLD_MPI*      f << "  " + s + "  // Load a signle record of integals                                                             " << endl;               
//*OLD_MPI*      f << "  " + s + "  const int &imo2 = loadbuf_ptr->i0;                                                              " << endl;                 
//*OLD_MPI*      f << "  " + s + "  const int &imo3 = loadbuf_ptr->i1;                                                              " << endl;            
//*OLD_MPI*      f << "  " + s + "  const int &imo4 = loadbuf_ptr->i2;                                                              " << endl;                 
//*OLD_MPI*      f << "  " + s + "  const double &v = loadbuf_ptr->v;                                                               " << endl;                   
//*OLD_MPI*      f << "  " + s + "  " + t.get_name() + "_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       " << endl;                      
//*OLD_MPI*      f << "  " + s + "}                                                                                                 " << endl;                        
//*OLD_MPI*      f << "  " + s + "const orz::DTensor " + t.get_name() + "_sym = orz::mr::sympack_int2(symblockinfo, i" + t.get_indices()[exth2_]->get_index() + ", s" + t.get_indices()[exth2_]->get_index() + ", " + t.get_name() + "); // V2=(IR-COV index) " << endl;
//*OLD_MPI*   }


//*NO_MPI*   // *********************************************************
//*NO_MPI*   // 
//*NO_MPI*   // *********************************************************
//*NO_MPI*   void SQreaktor::ReadERI(const string s, const SQtensor &t, ofstream &f)
//*NO_MPI*   {
//*NO_MPI*      f << "  " + s + "// Load ERIs from somewhere, e.g. disk, GA, etc..                                                 " << endl;  
//*NO_MPI*      f << "  " + s + t.get_name() + " <<= 0.0;                                                                          " << endl;      
//*NO_MPI*      f << "  " + s + "shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(i" + t.get_indices()[exth2_]->get_index() + ");" << endl;         
//*NO_MPI*      f << "  " + s + "for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {                                           " << endl;            
//*NO_MPI*      f << "  " + s + "  // Load a signle record of integals                                                             " << endl;               
//*NO_MPI*      f << "  " + s + "  const int &imo2 = loadbuf_ptr->i0;                                                              " << endl;                 
//*NO_MPI*      f << "  " + s + "  const int &imo3 = loadbuf_ptr->i1;                                                              " << endl;            
//*NO_MPI*      f << "  " + s + "  const int &imo4 = loadbuf_ptr->i2;                                                              " << endl;                 
//*NO_MPI*      f << "  " + s + "  const double &v = loadbuf_ptr->v;                                                               " << endl;                   
//*NO_MPI*      f << "  " + s + "  " + t.get_name() + "_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;                                       " << endl;                      
//*NO_MPI*      f << "  " + s + "}                                                                                                 " << endl;                        
//*NO_MPI*      f << "  " + s + "const orz::DTensor " + t.get_name() + "_sym = orz::mr::sympack_int2(symblockinfo, i" + t.get_indices()[exth2_]->get_index() + ", s" + t.get_indices()[exth2_]->get_index() + ", " + t.get_name() + "); // V2=(IR-COV index) " << endl;
//*NO_MPI*   }


//*SIMPLE_PTR*   // *********************************************************
//*SIMPLE_PTR*   // 
//*SIMPLE_PTR*   // *********************************************************
//*SIMPLE_PTR*   void SQreaktor::ReadD4(const string s, const SQtensor &t, ofstream &f)
//*SIMPLE_PTR*   {
//*SIMPLE_PTR*     SQindex* i1(t.get_indices()[0]);
//*SIMPLE_PTR*     SQindex* i2(t.get_indices()[1]);
//*SIMPLE_PTR* 
//*SIMPLE_PTR*     f << "  " + s + "// Load D4 from disk, or GA ....                                                     " << endl;
//*SIMPLE_PTR*     f << "  " + s + "int imoi = amo2imo[i" + i1->get_index() + "] - nclosed;                              " << endl;
//*SIMPLE_PTR*     f << "  " + s + "int imoj = amo2imo[i" + i2->get_index() + "] - nclosed;                              " << endl;
//*SIMPLE_PTR*     f << "  " + s + "                                                                                     " << endl;
//*SIMPLE_PTR*     f << "  " + s + "orz::DTensor rdm4_ij_sliced = rdm4(orz::Slice(imoi,imoi+1), orz::Slice(imoj,imoj+1), " << endl;
//*SIMPLE_PTR*     f << "  " + s + "                                   orz::Slice(),            orz::Slice(),            " << endl;
//*SIMPLE_PTR*     f << "  " + s + "                                   orz::Slice(),            orz::Slice(),            " << endl;
//*SIMPLE_PTR*     f << "  " + s + "                                   orz::Slice(),            orz::Slice()).copy();    " << endl;
//*SIMPLE_PTR*     f << "  " + s + "rdm4_sym = orz::mr::sympack_rdm4_2(symblockinfo, i" + i1->get_index() + ", s" + i1->get_index() + ", i" + i2->get_index() + ", s" + i2->get_index() + ", rdm4_ij_sliced);    " << endl;
//*SIMPLE_PTR*     f << "  " + s + "FC_FUNC(g_if_set_d4,G_IF_SET_D4)(s" + i1->get_index() + ", s" + i2->get_index() + ", i" + i1->get_index() + ", i" + i2->get_index() + ", rdm4_sym.cptr(), nir, nsym, psym);  " << endl;
//*SIMPLE_PTR*   }


  // *********************************************************
  // 
  // *********************************************************
  void SQreaktor::Cparallel(const string s, const SQindex &i, ofstream &f)
  { f << "  " + s + "if(hintmo.iproc_havingimo()[i" + i.get_index() + "] == myrank) {           " << endl; }


  // *********************************************************
  // 
  // *********************************************************
  void SQreaktor::ReadD4(const string s, const SQtensor &t, ofstream &f)
  {
    SQindex* i1(t.get_indices()[0]);
    SQindex* i2(t.get_indices()[1]);

    f << "  " + s + "// Load D4 from disk, or GA ....                                                     " << endl;
    f << "  " + s + "int imoi = amo2imo[i" + i1->get_index() + "] - nclosed;                              " << endl;
    f << "  " + s + "int imoj = amo2imo[i" + i2->get_index() + "] - nclosed;                              " << endl;
    f << "  " + s + "                                                                                     " << endl;
    f << "  " + s + "const double* rdm4_ij_sliced = rdm4.cptr() + (imoi*nocc+imoj)*nocc*nocc*nocc*nocc*nocc*nocc;" << endl; 
    f << "  " + s + "rdm4_sym = orz::mr::sympack_rdm4_2x(symblockinfo, i" + i1->get_index() + ", s" + i1->get_index() + ", i" + i2->get_index() + ", s" + i2->get_index() + ", rdm4_ij_sliced);    " << endl;
    f << "  " + s + "FC_FUNC(g_if_set_d4,G_IF_SET_D4)(s" + i1->get_index() + ", s" + i2->get_index() + ", i" + i1->get_index() + ", i" + i2->get_index() + ", rdm4_sym.cptr(), nir, nsym, psym);  " << endl;
  }


//*old*//   // *********************************************************
//*old*//   // 
//*old*//   // *********************************************************
//*old*//   void SQreaktor::ReadD4_Cumulant(const string s, const SQtensor &t, ofstream &f)
//*old*//   {
//*old*//     SQindex* i1(t.get_indices()[0]);
//*old*//     SQindex* i2(t.get_indices()[1]);
//*old*// 
//*old*//     f << endl;
//*old*//     f << "  " + s + "int imoi = amo2imo[i" + i1->get_index() + "] - nclosed;                              " << endl;
//*old*//     f << "  " + s + "int imoj = amo2imo[i" + i2->get_index() + "] - nclosed;                              " << endl;
//*old*//     f << "  " + s + "                                                                                     " << endl;
//*old*//     f << "  " + s + "// Generate D4 by cumulant expansion ....                                            " << endl;
//*old*//     f << "  " + s + "if(ctinp.use_d4cum_of()){"                                                            << endl;
//*old*//     f << "  " + s + "FC_FUNC(f_mrci_rdm4_cumulant_partial_opt,F_MRCI_RDM4_CUMULANT_PARTIAL_OPT)                  " << endl;        
//*old*//     f << "  " + s + "  (nocc, 0, rdmPack.rdm1().cptr(), rdmPack.rdm2().cptr(), rdmPack.cum2().cptr(), rdmPack.rdm3().cptr(),     " << endl;     
//*old*//     f << "  " + s + "   rdm4_ij_sliced.cptr(), imoi, imoj);                                               " << endl;        
//*old*//     f << "  " + s + "rdm4_sym = orz::mr::sympack_rdm4_2(symblockinfo, i" + i1->get_index() + ", s" + i1->get_index() + ", i" + i2->get_index() + ", s" + i2->get_index() + ", rdm4_ij_sliced);    " << endl;
//*old*//     f << "  " + s + "flops += nocc*nocc*nocc*nocc*nocc*nocc*70;" << endl;
//*old*//     f << "  " + s + "} // End if" << endl;
//*old*//     f << "  " + s + "// Slice the already existing 8-index 4-RDM ....                                            " << endl;
//*old*//     f << "  " + s + "else{" << endl;
//*old*//     f << "  " + s + "const double* rdm4_ij_sliced = rdm4.cptr() + (imoi*nocc+imoj)*nocc*nocc*nocc*nocc*nocc*nocc;" << endl; 
//*old*//     f << "  " + s + "rdm4_sym = orz::mr::sympack_rdm4_2x(symblockinfo, i" + i1->get_index() + ", s" + i1->get_index() + ", i" + i2->get_index() + ", s" + i2->get_index() + ", rdm4_ij_sliced);    " << endl;
//*old*//     f << "  " + s + "}" << endl;
//*old*//     f << "  " + s + "FC_FUNC(g_if_set_d4,G_IF_SET_D4)(s" + i1->get_index() + ", s" + i2->get_index() + ", i" + i1->get_index() + ", i" + i2->get_index() + ", rdm4_sym.cptr(), nir, nsym, psym);  " << endl;
//*old*//   }

  // *********************************************************
  // 
  // *********************************************************
  void SQreaktor::ReadD4_Cumulant(const string s, const SQtensor &t, ofstream &f)
  {
    SQindex* i1(t.get_indices()[0]);
    SQindex* i2(t.get_indices()[1]);

    f << endl;
    f << "  " + s + "int imoi = amo2imo[i" + i1->get_index() + "] - nclosed;                              " << endl;
    f << "  " + s + "int imoj = amo2imo[i" + i2->get_index() + "] - nclosed;                              " << endl;
    f << "  " + s + "                                                                                     " << endl;
    f << "  " + s + "// Generate D4 by cumulant expansion ....                                            " << endl;
    f << "  " + s + "if(ctinp.use_d4cum_of()){"                                                            << endl;
    f << "  " + s + "FC_FUNC(f_mr_rdm4_cumulant_partial_opt,F_MR_RDM4_CUMULANT_PARTIAL_OPT)                  " << endl;        
    f << "  " + s + "  (nocc, 0, rdmPack.rdm1().cptr(), rdmPack.rdm2().cptr(), rdmPack.cum2().cptr(), rdmPack.rdm3().cptr(),     " << endl;     
    f << "  " + s + "   rdm4_ij_sliced.cptr(), imoi, imoj);                                               " << endl;        
    f << "  " + s + "rdm4_sym = orz::mr::sympack_rdm4_2(symblockinfo, i" + i1->get_index() + ", s" + i1->get_index() + ", i" + i2->get_index() + ", s" + i2->get_index() + ", rdm4_ij_sliced);    " << endl;
    f << "  " + s + "flops += nocc*nocc*nocc*nocc*nocc*nocc*70;" << endl;
    f << "  " + s + "} // End if" << endl;
    f << "  " + s + "// Slice the already existing 8-index 4-RDM ....                                            " << endl;
    f << "  " + s + "else{" << endl;
    f << "  " + s + "const double* rdm4_ij_sliced = rdm4.cptr() + (imoi*nocc+imoj)*nocc*nocc*nocc*nocc*nocc*nocc;" << endl; 
    f << "  " + s + "rdm4_sym = orz::mr::sympack_rdm4_2x(symblockinfo, i" + i1->get_index() + ", s" + i1->get_index() + ", i" + i2->get_index() + ", s" + i2->get_index() + ", rdm4_ij_sliced);    " << endl;
    f << "  " + s + "}" << endl;
    f << "  " + s + "FC_FUNC(g_if_set_d4,G_IF_SET_D4)(s" + i1->get_index() + ", s" + i2->get_index() + ", i" + i1->get_index() + ", i" + i2->get_index() + ", rdm4_sym.cptr(), nir, nsym, psym);  " << endl;
  }

//*OLD_CUM*   // *********************************************************
//*OLD_CUM*   // 
//*OLD_CUM*   // *********************************************************
//*OLD_CUM*   void SQreaktor::ReadD4_Cumulant(const string s, const SQtensor &t, ofstream &f)
//*OLD_CUM*   {
//*OLD_CUM*     SQindex* i1(t.get_indices()[0]);
//*OLD_CUM*     SQindex* i2(t.get_indices()[1]);
//*OLD_CUM* 
//*OLD_CUM*     f << "  " + s + "// Generate D4 by cumulant expansion ....                                            " << endl;
//*OLD_CUM*     f << "  " + s + "int imoi = amo2imo[i" + i1->get_index() + "] - nclosed;                              " << endl;
//*OLD_CUM*     f << "  " + s + "int imoj = amo2imo[i" + i2->get_index() + "] - nclosed;                              " << endl;
//*OLD_CUM*     f << "  " + s + "                                                                                     " << endl;
//*OLD_CUM*     f << "  " + s + "FC_FUNC(f_mrci_rdm4_cumulant_partial, F_MRCI_RDM4_CUMULANT_PARTIAL)                  " << endl;        
//*OLD_CUM*     f << "  " + s + "  (nocc, 0, rdmPack.rdm1().cptr(), rdmPack.rdm2().cptr(), rdmPack.rdm3().cptr(),     " << endl;     
//*OLD_CUM*     f << "  " + s + "   rdm4_ij_sliced.cptr(), imoi, imoj);                                               " << endl;        
//*OLD_CUM*     f << "  " + s + "rdm4_sym = orz::mr::sympack_rdm4_2(symblockinfo, i" + i1->get_index() + ", s" + i1->get_index() + ", i" + i2->get_index() + ", s" + i2->get_index() + ", rdm4_ij_sliced);    " << endl;
//*OLD_CUM*     f << "  " + s + "FC_FUNC(g_if_set_d4,G_IF_SET_D4)(s" + i1->get_index() + ", s" + i2->get_index() + ", i" + i1->get_index() + ", i" + i2->get_index() + ", rdm4_sym.cptr(), nir, nsym, psym);  " << endl;
//*OLD_CUM*   }


  // *********************************************************
  // 
  // *********************************************************
  void SQreaktor::AccAmp(const string s, const SQtensor &t, ofstream &f)
  {
    f << "  " + s + "retval.acc_amp2(i" + t.get_indices()[extamp_]->get_index() + ", " + t.get_name() + "b);" << endl;
  }


  // *********************************************************
  // 
  // *********************************************************
  void SQreaktor::LoopEnd(string s, ofstream &f)
  {
    f << "  " + s + "}" << endl;
    f << "  " + s + "}" << endl;
  }


  // *********************************************************
  // 
  // *********************************************************
    void SQreaktor::LoopEnd_i(string s, const SQindex &i, ofstream &f)
  {
    f << "  " + s + "} // End i" << i.get_index() << endl;
    f << "  " + s + "} // End s" << i.get_index() << endl;
  }


}} // Femto::
