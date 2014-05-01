//
//  Factorize_cas.cc
//  
//
//  Created by Masaaki Saitow on 12/11/06.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <algorithm>
#include <fstream>
#include <cmath>
#include <string>
#include <boost/format.hpp>
#include <SQreaktor.hpp>

//#define _O1_DEBUG

namespace femto {

  // *********************************************************
  // 
  // *********************************************************
  void SQreaktor::factorize_cas()
  {
    if(LTensor_.get_indices().size() != 4 && isBareLHS_){
      cout << "Argument 0 cannnot be treated as a bareampack .... " << endl;
      abort();
    } // End if

    vector<SQindex*> Linds = LTensor_.get_indices();
    for(vector<SQindex*>::iterator i = Linds.begin();i != Linds.end();++i)
      if((*i)->get_isSummed()){
        cout << *i << " in 0th argument is a kind of dummy index." << endl;
        abort(); 
      } // End if

    for(vector<SQterm>::iterator t = inTerms_.begin();t != inTerms_.end();++t){
      if(!t->get_isInCanonical()){
        cout << "Term, " << *t << " is not in canonical form .... " << endl;
        abort();
      } // End if
    } // End t

    if(is_sfGen(LTensor_.get_name())){
      cout << "Factorize: 1st argument cannot be a spin-free unitary group generator ... " << endl;
      abort();
    }
    for(vector<SQterm>::iterator t = inTerms_.begin();t != inTerms_.end();++t){
      for(size_t num_t = 0;num_t < t->get_tensors().size();++num_t){
        if(is_sfGen(t->get_tensors()[num_t].get_name())){
          cout << "Factorize: A spin-free unitary group generator is detected in the 2nd argument .... " << endl;
          abort();
	} // End if
      } // End num_t
    } // End t

    if(!inTerms_.size()) return;

    // Header file for C++ 
    string CHname("c_" + title_ + ".h");
    string CHguard("c_" + title_ + "_h");
    // It's a kind of trick, though .... 
    transform(CHguard.begin(), CHguard.end(), CHguard.begin(), (int(*)(int))toupper); 
    ofstream CHfile( CHname.c_str() );
    CHfile << "\n" << endl;
    CHfile << "#ifndef " << CHguard << endl;
    CHfile << "#define " << CHguard << endl;
    CHfile << "" << endl;
    CHfile << "// #include <tensor/tensor.h>                  " << endl;  
    CHfile << "// #include <sci/hint/hintmo/hintmo.h>         " << endl;  
    CHfile << "// #include <sci/ctnew2/ctclass_input.h>       " << endl;   
    CHfile << "// #include <sci/ctnew2/ctclass_symblock.h>    " << endl;   
    CHfile << "// #include <sci/ctnew2/ctclass_rdmpack.h>     " << endl;   
    CHfile << "// #include <sci/ctnew2/ctclass_bareamppack.h> " << endl;       
    CHfile << "                                               " << endl;            
    CHfile << "extern \"C\"{                                  " << endl;      
    CHfile << "                                               " << endl;           
    CHfile << "                                               " << endl;            
    CHfile  << femto_logo("// ") << endl; 

    string CHend("");
    CHend += "      \n ";
    CHend += "}     \n ";       
    CHend += "      \n ";  
    CHend += "      \n ";      
    CHend += "#endif\n ";       
    CHend += "      \n ";        
    CHend += "      \n ";       

    // Output for tensorial contractions .... 
    string Fname = "f_" + title_ + ".F90";
    ofstream F90file( Fname.c_str() );
    F90file << "#include \"../f_ct.fh\"\n\n" << endl;    
    F90file << femto_logo("! ")              << endl;

    // Output for the body of tensorial contractions .... 
    bool Bareflag = false; // LHS has amplitude or not
    bool D4flag   = false; // LHS has 4-RDM or not
    for(vector<SQterm>::iterator t = inTerms_.begin();t != inTerms_.end();++t){
      for(size_t num_t = 0;num_t < t->get_tensors().size();++num_t){
        if(t->get_tensors()[num_t].get_name() == name_amp_) Bareflag = true;
        if(t->get_tensors()[num_t].get_name() == name_d4_)  D4flag   = true; 
      } // End t    
    } // End term

    // inTerms_ have constants or not
    bool Constsflag = false;
    vector<string> Consts;
    for(vector<SQterm>::iterator t = inTerms_.begin();t != inTerms_.end();++t){
      for(size_t num_c = 0;num_c < t->get_Consts().size();++num_c)
        if(find(Consts.begin(), Consts.end(), t->get_Consts()[num_c]) == Consts.end() && t->get_Consts()[num_c] != "")
          Consts.push_back(t->get_Consts()[num_c]);
    } // End t
    if(Consts.size()) Constsflag = true;

    string CPname = "c_" + title_ + ".cpp";
    ofstream CPfile( CPname.c_str() );
    CPfile << "                                                                                " << endl;      
    CPfile << "#include <orz/orz.h>                                                            " << endl;    
    CPfile << "#include <orz/openmp.h>                                                         " << endl;     
    CPfile << "#include <orz/cblas.h>                                                          " << endl;      
    CPfile << "#include <orz/clapack.h>                                                        " << endl;   
    CPfile << "#include <tensor/tensor.h>                                                      " << endl;     
    CPfile << "#include <sci/hint/para_disttools.h>                                            " << endl; 
    CPfile << "#include <sci/ctnew2/ct.h>                                                      " << endl; 
    CPfile << "#include <sci/ctnew2/ct_f.h>                                                    " << endl;  
    CPfile << "#include <sci/ctnew2/ctclass_input.h>                                           " << endl;     
    CPfile << "#include <sci/ctnew2/ctclass_symblock.h>                                        " << endl;     
    CPfile << "#include <sci/ctnew2/ctclass_hintmo.h>                                          " << endl;  
    CPfile << "#include <sci/ctnew2/ctclass_rdmpack.h>                                         " << endl;     
    CPfile << "#include <sci/ctnew2/ctclass_bareamppack.h>                                     " << endl;  
    CPfile << "#include <sci/ctnew2/ctclass_orthamppack.h>                                     " << endl;      
    CPfile << "#include <sci/ctnew2/diaghessian.h>                                             " << endl;      
    CPfile << "#include <sci/ctnew2/symamp2.h>                                                 " << endl;    
    CPfile << "#include <sci/ctnew2/mrci.h>                                                    " << endl;     
    CPfile << "#include <sci/ctnew2/" + CHname + ">                                            " << endl;
    CPfile << "                                                                                " << endl;        
    CPfile << "using std::cout;                                                                " << endl;      
    CPfile << "using std::endl;                                                                " << endl;  
    CPfile << "                                                                                " << endl;      
    CPfile << "#define FLOPCOUNT                                                               " << endl;  
    CPfile << "                                                                                " << endl;   
    if(do_timing_){
    CPfile << "//Timing object                                                                 " << endl;         
    CPfile << "extern std::vector<boost::tuple<std::string, double, double, double> > my_timer;" << endl;         
    CPfile << "extern double Fc0;                                                              " << endl << endl;         
    }
    CPfile << femto_logo("// ")                                                                  << endl;
    CPfile << "                                                                                " << endl;
    CPfile << "// ***************************************************************************  " << endl;           
    CPfile << "// orz::ct::mrci                                                                " << endl;      
    CPfile << "// ***************************************************************************  " << endl;          
    CPfile << "									    /*!        " << endl; 
    CPfile << "   @brief CT input                                                              " << endl;    
    CPfile << "									     */        " << endl;    
    CPfile << "                                                                                " << endl;     

    if(isBareLHS_)
    CPfile << "orz::ct::BareAmpPack orz::ct::" + title_ + "(const orz::ct::Input &ctinp,                                    " << endl;
    else if(!LTensor_.get_indices().size() && Bareflag)
    CPfile << "double orz::ct::" + title_ + "(const orz::ct::Input &ctinp,                                                  " << endl;
    else if(!isBareLHS_ && LTensor_.get_indices().size())
    CPfile << "orz::DTensor" + title_ +  "(const orz::ct::Input &ctinp,                                                  " << endl;
    else{
      cout << "Cannot detect the output type .... " << endl;
      abort();
    }
    CPfile << "                                  const orz::ct::SymBlockInfo &symblockinfo,                                " << endl;
    CPfile << "                                  const orz::ct::HintMO &hintmo,                                            " << endl;
    if(isBareLHS_)
    CPfile << "                                  const int alloc_type,                                            " << endl;
    if(D4flag){
    CPfile << "                                  const orz::ct::RdmPack &rdmPack,                                      " << endl;
    CPfile << "                                  const orz::DTensor &rdm4,                                                 " << endl;
    }
    if(!LTensor_.get_indices().size() && !isBareLHS_)
    CPfile << "                                  const double init_value,                                                  " << endl;
    if(Bareflag)
    CPfile << "                                  const orz::ct::BareAmpPack &" + name_amp_ + ",                             " << endl;
    CPfile << "                                  const int num_sigma";
    if(Consts.size()) CPfile << "," << endl;

    for(vector<string>::iterator c = Consts.begin();c != Consts.end();++c){
    CPfile << "                                  const double " + *c;
    if(*c != Consts.back()) CPfile << "," << endl;     
    } // End c
    CPfile << ") {" << endl;     

    CPfile << "                                                                                                                 " << endl;         
    CPfile << "                                                                                                                 " << endl;
    CPfile << "  // set up nmo nclosed, nocc                                                                                    " << endl;    
    CPfile << "  const FC_INT nclosed = ctinp.nclosed();                                                                        " << endl;         
    CPfile << "  const FC_INT nocc    = ctinp.nocc();                                                                           " << endl;                
    CPfile << "  const FC_INT nvir    = ctinp.nvir();                                                                           " << endl;                  
    CPfile << "  const FC_INT nmo     = nclosed + nocc + nvir;                                                                  " << endl;               
    CPfile << "  const FC_INT nir     = symblockinfo.nir();                                                                     " << endl;               
    CPfile << "  const FC_INT * const nsym    = symblockinfo.nsym().cptr();                                                     " << endl;                       
    CPfile << "  const FC_INT * const psym    = symblockinfo.psym().cptr();                                                     " << endl;                   
    CPfile << "  const FC_INT * const amo2imo = symblockinfo.amo2imo().cptr();                                                  " << endl;               
    CPfile << "                                                                                                                 " << endl;
    if(isBareLHS_){
    CPfile << "  std::ostringstream stm;                                                                                        " << endl;                       
    CPfile << "  stm << num_sigma;                                                                                              " << endl;                        
    CPfile << "  std::string name_of_sigma = \"" + LTensor_.get_name() +  "\" + stm.str() + \"]\"; // Name of the Sigma vector  " << endl;                              
    CPfile << "  orz::ct::BareAmpPack retval                                                                                    " << endl;                         
    CPfile << "    = orz::ct::BareAmpPack(ctinp, symblockinfo, name_of_sigma, alloc_type); // Sigma(a, a', e, e') tensor        " << endl;                 
    CPfile << "                                                                                                                 " << endl;                           
    CPfile << "  orz::DTensor " + LTensor_.get_name() + "b; // Container of S2_aae,[b] tensor                                   " << endl;                       
    CPfile << "                                                                                                                 " << endl;
    } // End if
    else if(!LTensor_.get_indices().size() && Bareflag)
    CPfile << "  double " + LTensor_.get_name() + " = init_value;                                                               " << endl;
    if(Bareflag)
    CPfile << "  orz::DTensor " + name_amp_ + "b; // Container of T2_aae,[b] tensor                                             " << endl;                       
    if(D4flag){
    CPfile << "  orz::DTensor rdm4_sym;                                                                                         " << endl;
    CPfile << "  orz::DTensor rdm4_ij_sliced(ctinp.use_d4cum_of() ? nocc*nocc*nocc*nocc*nocc*nocc : 0);                         " << endl;
    }

    // Print OpenMPI parameters
    CPfile << "  // set nproc, myrank                      " << endl;                           
    CPfile << "  const int nproc = orz::world().size();    " << endl;                    
    CPfile << "  const int myrank = orz::world().rank();   " << endl;                     
    CPfile << endl;
    // Setting up some integrals 
    CPfile << "  orz::DTensor moint1 = hintmo.int1(); // Setting up one-body integrals                                         " << endl;                        
    CPfile << "  const orz::DTensor moint1_sym = (myrank == 0) ? orz::ct::sympack_int1(symblockinfo, moint1) : orz::DTensor(); // moint1=(IR-COV index)" << endl;         
    CPfile << "  orz::DTensor " + name_h2_ + "(nmo,nmo,nmo);                                                                    " << endl;                             
    CPfile << "  double * const " + name_h2_ + "_ptr = " + name_h2_ + ".cptr();                                                  " << endl;                          
    CPfile << endl;

    // Print logos
    cout << endl;
    cout << femto_logo("");
    cout << endl;

    // Process the unlinked tensors to replace with core Fock matrix
    replace_Fock();

    // Rotate each index in inTerms_ and achieve the best matching of the external indices
    regulate_indices();

    // Decompose each group of terms into the groups of optimal binary comtractions
    if (type_terms_["noeri"].size()){
      binary_decomposition(type_terms_["noeri"]);
      generate_contract(CPfile, CHfile, F90file, binaries_, "noeri");
    } // End if
    if(type_terms_["eri_c"].size()){
      binary_decomposition(type_terms_["eri_c"]);
      generate_contract(CPfile, CHfile, F90file, binaries_, "eri_c");
    } // End if
    if(type_terms_["eri_o"].size()){
      binary_decomposition(type_terms_["eri_o"]);
      generate_contract(CPfile, CHfile, F90file, binaries_, "eri_o");
    } // End if
    if(type_terms_["eri_v"].size()){
      binary_decomposition(type_terms_["eri_v"]);
      generate_contract(CPfile, CHfile, F90file, binaries_, "eri_v");
    } // End if

    if(!LTensor_.get_indices().size() && Bareflag){
      CPfile << endl;
      CPfile << "  return  " << LTensor_.get_name()  << ";" << endl;
      CPfile << "} " << endl;
    } // End if
    else{
      if(isBareLHS_ && LTensor_.get_indices()[0]->get_char() == (char_state)1 && LTensor_.get_indices()[1]->get_char() == (char_state)1
	           && LTensor_.get_indices()[2]->get_char() == (char_state)2 && LTensor_.get_indices()[3]->get_char() == (char_state)2){

        CPfile << endl;
        CPfile << "  if(myrank == 0 || retval.alloc_type() != 2){                                                                                                          " << endl;
        CPfile << "  for(int ssig = 0;ssig < nir;++ssig){                                                                                      " << endl;        
        CPfile << "  for(int isig = symblockinfo.psym()(ssig,I_V,I_BEGIN);isig <= symblockinfo.psym()(ssig,I_V,I_END);++isig){                 " << endl;    
        CPfile << "    " + LTensor_.get_name() + "b = retval.get_amp2(isig);                                                                   " << endl;      
        CPfile << "    FC_FUNC(g_if_sigma_oovv_scale,G_IF_SIGMA_OOVV_SCALE) // to make hamiltonian symmetric (scale by some factor of overlap) " << endl;       
        CPfile << "      (ssig, isig, " + LTensor_.get_name() + "b.cptr(), nir, nsym, psym);                                                   " << endl;              
        CPfile << "    retval.put_amp2(isig, " + LTensor_.get_name() + "b); // S2ija, [b] <<-- Sb                                              " << endl;                   
        CPfile << "  } // End isig                                                                                                             " << endl;            
        CPfile << "  } // End ssig                                                                                                             " << endl;
        CPfile << "  } // End if                                                                                                               " << endl;     
	CPfile << "  orz::world().barrier();" << endl;
      } // End if
      if(isBareLHS_ && LTensor_.get_indices()[0]->get_char() == (char_state)0 && LTensor_.get_indices()[1]->get_char() == (char_state)0
	           && LTensor_.get_indices()[2]->get_char() == (char_state)2 && LTensor_.get_indices()[3]->get_char() == (char_state)2){

        CPfile << endl;
        CPfile << "  if(myrank == 0 || retval.alloc_type() != 2){                                                                                                          " << endl;
        CPfile << "  for(int ssig = 0;ssig < nir;++ssig){                                                                                      " << endl;        
        CPfile << "  for(int isig = symblockinfo.psym()(ssig,I_V,I_BEGIN);isig <= symblockinfo.psym()(ssig,I_V,I_END);++isig){                 " << endl;    
        CPfile << "    " + LTensor_.get_name() + "b = retval.get_amp2(isig);                                                                   " << endl;      
        CPfile << "    FC_FUNC(g_if_sigma_ccvv_scale,G_IF_SIGMA_CCVV_SCALE) // to make hamiltonian symmetric (scale by some factor of overlap) " << endl;       
        CPfile << "      (ssig, isig, " + LTensor_.get_name() + "b.cptr(), nir, nsym, psym);                                                   " << endl;              
        CPfile << "    retval.put_amp2(isig, " + LTensor_.get_name() + "b); // S2ija, [b] <<-- Sb                                              " << endl;                   
        CPfile << "  } // End isig                                                                                                             " << endl;            
        CPfile << "  } // End ssig                                                                                                             " << endl;     
        CPfile << "  } // End if                                                                                                               " << endl;
	CPfile << "  orz::world().barrier();" << endl;
      } // End if
      if(isBareLHS_ && LTensor_.get_indices()[0]->get_char() == (char_state)0 && LTensor_.get_indices()[1]->get_char() == (char_state)0
	           && LTensor_.get_indices()[2]->get_char() == (char_state)1 && LTensor_.get_indices()[3]->get_char() == (char_state)1){

        CPfile << endl;
        CPfile << "  if(myrank == 0 || retval.alloc_type() != 2){                                                                                                          " << endl;
        CPfile << "  for(int ssig = 0;ssig < nir;++ssig){                                                                                      " << endl;        
        CPfile << "  for(int isig = symblockinfo.psym()(ssig,I_O,I_BEGIN);isig <= symblockinfo.psym()(ssig,I_O,I_END);++isig){                 " << endl;    
        CPfile << "    " + LTensor_.get_name() + "b = retval.get_amp2(isig);                                                                   " << endl;      
        CPfile << "    FC_FUNC(g_if_sigma_ccoo_scale,G_IF_SIGMA_CCOO_SCALE) // to make hamiltonian symmetric (scale by some factor of overlap) " << endl;       
        CPfile << "      (ssig, isig, " + LTensor_.get_name() + "b.cptr(), nir, nsym, psym);                                                   " << endl;              
        CPfile << "    retval.put_amp2(isig, " + LTensor_.get_name() + "b); // S2ija, [b] <<-- Sb                                              " << endl;                   
        CPfile << "  } // End isig                                                                                                             " << endl;            
        CPfile << "  } // End ssig                                                                                                             " << endl;     
        CPfile << "  } // End if                                                                                                               " << endl;
	CPfile << "  orz::world().barrier();" << endl;
      } // End if

      CPfile << endl;
      CPfile << "  return retval; " << endl;
      CPfile << "} " << endl;
    } // End else

    CHfile << CHend;
    
  }                                               

} //End femto
