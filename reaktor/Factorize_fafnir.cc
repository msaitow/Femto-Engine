//
//  Factorize_new.cc
//  
//
//  Created by Masaaki Saitow on 12/10/24.
//  Copyright (c) 2012 Masaaki Saitow. All rights reserved.
//

#include <algorithm>
#include <fstream>
#include <cmath>
#include <string>
#include <boost/format.hpp>
#include <SQreaktor.hpp>

#define _DEBUG_OLD

using namespace std;
using namespace Femto;
using namespace Femto::Core;

namespace Femto { namespace Reaktor {

  // *********************************************************
  // 
  // *********************************************************
  string SQreaktor::factorize_fafnir()
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

    if(!binaries_.size()) return "";

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
    CHfile  << Femto_logo("// ") << endl; 

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
    F90file << "#include <sci/icmr/fsrc/f_mr.fh>\n\n" << endl;    
    F90file << Femto_logo("! ")              << endl;
    F90file << "!                                    Generated date : " + Femto_date() << endl;

    // Output for the body of tensorial contractions .... 
    bool Bareflag(false); // RHS has T2 amplitude or not
    bool DTenflag(false); // RHS has T1 amplitude or not
    bool D4flag  (false); // RHS has 4-RDM or not
    bool Fockflag(false); // RHS has CAS-Fock or not
    for(auto ts = binaries_.begin();ts != binaries_.end();++ts)
      for(auto t = ts->begin();t != ts->end();++t){
	auto temp(t->get_Rtensors());
	for(auto tt = temp.begin();tt != temp.end();++tt)
	  if     (tt->get_name() == name_amp_)   Bareflag = true;
	  else if(tt->get_name() == name_d4_)    D4flag   = true; 
	  else if(tt->get_name() == T1_name())   DTenflag = true; 
	  else if(tt->get_name() == Fock_name()) Fockflag = true; 
      } // End t

    // inTerms_ have constants or not
    bool Constsflag = false;
    vector<string> Consts;
    for(auto ts = binaries_.begin();ts != binaries_.end();++ts)
      for(auto t = ts->begin();t != ts->end();++t){
	auto temp(t->get_Consts());
	for(auto tt = temp.begin();tt != temp.end();++tt) if(find(Consts.begin(), Consts.end(), *tt) == Consts.end() && *tt != "Fc0") Consts.push_back(*tt); // The trace of the core Fock "Fc0" is assumed to be extern veriable
      } // End t
    if(Consts.size()) Constsflag = true;

    string CPname = "c_" + title_ + ".cpp";
    ofstream CPfile( CPname.c_str() );
    CPfile << "                                                                                " << endl;      
    CPfile << cont_header_;
    CPfile << "#include <sci/icmr/Femto/elems/" + CHname + ">                                  " << endl;
    CPfile << "                                                                                " << endl;        
    CPfile << "using std::cout;                                                                " << endl;      
    CPfile << "using std::endl;                                                                " << endl;  
    CPfile << "                                                                                " << endl;      
    CPfile << "#define FLOPCOUNT                                                               " << endl;  
    CPfile << "                                                                                " << endl;   
    CPfile << "using orz::mr::Femto::Fc0;                                                      " << endl;   
    CPfile << "using orz::mr::Femto::h1_int;                                                   " << endl;   
    CPfile << "using orz::mr::Femto::h6_int;                                                   " << endl;   
    CPfile << "                                                                                " << endl;   
    if(do_timing_){
//*old*//     CPfile << "//Timing object                                                                 " << endl;         
//*old*//     CPfile << "extern std::vector<boost::tuple<std::string, double, double> > my_timer;" << endl << endl;
//*old*//     CPfile << "// File stream object to write timing data                              " << endl;
//*old*//     CPfile << "extern std::ofstream file_timing;                                       " << endl << endl;
//*old*//     CPfile << "// Core integrals                                                               " << endl;                  
//*old*//     CPfile << "extern double Fc0;                                                              " << endl;         
//*old*//     CPfile << "extern double h1_int;                                                           " << endl << endl;         
//*old*//     CPfile << "// CAS-Fock matrix                                                              " << endl;                  
//*old*//     CPfile << "extern double h6_int;                                                           " << endl << endl;         
    }
    CPfile << Femto_logo("// ")                                                                  << endl;
    CPfile << "//                                   Generated date : " + Femto_date() << endl;
    CPfile << "                                                                                " << endl;
    CPfile << "// ***************************************************************************  " << endl;           
    CPfile << "// orz::mr::mrci                                                                " << endl;      
    CPfile << "// ***************************************************************************  " << endl;          
    CPfile << "									    /*!        " << endl; 
    CPfile << "   @brief CT input                                                              " << endl;    
    CPfile << "									     */        " << endl;    
    CPfile << "                                                                                " << endl;     

#ifdef _DEBUG_OLD
    /////////////////////// Toriaezu 2014/03/07 //////////////////
//no*good//    if(title_[0] == 's' && title_[1] == 'i' && title_[2] == 'g' && title_[3] == 'm' && title_[4] == 'a'){
//no*good//
//no*good//      int oLcount(0);
//no*good//      //for(auto i = title_.cbegin();i != title_.cend();++i) if(*i == 'o') ++ocount;
//no*good//      if(title_[ 6] == 'o') ++oLcount;
//no*good//      if(title_[ 7] == 'o') ++oLcount;
//no*good//      if(title_[ 8] == 'o') ++oLcount;
//no*good//      if(title_[ 9] == 'o') ++oLcount;
//no*good//
//no*good//      int oRcount(0);
//no*good//      if(title_[11] == 'o') ++oRcount;
//no*good//      if(title_[12] == 'o') ++oRcount;
//no*good//      if(title_[13] == 'o') ++oRcount;
//no*good//      if(title_[14] == 'o') ++oRcount;
//no*good//
//no*good//      if(oLcount >= 2 && oRcount >= 2) D4flag = true;
//no*good//    } // End if
    if(LTensor_.get_indices().size() == 4){
      bool hadD4(false);
      for(auto ts = binaries_.begin();ts != binaries_.end();++ts)
	for(auto t = ts->begin();t != ts->end();++t){
	  auto myTensors(t->get_Rtensors());
	  for(auto m = myTensors.begin();m != myTensors.end();++m)
	    if(m->get_name() == "D4" || m->get_name() == "C2" || m->get_name() == "C5" || m->get_name() == "C4" || m->get_name() == "C6") hadD4 = true;
	} // End t
      D4flag = hadD4;
    } // End if
    //////////////////////////////////////////////////////////////
#endif

    if(isBareLHS_)
    CPfile << "orz::mr::BareAmpPack orz::mr::Femto::" + title_ + "(const orz::mr::Input &ctinp,                                    " << endl;
    else if(!LTensor_.get_indices().size() && Bareflag)
    CPfile << "double orz::mr::Femto::" + title_ + "(const orz::mr::Input &ctinp,                                                  " << endl;
    else if(LTensor_.get_name() == D4C_nameL())
    CPfile << "void orz::mr::Femto::" + title_ +  "(const orz::mr::Input &ctinp,                                                     " << endl;
    else if(!isBareLHS_ && LTensor_.get_indices().size())
    CPfile << "orz::DTensor orz::mr::Femto::" + title_ +  "(const orz::mr::Input &ctinp,                                                     " << endl;
    else{
      cout << "Cannot detect the output type .... " << endl;
      abort();
    }
    CPfile << "                                  const orz::mr::SymBlockInfo &symblockinfo,                                 " << endl;
    CPfile << "                                  const orz::mr::HintMO &hintmo,                                             " << endl;
    if(Fockflag)
    CPfile << "                                  const orz::DTensor &CFock,                                                 " << endl;
    if(isBareLHS_)
    CPfile << "                                  const int alloc_type,                                                      " << endl;
    if(D4flag){
    CPfile << "                                  const orz::mr::RdmPack &rdmPack,                                           " << endl;
    CPfile << "                                  const orz::DTensor &rdm4,                                                  " << endl;
    }
    if(!LTensor_.get_indices().size() && !isBareLHS_)
    CPfile << "                                  const double init_value,                                                   " << endl;
    if(Bareflag)
    CPfile << "                                  const orz::mr::BareAmpPack &" + name_amp_ + ",                             " << endl;
    if(DTenflag)
    CPfile << "                                  const orz::DTensor &" + T1_name() + ",                                     " << endl;
    CPfile << "                                  const int num_sigma";
    if(Consts.size()) CPfile << "," << endl;

    for(vector<string>::iterator c = Consts.begin();c != Consts.end();++c){
    CPfile << "                                  const double " + *c;
    if(*c != Consts.back()) CPfile << "," << endl;     
    } // End c
    CPfile << ") {" << endl;     

    //////////////////////////////////////////////////////////////////////////////////////////////
    ostringstream retVal;
    { // Construction of the header file
      if(isBareLHS_)
      retVal << "orz::mr::BareAmpPack " + title_ + "(const orz::mr::Input &ctinp,                                    " << endl;
      else if(!LTensor_.get_indices().size() && Bareflag)
      retVal << "double " + title_ + "(const orz::mr::Input &ctinp,                                                  " << endl;
      else if(LTensor_.get_name() == D4C_nameL())
      retVal << "void " + title_ +  "(const orz::mr::Input &ctinp,                                                     " << endl;
      else if(!isBareLHS_ && LTensor_.get_indices().size())
      retVal << "orz::DTensor " + title_ +  "(const orz::mr::Input &ctinp,                                                     " << endl;
      else{
	cout << "Cannot detect the output type .... " << endl;
	abort();
      }
      retVal << "                                  const orz::mr::SymBlockInfo &symblockinfo,                                 " << endl;
      retVal << "                                  const orz::mr::HintMO &hintmo,                                             " << endl;
      if(Fockflag)
      retVal << "                                  const orz::DTensor &CFock,                                                 " << endl;
      if(isBareLHS_)
      retVal << "                                  const int alloc_type,                                                      " << endl;
      if(D4flag){
      retVal << "                                  const orz::mr::RdmPack &rdmPack,                                           " << endl;
      retVal << "                                  const orz::DTensor &rdm4,                                                  " << endl;
      }
      if(!LTensor_.get_indices().size() && !isBareLHS_)
      retVal << "                                  const double init_value,                                                   " << endl;
      if(Bareflag)
      retVal << "                                  const orz::mr::BareAmpPack &" + name_amp_ + ",                             " << endl;
      if(DTenflag)
      retVal << "                                  const orz::DTensor &" + T1_name() + ",                                     " << endl;
      retVal << "                                  const int num_sigma";
      if(Consts.size()) retVal << "," << endl;
      
      for(vector<string>::iterator c = Consts.begin();c != Consts.end();++c){
      retVal << "                                  const double " + *c;
      if(*c != Consts.back()) retVal << "," << endl;     
      } // End c
      retVal << ");" << endl << endl;           
    } // End scope
    //////////////////////////////////////////////////////////////////////////////////////////////

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
    if(LTensor_.get_name() != D4C_nameL()){
    if(isBareLHS_){
    CPfile << "  std::ostringstream stm;                                                                                        " << endl;                       
    CPfile << "  stm << num_sigma;                                                                                              " << endl;                        
    CPfile << "  std::string name_of_sigma = \"" + LTensor_.get_name() +  "\" + stm.str() + \"]\"; // Name of the Sigma vector  " << endl;                              
    CPfile << "  orz::mr::BareAmpPack retval                                                                                    " << endl;                         
    CPfile << "    = orz::mr::BareAmpPack(ctinp, symblockinfo, name_of_sigma, alloc_type); // Sigma(a, a', e, e') tensor        " << endl;                 
    CPfile << "                                                                                                                 " << endl;                           
    CPfile << "  orz::DTensor " + LTensor_.get_name() + "b; // Container of S2_aae,[b] tensor                                   " << endl;                       
    CPfile << "                                                                                                                 " << endl;
    } // End if
    else if(!LTensor_.get_indices().size() && Bareflag)
    CPfile << "  double " + LTensor_.get_name() + " = init_value;                                                               " << endl;
    else if( LTensor_.get_indices().size() && !Bareflag) {
    //CPfile << "  orz:DTensor " + LTensor_.get_name() + ";                                                                       " << endl;
    CPfile << "  orz::DTensor retval(";
    string Lind_name("orz::mr::sizeof_sympack_X");
    int Ccount(0);
    int Ocount(0);
    int Vcount(0);
    for(size_t num_i = 0;num_i < LTensor_.get_indices().size();++num_i) {
      if     (LTensor_.get_indices()[num_i]->get_char() == core) ++Ccount; 
      else if(LTensor_.get_indices()[num_i]->get_char() == act ) ++Ocount; 
      else if(LTensor_.get_indices()[num_i]->get_char() == virt) ++Vcount;
    } // End if
    for(int c = 0;c < Ccount;++c) Lind_name += "c";
    for(int o = 0;o < Ocount;++o) Lind_name += "a";
    for(int v = 0;v < Vcount;++v) Lind_name += "v";
    Lind_name += "(symblockinfo, 0)";    

    CPfile << Lind_name << ");" << endl;
    CPfile << endl;
    }
    if(Bareflag)
    CPfile << "  orz::DTensor " + name_amp_ + "b; // Container of T2_aae,[b] tensor                                             " << endl;                       
    if(D4flag){
    CPfile << "  orz::DTensor rdm4_sym;                                                                                         " << endl;
    CPfile << "  orz::DTensor rdm4_ij_sliced(ctinp.use_d4cum_of() ? nocc*nocc*nocc*nocc*nocc*nocc : 0);                         " << endl;
    }
    } // End if (D4C)

    // Print OpenMPI parameters
    CPfile << "  // set nproc, myrank                      " << endl;                           
    CPfile << "  const int nproc = orz::world().size();    " << endl;                    
    CPfile << "  const int myrank = orz::world().rank();   " << endl;                     
    CPfile << endl;
    // Setting up some integrals 
    CPfile << "  orz::DTensor moint1 = hintmo.int1(); // Setting up one-body integrals                                         " << endl;                        
    CPfile << "  const orz::DTensor moint1_sym = (myrank == 0) ? orz::mr::sympack_int1(symblockinfo, moint1) : orz::DTensor(); // moint1=(IR-COV index)" << endl;         
    CPfile << "  orz::DTensor " + name_h2_ + "(nmo,nmo,nmo);                                                                    " << endl;                             
    CPfile << "  double * const " + name_h2_ + "_ptr = " + name_h2_ + ".cptr();                                                  " << endl;                          
    CPfile << endl;

    // Preparation for the timing object
    if(do_timing_){
      CPfile << "  // Timing object"                     << endl;
      CPfile << "  orz::ProgressTimer time_sigma(false);" << endl << endl;
    } // End if

    // Time each term or contraction for verbose purpose
#ifdef _VERBOSE_MODE
    CPfile << "#ifdef _VERBOSE" << endl;
    CPfile << "  // Timing to write (VERBOSE_MODE)" << endl;
    CPfile << "  std::ostringstream m_stm;" << endl;
    CPfile << "  m_stm << myrank;" << endl;
    CPfile << "  std::string file_name(\"time_" + title_ + "[\" + m_stm.str() + \"].out\");" << endl;
    CPfile << "  std::ofstream outfile(file_name.c_str());" << endl;
    CPfile << "#endif" << endl << endl;
    //CPfile << "std::ofstream outfile(\"time_" + title_ + "[].out\");" << endl << endl;
#endif

    // Print logos
    cout << endl;
    cout << Femto_logo("");
    cout << endl;

    // Print guard if specific orbital type is empty
    int is_guarded(0);
    if(guard_core_)
      { CPfile << "  if(nclosed){" << endl; ++is_guarded; }
    if(guard_act_)
      { CPfile << "  if(nocc){" << endl; ++is_guarded; }

    cout << "LTensor_ : " << LTensor_ << endl;
    for(auto ts = binaries_.begin();ts != binaries_.end();++ts){
      SQcont<SQbinary> temp(*ts);
      cout << temp << endl << endl;
    }

    // Rotate each index in inTerms_ and set the external flag
    set_and_regulate_indices(binaries_);

#ifndef _GEN_DEBUG
    // Decompose each group of terms into the groups of optimal binary comtractions
    if (FafBinaries_["noeri"].size()){
      generate_contract_new(CPfile, CHfile, F90file, FafBinaries_["noeri"], "noeri");
    } // End if

    ///////////////////////////////////////////////////////////////////////////////////
    CPfile << "//-@ERI.contractions(begin)" << endl << endl;
    if(FafBinaries_["eri_c"].size()){
      generate_contract_new(CPfile, CHfile, F90file, FafBinaries_["eri_c"], "eri_c");
    } // End if
    if(FafBinaries_["eri_o"].size()){
      generate_contract_new(CPfile, CHfile, F90file, FafBinaries_["eri_o"], "eri_o");
    } // End if
    if(FafBinaries_["eri_v"].size()){
      generate_contract_new(CPfile, CHfile, F90file, FafBinaries_["eri_v"], "eri_v");
    } // End if
    CPfile << "//-@ERI.contractions(end)" << endl << endl;
    ///////////////////////////////////////////////////////////////////////////////////
    CPfile << "//-@D4C.contractions(begin)" << endl << endl;
    if(FafBinaries_["d4c_c"].size()){
      generate_contract_new(CPfile, CHfile, F90file, FafBinaries_["d4c_c"], "d4c_c");
    } // End if
    if(FafBinaries_["d4c_o"].size()){
      generate_contract_new(CPfile, CHfile, F90file, FafBinaries_["d4c_o"], "d4c_o");
    } // End if
    if(FafBinaries_["d4c_v"].size()){
      generate_contract_new(CPfile, CHfile, F90file, FafBinaries_["d4c_v"], "d4c_v");
    } // End if
    CPfile << "//-@D4C.contractions(end)"   << endl << endl;
#endif

    for(int i = 0;i < is_guarded;++i)
      //if(is_guarded)
      CPfile << "  } // Guard" << endl;

    if(do_timing_){
     CPfile << "  // Do timing!" << endl;
     CPfile << "  orz::mr::Femto::my_timer.push_back(boost::make_tuple(" << "\"" << title_ << "\", time_sigma.elapsed_cputime(), time_sigma.elapsed_wallclocktime()));" << endl;
     CPfile << "  orz::mr::Femto::file_timing << \"* \" << boost::format(\"%20s : %10.7f %10.7f \") % orz::mr::Femto::my_timer.back().get<0>() % orz::mr::Femto::my_timer.back().get<1>() % orz::mr::Femto::my_timer.back().get<2>() << endl;" << endl;
     CPfile << "  flush(orz::mr::Femto::file_timing);" << endl;
//*old*//       CPfile << "  // Do timing!" << endl;
//*old*//       CPfile << "  my_timer.push_back(boost::make_tuple(" << "\"" << title_ << "\", time_sigma.elapsed_cputime(), time_sigma.elapsed_wallclocktime()));" << endl;
//*old*//       CPfile << "  file_timing << \"* \" << boost::format(\"%20s : %10.7f %10.7f \") % my_timer.back().get<0>() % my_timer.back().get<1>() % my_timer.back().get<2>() << endl;" << endl;
//*old*//       CPfile << "  flush(file_timing);" << endl;
    } // End if
    CPfile << endl;

    if(LTensor_.get_name() != D4C_nameL()){
      if(!LTensor_.get_indices().size() && Bareflag){
	
	CPfile << "  double sum_" +  LTensor_.get_name() << " = 0.0;" << endl;
	CPfile << "  boost::mpi::reduce(orz::world()," + LTensor_.get_name() + ", sum_" + LTensor_.get_name() + ", std::plus<double>(), 0);" << endl;
	CPfile << "  " << LTensor_.get_name()  << " = sum_" + LTensor_.get_name() << ";"<< endl;        
	
	CPfile << endl;
	CPfile << "  return  " << LTensor_.get_name()  << ";" << endl;
	CPfile << "} " << endl;
      } // End if
      else{
	
	if(LTensor_.get_indices().size() && !isBareLHS_){
	  CPfile << "  orz::DTensor retval_(orz::tensor::mpi::reduce_plus(retval, 0));" << endl;
	  CPfile << "  return retval_; " << endl;
	} // End if      
	else {
	  CPfile << "  return retval; " << endl;
	} // Else
      CPfile << "} " << endl;
      } // End else
    } // End if
    else{
      CPfile << "} " << endl;
    } // End else


    CHfile << CHend;
    
    return retVal.str();
  }                                               

}} //End Femto
