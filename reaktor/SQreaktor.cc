//
//  SQreaktor.cc
//  
//
//  Created by Masaaki Saitow on 12/09/04.
//  Copyright (c) 2012 Masaaki Saitow. All rights reserved.
//

#include <cmath>
#include <iostream>
#include <SQreaktor.hpp>

using namespace std;
using namespace Femto;
using namespace Femto::Core;

namespace Femto { namespace Reaktor { 
  

  // *********************************************************
  // 
  // *********************************************************
  SQreaktor::SQreaktor(const SQtensor &LTensor, const vector<SQterm> &inTerms, const string title, 
                        const bool isBareLHS, const string name_h2, const string name_amp)
    : LTensor_(LTensor), 
      inTerms_(inTerms),
      title_(title),
      isBareLHS_(isBareLHS), // Whether LHS is represented by BareAmpPack object
      name_h2_(name_h2),     // Neme of ERI 
      name_amp_(name_amp),   // Name of *T2* amplitude
      exth2_(0),             // External indices in ERI
      extamp_(3),            // External indices in BareAmp
      extd4c_(5)             // External indices in D4C
  {
    cout << "Setting up parameters as default .... " << endl;
    ///////// Default parameters ///////////
    use_oldgemm_  = false;
    is_priority_  = polynomial;
    h2_prior_     = true;
    isOF_         = false;
    numCore_      = 30;
    numOcc_       = 30;
    numVirt_      = 500;
    use_cumulant_ = false;
    use_omp_      = false;
    use_gemm_     = true;
    Set_V_D_      = true;
    Set_T_D_      = false; // <-- Newly added 2013/01/08
    Set_D4C_      = false; // <-- NEW
    do_timing_    = true;
    guard_core_   = false;
    guard_act_    = false;
    MaxSize_      = pow(numOcc_, 6);
    name_h1_      = "h";
    name_d4_      = RDM_name() + "4";
    ///////// External indices ////////////
    //exth2_        = 0;
    //extamp_       = 3;
    ///////////////////////////////////////
  }

  // *********************************************************
  // 
  // *********************************************************
  SQreaktor::SQreaktor(const vector<vector<Femto::Core::SQbinary> > &inContras, const string title, 
                        const bool isBareLHS, const string name_h2, const string name_amp)
    : binaries_(inContras),
      LTensor_(inContras[0].back().get_Ltensor()),
      title_(title),
      isBareLHS_(isBareLHS), // Whether LHS is represented by BareAmpPack object
      name_h2_(name_h2),     // Neme of ERI 
      name_amp_(name_amp),   // Name of *T2* amplitude
      exth2_(0),             // External indices in ERI
      extamp_(3),            // External indices in BareAmp
      extd4c_(5)             // External indices in D4C
  {
    cout << "Setting up parameters as default .... " << endl;
    ///////// Default parameters ///////////
    use_oldgemm_  = false;
    is_priority_  = polynomial;
    h2_prior_     = true;
    isOF_         = false;
    numCore_      = 30;
    numOcc_       = 30;
    numVirt_      = 500;
    use_cumulant_ = false;
    use_omp_      = false;
    use_gemm_     = true;
    Set_V_D_      = true;
    Set_T_D_      = false; // <-- Newly added 2013/01/08
    Set_D4C_      = false; // <-- NEW
    do_timing_    = true;
    guard_core_   = false;
    guard_act_    = false;
    MaxSize_      = pow(numOcc_, 6);
    name_h1_      = "h";
    name_d4_      = RDM_name() + "4";
    ///////// External indices ////////////
    //exth2_        = 0;
    //extamp_       = 3;
    ///////////////////////////////////////
  }

  // *********************************************************
  // 
  // *********************************************************
  //SQreaktor::~SQreaktor(){ }

  // *********************************************************
  // 
  // *********************************************************
  std::string
  SQreaktor::generate(const Codes code_name, const Modes mode_name)
  {

//     cout << "LTensor_ " << LTensor_ << endl;
//     for(auto ts = binaries_.begin();ts != binaries_.end();++ts){
//       SQcont<SQbinary> temp(*ts);
//       cout << temp << endl << endl;
//     }

    if(code_name == Orz){
      if      (mode_name == SimpleLoops) return simpleloops();
      //else if (mode_name == Factorize)   factorize();
      else if (mode_name == Factorize)   return factorize_new();
      //else if (mode_name == Cas)         factorize_cas();
      else if (mode_name == Fafnir)      return factorize_fafnir();
      else {
        cout << "Mode can be specified by 0=simpleloops, 1=factorize and 2=cas." << endl;
        abort();
      } // End if
    } // End if
    //else if(code_name == Psi)
    else{
      cout << "Not yet implemented:" << endl;
      abort();
    } // End else
    //~SQreaktor();
  }

  // *********************************************************
  // 
  // *********************************************************
  void
  SQreaktor::set_oldgemm(const bool flag)
  {  use_oldgemm_ = flag; }

  // *********************************************************
  // 
  // *********************************************************
  void
  SQreaktor::set_priority(const priority flag)
  { is_priority_ = flag; }

  // *********************************************************
  // 
  // *********************************************************
  void
  SQreaktor::set_h2_prior(const bool flag)
  { h2_prior_ = flag; }

  // *********************************************************
  // 
  // *********************************************************
  void
  SQreaktor::set_isOF(const bool flag)
  { isOF_ = flag; }

  // *********************************************************
  // 
  // *********************************************************
  void
  SQreaktor::set_use_cumulant(const bool flag)
  { use_cumulant_ = flag; }

  // *********************************************************
  // 
  // *********************************************************
  void
  SQreaktor::set_use_omp(const bool flag)
  { use_omp_ = flag; }

  // *********************************************************
  // 
  // *********************************************************
  void
  SQreaktor::set_use_gemm(const bool flag)
  { use_gemm_ = flag; }

  // *********************************************************
  // 
  // *********************************************************
  void
  SQreaktor::set_V_D(const bool flag)
  { Set_V_D_ = flag; }

  // *********************************************************
  // 
  // *********************************************************
  void
  SQreaktor::set_T_D(const bool flag)
  { Set_T_D_ = flag; }

  // *********************************************************
  // 
  // *********************************************************
  void
  SQreaktor::set_D4C(const bool flag)
  { Set_D4C_ = flag; }

  // *********************************************************
  // 
  // *********************************************************
  void
  SQreaktor::set_do_timing(const bool flag)
  { do_timing_ = flag; }

  // *********************************************************
  // 
  // *********************************************************
  void
  SQreaktor::set_guard_core(const bool flag)
  { guard_core_ = flag; }

  // *********************************************************
  // 
  // *********************************************************
  void
  SQreaktor::set_guard_act(const bool flag)
  { guard_act_ = flag; }

  // *********************************************************
  // 
  // *********************************************************
  void
  SQreaktor::set_name_h1(const string &name)
  { name_h1_ = name; }

  // *********************************************************
  // 
  // *********************************************************
  void
  SQreaktor::set_numCore(const int num)
  { numCore_ = num; }

  // *********************************************************
  // 
  // *********************************************************
  void
  SQreaktor::set_numOcc(const int num)
  { numOcc_ = num; }

  // *********************************************************
  // 
  // *********************************************************
  void
  SQreaktor::set_numVirt(const int num)
  { numVirt_ = num; }

  // *********************************************************
  // 
  // *********************************************************
  vector<SQterm> SQreaktor::get_inTerms() const
  { return inTerms_; }

  // *********************************************************
  // 
  // *********************************************************
  size_t SQreaktor::num_inTerms() const
  { return inTerms_.size(); }

  // *********************************************************
  // 
  // *********************************************************
  vector<SQterm> SQreaktor::get_terms_d4() const
  {
    vector<SQterm> term_d4; term_d4.reserve(inTerms_.size());
    for(vector<SQterm>::const_iterator t = inTerms_.begin();t != inTerms_.end();++t){
      vector<SQtensor> ten(t->get_tensors());
      for(vector<SQtensor>::const_iterator tt = ten.begin();tt != ten.end();++tt)
	if(tt->get_name() == name_d4_ && find(term_d4.begin(), term_d4.end(), *t) == term_d4.end()) term_d4.push_back(*t);
    } // End t
    return term_d4;
  }

  // *********************************************************
  // 
  // *********************************************************
  void SQreaktor::init_header(string file_name)
  {
    ifstream ifs(file_name.c_str());
    if(ifs.fail()){
      cout << "init_header: Can't open file that contains the information of header files" << endl;
      abort();
    } // End if

    string line;
    while(!ifs.eof()){
      ifs >> line;
      // * In case "#include<**.h>" 
      if     (line.at(0) == '#' && *(line.end()-1) == '>') { line += "\n"; cont_header_ += line; }
      // * In case "#include `space'" 
      else if(line.at(0) == '#' && *(line.end()-1) != '>') { line += " ";  cont_header_ += line; }
      // * In case `space' <**.h>
      else if(line.at(0) != '#' && *(line.end()-1) == '>') { line += "\n"; cont_header_ += line; }
//*NO_GOOD*       else { 
//*NO_GOOD* 	cout << "init_header: Maybe, there's inappropriate information detected in " << file_name << endl; 
//*NO_GOOD* 	cout << "      >> " << line << endl;
//*NO_GOOD* 	abort(); 
//*NO_GOOD*       }
    } // End while

  }

  // *********************************************************
  // 
  // *********************************************************
  ostream&
  operator<<(std::ostream &os, const SQreaktor &i)
  {
    os << "::::: Parameters used in the scheme :::::"                                     << endl;
    os << "* Priority paradigm                    : " << 
      (i.is_priority_ == memory ? "memory" : "polynomial")                               << endl;
    os << "* Intermediate is processed on-the-fly : " << (i.isOF_ ? "Yes" : "No")         << endl;
    os << "* Use cumulant in 4-RDM                : " << (i.use_cumulant_ ? "Yes" : "No") << endl;
    os << "* Use OpenMP                           : " << (i.use_omp_ ? "Yes" : "No")      << endl;
    os << "* Use BLAS in the contraction          : " << (i.use_gemm_ ? "Yes" : "No")     << endl;
    os << "* Do timing for each contraction       : " << (i.do_timing_ ? "Yes" : "No")    << endl << endl;
    os << "::::::::: Size of target system :::::::::" << endl;
    os << "* Number of core    orbitals : " << i.numCore_ << endl;
    os << "* Number of active  orbitals : " << i.numOcc_  << endl;
    os << "* Number of virtual orbitals : " << i.numVirt_ << endl;
    os << ":::::::::::::::::::::::::::::::::::::::::" << endl;
    return os;
  }

}} //Femto::
