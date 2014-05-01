//
//  SQinput.cc
//  
//
//  Created by Masaaki Saitow on 12/11/16.
//  Copyright (c) 2012 Masaaki Saitow. All rights reserved.
//

#include <cmath>
#include <iostream>
#include <SQinput.hpp>

using namespace std;
using namespace Femto;
//using namespace Femto::Core;

namespace Femto { namespace Reaktor2 { 
  
  // *********************************************************
  // 
  // *********************************************************
  SQinput::SQinput(const string title, const bool isBareLHS)
    : title_(title),
      isBareLHS_(isBareLHS) // Whether LHS is represented by BareAmpPack object
  {
    cout << "Setting up parameters as default .... " << endl;
    ///////// Default parameters ///////////
    is_priority_  = polynomial;
    numCore_      = 30;
    numOcc_       = 30;
    numVirt_      = 500;
    use_cumulant_ = false;
    use_omp_      = false;
    use_gemm_     = true;
    set_V_D_      = true;
    do_timing_    = true;
    MaxSize_      = pow(numOcc_, 6);
    name_h1_      = "h"; 
    name_h2_      = "V2";
    name_t1_      = "T1";
    name_t2_      = "T2";
    name_d4_      = RDM_name() + "4";

    ///////// External indices ////////////
    exth2_        = 0;
    extamp_       = 3;
    ///////////////////////////////////////
  }

  // *********************************************************
  // 
  // *********************************************************
  SQinput::~SQinput(){ }

  // *********************************************************
  // 
  // *********************************************************
  void
  SQinput::set_priority(const priority flag)
  { is_priority_ = flag; }

  // *********************************************************
  // 
  // *********************************************************
  void
  SQinput::set_use_cumulant(const bool flag)
  { use_cumulant_ = flag; }

  // *********************************************************
  // 
  // *********************************************************
  void
  SQinput::set_use_omp(const bool flag)
  { use_omp_ = flag; }

  // *********************************************************
  // 
  // *********************************************************
  void
  SQinput::set_use_gemm(const bool flag)
  { use_gemm_ = flag; }

  // *********************************************************
  // 
  // *********************************************************
  void
  SQinput::set_V_D(const bool flag)
  { set_V_D_ = flag; }

  // *********************************************************
  // 
  // *********************************************************
  void
  SQinput::set_do_timing(const bool flag)
  { do_timing_ = flag; }

  // *********************************************************
  // 
  // *********************************************************
  void
  SQinput::set_name_h1(const string &name)
  { name_h1_ = name; }

  // *********************************************************
  // 
  // *********************************************************
  void
  SQinput::set_name_h2(const string &name)
  { name_h2_ = name; }

  // *********************************************************
  // 
  // *********************************************************
  void
  SQinput::set_name_t1(const string &name)
  { name_t1_ = name; }

  // *********************************************************
  // 
  // *********************************************************
  void
  SQinput::set_name_t2(const string &name)
  { name_t2_ = name; }

  // *********************************************************
  // 
  // *********************************************************
  void
  SQinput::set_numCore(const int num)
  { numCore_ = num; }

  // *********************************************************
  // 
  // *********************************************************
  void
  SQinput::set_numOcc(const int num)
  { numOcc_ = num; }

  // *********************************************************
  // 
  // *********************************************************
  void
  SQinput::set_numVirt(const int num)
  { numVirt_ = num; }

  // *********************************************************
  // 
  // *********************************************************
  void 
  SQinput::set_exth2(const int num){
    if(num < 0 || num > 3){
      cout << "SQinput: exth2 cannot exceed 0 - 3." << endl;
      abort();
    } // End if
    exth2_ = num;
  }

  // *********************************************************
  // 
  // *********************************************************
  void 
  SQinput::set_extamp(const int num){
    if(num < 0 || num > 3){
      cout << "SQinput: extamp cannot exceed 0 - 3." << endl;
      abort();
    } // End if
    extamp_ = num;
  }

  // *********************************************************
  // 
  // *********************************************************
  priority
  SQinput::get_priority() const
  { return is_priority_; }

  // *********************************************************
  // 
  // *********************************************************
  bool 
  SQinput::get_use_cumulant() const
  { return use_cumulant_; }

  // *********************************************************
  // 
  // *********************************************************
  bool 
  SQinput::get_use_omp() const
  { return use_omp_; }

  // *********************************************************
  // 
  // *********************************************************
  bool 
  SQinput::get_use_gemm() const
  { return use_gemm_; }

  // *********************************************************
  // 
  // *********************************************************
  bool 
  SQinput::get_do_timing() const
  { return do_timing_; }

  // *********************************************************
  // 
  // *********************************************************
  bool 
  SQinput::get_V_D() const
  { return set_V_D_; }

  // *********************************************************
  // 
  // *********************************************************
  ostream&
  operator<<(std::ostream &os, const SQinput &i)
  {

    os << ">>======================================================<<" << endl;
    os << ">> SQinput :: An Input Class for SQreaktor              <<" << endl;
    os << ">>======================================================<<" << endl << endl; 

    os << ">> Priority is set to ..................... " << (i.get_priority() == polynomial ? "Polynomial" : "memory") << endl;
    os << ">> Use of cumulant 4-RDM .................. " << (i.get_use_cumulant() ? "Y" : "N") << endl;
    os << ">> Use of OpenMP .......................... " << (i.get_use_omp()      ? "Y" : "N") << endl;
    os << ">> Use of dgemm ........................... " << (i.get_use_gemm()     ? "Y" : "N") << endl;
    os << ">> Do timing .............................. " << (i.get_do_timing()    ? "Y" : "N") << endl;
    os << ">> Set loading indices of ERI and 4-RDM ... " << (i.get_V_D()          ? "Y" : "N") << endl;
    os << ">> Name of one-dody integral .............. " << i.get_name_h1()                    << endl; 
    os << ">> Name of two-dody integral .............. " << i.get_name_h2()                    << endl; 
    os << ">> Name of T1-amplitude ................... " << i.get_name_t1()                    << endl; 
    os << ">> Name of T2-amplitude ................... " << i.get_name_t2()                    << endl; 
    os << ">> Loading index of BareAmpPack ........... " << i.get_extamp() << "-th"            << endl;
    os << ">> Loading index of two-body integral ..... " << i.get_exth2()  << "-th"            << endl;
    os << ">> Number of core orbitals[target] ........ " << i.get_numCore()                    << endl; 
    os << ">> Number of active orbitals[target] ...... " << i.get_numOcc()                     << endl;
    os << ">> Number of virtual orbitals[target] ..... " << i.get_numVirt()                    << endl << endl;
  
    os << ">>======================================================<<" << endl << endl; 

    return os;
  }

}} //Femto::
