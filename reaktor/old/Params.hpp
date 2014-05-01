//
//  Params.hpp
//
//  Some auxiliary parameters are defined in this file  
//
//
//  Created by Masaaki Saitow on 12/03/26.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#pragma once

#include <cmath>
//#include <femto.hpp>
//#include <SQterm.hpp>
//#include <SQtensor.hpp>
//#include <SQindex.hpp>
//#include <Orz.hpp>

//using namespace std;

namespace femto { namespace Reaktor {

  // Priority flag
  enum priority {polynomial=0, memory=1};

  inline priority is_priority()
  { return polynomial; } 

  // Whether ERI is processed with the top priority
  inline bool h2_prior()
  { return  true; }

  // Whether intermediate is processed in ad hoc fashoion in every situation
  inline bool isOF()
  { return false; } // Currently only option allowed is true!!!

  // Number of core orbitals in the model systems 
  inline int numCore()
  { return 30; }

  // Number of active  orbitals in the model systems 
  inline int numOcc()
  { return 30; }

  // Number of virtual orbitals in the model systems 
  inline int numVirt()
  { return 500; }

  // Whether cumulant decomposition is turned on
  inline bool use_cumulant()
  { return false; }

  // Whether openmp is turned on
  inline bool use_omp()
  { return false; }

  // Whether dgemm is utilized
  inline bool use_gemm()
  { return true; }

  // Whether time the tensor contraction
  inline bool do_timing()
  { return true; }

  // Maximum length of the tensor that can be processed in ad hoc fashion
  inline unsigned long int MaxSize()
  { return pow(numOcc(),6); }

  // External index of ERI
  inline int exth2()
  { return 0; }

  // Name of the one-body integral
  inline string name_h1()
  { return "h"; }

  // External index of BareAmp
  inline int extamp()
  { return 3; }

}}

