//
//  SQdf_aux.cc
//  
//
//  Created by Masaaki Saitow on 12/11/27.
//  Copyright (c) 2012 Masaaki Saitow. All rights reserved.
//

#include <iostream>
#include <SQdf_aux.hpp>

using namespace std;

namespace Femto { namespace Core {
  

  // *********************************************************
  // 
  // *********************************************************
  SQdf_aux::SQdf_aux(const string name, const bool extFlag)
  {
    index_     = name;
    charactar_ = aux;  // must be this
    isSummed_  = true;
    isExt_     = extFlag;
  }


}} //Femto::core
