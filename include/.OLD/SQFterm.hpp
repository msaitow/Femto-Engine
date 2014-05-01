//
//  SQFtensor.hpp
//  
//  Class that represents the factorized term object, composed of factor, numerical constants, 
//  a vector of tensors and a tensor. Structure of this term is like,
//
//  Const [ c1 Term1 + c2 Term2 + ... ] Terms
//
//  Created by Masaaki Saitow on 12/03/26.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//


#pragma once

#include <cstring>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <utility>
#include <femto.hpp>
#include <SQtensor.hpp>
#include <SQindex.hpp>

using namespace std;

namespace femto {

  class SQFterm{

  public:
    SQFterm();
    SQFterm(const SQFterm &obj);
    SQFterm(const vector<double>   numConsts, const vector<string> Consts, 
            const vector<SQtensor> Tensors1,  const vector<SQtensor> Tensors2);

    friend std::ostream &operator<<(std::ostream &os, const SQterm &t);

    vector<SQtensor> get_tensors1() const;
    vector<SQtensor*> get_tensors1_ptr();
    vector<SQtensor> get_tensors2() const;
    vector<SQtensor*> get_tensors2_ptr();

    vector<string> get_Consts() const;
    
    void set_tensors1(const vector<SQtensor> &Tensors);
    void set_tensors2(const vector<SQtensor> &Tensors);
    void set_numConst(const vector<double> nums);
    void set_Consts(const vector<string> consts);

  private:

    vector<string> Consts_;        // Constants of this term
    vector<double> numConsts_;     // Numerical coefficients of tensors1 
    vector<SQtensor> Tensors1_;    // Tensors of this term
    vector<SQtensor> Tensors2_;    // Tensors of this term
  }; 


} // femto::term
