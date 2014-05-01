//
//  SQbinary.hpp
//  
//  Class that represents the binary contraction, which is composed of a tensor on LHS 
//  and at least one tensor on RHS. Unlikely to SQcontract, SQbinary contains up to only
//  two tensors on RHS 
//
//  Created by Masaaki Saitow on 12/11/12.
//  Copyright (c) 2012 Masaaki Saitow. All rights reserved.
//


#pragma once

#include <SQcontract.hpp>


namespace Femto { namespace Core {

  // Class that represents the binary contraction
  class SQbinary : public SQcontract{
  public:
    SQbinary();
    //SQbinary(const SQbinary &obj);
    SQbinary(const double numConst, const std::vector<std::string> Consts,
             const SQtensor &Ltensor, const std::vector<SQtensor> &Rtensors);
    
  private:
    // Some stuffs for serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int ver){
      ar & boost::serialization::base_object<SQtensor>(*this);
    }     
  };

}} // Femto::SQbinary

