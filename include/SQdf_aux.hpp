//
//  SQdf_aux.hpp
//  
//  Class for the orbital symbol, which is used as an index of the SQtensor object.
//  Unlikely to SQindex, this class is intended to represent none of the generic orbital symbols,
//  but the auxiliary basis used in the DF/RI like decomposition of the ERI.
//
//  Created by Masaaki Saitow on 12/11/27.
//  Copyright (c) 2012 Masaaki Saitow. All rights reserved.
//

#pragma once

#include <cstring>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/nvp.hpp>
#include <SQindex.hpp>
#include <Femto.hpp>


namespace Femto { namespace Core {
  
  class SQdf_aux : public SQindex{
      
  public:
    SQdf_aux(const std::string name, const bool extdFlag=false);
      
  private:    
    // Some stuffs for serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int ver){
      ar & boost::serialization::make_nvp("index",     index_);
      ar & boost::serialization::make_nvp("charactar", charactar_);
      ar & boost::serialization::make_nvp("isSummed",  isSummed_);
      ar & boost::serialization::make_nvp("isExt",     isExt_);
    } 

  };
  
}} // Femto::core

