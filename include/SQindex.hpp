//
//  SQindex.hpp
//  
//  Class for the orbital symbol, which is used as an index of the SQtensor object
//
//  Created by Masaaki Saitow on 12/03/26.
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
#include <boost/shared_ptr.hpp>
#include <Femto.hpp>

namespace Femto { namespace Core {
  
  class SQindex{
      
  public:
    SQindex(){};
    SQindex(const SQindex &obj);
    SQindex(const std::string name, const enum char_state nature, 
            const bool summedFlag=false, const bool extdFlag=false);
      
    std::string get_index() const;
    enum char_state get_char() const;
    bool get_isSummed() const;
    bool get_isExt() const;
    bool operator>(const SQindex &other) const;  // Compare except isExt_
    bool operator<(const SQindex &other) const;  // Compare except isExt_
    bool operator==(const SQindex &other) const; // Compare except isExt_
    bool operator!=(const SQindex &other) const; // Compare except isExt_
    bool operator%=(const SQindex &other) const; // Compare all the members
    void switch_isExt(const bool extFlag);
    void put_index(const std::string name);

    void copy(const SQindex &other); // Copy the content of other to this. If DF/RI is not utilized, or in the expression generating step,
                                     // use of this method is unrecommendable. 
    
    friend std::ostream &operator<<(std::ostream &os, const SQindex &i);
    SQindex &operator=(const SQindex &obj);
      
  protected:
    std::string index_;         // Name of the index
    enum char_state charactar_; // Orbital type
    bool isSummed_;             // Whether dummy or not
    bool isExt_;                // whether external or not
                                // Here, the external index means the one whose the associated loop is declared at C++ level.
                                // In usual the ERI, T2-amplitude, or 4-RDM is loaded according to the external indices.   

  private:    
    // Some stuffs for serialization
    // 2013/01/07 :: This serialization function is verified to work correctly (MS)
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

typedef boost::shared_ptr<Femto::Core::SQindex>  SharedIndex;

