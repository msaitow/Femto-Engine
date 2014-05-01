//
//  SQterm.hpp
//  
//  Class that represents the term object, composed of factor, numerical constants, tensors
//  and the entity of the indices commomly shared by the tensors
//
//  Created by Masaaki Saitow on 12/03/26.
//  Copyright (c) 2012 Masaaki Saitow. All rights reserved.
//


#pragma once

#include <cstring>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <utility>
#include <Femto.hpp>
#include <SQtensor.hpp>
#include <SQindex.hpp>

//#define DEBUG_ false

namespace Femto { namespace Core { 

  class SQterm{

  public:
    SQterm();
    ~SQterm();
    SQterm(const SQterm &obj);
    SQterm(const double numConst, const std::vector<std::string> Consts, 
           const std::vector<SQtensor> Tensors, const bool isInCanonical=false);

    SQterm  operator=(const SQterm &obj);
    SQterm  operator*(const SQterm &obj);
    bool    operator==(const SQterm &obj);

    friend std::ostream &operator<<(std::ostream &os, const SQterm &t);

    std::vector<SQtensor> get_tensors() const;
    std::vector<SQtensor*> get_tensors_ptr();
    std::vector<SQindex*> get_summedBody();
    double get_numConst() const;
    std::vector<std::string> get_Consts() const;
    bool get_isInCanonical() const;
    
    void set_isInCanonical(const bool isInCanonical);
    void set_tensors(const std::vector<SQtensor> &Tensors);
    void set_numConst(const double num);
    void set_Consts(const std::vector<std::string> consts);
    void print_summedBody();
    std::string convert2LaTeX() const;
    void masquerade();

    void transform2RDM();    // Convert sfGens to the conrresponding RDMs, only if isInCanonical = true
    void contractkDeltas();  // Contract Kronecker's deltas
    void decomposeRDMVirt(); // Decompose RDM, or sfGen if it contains virtual index

  private:
    void set_summedBody();
    void erase_duplication();

    bool isInCanonical_;                 // Whether this term is in canonical form
    double numConst_;                    // Numerical factor of this term
    std::vector<std::string> Consts_;    // Constants of this term
    std::vector<SQtensor> Tensors_;      // Tensors of this term
    std::vector<SQindex> summedIndices_; // Body of summed indices, of which shared_ptrs are shared
                                         // among the Tensors_ as vector of the SQindex(not as pointer)   

    // Some stuffs for serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int ver){
      ar & boost::serialization::make_nvp("isInCanonical", isInCanonical_);
      ar & boost::serialization::make_nvp("numConst",      numConst_);
      ar & boost::serialization::make_nvp("Consts",        Consts_);
      ar & boost::serialization::make_nvp("Tensors",       Tensors_);
      ar & boost::serialization::make_nvp("summedIndices", summedIndices_);
    }
   
  }; 

  // Returns a term composed of term a multiplied by term b 
  SQterm  times(SQterm &a, SQterm &b);
  // Whether two terms are factorizable or not (in other words, same except numConst and Consts)
  bool isFactorizable(SQterm &a, SQterm &b);
  // Whether two terms are factorizable or not (in other words, same except numConst and Consts)
  bool isFactorizable2(SQterm &a, SQterm &b);
  // Whether two terms givens are additive or not
  bool isAdditive(SQterm &a, SQterm &b);
  // Delete all the terms with zero facto
  void screenTerms(const std::vector<SQterm> &inTerms, std::vector<SQterm> *outTerms, const double crit=0.0);
  // Delete all the terms with zero factor
  std::vector<SQterm> screenTerms(const std::vector<SQterm> &inTerms, const double crit=0.0);

  // Normal order inTerm, which is assumed to be composed of the usual commutable tensors like T2 and,
  // the uncommutable spin-free unitary group generators 
  void normalOrder(std::vector<SQterm> *inTerms, const bool OF_flag=true);
  // In case only spin-free generator exist in given term
  void normalOrder_sf(std::vector<SQterm> *inTerms, const bool OF_flag=true);

  // Normal order inTerm in terms of the commutator manner, in which all the spin-free generators
  // are supposed to be composed of multiple commutators from right most pairs.
  // For example, if inTerm is given like this, {T2, Ea, Eb, Ec}, which is 
  // interpreted as normal ordering of T2 [Ea,[Eb,Ec]
  void normalOrderComm(std::vector<SQterm> *inTerms, const bool OF_flag=true);  

  // Combine terms with isAdditive() == true
  std::vector<SQterm> combineTerms(std::vector<SQterm> &inTerms);
  // Combine terms with isAdditive() == true
  void combineTerms(std::vector<SQterm> &inTerms, std::vector<SQterm> *outTerms);
  // Combine terms with isAdditive() == true
  void combineTerms2(std::vector<SQterm> &inTerms, std::vector<SQterm> *outTerms);
  // Combine terms with isAdditive() == true
  void combineTerms_new(std::vector<SQterm> &inTerms);

  // contracts core indices from the RDMs
  void decomposeRDMCore(SQterm &inTerm, std::vector<SQterm> *outTerms);

  // Function to construct the core Fock matrix
  void makeFock(std::vector<SQterm> &inTerm);

}} // Femto::core
