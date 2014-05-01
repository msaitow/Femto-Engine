//
//  SQcontract.hpp
//  
//  Class that represents the tensor contraction, which is composed of a tensor on LHS 
//  and at least one tensor on RHS. And the associated indices of tensors on both LHS and
//  RHS are linked and stored as summedIndices_
//
//  Created by Masaaki Saitow on 12/11/12.
//  Copyright (c) 2012 Masaaki Saitow. All rights reserved.
//


#pragma once

#include <cstring>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <utility>
#include <SQindex.hpp>
#include <SQtensor.hpp>
#include <Femto.hpp>


namespace Femto { namespace Core {

  // Class that represents the tensor contraction
  class SQcontract{
  
  public:
    SQcontract();
    SQcontract(const SQcontract &obj);
    SQcontract(const double numConst, const std::vector<std::string> Consts,
               const SQtensor &Ltensor, const std::vector<SQtensor> &Rtensors);

    SQcontract operator=(const SQcontract &obj);  
    friend std::ostream &operator<<(std::ostream &os, const SQcontract &t);

    SQtensor get_Ltensor() const;
    SQtensor* get_Ltensor_ptr();
    std::vector<SQtensor> get_Rtensors() const;
    std::vector<SQtensor*> get_Rtensors_ptr();
    std::vector<SQindex*> get_summedBody();
    double get_numConst() const;
    std::vector<std::string> get_Consts() const;

    // These three functions are for specifications on the internal indices
    // that defines the body of the Ltensor as a tensor quantities at C++/FORTRAN levels
    void set_Lindices(bool val);          
    void set_Rindices(int num, bool val);
    bool get_Lindices();          
    std::vector<bool> get_Rindices();

    std::vector<SQindex*> get_RInnerIndices(int num);
    std::vector<SQindex*> get_LInnerIndices();
    std::vector<SQindex*> get_ROuterIndices(int num);
    std::vector<SQindex*> get_LOuterIndices();
    void set_ROuterIndices(int num, Femto::SQcont<Femto::Core::SQindex> inds);
    void set_LOuterIndices(Femto::SQcont<Femto::Core::SQindex> inds);
    void clear_ROuterIndices(int num);
    void clear_LOuterIndices();

    //*-- --*// std::vector<SQindex*> get_Lindices();
    //*-- --*// void set_Lindices(std::vector<std::string> &name_list); 
    //*-- --*// void clear_Lindices();

    void set_Ltensor(const SQtensor &Tensor);
    void set_Rtensors(const std::vector<SQtensor> &Tensors);
    void set_numConst(const double num);
    void set_Consts(const std::vector<std::string> consts);
    void masquerade();
    void print_summedBody();
    void contractkDeltas();

  protected:
    void set_summedBody();

    double   numConst_;               // Numerical factor
    std::vector<std::string> Consts_; // Constants of this contraction
    //*-- --*// std::vector<std::string> InternalNames_;   // Names of indices of the FORTRAN level representation of Ltensor_,
    //*-- --*//                                  // if nothing is set to this, internal indices are used

    bool Lindices_;                  // If this is set true, all the indices of Ltensor_ becomes the indices used at the FORTRAN level. 
                                     // By default this is set to false, which means only internal indices are used to allocate Ltensor_ at the FORTRAN level
    std::vector<bool> Rindices_;     // If this is set true, all the indices of Rtensors_ becomes the indices used at the FORTRAN level. 
                                     // By default this is set to false, which means only internal indices are used to allocate Rtensors_ at the FORTRAN level

    std::vector<Femto::Core::SQindex*> LInnerInds_;
    std::vector<std::vector<Femto::Core::SQindex*> > RInnerInds_;

    SQtensor Ltensor_;                    // Tensor on the left-hand side
    std::vector<SQtensor> Rtensors_;      // Tensors on the right-hand side
    std::vector<SQindex> summedIndices_;  // Body of summed indices, of which shared_ptrs are shared
                                          // among the Tensors_ as vector of the SQindex(not as pointer)
  private:
    // Some stuffs for serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int ver){
      ar & boost::serialization::make_nvp("numConst",       numConst_);
      ar & boost::serialization::make_nvp("Consts",         Consts_);
      ar & boost::serialization::make_nvp("Ltensor",        Ltensor_);
      ar & boost::serialization::make_nvp("Rtensors",       Rtensors_);
      ar & boost::serialization::make_nvp("Lindices_",      Lindices_);
      ar & boost::serialization::make_nvp("Rindices_",      Rindices_);
      ar & boost::serialization::make_nvp("Rtensors",       Rtensors_);
      ar & boost::serialization::make_nvp("summedIndices_", summedIndices_);
    }

  };
  
  // Whether two tensor contractions are factorizable or not
  bool isFactorizable(SQcontract &a, SQcontract &b);

}} // Femto::SQcontract

