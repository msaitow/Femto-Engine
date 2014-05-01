//
//  SQFterm.cc
//  
//
//  Created by Masaaki Saitow on 12/06/30.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <cmath>
#include <vector>
#include <algorithm>
#include <boost/format.hpp>
#include <femto.hpp>
#include <SQterm.hpp>
#include <SQtensor.hpp>
#include <SQindex.hpp>

namespace femto {

  // *********************************************************
  // 
  // *********************************************************
  SQFterm::SQFterm()
  {}

  // *********************************************************
  // 
  // *********************************************************
  SQFterm::SQFterm(const vector<double> numConsts,  const vector<string> Consts, 
                   const vector<SQtensor> Tensors1, const vector<SQtensor> Tensors2)
    : numConsts_(numConsts),
      Consts_(Consts),
  {
    vector<SQtensor> Tnon_commutes;
    for(size_t I = 0;I < Tensors1.size();++I){
      if(Tensors1[I].isCommutable()) Tensors1_.push_back(Tensors1[I]);
      else                           Tnon_commutes.push_back(Tensors1[I]);
    }
    if(Tnon_commutes.size()){
      cout << "Non-commutable tensor cannot be treated in this level .... " << endl;
      abort();
    }
    Tensors1_(Tensors1);
    for(size_t I = 0;I < Tensors2.size();++I){
      if(Tensors2[I].isCommutable()) Tensors2_.push_back(Tensors2[I]);
      else                           Tnon_commutes.push_back(Tensors2[I]);
    }
    if(Tnon_commutes.size()){
      cout << "Non-commutable tensor cannot be treated in this level .... " << endl;
      abort();
    }    
    Tensors2_(Tensors2);
  }

  // *********************************************************
  // 
  // *********************************************************
  SQFterm::SQFterm(const SQterm &obj)
    : numConsts_(obj.numConsts_),
      Consts_(obj.Consts_),
      Tensors1_(obj.Tensors1_),
      Tensors2_(obj.Tensors2_)
  {}

  // *********************************************************
  // 
  // *********************************************************
  vector<SQtensor> SQterm::get_tensors1() const
  { return Tensors1_; }

  // *********************************************************
  // 
  // *********************************************************
  vector<SQtensor*> SQterm::get_tensors1_ptr()
  {
    vector<SQtensor*> retval;
    for(size_t i = 0;i < Tensors1_.size();++i) retval.push_back(&(Tensors1_[i])); 
    return retval; 
  }

  // *********************************************************
  // 
  // *********************************************************
  vector<SQtensor> SQterm::get_tensors2() const
  { return Tensors2_; }

  // *********************************************************
  // 
  // *********************************************************
  vector<SQtensor*> SQterm::get_tensors2_ptr()
  {
    vector<SQtensor*> retval;
    for(size_t i = 0;i < Tensors2_.size();++i) retval.push_back(&(Tensors2_[i])); 
    return retval; 
  }

  // *********************************************************
  // 
  // *********************************************************
  vector<double> SQterm::get_numConsts() const
  { return numConsts_; }

  // *********************************************************
  // 
  // *********************************************************
  vector<string> SQterm::get_Consts() const
  { return Consts_; }

  // *********************************************************
  // 
  // *********************************************************
  void SQterm::set_tensors1(const vector<SQtensor> &Tensors)
  {
    vector<SQtensor> Tcommutes, Tnon_commutes;
    for(size_t I = 0;I < Tensors.size();++I){
      if(Tensors[I].isCommutable()) Tcommutes.push_back(Tensors[I]);
      else                          Tnon_commutes.push_back(Tensors[I]);
    }
    if(Tnon_commutes.size()){
      cout << "Non-commutable tensor cannnot be treat at this level .... " << endl;
      abort();
    }
    Tensors1_(Tensors);
  }

  // *********************************************************
  // 
  // *********************************************************
  void SQterm::set_tensors2(const vector<SQtensor> &Tensors)
  {
    vector<SQtensor> Tcommutes, Tnon_commutes;
    for(size_t I = 0;I < Tensors.size();++I){
      if(Tensors[I].isCommutable()) Tcommutes.push_back(Tensors[I]);
      else                          Tnon_commutes.push_back(Tensors[I]);
    }
    if(Tnon_commutes.size()){
      cout << "Non-commutable tensor cannnot be treat at this level .... " << endl;
      abort();
    }
    Tensors2_(Tensors);
  }

  // *********************************************************
  // 
  // *********************************************************
  SQterm SQterm::operator=(const SQterm &obj){
    isInCanonical_ = obj.isInCanonical_;
    numConst_      = obj.numConst_;     
    Consts_        = obj.Consts_;       
    Tensors_       = obj.Tensors_;      

    set_summedBody();

    return (*this);
  }

  // *********************************************************
  // 
  // *********************************************************
  ostream& operator <<(std::ostream &os, const SQterm &t)
  {
    os << boost::format("(%14.8f)") % t.get_numConst() << " ";

    vector<string> strs(t.get_Consts());
    for(size_t i = 0;i < strs.size();++i) os << strs[i] << " ";

    vector<SQtensor> tensors(t.get_tensors());    
    for(size_t i = 0;i < tensors.size();++i) os << tensors[i];
    return os;
  }

} // femto::

