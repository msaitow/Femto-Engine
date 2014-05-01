//
//  SQcontract.cc
//  
//
//  Created by Masaaki Saitow on 12/11/12.
//  Copyright (c) 2012 Masaaki Saitow. All rights reserved.
//

#include <cmath>
#include <vector>
#include <algorithm>
#include <boost/format.hpp>
#include <Femto.hpp>
#include <SQcontract.hpp>
#include <SQbinary.hpp>
#include <SQtensor.hpp>

using namespace std;

namespace Femto { namespace Core {

  // *********************************************************
  // 
  // *********************************************************
  SQcontract::SQcontract()
  {}

  // *********************************************************
  // 
  // *********************************************************
  SQcontract::SQcontract(const double numConst, const vector<string> Consts, 
		         const SQtensor &Ltensor, const vector<SQtensor> &Rtensors)
    : numConst_(numConst),
      Consts_(Consts),
      Ltensor_(Ltensor),
      Rtensors_(Rtensors),
      Lindices_(false),
      Rindices_(Rtensors.size(),false),
      RInnerInds_(Rtensors.size())
  {   
    if(is_RDM(Ltensor_.get_name()) || is_sfGen(Ltensor_.get_name()) || Ltensor_.get_name() == kDelta_name())
      { cout << "SQcontract: Ltensor cannot be either of sfGen, RDM or kDelta" << endl; abort(); }
    for(vector<SQtensor>::iterator t = Rtensors_.begin();t != Rtensors_.end();++t){
      if(is_sfGen(t->get_name())){
        cout << "Spin-free unitary group generator is detected: " << endl;
        abort();
      } // End if
    } // End t
    set_summedBody();
  }

  // *********************************************************
  // 
  // *********************************************************
  SQcontract::SQcontract(const SQcontract &obj)
    : Ltensor_(obj.Ltensor_),
      numConst_(obj.numConst_),
      Consts_(obj.Consts_),
      Lindices_(obj.Lindices_),
      Rindices_(obj.Rindices_),
      Rtensors_(obj.Rtensors_),
      LInnerInds_(obj.LInnerInds_),
      RInnerInds_(obj.RInnerInds_)
  { set_summedBody(); }

  // *********************************************************
  // 
  // *********************************************************
  SQcontract SQcontract::operator=(const SQcontract &obj){
    numConst_      = obj.numConst_;
    Consts_        = obj.Consts_;
    Ltensor_       = obj.Ltensor_;
    Rtensors_      = obj.Rtensors_;
    Lindices_      = obj.Lindices_;    
    Rindices_      = obj.Rindices_;    

    set_summedBody();
  }

  // *********************************************************
  // 
  // *********************************************************
  void SQcontract::set_summedBody()
  {
    // Firstly, assign all the &summedIndices_[:] to &tempIndices[:]
    vector<SQindex> tempIndices(summedIndices_);
    for(size_t i = 0;i < tempIndices.size();++i){
      for(size_t I = 0;I < Rtensors_.size();++I){
        vector<SQindex*> temp1(Rtensors_[I].get_indices());
        for(size_t j = 0;j < temp1.size();++j){
          if(*temp1[j]==tempIndices[i]) Rtensors_[I].put_indices(j, &tempIndices[i]);
	} // End j
      } // End I
      vector<SQindex*> temp1(Ltensor_.get_indices());
      for(size_t j = 0;j < temp1.size();++j){
        if(*temp1[j]==tempIndices[i]) Ltensor_.put_indices(j, &tempIndices[i]);
      } // End j
    } // End i
    summedIndices_.clear();

    // Search for the dummy indices .... 
    // Maybe it's nice for the sigma construction, but not for the diagonal preconditionor (2012/10/27)
    for(size_t I = 0;I < Rtensors_.size();++I){
      vector<SQindex*> temp1(Rtensors_[I].get_indices());
      vector<SQindex> indices;
      for(size_t i = 0;i < temp1.size();++i) indices.push_back(*temp1[i]);
      for(size_t j = 0;j < indices.size();++j){
        if(find(summedIndices_.begin(), summedIndices_.end(), indices[j])==summedIndices_.end()) {
          summedIndices_.push_back(*temp1[j]);
        } // End if
      } // End j
    } // End I
    vector<SQindex*> temp1(Ltensor_.get_indices());
    vector<SQindex> indices;
    for(size_t i = 0;i < temp1.size();++i) indices.push_back(*temp1[i]);
    for(size_t j = 0;j < indices.size();++j){
      if(find(summedIndices_.begin(), summedIndices_.end(), indices[j])==summedIndices_.end()) {
        summedIndices_.push_back(*temp1[j]);
      } // End if
    } // End j
    
    // Replace all the dummy indices, which are shared among all the terms, with those
    // copied in the summedIndices
    for(size_t i = 0;i < summedIndices_.size();++i){
      SQindex temp_i(summedIndices_[i]);
      for(size_t I = 0;I < Rtensors_.size();++I){
        vector<SQindex*> temp1 = Rtensors_[I].get_indices();
        vector<SQindex> indices;        
        for(size_t j = 0;j < temp1.size();++j) indices.push_back(*temp1[j]);
        for(size_t j = 0;j < indices.size();++j){
          if(indices[j] == temp_i) Rtensors_[I].put_indices(j, &summedIndices_[i]);
	}
      } // End I
      for(size_t num_i = 0;num_i < Ltensor_.get_indices().size();++num_i)
        if(*(Ltensor_.get_indices()[num_i]) == temp_i) Ltensor_.put_indices(num_i, &summedIndices_[i]);
    } // End i

    //masquerade();
  }

  // *********************************************************
  // 
  // *********************************************************
  void SQcontract::masquerade()
  {
     // Set the dummy names for all the dummy indices
     int Ccount(0);
     int Ocount(0);
     int Vcount(0);
     for(size_t i = 0;i < summedIndices_.size();++i){
       // In case of core
       if(summedIndices_[i].get_char() == 0 && summedIndices_[i].get_isSummed()){
         ++Ccount;
         ostringstream stm;
         stm << Ccount;
         summedIndices_[i].put_index("c" + stm.str());
       }
       else
       // In case of active
 	if(summedIndices_[i].get_char() == 1 && summedIndices_[i].get_isSummed()){
         ++Ocount;
         ostringstream stm;
         stm << Ocount;
         summedIndices_[i].put_index("o" + stm.str());
       }
       else
       // In case of virtual
 	if(summedIndices_[i].get_char() == 2 && summedIndices_[i].get_isSummed()){
         ++Vcount;
         ostringstream stm;
         stm << Vcount;
         summedIndices_[i].put_index("v" + stm.str());
       }     
     } // End i

  }

  // *********************************************************
  // 
  // *********************************************************
  vector<SQtensor> SQcontract::get_Rtensors() const
  { return Rtensors_; }

  // *********************************************************
  // 
  // *********************************************************
  SQtensor SQcontract::get_Ltensor() const
  { return Ltensor_; }

  // *********************************************************
  // 
  // *********************************************************
  vector<SQtensor*> SQcontract::get_Rtensors_ptr()
  {
    vector<SQtensor*> retval; retval.reserve(Rtensors_.size());
    for(size_t i = 0;i < Rtensors_.size();++i) retval.push_back(&(Rtensors_[i])); 
    return retval; 
  }

  // *********************************************************
  // 
  // *********************************************************
  SQtensor* SQcontract::get_Ltensor_ptr()
  { return &Ltensor_; }

  // *********************************************************
  // 
  // *********************************************************
  double SQcontract::get_numConst() const
  { return numConst_; }

  // *********************************************************
  // 
  // *********************************************************
  vector<string> SQcontract::get_Consts() const
  { return Consts_; }

  // *********************************************************
  // 
  // *********************************************************
  void SQcontract::set_Rtensors(const vector<SQtensor> &Tensors)
  {
    Rtensors_ = Tensors;
    for(size_t i = 0;i < Rtensors_.size();++i) Rindices_[i] = false;
    set_summedBody();    
  }

  //*-- --*// // *********************************************************
  //*-- --*// // 
  //*-- --*// // *********************************************************
  //*-- --*// vector<SQindex*> SQcontract::get_Lindices()
  //*-- --*// {
  //*-- --*//   vector<SQindex*> retval;
  //*-- --*//   // If InternalNames are not set, just return the internal indices.
  //*-- --*//   if(!InternalNames_.size()){
  //*-- --*//     for(size_t num_i = 0;num_i < summedIndices_.size();++num_i)
  //*-- --*//       if(!summedIndices_[num_i].get_isExt()) retval.push_back(&summedIndices_[num_i]);
  //*-- --*//   }
  //*-- --*//   // If InternalNames are given specifically, return only these indices.
  //*-- --*//   {
  //*-- --*//     for(size_t num_i = 0;num_i < summedIndices_.size();++num_i)
  //*-- --*//       if(find(InternalNames_.begin(), InternalNames_.end(), summedIndices_[num_i].get_index()) != InternalNames_.end()) 
  //*-- --*// 	  retval.push_back(&summedIndices_[num_i]);      
  //*-- --*//   }
  //*-- --*//   return retval;
  //*-- --*// }

  //*-- --*// // *********************************************************
  //*-- --*// // 
  //*-- --*// // *********************************************************
  //*-- --*// void SQcontract::set_Lindices(vector<string> &name_list)
  //*-- --*// { InternalNames_ = name_list; }
  //*-- --*// 
  //*-- --*// // *********************************************************
  //*-- --*// // 
  //*-- --*// // *********************************************************
  //*-- --*// void SQcontract::clear_Lindices()
  //*-- --*// { InternalNames_.clear(); }

  // *********************************************************
  // 
  // *********************************************************
  // NOTE :: 
  //   Basically, Ltensor_ can be sigma vector, preconditionor or the intermediate tensor.
  //   Lindices_ are intended to be used to represent how the Ltensor_ is treated in 
  //   the code generation step, which means, in case that Ltensor_ is an intermediate type
  //   tensor associated with the loading index of ERI, it's sufficient to use only the internal
  //   indices as the actual indices of DTensor, or symblock. But in the other case, in which, 
  //   Ltensor_ doesn't have the loading index, all the indices_ may have to be treated as the 
  //   actual indices of the DTensor representation of the Ltensor_. Lindices_ can stand for which types
  //   of the treatment is appropriate for LTensor_. In case of false (by default), only the internal indices
  //   are treated as the actual indices of DTensor representation. But if it's true, all the indices are used.  
  void SQcontract::set_Lindices(bool val)
  { Lindices_ = val; }

  // *********************************************************
  // 
  // *********************************************************
  // NOTE ::
  //   Meaning of Rindices_ is same to Lindices_ for LTensor_. For n-th member of Rtensors_, 
  //   Rindices_[n] stands for how the associated indices are treated in the code generating steps.
  void SQcontract::set_Rindices(int num, bool val)
  { Rindices_[num] = val; }

  // *********************************************************
  // 
  // *********************************************************
  vector<SQindex*> SQcontract::get_RInnerIndices(int num)
  {
    if     ( Rindices_[num] && !RInnerInds_[num].size())
      return Rtensors_[num].get_indices();
    else if( Rindices_[num] &&  RInnerInds_[num].size())
      return RInnerInds_[num];
    else if(!Rindices_[num] &&  RInnerInds_[num].size()) //NEW//
      return RInnerInds_[num];
    else {     
      vector<SQindex*> temp;
      for(size_t i = 0;i < Rtensors_[num].get_indices().size();++i)
	if(!Rtensors_[num].get_indices()[i]->get_isExt()) temp.push_back(Rtensors_[num].get_indices()[i]);
      return temp;
    } // End else
  }

  // *********************************************************
  // 
  // *********************************************************
  vector<SQindex*> SQcontract::get_LInnerIndices()
  {
    if     ( Lindices_ && !LInnerInds_.size())
      return Ltensor_.get_indices();
    else if( Lindices_ &&  LInnerInds_.size())
      return LInnerInds_;
    else if(!Lindices_ &&  LInnerInds_.size()) //NEW//
      return LInnerInds_;
    else {
      vector<SQindex*> temp;
      for(size_t i = 0;i < Ltensor_.get_indices().size();++i)
	if(!Ltensor_.get_indices()[i]->get_isExt()) temp.push_back(Ltensor_.get_indices()[i]);
      return temp;
    } // End else
  }

  // *********************************************************
  // 
  // *********************************************************
  vector<SQindex*> SQcontract::get_ROuterIndices(int num)
  {
    vector<SQindex*> tempInner(this->get_RInnerIndices(num));
    vector<SQindex*> tempAll(Rtensors_[num].get_indices());
    vector<SQindex*> retVal;
    for(auto i = tempAll.begin();i != tempAll.end();++i)
      if(find(tempInner.begin(), tempInner.end(), *i) == tempInner.end()) retVal.push_back(*i);
    return retVal;
  }

  // *********************************************************
  // 
  // *********************************************************
  vector<SQindex*> SQcontract::get_LOuterIndices()
  {
    vector<SQindex*> tempInner(this->get_LInnerIndices());
    vector<SQindex*> tempAll(Ltensor_.get_indices());
    vector<SQindex*> retVal;
    for(auto i = tempAll.begin();i != tempAll.end();++i)
      if(find(tempInner.begin(), tempInner.end(), *i) == tempInner.end()) retVal.push_back(*i);
    return retVal;
  }

  // *********************************************************
  // 
  // *********************************************************
  void SQcontract::set_ROuterIndices(int num, SQcont<SQindex> inds)
  {
    RInnerInds_[num].clear();
    vector<SQindex*> temp(Rtensors_[num].get_indices());
    for(auto t = temp.begin();t != temp.end();++t) if(inds.count(**t) == 0) RInnerInds_[num].push_back(*t);
//*//    //  old impl 
//*//    SQcont<SQindex*> oldOuters(this->get_ROuterIndices(num));
//*//    vector<SQindex*> temp(Rtensors_[num].get_indices());
//*//    RInnerInds_[num].clear();
//*//    for(auto i = oldOuters.begin();i != oldOuters.end();++i) if(inds.count(**i) == 0) inds <= **i;
//*//    for(auto t = temp.begin();t != temp.end();++t)
//*//      if(inds.count(**t) == 0) RInnerInds_[num].push_back(*t);
  }

  // *********************************************************
  // 
  // *********************************************************
  void SQcontract::clear_ROuterIndices(int num)
  { RInnerInds_[num].clear(); }

  // *********************************************************
  // 
  // *********************************************************
  void SQcontract::set_LOuterIndices(SQcont<SQindex> inds)
  {
    LInnerInds_.clear();
    vector<SQindex*> temp(Ltensor_.get_indices());
    for(auto i = temp.begin();i != temp.end();++i) if(inds.count(**i) == 0) LInnerInds_.push_back(*i);
//*//    //  old impl 
//*//    SQcont<SQindex*> oldOuters(this->get_LOuterIndices());
//*//    vector<SQindex*> temp(Ltensor_.get_indices());
//*//    LInnerInds_.clear();
//*//    for(auto i = oldOuters.begin();i != oldOuters.end();++i) if(inds.count(**i) == 0) inds <= **i;
//*//    for(auto t = temp.begin();t != temp.end();++t)
//*//      if(inds.count(**t) == 0) LInnerInds_.push_back(*t);
  }

  // *********************************************************
  // 
  // *********************************************************
  void SQcontract::clear_LOuterIndices()
  { LInnerInds_.clear(); }

  // *********************************************************
  // 
  // *********************************************************
  bool SQcontract::get_Lindices()
  { return Lindices_; } 

  // *********************************************************
  // 
  // *********************************************************
  vector<bool> SQcontract::get_Rindices()
  { return Rindices_; }  

  // *********************************************************
  // 
  // *********************************************************
  void SQcontract::set_Ltensor(const SQtensor &Tensor)
  {
    Ltensor_ = Tensor;
    set_summedBody();    
  }

  // *********************************************************
  // 
  // *********************************************************
  void SQcontract::set_numConst(const double num)
  { numConst_ = num; }

  // *********************************************************
  // 
  // *********************************************************
  void SQcontract::set_Consts(const vector<string> consts)
  { Consts_ = consts; }

  // *********************************************************
  // 
  // *********************************************************
  ostream& operator <<(std::ostream &os, const SQcontract &t)
  {
    os << t.get_Ltensor() << "+= ";
    os << boost::format("(%14.8f)") % t.get_numConst() << " ";

    vector<string> strs(t.get_Consts());
    for(size_t i = 0;i < strs.size();++i) os << strs[i] << " ";

    vector<SQtensor> tensors(t.get_Rtensors());    
    for(size_t i = 0;i < tensors.size();++i) os << tensors[i];
    return os;
  }

  // *********************************************************
  // 
  // *********************************************************
  vector<SQindex*> SQcontract::get_summedBody()
  {
    vector<SQindex*> retval;
    vector<SQindex>::iterator ind = summedIndices_.begin();
    for(;ind != summedIndices_.end();++ind) retval.push_back(&(*ind));
    return retval;
  }

  // *********************************************************
  // Contract Kronecker deltas
  // *********************************************************
  void SQcontract::contractkDeltas()
  {

    // *******************************************************************
    // * The case that Kronecker's delta has at least one dummy index is *
    // * not considered because this class is not intended for the deri- *
    // * vation of the many-body equation                                *
    // *******************************************************************
    for(vector<SQtensor>::iterator t = Rtensors_.begin(); t != Rtensors_.end();){
      // In case of delta_{p}^{p} ....
      if((t->get_name()==kDelta_name()) && (*(t->get_indices()[0]) == *(t->get_indices()[1])))
        t = Rtensors_.erase(t);
//*- -*//      else if(t->get_name()==kDelta_name() && t->get_indices()[0]->get_isSummed() || t->get_indices()[1]->get_isSummed()){
//*- -*//	cout << "SQcontract: kDelta with dummy index is detected. It's not the one to be handled at code-generation process" << endl;
//*- -*//	abort();
//*- -*//      } // End if
      // In case of both indices of kDelta are not dummy ....
      else if((t->get_name()==kDelta_name()) && !(t->get_indices()[0]->get_isSummed()) && !(t->get_indices()[1]->get_isSummed())){
        if(t->get_indices()[0]->get_char() != t->get_indices()[1]->get_char()){
	  numConst_ = 0; // This term shoud be zero
	  return;
	} // End if
	// If kdelta is composed if purely, non-dummy indices (not implemented in SQterm) ... 
	else{
	  SQindex* killed;
	  SQindex* killer;
	  if(*(t->get_indices()[0]) > *(t->get_indices()[1])) { killer = t->get_indices()[1]; killed = t->get_indices()[0]; }
	  else                                                { killer = t->get_indices()[0]; killed = t->get_indices()[1]; }
	  
	  // Kill the index to be killed 
	  vector<SQindex>::iterator killed_ptr = find(summedIndices_.begin(), summedIndices_.end(), *killed);
	  if(killed_ptr == summedIndices_.end()){
	    cout << "SQcontract: Algorithmic Error ..... " << endl;
	    abort();
	  } // End if
	  else{
	    *killed_ptr = *killer; 
	  } // End else
	  t = Rtensors_.erase(t);
	} // End else
	++t;
      } // End else
      // In case of either of indices in kDelta is dummy .... 
      else if((t->get_name()==kDelta_name()) && (t->get_indices()[0]->get_isSummed() || t->get_indices()[1]->get_isSummed())){
        if(t->get_indices()[0]->get_char() != t->get_indices()[1]->get_char()){
          numConst_ = 0;
	  return;
	} // End if
	else{
          SQindex *killed;
          SQindex *killer;
          if(t->get_indices()[1]->get_isSummed()) { killed = t->get_indices()[1]; killer = t->get_indices()[0];}
          else                                    { killed = t->get_indices()[0]; killer = t->get_indices()[1];}
 	  // Kill the index to be killed 
          vector<SQindex>::iterator killed_ptr = find(summedIndices_.begin(), summedIndices_.end(), *killed);
          if(killed_ptr == summedIndices_.end()){
            cout << "SQcontract: Algorithmic Error ..... " << endl;
            abort();
	  }
	  else{
            *killed_ptr = *killer; 
	  }
 	  t = Rtensors_.erase(t);
	} // End else
      } // End else if
      else ++t;

    } // End for

  }

  // *********************************************************
  // 
  // *********************************************************
  void SQcontract::print_summedBody()
  {
    cout << ">> summedBody <<" << endl;
    vector<SQindex>::iterator i = summedIndices_.begin();
    int count = 0;
    for(;i != summedIndices_.end();++i)
      cout << (boost::format("[%10d] ") % count++) << *i << "(" << &(*i) << ") " <<  (i->get_isExt() ? "E" : "I") << endl; 
  }

  // *********************************************************
  // Return whether two SQcontract are factorizable or not
  // *********************************************************
  bool isFactorizable(SQcontract &a, SQcontract &b)
  {
    // Compare constants
    vector<string> Const_a = a.get_Consts();
    vector<string> Const_b = b.get_Consts();
    if(Const_a != Const_b) return false;

    // Compare number of indices and names of (L)tensors
    if(a.get_Ltensor().get_indices().size() != b.get_Ltensor().get_indices().size() || 
       a.get_Ltensor().get_name() != b.get_Ltensor().get_name()) return false;

    // Compare names of (R)tensors and associated indices
    vector<SQtensor> ten_a = a.get_Rtensors();
    vector<SQtensor> ten_b = b.get_Rtensors();
    if(ten_a.size() != ten_b.size()) return false;
    sort(ten_a.begin(), ten_a.end());
    sort(ten_b.begin(), ten_b.end());
    for(size_t i = 0;i < ten_a.size();++i){
      if(ten_a[i].get_name()           != ten_b[i].get_name())           
        return false;
      if(ten_a[i].get_indices().size() != ten_b[i].get_indices().size())
        return false;
    } // End i

    // Count all the dummies in each orbital group
    vector<SQindex*> a_indices = a.get_summedBody();
    vector<SQindex*> c_ptr;
    vector<SQindex*> o_ptr;
    vector<SQindex*> v_ptr;
    for(vector<SQindex*>::iterator i = a_indices.begin();i != a_indices.end();++i){
      if     ((*i)->get_char()==(char_state)0 && (*i)->get_isSummed()) c_ptr.push_back(*i); 
      else if((*i)->get_char()==(char_state)1 && (*i)->get_isSummed()) o_ptr.push_back(*i); 
      else if((*i)->get_char()==(char_state)2 && (*i)->get_isSummed()) v_ptr.push_back(*i);
    } // End i

    // Before comparing the tensors, masquerade the other one!
    //b.masquerade();

    // Permute pairs of indices in each orbital group
    IIvector c_perms = makePermutations((int)c_ptr.size());
    IIvector o_perms = makePermutations((int)o_ptr.size());
    IIvector v_perms = makePermutations((int)v_ptr.size());
    int c_max = (c_perms.size() ? c_perms.size() : 1);
    int o_max = (o_perms.size() ? o_perms.size() : 1);
    int v_max = (v_perms.size() ? v_perms.size() : 1);
    for(int ic = 0;ic < c_max;++ic){
      for(int io = 0;io < o_max;++io){
        for(int iv = 0;iv < v_max;++iv){
          // Replace each name of index
          if(c_perms.size())
            //cout << c_perms.size() << endl;
            for(size_t i = 0;i < c_ptr.size();++i){
              stringstream num; num << c_perms[ic][i]+1;
              c_ptr[i]->put_index("c"+num.str());
            }
	  if(o_perms.size())
            for(size_t i = 0;i < o_ptr.size();++i){
              stringstream num; num << o_perms[io][i]+1;
              o_ptr[i]->put_index("o"+num.str());
            }
          if(v_perms.size())
            for(size_t i = 0;i < v_ptr.size();++i){
              stringstream num; num << v_perms[iv][i]+1;
              v_ptr[i]->put_index("v"+num.str()); 
            }
	  vector<SQtensor> a_tensors(a.get_Rtensors()); sort(a_tensors.begin(), a_tensors.end());
	  vector<SQtensor> b_tensors(b.get_Rtensors()); sort(b_tensors.begin(), b_tensors.end());
	  if(a_tensors == b_tensors) return true;
        } // End iv
      } // End io
    } // End ic

    return false;
  }

  // *********************************************************
  // 
  // *********************************************************
  SQbinary::SQbinary()
  {}

  // *********************************************************
  // 
  // *********************************************************
  SQbinary::SQbinary(const double numConst, const vector<string> Consts, 
		     const SQtensor &Ltensor, const vector<SQtensor> &Rtensors)
//     : numConst_(numConst),
//       Consts_(Consts),
//       Ltensor_(Ltensor),
//       Rtensors_(Rtensors),
//       Lindices_(false),
//       Rindices_(2,false)
  {

    numConst_ = numConst;
    Consts_   = Consts;
    Ltensor_  = Ltensor;
    Rtensors_ = Rtensors;
    Lindices_ = false;
    Rindices_;
    Rindices_.push_back(false);
    Rindices_.push_back(false);
    LInnerInds_;
    RInnerInds_.resize(2);

    if(is_RDM(Ltensor_.get_name()) || is_sfGen(Ltensor_.get_name()) || Ltensor_.get_name() == kDelta_name())
      { cout << "SQcontract: Ltensor cannot be either of sfGen, RDM or kDelta" << endl; abort(); }
   
    int num_non_kdeltas = 0;
    for(vector<SQtensor>::iterator t = Rtensors_.begin();t != Rtensors_.end();++t)
      if(t->get_name() != kDelta_name()) ++num_non_kdeltas;
    if(num_non_kdeltas >= 3){
      cout << "Number of non Kronecker's delta type tensors has to be less than 3." << endl;
      for(vector<SQtensor>::iterator t = Rtensors_.begin();t != Rtensors_.end();++t) cout << *t;
      cout << endl;
      abort();
    } // End if
    for(vector<SQtensor>::iterator t = Rtensors_.begin();t != Rtensors_.end();++t){
      if(is_sfGen(t->get_name())){
        cout << "Spin-free unitary group generator is detected: " << endl;
        abort();
      } // End if
    } // End t
    set_summedBody();
  }

}} // Femto::

