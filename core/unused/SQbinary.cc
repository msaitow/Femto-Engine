//
//  SQbinary.cc
//  
//
//  Created by Masaaki Saitow on 12/09/09.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <cmath>
#include <vector>
#include <algorithm>
#include <boost/format.hpp>
#include <femto.hpp>
#include <SQbinary.hpp>
#include <SQtensor.hpp>

namespace femto {

  // *********************************************************
  // 
  // *********************************************************
  SQbinary::SQbinary()
  {}

  // *********************************************************
  // 
  // *********************************************************
  SQbinary::SQbinary(const double numConst, const vector<string> Consts, 
		     const SQtensor &Ltensor, const vector<SQtensor> Rtensors)
    : numConst_(numConst),
      Consts_(Consts),
      Ltensor_(Ltensor),
      Rtensors_(Rtensors),
      Lindices_(false),
      Rindices_(2,false)
  {
   
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

  // *********************************************************
  // 
  // *********************************************************
  SQbinary::SQbinary(const SQbinary &obj)
    : Ltensor_(obj.Ltensor_),
      numConst_(obj.numConst_),
      Consts_(obj.Consts_),
      //InternalNames_(obj.InternalNames_),
      Lindices_(obj.Lindices_),
      Rindices_(obj.Rindices_),
      Rtensors_(obj.Rtensors_)
  { set_summedBody(); }

  // *********************************************************
  // 
  // *********************************************************
  SQbinary SQbinary::operator=(const SQbinary &obj){
    numConst_      = obj.numConst_;
    Consts_        = obj.Consts_;
    //InternalNames_ = obj.InternalNames_;
    Ltensor_       = obj.Ltensor_;
    Rtensors_      = obj.Rtensors_;
    Lindices_      = obj.Lindices_;    
    Rindices_      = obj.Rindices_;    

    set_summedBody();
  }

  // *********************************************************
  // 
  // *********************************************************
  void SQbinary::set_summedBody()
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
    // Maybe nice the sigma construction, but not the diagonal preconditionor (2012/10/27)
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

//*-- OLD --*//   // *********************************************************
//*-- OLD --*//   // 
//*-- OLD --*//   // *********************************************************
//*-- OLD --*//   void SQbinary::set_summedBody()
//*-- OLD --*//   {
//*-- OLD --*//     // Firstly, assign all the &summedIndices_[:] to &tempIndices[:]
//*-- OLD --*//     vector<SQindex> tempIndices(summedIndices_);
//*-- OLD --*//     for(size_t i = 0;i < tempIndices.size();++i){
//*-- OLD --*//       for(size_t I = 0;I < Rtensors_.size();++I){
//*-- OLD --*//         vector<SQindex*> temp1(Rtensors_[I].get_indices());
//*-- OLD --*//         for(size_t j = 0;j < temp1.size();++j){
//*-- OLD --*//           if(*temp1[j]==tempIndices[i]) Rtensors_[I].put_indices(j, &tempIndices[i]);
//*-- OLD --*// 	} // End j
//*-- OLD --*//       } // End I
//*-- OLD --*//       vector<SQindex*> temp1(Ltensor_.get_indices());
//*-- OLD --*//       for(size_t j = 0;j < temp1.size();++j){
//*-- OLD --*//         if(*temp1[j]==tempIndices[i]) Ltensor_.put_indices(j, &tempIndices[i]);
//*-- OLD --*//       } // End j
//*-- OLD --*//     } // End i
//*-- OLD --*//     summedIndices_.clear();
//*-- OLD --*// 
//*-- OLD --*//     // Search for the dummy indices .... 
//*-- OLD --*//     // Maybe nice the sigma construction, but not the diagonal preconditionor (2012/10/27)
//*-- OLD --*//     for(size_t I = 0;I < Rtensors_.size();++I){
//*-- OLD --*//       vector<SQindex*> temp1(Rtensors_[I].get_indices());
//*-- OLD --*//       vector<SQindex> indices;
//*-- OLD --*//       for(size_t i = 0;i < temp1.size();++i) indices.push_back(*temp1[i]);
//*-- OLD --*//       for(size_t j = 0;j < indices.size();++j){
//*-- OLD --*//         if(find(summedIndices_.begin(), summedIndices_.end(), indices[j])==summedIndices_.end()) {
//*-- OLD --*//           summedIndices_.push_back(*temp1[j]);
//*-- OLD --*//         } // End if
//*-- OLD --*//       } // End j
//*-- OLD --*//     } // End I
//*-- OLD --*//     vector<SQindex*> temp1(Ltensor_.get_indices());
//*-- OLD --*//     vector<SQindex> indices;
//*-- OLD --*//     for(size_t i = 0;i < temp1.size();++i) indices.push_back(*temp1[i]);
//*-- OLD --*//     for(size_t j = 0;j < indices.size();++j){
//*-- OLD --*//       if(find(summedIndices_.begin(), summedIndices_.end(), indices[j])==summedIndices_.end()) {
//*-- OLD --*//         summedIndices_.push_back(*temp1[j]);
//*-- OLD --*//       } // End if
//*-- OLD --*//     } // End j
//*-- OLD --*//     
//*-- OLD --*//     // Replace all the dummy indices, which are shared among all the terms, with those
//*-- OLD --*//     // copied in the summedIndices
//*-- OLD --*//     for(size_t i = 0;i < summedIndices_.size();++i){
//*-- OLD --*//       SQindex temp_i(summedIndices_[i]);
//*-- OLD --*//       for(size_t I = 0;I < Rtensors_.size();++I){
//*-- OLD --*//         vector<SQindex*> temp1 = Rtensors_[I].get_indices();
//*-- OLD --*//         vector<SQindex> indices;        
//*-- OLD --*//         for(size_t j = 0;j < temp1.size();++j) indices.push_back(*temp1[j]);
//*-- OLD --*//         for(size_t j = 0;j < indices.size();++j){
//*-- OLD --*//           if(indices[j] == temp_i) Rtensors_[I].put_indices(j, &summedIndices_[i]);
//*-- OLD --*// 	}
//*-- OLD --*//       } // End I
//*-- OLD --*//       for(size_t num_i = 0;num_i < Ltensor_.get_indices().size();++num_i)
//*-- OLD --*//         if(*(Ltensor_.get_indices()[num_i]) == temp_i) Ltensor_.put_indices(num_i, &summedIndices_[i]);
//*-- OLD --*//     } // End i
//*-- OLD --*// 
//*-- OLD --*//     //masquerade();
//*-- OLD --*//   }

  // *********************************************************
  // 
  // *********************************************************
  void SQbinary::masquerade()
  {
     // Set the dummy names for all the dummy indices
     int Ccount = 0;
     int Ocount = 0;
     int Vcount = 0;
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
  vector<SQtensor> SQbinary::get_Rtensors() const
  { return Rtensors_; }

  // *********************************************************
  // 
  // *********************************************************
  SQtensor SQbinary::get_Ltensor() const
  { return Ltensor_; }

  // *********************************************************
  // 
  // *********************************************************
  vector<SQtensor*> SQbinary::get_Rtensors_ptr()
  {
    vector<SQtensor*> retval; retval.reserve(Rtensors_.size());
    for(size_t i = 0;i < Rtensors_.size();++i) retval.push_back(&(Rtensors_[i])); 
    return retval; 
  }

  // *********************************************************
  // 
  // *********************************************************
  SQtensor* SQbinary::get_Ltensor_ptr()
  { return &Ltensor_; }

  // *********************************************************
  // 
  // *********************************************************
  double SQbinary::get_numConst() const
  { return numConst_; }

  // *********************************************************
  // 
  // *********************************************************
  vector<string> SQbinary::get_Consts() const
  { return Consts_; }

  // *********************************************************
  // 
  // *********************************************************
  void SQbinary::set_Rtensors(const vector<SQtensor> &Tensors)
  {
    Rtensors_ = Tensors;
    //clear_Lindices();
    Rindices_[0] = false;
    Rindices_[1] = false;
    set_summedBody();    
  }

  //*-- --*// // *********************************************************
  //*-- --*// // 
  //*-- --*// // *********************************************************
  //*-- --*// vector<SQindex*> SQbinary::get_Lindices()
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
  //*-- --*// void SQbinary::set_Lindices(vector<string> &name_list)
  //*-- --*// { InternalNames_ = name_list; }
  //*-- --*// 
  //*-- --*// // *********************************************************
  //*-- --*// // 
  //*-- --*// // *********************************************************
  //*-- --*// void SQbinary::clear_Lindices()
  //*-- --*// { InternalNames_.clear(); }

  // *********************************************************
  // 
  // *********************************************************
  void SQbinary::set_Lindices(bool val)
  { Lindices_ = val; }

  // *********************************************************
  // 
  // *********************************************************
  void SQbinary::set_Rindices(int num, bool val)
  { Rindices_[num] = val; }

  // *********************************************************
  // 
  // *********************************************************
  bool SQbinary::get_Lindices()
  { return Lindices_; } 

  // *********************************************************
  // 
  // *********************************************************
  vector<bool> SQbinary::get_Rindices()
  { return Rindices_; }  

  // *********************************************************
  // 
  // *********************************************************
  void SQbinary::set_Ltensor(const SQtensor &Tensor)
  {
    Ltensor_ = Tensor;
    set_summedBody();    
  }

  // *********************************************************
  // 
  // *********************************************************
  void SQbinary::set_numConst(const double num)
  { numConst_ = num; }

  // *********************************************************
  // 
  // *********************************************************
  void SQbinary::set_Consts(const vector<string> consts)
  { Consts_ = consts; }

  // *********************************************************
  // 
  // *********************************************************
  ostream& operator <<(std::ostream &os, const SQbinary &t)
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
  vector<SQindex*> SQbinary::get_summedBody()
  {
    vector<SQindex*> retval;
    vector<SQindex>::iterator ind = summedIndices_.begin();
    for(;ind != summedIndices_.end();++ind) retval.push_back(&(*ind));
    return retval;
  }

  // *********************************************************
  // Contract Kronecker deltas
  // *********************************************************
  void SQbinary::contractkDeltas()
  {

    // *******************************************************************
    // * In case that LTensor is a Kronecker's delta is *NOT* considered *
    // *******************************************************************
    for(vector<SQtensor>::iterator t = Rtensors_.begin(); t != Rtensors_.end();){
      // In case of delta_{p}^{p} ....
      if((t->get_name()==kDelta_name()) && (*(t->get_indices()[0]) == *(t->get_indices()[1])))
        t = Rtensors_.erase(t);
      // In case of both indices of kDelta are dummy ....
      else if((t->get_name()==kDelta_name()) && !(t->get_indices()[0]->get_isSummed()) && !(t->get_indices()[1]->get_isSummed())){
        if(t->get_indices()[0]->get_char() != t->get_indices()[1]->get_char()){
	  numConst_ = 0; // This term shoud be zero
	  return;
	} // End if
        // If kdelta is composed if purely, non-dummy indices (not implemented in SQterm) ... 
        else{
          SQindex* killed;
          SQindex* killer;
          if(*(t->get_indices()[0]) > *(t->get_indices()[1])) { killed = t->get_indices()[1]; killer = t->get_indices()[0]; }
          else                                                { killed = t->get_indices()[0]; killer = t->get_indices()[1]; }

 	  // Kill the index to be killed 
          vector<SQindex>::iterator killed_ptr = find(summedIndices_.begin(), summedIndices_.end(), *killed);
          if(killed_ptr == summedIndices_.end()){
            cout << "Algorithmic Error ..... " << endl;
            abort();
	  }
	  else{
            *killed_ptr = *killer; 
	  }
 	  t = Rtensors_.erase(t);

	} // End else
	++t;
      } // End else if
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
            cout << "Algorithmic Error ..... " << endl;
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
  void SQbinary::print_summedBody()
  {
    cout << ">> summedBody <<" << endl;
    vector<SQindex>::iterator i = summedIndices_.begin();
    int count = 0;
    for(;i != summedIndices_.end();++i)
      cout << (boost::format("[%10d] ") % count++) << *i << "(" << &(*i) << ") " <<  (i->get_isExt() ? "E" : "I") << endl; 
  }

  // *********************************************************
  // Return whether two SQbinary are factorizable or not
  // *********************************************************
  bool isFactorizable(SQbinary &a, SQbinary &b)
  {
    // Compare constants
    vector<string> Const_a = a.get_Consts();
    vector<string> Const_b = b.get_Consts();
    if(Const_a != Const_b) return false;

    // Compare number of indices of (L)tensors
    if(a.get_Ltensor().get_indices().size() != b.get_Ltensor().get_indices().size()) return false;

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


} // femto::

