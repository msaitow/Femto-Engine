//
//  SQtensor.cc
//  
//
//  Created by Masaaki Saitow on 12/03/26.
//  Copyright (c) 2012 Masaaki Saitow. All rights reserved.
//

#include <cassert>
#include <algorithm>
#include <boost/format.hpp>
#include <Femto.hpp>
#include <SQtensor.hpp>
#include <SQindex.hpp>

// Now permutational symmetry is *partially* available.

using namespace std;

namespace Femto { namespace Core {

  // Some auxiliary function
  struct PtrGreat { bool operator()(const SQindex* i, const SQindex* j) const 
    { return (*i) > (*j); } };

  // *********************************************************
  // 
  // *********************************************************
  SQtensor::SQtensor()
  {}

  // *********************************************************
  // 
  // *********************************************************
  SQtensor::SQtensor(const SQtensor &obj)
    : name_(obj.name_),
      //sign_(obj.sign_),
      indices_(obj.indices_),
      isCommutable_(obj.isCommutable_),
      not_(obj.not_)
  { 
    symm_         = obj.symm_;
    permutations_ = obj.permutations_; 
    factors_      = obj.factors_;
  }

  // *********************************************************
  // 
  // *********************************************************
  SQtensor::SQtensor(const string &name, const vector<SQindex*> indices, 
                     const Symmetry &symm, const bool isCommutable, const notation nota)
    : name_(name),
      indices_(indices),
      //sign_(sign),
      isCommutable_(isCommutable)      ,
      not_(nota)
      //  { symm_.first = symm.first; symm_.second = symm.second; }
  { 

    //
    // Symmetry = < < .... >, < .... >, .... >, < ..... >
    //              --------------------------  ---------
    //                       Patterns            Factors
    //   

    // assert doesn't work. So why????

    symm_ = symm;
    int len_s = symm.first.size();
    if(len_s != symm.second.size()){
      cout << "Inconsisitency occured in symm and factors" << endl;
      abort();   
    }
    size_t len_t = indices.size();
    vector<int> temp_i;
    for(size_t i = 0;i < len_t;++i) temp_i.push_back(int(i));

    if(not_ != (notation)0 && not_ != (notation)1){
      cout << "Notation has to be given in either Dirac or Mulliken." << endl;
      abort();
    }

    IIvector::iterator it = symm_.first.begin();
    while(it != symm_.first.end()){
      size_t len_s = (*it).size();
      if (len_t != len_s && len_s != 0){
        cout << "Inconsistency occured in length of symm and indices" << endl;
        abort();
      }
      Ivector sorted(*(it));
      sort(sorted.begin(), sorted.end());
      if(sorted != temp_i){
        cout << "Symm  has to be given as a vector of contiguous numbers" << endl;
        abort();
      }

      ++it;
    }

    makeAllPerms();
    //sortIndices(); //*Not_ALWAYS_GOOD*
  }

  // *********************************************************
  // 
  // *********************************************************
  string SQtensor::get_name() const
  { return name_; }

  // *********************************************************
  // 
  // *********************************************************
  string SQtensor::get_notation_str() const
  { return (not_ ? "Dirac" : "Mulliken"); }

  // *********************************************************
  // 
  // *********************************************************
  notation SQtensor::get_notation() const
  { return not_; }

  // *********************************************************
  // 
  // *********************************************************
  bool SQtensor::isCommutable() const
  { return isCommutable_; }

  // *********************************************************
  // 
  // *********************************************************
  SQtensor& SQtensor::operator=(const SQtensor &obj){
    name_         = obj.name_;
    isCommutable_ = obj.isCommutable_;
    indices_      = obj.indices_;
    symm_         = obj.symm_;
    permutations_ = obj.permutations_;
    factors_      = obj.factors_;
    not_          = obj.not_;
    //sign_         = obj.sign_;
    
    return *this;
  }

  // *********************************************************
  // Added 2012/10/22
  // *********************************************************
  bool SQtensor::operator==(const SQtensor &obj) const{
    
    if(name_         != obj.name_)         return false;
    if(isCommutable_ != obj.isCommutable_) return false;
    //    if(symm_         != obj.symm_)         return false;
    if(indices_.size() != obj.indices_.size()) return false;
    if(permutations_ != obj.permutations_) return false;    
    if(factors_      != obj.factors_)      return false;
    if(not_          != obj.not_)          return false;
    //if(sign_         != obj.sign_)         return false;

    // Compare the non-dummy indices first
    vector<SQindex*> Is; Is.reserve(indices_.size()); 
    vector<SQindex*> Js; Js.reserve(obj.indices_.size());
    for(size_t num_i = 0;num_i < indices_.size();++num_i) 
      if(!indices_[num_i]->get_isSummed()) Is.push_back(indices_[num_i]); 
    for(size_t num_j = 0;num_j < obj.indices_.size();++num_j) 
      if(!obj.indices_[num_j]->get_isSummed()) Js.push_back(obj.indices_[num_j]); 

    if(Is.size() != Js.size()) return false;

    sort(Is.begin(), Is.end(), PtrGreat());
    sort(Js.begin(), Js.end(), PtrGreat());
    for(size_t num_i = 0;num_i < Is.size();++num_i)
      if(*Is[num_i] != *Js[num_i]) return false;

    // If matched, then compare all the possibilities
    // even though this algorithm may be improved
    bool matched = false;
    for(size_t i = 0;i < permutations_.size();++i){
      bool allOK = true;
      for(size_t k = 0;k < indices_.size();++k){
        if(!(*(indices_[permutations_[i][k]]) == *(obj.indices_[k]))) 
          { allOK = false; break; }
      }
      if(allOK) { matched = true; break;}
    } // End i

    return matched;
  }

//*OLD_IMPL*   // *********************************************************
//*OLD_IMPL*   // 
//*OLD_IMPL*   // *********************************************************
//*OLD_IMPL*   bool SQtensor::operator==(const SQtensor &obj) const{
//*OLD_IMPL*     
//*OLD_IMPL*     if(name_         != obj.name_)         return false;
//*OLD_IMPL*     if(isCommutable_ != obj.isCommutable_) return false;
//*OLD_IMPL*     //    if(symm_         != obj.symm_)         return false;
//*OLD_IMPL*     if(indices_.size() != obj.indices_.size()) return false;
//*OLD_IMPL*     if(permutations_ != obj.permutations_) return false;    
//*OLD_IMPL*     if(factors_      != obj.factors_)      return false;
//*OLD_IMPL*     if(not_          != obj.not_)          return false;
//*OLD_IMPL*     //if(sign_         != obj.sign_)         return false;
//*OLD_IMPL* 
//*OLD_IMPL*     bool matched = false;
//*OLD_IMPL*     for(size_t i = 0;i < permutations_.size();++i){
//*OLD_IMPL*       bool allOK = true;
//*OLD_IMPL*       for(size_t k = 0;k < indices_.size();++k){
//*OLD_IMPL*         if(!(*(indices_[permutations_[i][k]]) == *(obj.indices_[k]))) 
//*OLD_IMPL*           { allOK = false; break; }
//*OLD_IMPL*       }
//*OLD_IMPL*       if(allOK) { matched = true; break;}
//*OLD_IMPL*     } // End i
//*OLD_IMPL* 
//*OLD_IMPL*     return matched;
//*OLD_IMPL*   }

  // *********************************************************
  //  *NOTE* Make sure this order is not always unique if 
  //         there are several similar tensors.
  // *********************************************************
  bool SQtensor::operator<(const SQtensor &obj) const
  {
    bool retval;
    if      (name_  < obj.name_) retval = true;
    else if (name_  > obj.name_) retval = false;
    else if (name_ == obj.name_){
      //if     (sign_ < obj.sign_) retval = true;
      //else if(sign_ > obj.sign_) retval = false;
      //else if(sign_ == obj.sign_){
      if     (not_  < obj.not_) retval = true;
      else if(not_  > obj.not_) retval = false;
      else if(not_ == obj.not_){ 
        if     (isCommutable_  < obj.isCommutable_) retval = true;
        else if(isCommutable_  > obj.isCommutable_) retval = false;
        else if(isCommutable_ == obj.isCommutable_){
          if      (permutations_  < obj.permutations_) retval = true;
          else if (permutations_  > obj.permutations_) retval = false;
          else if (permutations_ == obj.permutations_){
            if      (factors_  < obj.factors_) retval = true;
            else if (factors_  > obj.factors_) retval = false;
            else if (factors_ == obj.factors_){
              if      (indices_.size()  < obj.indices_.size()) retval = true;
              else if (indices_.size()  > obj.indices_.size()) retval = false;
              else if (indices_.size() == obj.indices_.size()){

                // Compare the number of                 
                vector<SQindex*> my_indices;  my_indices.reserve(indices_.size());               
		for(auto i = indices_.begin();i != indices_.end();++i) if(!(*i)->get_isSummed()) my_indices.push_back(*i);               
		sort(my_indices.begin(), my_indices.end());
                vector<SQindex*> obj_indices; obj_indices.reserve(indices_.size());               
		for(auto i = obj.indices_.begin();i != obj.indices_.end();++i) if(!(*i)->get_isSummed()) obj_indices.push_back(*i);               
		sort(obj_indices.begin(), obj_indices.end());

		if     (my_indices  < obj_indices) retval = true;
		else if(my_indices  > obj_indices) retval = false;
		else if(my_indices == obj_indices){

		  // Count the scores
		  int ncore_(0), obj_ncore_(0);
		  int nact_ (0), obj_nact_ (0);
		  int nvir_ (0), obj_nvir_ (0);
		  for(size_t i = 0;i < indices_.size();++i){
		    if      (    indices_[i]->get_char() == 0) ++ncore_;
		    else if (    indices_[i]->get_char() == 1) ++nact_;
		    else if (    indices_[i]->get_char() == 2) ++nvir_;
		    if      (obj.indices_[i]->get_char() == 0) ++obj_ncore_;
		    else if (obj.indices_[i]->get_char() == 1) ++obj_nact_;
		    else if (obj.indices_[i]->get_char() == 2) ++obj_nvir_;
		  } // End i
		  
		  if      (ncore_  < obj_ncore_) retval = true;
		  else if (ncore_  > obj_ncore_) retval = false;
		  else if (ncore_ == obj_ncore_) {
		    if      (nact_   < obj_nact_ ) retval = true;
		    else if (nact_   > obj_nact_ ) retval = false;
		    else if (nact_  == obj_nact_ ) {
		      if      (nvir_   < obj_nvir_ ) retval = true;
		      else if (nvir_   > obj_nvir_ ) retval = false;
		      else if (nvir_  == obj_nvir_ ) {
			// Count number of non-dummy indices
			int count_ = 0, obj_count_ = 0;
			for(size_t i = 0;i < indices_.size();++i){
			  if(!    indices_[i]->get_isSummed()) ++count_;
			  if(!obj.indices_[i]->get_isSummed()) ++obj_count_;
			} // End i
			if      (count_  < obj_count_) retval = true;
			else if (count_  > obj_count_) retval = false;
			else if (count_ == obj_count_) retval = false;
		      } // End if(ncore_)
		    } // End if(nact_ )
		  } // End if(nvirt_)
		} // End if(my_indices)
  	      } // End if(indices_.size())
  	    } // End if(factors_)
  	  } // End if(permutations_)
        } // End if(isCommutable_)
      } // End if(not_)
	//}
    } // End if(name_)

    return retval;    
  }

//*OBSOLETE*   // *********************************************************
//*OBSOLETE*   //  *NOTE* Make sure this order is not always unique if 
//*OBSOLETE*   //         there are several similar tensors.
//*OBSOLETE*   // *********************************************************
//*OBSOLETE*   bool SQtensor::operator<(const SQtensor &obj) const
//*OBSOLETE*   {
//*OBSOLETE*     bool retval;
//*OBSOLETE*     if      (name_  < obj.name_) retval = true;
//*OBSOLETE*     else if (name_  > obj.name_) retval = false;
//*OBSOLETE*     else if (name_ == obj.name_){
//*OBSOLETE*       //if     (sign_ < obj.sign_) retval = true;
//*OBSOLETE*       //else if(sign_ > obj.sign_) retval = false;
//*OBSOLETE*       //else if(sign_ == obj.sign_){
//*OBSOLETE*       if     (not_  < obj.not_) retval = true;
//*OBSOLETE*       else if(not_  > obj.not_) retval = false;
//*OBSOLETE*       else if(not_ == obj.not_){ 
//*OBSOLETE*         if     (isCommutable_  < obj.isCommutable_) retval = true;
//*OBSOLETE*         else if(isCommutable_  > obj.isCommutable_) retval = false;
//*OBSOLETE*         else if(isCommutable_ == obj.isCommutable_){
//*OBSOLETE*           if      (permutations_  < obj.permutations_) retval = true;
//*OBSOLETE*           else if (permutations_  > obj.permutations_) retval = false;
//*OBSOLETE*           else if (permutations_ == obj.permutations_){
//*OBSOLETE*             if      (factors_  < obj.factors_) retval = true;
//*OBSOLETE*             else if (factors_  > obj.factors_) retval = false;
//*OBSOLETE*             else if (factors_ == obj.factors_){
//*OBSOLETE*               if      (indices_.size()  < obj.indices_.size()) retval = true;
//*OBSOLETE*               else if (indices_.size()  > obj.indices_.size()) retval = false;
//*OBSOLETE*               else if (indices_.size() == obj.indices_.size()){
//*OBSOLETE*                 // Count the scores
//*OBSOLETE*                 int ncore_ = 0, obj_ncore_ = 0;
//*OBSOLETE*                 int nact_  = 0, obj_nact_  = 0;
//*OBSOLETE*                 int nvir_  = 0, obj_nvir_  = 0;
//*OBSOLETE*                 for(size_t i = 0;i < indices_.size();++i){
//*OBSOLETE*                   if      (    indices_[i]->get_char() == 0) ++ncore_;
//*OBSOLETE*                   else if (    indices_[i]->get_char() == 1) ++nact_;
//*OBSOLETE*                   else if (    indices_[i]->get_char() == 2) ++nvir_;
//*OBSOLETE*                   if      (obj.indices_[i]->get_char() == 0) ++obj_ncore_;
//*OBSOLETE*                   else if (obj.indices_[i]->get_char() == 1) ++obj_nact_;
//*OBSOLETE*                   else if (obj.indices_[i]->get_char() == 2) ++obj_nvir_;
//*OBSOLETE*                 } // End i
//*OBSOLETE*                 int score_     =     nact_ + 2 *     nvir_;
//*OBSOLETE*                 int obj_score_ = obj_nact_ + 2 * obj_nvir_;
//*OBSOLETE*   
//*OBSOLETE*                 if      (score_  < obj_score_) retval = true;
//*OBSOLETE*                 else if (score_  > obj_score_) retval = false;
//*OBSOLETE*                 else if (score_ == obj_score_) {
//*OBSOLETE*                   // Count number of non-dummy indices
//*OBSOLETE*                   int count_ = 0, obj_count_ = 0;
//*OBSOLETE*                   for(size_t i = 0;i < indices_.size();++i){
//*OBSOLETE*                     if(!    indices_[i]->get_isSummed()) ++count_;
//*OBSOLETE*                     if(!obj.indices_[i]->get_isSummed()) ++obj_count_;
//*OBSOLETE*     	          } // End i
//*OBSOLETE*                   if      (count_  < obj_count_) retval = true;
//*OBSOLETE*                   else if (count_  > obj_count_) retval = false;
//*OBSOLETE*                   else if (count_ == obj_count_) retval = false;
//*OBSOLETE*   	        } // End if(score_)
//*OBSOLETE*   	      } // End if(indices_.size())
//*OBSOLETE*   	    } // End if(factors_)
//*OBSOLETE*   	  } // End if(permutations_)
//*OBSOLETE*         } // End if(isCommutable_)
//*OBSOLETE*       } // End if(not_)
//*OBSOLETE* 	//}
//*OBSOLETE*     } // End if(name_)
//*OBSOLETE* 
//*OBSOLETE*     return retval;    
//*OBSOLETE*   }

  // *********************************************************
  // 
  // *********************************************************
  bool SQtensor::operator>(const SQtensor &obj) const
  {
    bool retval = true;
    if(*this == obj) retval = false;
    else if(*this<obj) retval = false;
    return retval;
  }

  // *********************************************************
  // 
  // *********************************************************
  vector<SQindex*> SQtensor::get_indices() const
  { return indices_; }

  // *********************************************************
  // 
  // *********************************************************
  IIvector SQtensor::get_perms() const
  { return permutations_; }

  // *********************************************************
  // 
  // *********************************************************
  void SQtensor::put_indices(size_t i, SQindex *v)
  { 
    vector<SQindex*>::iterator it = indices_.begin(); it += i;
    it = indices_.erase(it);
    indices_.insert(it, v);
    //sortIndices(); //*TEST* 
  }

  // *********************************************************
  // 
  // *********************************************************
  ostream& operator <<(std::ostream &os, const SQtensor &s)
  {
    std::string IndLabel;
    vector<SQindex*> indices(s.get_indices());
    vector<SQindex*>::iterator it = indices.begin();
    while(it != indices.end()){
      IndLabel += (*it)->get_index();
      if (it != indices.end()-1) IndLabel += ",";
      ++it;
    }
    os << s.get_name() << "(" << IndLabel << ") ";
    return os;
  }

  // *********************************************************
  //  Taken from symPermutes in SQA
  // *********************************************************
  void SQtensor::makeAllPerms()
  {
    IIvector tuples(1);
    for(size_t i = 0;i < indices_.size();++i) tuples[0].push_back(int(i)); 
    Ivector factors;
    factors.push_back(1);
    bool allFound = false;
    while(!allFound){
      allFound = true;
      IIvector newTuples;
      Ivector newFactors;
      for(size_t j = 0;j < tuples.size();++j){
        for(size_t num_f = 0;num_f < symm_.first.size();++num_f){
          Ivector symm = symm_.first[num_f];
          newFactors.push_back(symm_.second[num_f]*factors[j]);
          Ivector temp;
          for(size_t num_s = 0;num_s < symm.size();++num_s){
            int i = symm[num_s];
            temp.push_back(tuples[j][i]);
          } // End num_s
          newTuples.push_back(temp);
	} // End num_f
      } // End j
      while(newTuples.size()){
        bool isNew = true;
        for(size_t num_tup = 0;num_tup < tuples.size();++num_tup){
          Ivector tup(tuples[num_tup]);
          if(newTuples[0] == tup){
            isNew = false;
            break;
	  } // End if
	} // End num_tup
        if(isNew){
          allFound = false;
          tuples.push_back(newTuples[0]);
          factors.push_back(newFactors[0]);
	} // End if
        newTuples.erase(newTuples.begin());
        newFactors.erase(newFactors.begin());
      } // End newTuples.size()
      //      cout << "MADE .... " << endl;
    } // End allFound
    permutations_ = tuples;
    factors_      = factors;
  }

  // *********************************************************
  //  Taken from sortIndices in SQA
  // *********************************************************
  int SQtensor::sortIndices()
  {
    IIvector tup, tuples;
    Ivector fac, factors;
    //    makeAllPerms();

    tup = permutations_;
    fac = factors_;

    tuples = tup;
    factors = fac;

    Ivector scores; 
    for(size_t num_tup = 0;num_tup < tuples.size();++num_tup) scores.push_back(0);
    for(size_t i = 0;i < indices_.size()-1;++i){
      for(size_t j = i+1;j < indices_.size();++j){
        for(size_t k = 0;k < tuples.size();++k){
          if(*indices_[tuples[k][i]] < *indices_[tuples[k][j]]) ++scores[k];
	} // End k
      } // End j
//*OLD*       int maxScore = 0;
//*OLD*       for(size_t I = 0;I < scores.size();++I) maxScore = (scores[I] > maxScore ? scores[I] : maxScore);
      int maxScore = *max_element(scores.begin(), scores.end()); //*NEW* 
      int I = 0;
      while(I < scores.size()){
        if(scores[I] < maxScore) {
          scores.erase(scores.begin()+I);
          tuples.erase(tuples.begin()+I);
          factors.erase(factors.begin()+I);
        }
        else ++I;
      } // End while
      if(!scores.size()) break;
    } // End i

    if(scores.size() > 1){
      bool unique = true;
      for(size_t i = 1;i < scores.size();++i){
        for(size_t j = 0;j < tuples[i].size();++j){
          if(indices_[tuples[i][j]] != indices_[tuples[0][j]]) {
            unique = false;
            break;
	  }
	} // End j
        if(!unique) break;
      } // End i
      if(!unique)
        throw, "No unique winner produced when sorting the indices";
    } // End if

    vector<SQindex*> newIndices;
    for(size_t i = 0;i < tuples[0].size();++i)
      newIndices.push_back(indices_[tuples[0][i]]);

    indices_ = newIndices;
    return factors[0];
  }

  // *********************************************************
  // 
  // *********************************************************
  bool SQtensor::hasIndex(const SQindex *i) const
  {
    bool retval = false;
    for(size_t num = 0;num < indices_.size();++num)
      if(*indices_[num] == *i){
	//        cout << "This is it >> " << *indices_[num] << endl; // TEST
        retval = true;
      }
//*BUG*     vector<SharedIndex>::iterator it = find(indices_.begin(), indices_.end(), const_cast<SharedIndex>(i));
//*BUG*     if(it != indices_.end()) retval = true;
    return retval; 
  }


  // *********************************************************
  // Rotate indices according to i-th permutations
  // *********************************************************
  void SQtensor::rotateIndices(const size_t i)
  {
    vector<SQindex*> newIndices;
    newIndices.reserve(indices_.size());
    for(size_t j = 0;j < permutations_[i].size();++j){
      newIndices.push_back(indices_[permutations_[i][j]]);
    } // End j
    indices_ = newIndices;
    // Consider sign if, sign_ is enabled.
  }


  // *********************************************************
  // 
  // *********************************************************
  void SQtensor::print_symm() const
  {
    size_t dim = indices_.size();
    cout << endl;
    cout << boost::format("@@@ Symm of %s") % name_ << endl;
    int count = 0;
    for(size_t i = 0;i < permutations_.size();++i){
      cout << boost::format("[%5d] ") % count << "[";
      for(int j = 0;j < permutations_[i].size();++j){
        cout << permutations_[i][j] << (j != dim-1 ? "," : "], ");
      }
      cout << factors_[i] << endl;
      ++count;
    }

  }


  // *********************************************************
  // Convert from Dirac 2 Mulliken
  // *********************************************************
  void SQtensor::convertD2M()
  {
    if(not_ == Mulliken) return;
    size_t order = indices_.size() / 2;
    // Replace indices in Mulliken order ....
    Ivector temp_i; 
    vector<SQindex*> newIndices; newIndices.reserve(2*order);
    for(size_t i = 0;i < order;++i){
      newIndices.push_back(indices_[i      ]);
      newIndices.push_back(indices_[i+order]);
      temp_i.push_back(i      );
      temp_i.push_back(i+order);
    } // End i
    indices_ = newIndices;

    // Prepare replacement array
    typedef pair<int, int> p_int;
    vector<p_int> x;
    for(size_t i = 0;i < 2*order;++i){
      p_int temp;
      temp.first  = i;
      temp.second = temp_i[i];
      x.push_back(temp);
    } // End i
    sort(x.begin(), x.end(), SecGreat());
    
//    for(size_t i = 0;i < x.size();++i) cout << x[i].first << ", " << x[i].second << " : ";
//    cout << endl;

    IIvector newPerm; newPerm.reserve(permutations_.size());
    IIvector::iterator p = permutations_.begin();
    for(;p != permutations_.end();++p){
      vector<int> q;
      newPerm.push_back(q);
      for(size_t i = 0;i < order;++i){
	IIvector::iterator it = newPerm.end(); --it;
        it->push_back(x[(*p)[i      ]].first);
        it->push_back(x[(*p)[i+order]].first);
      } // End i
    } // End p
    permutations_ = newPerm;
    sortIndices();
    not_ = Mulliken; // Notation is changed to Mulliken
  }


  // *********************************************************
  // Convert from Mulliken 2 Dirac (not tested yet!!!)
  // *********************************************************
  void SQtensor::convertM2D()
  {
//*OLD*     cout << "SQtensor::convertM2D() is not yet implemented, sorry" << endl;
//*OLD*     abort();
    if(not_ == Dirac) return;
    size_t order = indices_.size() / 2;
    // Replace indices in Dirac order ....
    Ivector temp_i;
    vector<SQindex*> newIndices; newIndices.reserve(2*order);
    for(size_t i = 0;i < order;++i){
      newIndices.push_back(indices_[2*i  ]);
      temp_i.push_back(2*i  );
    } // End i
    for(size_t i = 0;i < order;++i){
      newIndices.push_back(indices_[2*i+1]);
      temp_i.push_back(2*i+1);
    } // End i

    // Prepare replacement array
    typedef pair<int, int> p_int;
    vector<p_int> x;
    for(size_t i = 0;i < 2*order;++i){
      p_int temp;
      temp.first  = i;
      temp.second = temp_i[i];
      x.push_back(temp);
    } // End i
    sort(x.begin(), x.end(), SecGreat());

    IIvector newPerm; newPerm.reserve(permutations_.size());
    IIvector::iterator p = permutations_.begin();
    for(;p != permutations_.end();++p){
      vector<int> q;
      newPerm.push_back(q);
      for(size_t i = 0;i < order;++i){
	IIvector::iterator it = newPerm.end(); --it;
        it->push_back(x[(*p)[i      ]].first);
        it->push_back(x[(*p)[i+order]].first);
      } // End i
    } // End p
    permutations_ = newPerm;
    sortIndices();
    not_ = Dirac; // Notation is changed to Dirac
       
  }

  // *********************************************************
  // Now available for tensors in both Dirac and Mulliken notation.
  // *********************************************************
  string SQtensor::convert2LaTeX() const
  {
    string ten("");
    if(name_ == kDelta_name()) ten += "\\delta";
    else if(is_RDM(name_)) ten += RDM_name();
    else if(is_sfGen(name_)) ten += sfGen_name();
    else ten += name_;

    if(indices_.size()){
      if(not_ == Dirac){
        ten += "_{";
        for(size_t i = 0;i < indices_.size();++i){
          string num(indices_[i]->get_index());
          num.erase(num.begin());
          //cout << "num: " << num << endl;
          if(atoi(num.c_str())){
            ten += indices_[i]->get_index().at(0);
            ten += "_{";
            ten += num;
            ten += "}";
          } // End if
          else ten += indices_[i]->get_index();
          if     (i == indices_.size()/2-1) ten += "}^{";
          else if(i == indices_.size()  -1) ten += "} ";       
        } // End i
      } // End if
      else if(not_ == Mulliken){
        ten += "_{";
        for(size_t i = 0;i < indices_.size()/2;++i){
          string num(indices_[2*i]->get_index());
          num.erase(num.begin());
          //cout << "num: " << num << endl;
          if(atoi(num.c_str())){
            ten += indices_[2*i]->get_index().at(0);
            ten += "_{";
            ten += num;
            ten += "}";
          } // End if
          else ten += indices_[2*i]->get_index();         
	} // End i
        ten += "}^{";
        for(size_t i = 0;i < indices_.size()/2;++i){
          string num(indices_[2*i+1]->get_index());
          num.erase(num.begin());
          //cout << "num: " << num << endl;
          if(atoi(num.c_str())){
            ten += indices_[2*i+1]->get_index().at(0);
            ten += "_{";
            ten += num;
            ten += "}";
          } // End if
          else ten += indices_[2*i+1]->get_index();         
	} // End i
        ten += "}";
      } // End if
      else if(not_ == None){
        if     (name_ == aCre_name()) ten += "^{";
	else if(name_ == aDes_name()) ten += "_{";
	else { 
	  cout << "SQtensor.cc :: I don't know how to handle this tensor in convert2LaTeX(" << *this << ")" << endl;
	  abort();
	} // End else
	string num(indices_[0]->get_index());
	num.erase(num.begin());
	//cout << "num: " << num << endl;
	if(atoi(num.c_str())){
	  ten += indices_[0]->get_index().at(0);
	  ten += "_{";
	  ten += num;
	  ten += "}";
	} // End if
	else ten += indices_[0]->get_index();         
	if(indices_.size() != 1){
	  cout << "SQtensor.cc :: I don't know how to handle this tensor in convert2LaTeX(" << *this << ")" << endl;
	  abort();
	} // End if
      } // Else if
    } // End if
      
    return ten;
  }


  // *********************************************************
  // Kronecker's delta class
  // *********************************************************
  kDelta::kDelta()
  {}

  // *********************************************************
  // Kronecker's delta class
  // *********************************************************
  kDelta::kDelta(const vector<SQindex*> indices)
  { 
    // *VERY*, *VERY* bitchy! But, works!
    name_         = kDelta_name();
    indices_      = indices;
    isCommutable_ = true;
    not_          = Dirac;

    if(indices_.size() != 2){
      cout << "Kronecker's delta should contain two indices" << endl;
      abort();
    }

    IIvector perm;
    Ivector perm1;
    perm1.push_back(0);
    perm1.push_back(1);
    Ivector perm2;
    perm2.push_back(1);
    perm2.push_back(0);

    perm.push_back(perm1);
    perm.push_back(perm2);

    Ivector fac;
    fac.push_back(1);
    fac.push_back(1);
  
    symm_.first = perm;
    symm_.second = fac;

    makeAllPerms();
    sortIndices();
  }

  // *********************************************************
  // Spin-free unitary group generator class
  // *********************************************************

  // Definition:
  //
  // sfGen(p,q,r,s) = sum_{tau,sigma} Cre(p,tau) Cre(q,sigma) Des(s,sigma) Des(r,tau)
  //
  // in other words, sfGen(p,q,r,s) = E^{pq}_{rs}

  sfGen::sfGen(const vector<SQindex*> indices)
  { 
    name_         = sfGen_name();
    indices_      = indices;
    isCommutable_ = false;
    not_          = Dirac;

    if(indices_.size()%2){
      cout << "A spin-free unitary group generator must have even-number of indices" << endl;
      abort();
    }

    if(!indices_.size()){
      cout << "A spin-free unitary group generator must have at least 2 indices" << endl;
      abort();
    }
    
    int order_ = indices_.size() / 2; // This should become a member of this class
    ostringstream stm; stm << order_;
    name_ += stm.str();

//*OLD*     // The ordinary permutation order .... 
//*OLD*     Ivector temp;
//*OLD*     for(size_t i = 0;i < indices_.size();++i) temp.push_back((int)i);
//*OLD*     symm_.first.push_back(temp);
//*OLD*     symm_.second.push_back(1);

//*OLD*     for(size_t j = 0;j < temp.size();++j) cout << temp[j]; //*TEST*
//*OLD*     cout << endl;                                          //*TEST*
//*OLD* 
//*OLD*     cout << "made" << endl;
//*OLD*     cout << order_-1 << endl;
//*OLD* 
//*OLD*     // The pair-wise replaced orders ..... 
//*OLD*     for(int i = 0;i < (order_-1);++i){
//*OLD*     cout << "made" << endl;
//*OLD*       Ivector temp;
//*OLD*       if(i == 0) temp.push_back(1);
//*OLD*       else       temp.push_back(0);
//*OLD*       for(int j = 1;j < 2*order_;++j){
//*OLD*         if(j==i)               temp.push_back(i+1);
//*OLD*         else if(j==i+1)        temp.push_back(i);
//*OLD*         else if(j==i+order_)   temp.push_back(i+1+order_);
//*OLD*         else if(j==i+1+order_) temp.push_back(i  +order_);
//*OLD*         else                   temp.push_back(j);
//*OLD*       } // End j
//*OLD* 
//*OLD*       for(size_t j = 0;j < temp.size();++j) cout << temp[j]; //*TEST*
//*OLD*       cout << endl;                                          //*TEST*
//*OLD* 
//*OLD*       symm_.first.push_back(temp);
//*OLD*       symm_.second.push_back(1);
//*OLD*     } // End i

    IIvector perm(makePermutations(order_));
    for(size_t i = 0;i < perm.size();++i){
      Ivector p = perm[i];
      for(size_t j = 0;j < (size_t)order_;++j){
        int num = p[j] + order_;
        p.push_back(num);
      } //End j
      permutations_.push_back(p);
      factors_.push_back(1);
    } // End i

    sortIndices();

//*DEBUG*     for(size_t i = 0;i < permutations_.size();++i){
//*DEBUG*       cout << "[ ";
//*DEBUG*       for(size_t j = 0;j < permutations_[i].size();++j){
//*DEBUG*         cout << permutations_[i][j] << " ";
//*DEBUG*       }
//*DEBUG*       cout << " ], " << factors_[i] << endl;
//*DEBUG*     }

  }

//*OLD*   // *********************************************************
//*OLD*   // 
//*OLD*   // *********************************************************
//*OLD*   RDM::RDM()
//*OLD*   {}

  // *********************************************************
  // Reduced-Density Matrix (RDM) class
  // *********************************************************
  RDM::RDM(const vector<SQindex*> indices)
  { 
    name_         = RDM_name();
    indices_      = indices;
    isCommutable_ = true;
    not_          = Dirac;

    if(indices_.size()%2){
      cout << "A spin-free unitary group generator must have even-number of indices" << endl;
      abort();
    }
    else if(!indices_.size()){
      cout << "Spin-free unitary group generator must contain at least 2 indices" << endl;
      abort();
    }
    
    int order_ = indices_.size() / 2; // This should become a member of this class
    ostringstream stm; stm << order_;
    name_ += stm.str();

    Ivector temp;
    for(size_t i = 0;i < indices_.size();++i) temp.push_back((int)i);
    //    IIvector perm = makeTuples1(order_, temp);
    IIvector perm(makePermutations(order_));

//*DEBUG*     int count = 0;
//*DEBUG*     for(int i = 0;i < perm.size();++i){
//*DEBUG*       cout << boost::format("[%d] ") % count;
//*DEBUG*       for(int j = 0;j < perm[i].size();++j) cout << perm[i][j] << " ";
//*DEBUG*       cout << endl;
//*DEBUG*       ++count;
//*DEBUG*     }

    for(size_t I = 0;I < perm.size();++I){
      Ivector p(perm[I]);
      for(size_t i = 0;i < (size_t)order_;++i){
        int num = p[i] + order_;
        p.push_back(num);
      } // End p
      permutations_.push_back(p);
      factors_.push_back(1);

      Ivector q(perm[I]);
      for(size_t i = 0;i < (size_t)order_;++i) q[i] += order_;
      for(size_t i = 0;i < (size_t)order_;++i){
        int num = q[i] - order_;
        q.push_back(num);
      } // End q
      permutations_.push_back(q);
      factors_.push_back(1);      
    } // End I

    sortIndices();

//*DEBUG*     for(size_t i = 0;i < permutations_.size();++i){
//*DEBUG*       cout << "[ ";
//*DEBUG*       for(size_t j = 0;j < permutations_[i].size();++j){
//*DEBUG*         cout << permutations_[i][j] << " ";
//*DEBUG*       }
//*DEBUG*       cout << " ], " << factors_[i] << endl;
//*DEBUG*     }

  }

  // *********************************************************
  // Creation operator class
  // *********************************************************
  aCre::aCre(SQindex* index)
  { 
    vector<SQindex*> indices;
    indices.push_back(index);
    name_         = aCre_name();
    indices_      = indices;
    isCommutable_ = false;
    not_          = None;

    IIvector perm;
    Ivector perm1;
    perm1.push_back(0);
    perm.push_back(perm1);

    Ivector fac;
    fac.push_back(1);
  
    symm_.first  = perm;
    symm_.second = fac;
  }

  // *********************************************************
  // Destruction operator class
  // *********************************************************
  aDes::aDes(SQindex* index)
  { 
    vector<SQindex*> indices;
    indices.push_back(index);
    name_         = aDes_name();
    indices_      = indices;
    isCommutable_ = false;
    not_          = None;

    IIvector perm;
    Ivector perm1;
    perm1.push_back(0);
    perm.push_back(perm1);

    Ivector fac;
    fac.push_back(1);
  
    symm_.first  = perm;
    symm_.second = fac;
  }

}} // Femto::core


