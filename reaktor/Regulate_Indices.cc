//
//  Regulate_Indices.cc
//  
//
//  Created by Masaaki Saitow on 12/10/25.
//  Copyright (c) 2012 Masaaki Saitow. All rights reserved.
//

#include <algorithm>
#include <fstream>
#include <cmath>
#include <string>
#include <boost/format.hpp>
#include <SQreaktor.hpp>

using namespace std;
using namespace Femto;
using namespace Femto::Core;

namespace Femto { namespace Reaktor { 

  // *********************************************************
  // 
  // *********************************************************
  void SQreaktor::regulate_indices()
  {
    if(LTensor_.get_indices().size() != 4 && isBareLHS_){
      cout << "Argument 0 cannnot be treated as a bareampack .... " << endl;
      abort();
    } // End if

    vector<SQindex*> Linds = LTensor_.get_indices();
    for(vector<SQindex*>::iterator i = Linds.begin();i != Linds.end();++i)
      if((*i)->get_isSummed()){
        cout << *i << " in 0th argument is a kind of dummy index." << endl;
        abort(); 
      } // End if

    for(vector<SQterm>::iterator t = inTerms_.begin();t != inTerms_.end();++t){
      if(!t->get_isInCanonical()){
        cout << "Term, " << *t << " is not in canonical form .... " << endl;
        abort();
      } // End if
    } // End t

    if(is_sfGen(LTensor_.get_name())){
      cout << "Factorize: 1st argument cannot be a spin-free unitary group generator ... " << endl;
      abort();
    }
    for(vector<SQterm>::iterator t = inTerms_.begin();t != inTerms_.end();++t){
      for(size_t num_t = 0;num_t < t->get_tensors().size();++num_t){
        if(is_sfGen(t->get_tensors()[num_t].get_name())){
          cout << "Factorize: A spin-free unitary group generator is detected in the 2nd argument .... " << endl;
          abort();
	} // End if
      } // End num_t
    } // End t

    if(LTensors_.size() != inTerms_.size()){
      cout << "Regulate_Indices: LTensors_ are not set. Do process_kDeltas first. if already done, something is wrong." << endl;
      abort();
    } // End if

    if(!inTerms_.size()) return;

    // If ERI and RDM are not of Mulliken form, transform them
    for(vector<SQterm>::iterator t = inTerms_.begin();t != inTerms_.end();++t){
      for(size_t num_t = 0;num_t < t->get_tensors().size();++num_t){
        if(t->get_tensors()[num_t].get_name() == name_h2_){
          if(t->get_tensors()[num_t].get_notation() == (notation)Dirac) 
	    t->get_tensors_ptr()[num_t]->convertD2M(); // Convert Dirac->Mulliken
	} // End if
        else if(is_RDM(t->get_tensors()[num_t].get_name())){
          if(t->get_tensors()[num_t].get_notation() == (notation)Dirac) 
	    t->get_tensors_ptr()[num_t]->convertD2M(); // Convert Dirac->Mulliken         
	} // End if
      } // End num_t
    } // End t

    // Convert term by term ....
    int numTerm(0); 
    for(vector<SQterm>::iterator thisTerm = inTerms_.begin();thisTerm != inTerms_.end();++thisTerm, ++numTerm){
      //* --- * // Starts from L.1646 in sqaConvert.py ... 
      //* --- * cout << "! No."  << numTerm << endl;
      //* --- * cout << "! " << LTensor_ << " <-- " << endl;
      //* --- * cout << "! " << *thisTerm << endl;

      //* --- * // Retain the consistency between LTensor_ and thisTerm
      //* --- * vector<SQindex*> Linds(LTensor_.get_indices());
      //* --- * vector<SQindex*> thisBoddies = thisTerm->get_summedBody();
      //* --- * for(size_t num_i = 0;num_i < Linds.size();++num_i){
      //* --- *   SQindex* Li = Linds[num_i];
      //* --- *   for(vector<SQindex*>::iterator Bi = thisBoddies.begin();Bi != thisBoddies.end();++Bi){
      //* --- *     if(*Li == *(*Bi)) LTensor_.put_indices(num_i, *Bi);
      //* --- * 	} // End Bi
      //* --- * } // End Li

      //* --- * // Break Kronecker's delta
      //* --- * vector<SQtensor> new_t;
      //* --- * vector<SQindex> newLTind;
      //* --- * for(size_t i = 0;i < thisTerm->get_tensors().size();++i){
      //* --- *   SQtensor t = thisTerm->get_tensors()[i];
      //* --- *   if(t.get_name() != kDelta_name()) new_t.push_back(t);
      //* --- *   else{
      //* --- *     SQindex* i1 = t.get_indices()[0];
      //* --- *     SQindex* i2 = t.get_indices()[1];
      //* --- *     *i2 = *i1;
      //* --- *     newLTind.clear();
      //* --- *     for(size_t j = 0;j < LTensor_.get_indices().size();++j) newLTind.push_back(*LTensor_.get_indices()[j]);
      //* --- * 	}
      //* --- * } // End i
      //* --- * if(new_t.size() != thisTerm->get_tensors().size()){
      //* --- *   thisTerm->set_tensors(new_t);
      //* --- * 
      //* --- *   // If so, consistency between LTensor_ and thisTerm has to be retain again
      //* --- *   vector<SQindex*> thisBoddies = thisTerm->get_summedBody();
      //* --- *   for(size_t num_i = 0; num_i = newLTind.size();++num_i){
      //* --- *     SQindex Li = newLTind[num_i];
      //* --- *     for(vector<SQindex*>::iterator Bi = thisBoddies.begin();Bi != thisBoddies.end();++Bi){
      //* --- *       if(Li == *(*Bi)) LTensor_.put_indices(num_i, *Bi);
      //* --- * 	} // End Bi
      //* --- *   } // End Li  
      //* --- * } // End if

      // Achieve the best matching
      typedef pair<bool, SQtensor> conTen;
      conTen h2f, ampf, d4f;
      h2f.first  = false;
      ampf.first = false;
      d4f.first  = false;
      for(size_t num_t = 0;num_t < thisTerm->get_tensors().size();++num_t){
        if     (thisTerm->get_tensors()[num_t].get_name() == name_h2_){
          h2f.first  = true;
          h2f.second = (thisTerm->get_tensors()[num_t]);
	} // End if
        else if(thisTerm->get_tensors()[num_t].get_name() == name_amp_){
          ampf.first  = true;
          ampf.second = (thisTerm->get_tensors()[num_t]); 
          //cout << thisTerm->get_tensors()[num_t].get_indices()[extamp_];
	} // End if
        else if(thisTerm->get_tensors()[num_t].get_name() == name_d4_){
          d4f.first  = true;
          d4f.second = (thisTerm->get_tensors()[num_t]);
	} // End if
      } // End num_t

      // Match the external indices of the bareamp and eri, if possible.
      // If impossible, consider such the possibility of d4 and eri.
      bool SetT_V(false);
      bool SetT_D(false);
      bool SetL_V(false);
      bool SetL_T(false);
      bool SetV_D(false);
      bool SetV_v(false);
      // LHS v.s. T2
      if(ampf.first && isBareLHS_){
        vector<SQindex*> ampind(ampf.second.get_indices());
        vector<SQindex*> LTind(LTensors_[numTerm].get_tensors()[0].get_indices());
        for(vector<SQindex*>::iterator i = ampind.begin();i != ampind.end();++i){
          for(vector<SQindex*>::iterator j = LTind.begin();j != LTind.end();++j){
            if(**i == **j){
              IIvector ampperm(ampf.second.get_perms());
              IIvector LTperm(LTensors_[numTerm].get_tensors()[0].get_perms());
              typedef pair<bool, size_t> found;
              found amp, LT;
              amp.first = false;
              LT.first = false;
              amp.second = 0;
              LT.second = 0;
              for(size_t k = 0;k < ampperm.size();++k){
                if(*ampind[ampperm[k][extamp_]] == **i) { amp.first = true; amp.second = k; break; } 
                if(amp.first) break;
	      } //End k
              for(size_t k = 0;k < LTperm.size();++k){
                if(*LTind[LTperm[k][extamp_]] == **j) { LT.first = true; LT.second = k; break; } 
                if(LT.first) break;
	      } //End k

              if(amp.first && LT.first) {
		vector<SQtensor*>::iterator amp_t;
                vector<SQtensor*> ten_ptr(thisTerm->get_tensors_ptr()); 
                for(vector<SQtensor*>::iterator t = ten_ptr.begin();t != ten_ptr.end();++t){
                  if     ((*t)->get_name() == name_amp_) 
		    { (*t)->rotateIndices(amp.second); amp_t = t; break; }
	        } // End t
                //SQtensor* Lptr(LTensors_[numTerm].get_tensors()[0].rotateIndices(LT.second));
                SQtensor* Lptr(LTensors_[numTerm].get_tensors_ptr()[0]);
		Lptr->rotateIndices(LT.second);
                cout << ">> Indices of BareAmp are rotated to match with LHS.           -> " << numTerm << endl;
		cout << "   |--> ( " << LTensors_[numTerm].get_tensors()[0] << "/ " << **amp_t << ")" << endl;
                SetL_T = true;
	      } // End if
              if(SetL_T) break;
	    } // End if
	  } // End j
          if(SetL_T) break;
	} // End i
      } // End if

      // V2 v.s. D4
      if(Set_V_D_)
	//if(h2f.first && d4f.first && !SetT_V && !SetL_V){ // OLD This may be better (2013/01/08)
      if(h2f.first && d4f.first){ // <-- NEW 2013/01/08
        vector<SQindex*> h2ind(h2f.second.get_indices());
        vector<SQindex*> d4ind(d4f.second.get_indices());
        //SQindex* the_first = h2ind[exth2_];
        for(vector<SQindex*>::iterator i = h2ind.begin();i != h2ind.end();++i){
          for(vector<SQindex*>::iterator j = d4ind.begin();j != d4ind.end();++j){
            if(**i == **j){
              IIvector h2perm(h2f.second.get_perms());
              IIvector d4perm(d4f.second.get_perms());
              typedef pair<bool, size_t> found;
              found h2, d4;
              h2.first = false;
              d4.first = false;
              h2.second = 0;
              d4.second = 0;
              for(size_t k = 0;k < h2perm.size();++k){
                if(*h2ind[h2perm[k][exth2_]] == **i) { h2.first = true; h2.second = k; break; } 
                if(h2.first) break;
	      } //End k
              for(size_t k = 0;k < d4perm.size();++k){
                if(*d4ind[d4perm[k][0]] == **j) { d4.first = true; d4.second = k; break; } 
                if(d4.first) break;
	      } //End k

              if(h2.first && d4.first) {
                vector<SQtensor*> ten_ptr(thisTerm->get_tensors_ptr()); 
                for(vector<SQtensor*>::iterator t = ten_ptr.begin();t != ten_ptr.end();++t){
                  if     ((*t)->get_name() == name_h2_) 
                    (*t)->rotateIndices(h2.second);
                  else if((*t)->get_name() == name_d4_) 
                    (*t)->rotateIndices(d4.second);
	        } // End t
                cout << ">> Indices of ERI and D4 are rotated to match with each other. -> " << numTerm << endl;
                //cout << "H2: " << h2.second << " D4: " << d4.second << endl;
                //cout << "*TEST* " << *thisTerm << endl; 
                SetV_D = true;
	      } // End if
	    } // End if
            if(SetV_D) break;
	  } // End j
          if(SetV_D) break;
	} // End i
      } // End if

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // T2 v.s. D4 (Newly added 2013/01/08)
      if(Set_T_D_)
      if(ampf.first && d4f.first && !SetV_D){
        vector<SQindex*> t2ind(ampf.second.get_indices());
        vector<SQindex*> d4ind(d4f.second.get_indices());
        for(vector<SQindex*>::iterator i = t2ind.begin();i != t2ind.end();++i){
          for(vector<SQindex*>::iterator j = d4ind.begin();j != d4ind.end();++j){
            if(**i == **j){
              IIvector t2perm(ampf.second.get_perms());
              IIvector d4perm(d4f.second.get_perms());
              typedef pair<bool, size_t> found;
              found t2, d4;
              t2.first = false;
              d4.first = false;
              t2.second = 0;
              d4.second = 0;
              for(size_t k = 0;k < t2perm.size();++k){
                if(*t2ind[t2perm[k][extamp_]] == **i) { t2.first = true; t2.second = k; break; } 
                if(t2.first) break;
	      } //End k
              for(size_t k = 0;k < d4perm.size();++k){
                if(*d4ind[d4perm[k][0]] == **j) { d4.first = true; d4.second = k; break; } 
                if(d4.first) break;
	      } //End k

              if(t2.first && d4.first) {
                vector<SQtensor*> ten_ptr(thisTerm->get_tensors_ptr()); 
                for(vector<SQtensor*>::iterator t = ten_ptr.begin();t != ten_ptr.end();++t){
                  if     ((*t)->get_name() == name_amp_) 
                    (*t)->rotateIndices(t2.second);
                  else if((*t)->get_name() == name_d4_) 
                    (*t)->rotateIndices(d4.second);
	        } // End t
                cout << ">> Indices of T2 and D4 are rotated to match with each other. -> " << numTerm << endl;
                //cout << "H2: " << h2.second << " D4: " << d4.second << endl;
                //cout << "*TEST* " << *thisTerm << endl; 
                SetT_D = true;
	      } // End if
	    } // End if
            if(SetT_D) break;
	  } // End j
          if(SetT_D) break;
	} // End i
      } // End if
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      // LHS v.s. V2
      if(h2f.first && isBareLHS_ && !SetV_D && !SetL_T){
        vector<SQindex*> h2ind(h2f.second.get_indices());
        vector<SQindex*> LTind(LTensors_[numTerm].get_tensors()[0].get_indices());
        for(vector<SQindex*>::iterator i = h2ind.begin();i != h2ind.end();++i){
          for(vector<SQindex*>::iterator j = LTind.begin();j != LTind.end();++j){
            if(**i == **j){
              IIvector h2perm(h2f.second.get_perms());
              IIvector LTperm(LTensors_[numTerm].get_tensors()[0].get_perms());
              typedef pair<bool, size_t> found;
              found h2, LT;
              h2.first = false;
              LT.first = false;
              h2.second = 0;
              LT.second = 0;
              for(size_t k = 0;k < h2perm.size();++k){
                if(*h2ind[h2perm[k][exth2_]] == **i) { h2.first = true; h2.second = k; break; } 
                if(h2.first) break;
	      } //End k
              for(size_t k = 0;k < LTperm.size();++k){
                if(*LTind[LTperm[k][extamp_]] == **j) { LT.first = true; LT.second = k; break; } 
                if(LT.first) break;
	      } //End k

	      SQtensor* ERI(NULL);
              if(h2.first && LT.first) {
                vector<SQtensor*> ten_ptr(thisTerm->get_tensors_ptr()); 
                for(vector<SQtensor*>::iterator t = ten_ptr.begin();t != ten_ptr.end();++t){
                  if     ((*t)->get_name() == name_h2_) 
                    { (*t)->rotateIndices(h2.second); ERI = (*t); }
	        } // End t
                SQtensor* Lptr(LTensors_[numTerm].get_tensors_ptr()[0]);
		Lptr->rotateIndices(LT.second);
                cout << ">> Indices of ERI are rotated to match with LHS.               -> " << numTerm << endl;
		cout << "   |--> ( " << *Lptr << "/ " << *ERI << ")" << endl;
                SetL_V = true;
	      } // End if
              if(SetL_V) break;
	    } // End if
	  } // End j
          if(SetL_V) break;
	} // End i
      } // End if

      // T2 v.s. V2
      if(h2f.first && ampf.first && !SetV_D && !SetL_V || SetL_T){
        vector<SQindex*> h2ind(h2f.second.get_indices());
        vector<SQindex*> ampind(ampf.second.get_indices());
        for(vector<SQindex*>::iterator i = h2ind.begin();i != h2ind.end();++i){
          for(vector<SQindex*>::iterator j = ampind.begin();j != ampind.end();++j){
            if(**i == **j){
              IIvector h2perm(h2f.second.get_perms());
              IIvector ampperm(ampf.second.get_perms());
              typedef pair<bool, size_t> found;
              found h2, amp;
              h2.first = false;
              amp.first = false;
              h2.second = 0;
              amp.second = 0;
              for(size_t k = 0;k < h2perm.size();++k){
                if(*h2ind[h2perm[k][exth2_]] == **i) { h2.first = true; h2.second = k; break; } 
                if(h2.first) break;
	      } //End k
              for(size_t k = 0;k < ampperm.size();++k){
                if(*ampind[ampperm[k][extamp_]] == **j) { amp.first = true; amp.second = k; break; } 
                if(amp.first) break;
	      } //End k

              if(h2.first && amp.first) {
                vector<SQtensor*> ten_ptr(thisTerm->get_tensors_ptr()); 
                for(vector<SQtensor*>::iterator t = ten_ptr.begin();t != ten_ptr.end();++t){
                  if     ((*t)->get_name() == name_h2_) 
                    (*t)->rotateIndices(h2.second);
                  else if((*t)->get_name() == name_amp_) 
                    (*t)->rotateIndices(amp.second);
	        } // End t
                cout << ">> Indices of ERI are rotated to match with Bareamp.           -> " << numTerm << endl;
                SetT_V = true;
	      } // End if
              if(SetT_V) break;
	    } // End if
	  } // End j
          if(SetT_V) break;
	} // End i
      } // End if

      // Set loading index of ERI to become virtual if possible
      //if(!SetT_V && !SetL_V && /*!SetV_D &&*/ h2f.first){ // OLD setting (2013/01/08). This may be better
      if(!SetT_V && !SetL_V && !SetV_D && h2f.first){ // <-- NEW 2013/01/08
        vector<SQindex*> h2ind(h2f.second.get_indices());
        for(size_t num_i = 0;num_i < h2f.second.get_indices().size();++num_i){
          SQindex* i = h2f.second.get_indices()[num_i];
          if(i->get_char() == (char_state)2){
            IIvector h2perm(h2f.second.get_perms());
            typedef pair<bool, size_t> found;
            found h2;
            h2.first = false;
            h2.second = 0;
            for(size_t k = 0;k < h2perm.size();++k){
              if(*h2ind[h2perm[k][exth2_]] == *i) { h2.first = true; h2.second = k; break; } 
	    } // End k

            if(h2.first){
              vector<SQtensor*> ten_ptr(thisTerm->get_tensors_ptr()); 
              for(vector<SQtensor*>::iterator t = ten_ptr.begin();t != ten_ptr.end();++t){
                if((*t)->get_name() == name_h2_) (*t)->rotateIndices(h2.second);
	      } // End t
              cout << ">> The indices of ERI are rotated to became virtual.           -> " << numTerm << endl;
              SetV_v = true;              
	    } // End if
	  } // End if
          if(SetV_v) break;
	} // End num_i
      } // End if
    
      //* --- * // Set all the external indices 
      //* --- * vector<SQtensor> tensors = thisTerm->get_tensors();
      //* --- * for(vector<SQtensor>::iterator t = tensors.begin();t != tensors.end();++t){
      //* --- *   // In case of ERI
      //* --- *   if(t->get_name() == name_h2_) {
      //* --- *     vector<SQindex*> inds = t->get_indices();
      //* --- *     inds[exth2_]->switch_isExt(true);
      //* --- * 	} // End if
      //* --- *   // In case of Bareamp
      //* --- *   else if(t->get_name() == name_amp_) {
      //* --- *     vector<SQindex*> inds = t->get_indices();
      //* --- *     inds[extamp_]->switch_isExt(true);
      //* --- * 	} // End if
      //* --- *   else if(t->get_name() == name_d4_) {
      //* --- *     vector<SQindex*> inds = t->get_indices();
      //* --- *     inds[0]->switch_isExt(true);        
      //* --- *     inds[1]->switch_isExt(true);        
      //* --- * 	}
      //* --- * } // End t
      //* --- * if(isBareLHS_){
      //* --- *   vector<SQindex*> inds = LTensor_.get_indices();
      //* --- *   inds[extamp_]->switch_isExt(true);
      //* --- * } // End if

    } // End thisTerm

    int count = 0;
    cout << endl;
    cout << "++ Indices are rotated to achive the best matching :: " << endl;
    for(vector<SQterm>::iterator t = inTerms_.begin();t != inTerms_.end();++t, ++count){
      cout << count << " : " << LTensors_[count].get_tensors()[0] << " += " << *t << endl;
    } // End t
    cout << endl;

    // Divide each term in inTerms_ into type_terms_
    vector<SQterm> noeri; noeri.reserve(inTerms_.size());
    vector<SQterm> eri_c; eri_c.reserve(inTerms_.size());
    vector<SQterm> eri_o; eri_o.reserve(inTerms_.size());
    vector<SQterm> eri_v; eri_v.reserve(inTerms_.size());
    vector<SQterm> d4c_c; eri_c.reserve(inTerms_.size());
    vector<SQterm> d4c_o; eri_o.reserve(inTerms_.size());
    vector<SQterm> d4c_v; eri_v.reserve(inTerms_.size());


    // Divide each term in LTensors_ into type_LTensors_
    vector<SQterm> LT_noeri; LT_noeri.reserve(inTerms_.size());
    vector<SQterm> LT_eri_c; LT_eri_c.reserve(inTerms_.size());
    vector<SQterm> LT_eri_o; LT_eri_o.reserve(inTerms_.size());
    vector<SQterm> LT_eri_v; LT_eri_v.reserve(inTerms_.size());
    vector<SQterm> LT_d4c_c; LT_eri_c.reserve(inTerms_.size());
    vector<SQterm> LT_d4c_o; LT_eri_o.reserve(inTerms_.size());
    vector<SQterm> LT_d4c_v; LT_eri_v.reserve(inTerms_.size());


    int num(0);
    for(vector<SQterm>::iterator t = inTerms_.begin();t != inTerms_.end();++t,++num){
      int flag = 0;
      // 0 : No ERI
      // 1 : ERI core
      // 2 : ERI active
      // 3 : ERI virtual
      // 4 : D4C core
      // 5 : D4C active
      // 6 : D4C virtual
      for(size_t num_t = 0;num_t < t->get_tensors().size();++num_t){
        if     (t->get_tensors()[num_t].get_name() == name_h2_){
          if     (t->get_tensors()[num_t].get_indices()[exth2_]->get_char() == core)
            flag = 1;
          else if(t->get_tensors()[num_t].get_indices()[exth2_]->get_char() == act)
            flag = 2;
          else if(t->get_tensors()[num_t].get_indices()[exth2_]->get_char() == virt)
            flag = 3;
	} // End if
	else if(is_D4C(t->get_tensors()[num_t].get_name())){
          if     (t->get_tensors()[num_t].get_indices()[extd4c_]->get_char() == core)
            flag = 4;
          else if(t->get_tensors()[num_t].get_indices()[extd4c_]->get_char() == act)
            flag = 5;
          else if(t->get_tensors()[num_t].get_indices()[extd4c_]->get_char() == virt)
            flag = 6;
	} // End if
      } // End num_t

      if     (flag == 0){
        noeri.push_back(*t);
        LT_noeri.push_back(LTensors_[num]);
      } // End if
      else if(flag == 1){
        eri_c.push_back(*t);
        LT_eri_c.push_back(LTensors_[num]);
      } // End if
      else if(flag == 2){
        eri_o.push_back(*t);
        LT_eri_o.push_back(LTensors_[num]);
      } // End if
      else if(flag == 3){
        eri_v.push_back(*t);
        LT_eri_v.push_back(LTensors_[num]);
      } // End if
      else if(flag == 4){
        d4c_c.push_back(*t);
        LT_d4c_c.push_back(LTensors_[num]);
      } // End if
      else if(flag == 5){
        d4c_o.push_back(*t);
        LT_d4c_o.push_back(LTensors_[num]);
      } // End if
      else if(flag == 6){
        d4c_v.push_back(*t);
        LT_d4c_v.push_back(LTensors_[num]);
      } // End if

    } // End t

    cout << "++ Terms are classified into these groups :: " << endl;
    if(noeri.size()){
      cout << ">> No ERIs <<" << endl;
      int count = 0;
      for(vector<SQterm>::iterator t = noeri.begin();t != noeri.end();++t,++count){
	cout << count << " : " << LT_noeri[count] << " += " << *t << endl;      
      } // End t
      cout << endl;
    } // End if

    if(eri_c.size()){
      cout << ">> ERIs loaded by core indices <<" << endl;
      int count = 0;
      for(vector<SQterm>::iterator t = eri_c.begin();t != eri_c.end();++t,++count){
	cout << count << " : " << LT_eri_c[count] << " += " << *t << endl;      
      } // End t
      cout << endl;
    } // End if

    if(eri_o.size()){
      cout << ">> ERIs loaded by active indices <<" << endl;
      int count = 0;
      for(vector<SQterm>::iterator t = eri_o.begin();t != eri_o.end();++t,++count){
	cout << count << " : " << LT_eri_o[count] << " += " << *t << endl;      
      } // End t
      cout << endl;
    } // End if

    if(eri_v.size()){
      cout << ">> ERIs loaded by virtual indices <<" << endl;
      int count = 0;
      for(vector<SQterm>::iterator t = eri_v.begin();t != eri_v.end();++t,++count){
	cout << count << " : " << LT_eri_v[count] << " += " << *t << endl;      
      } // End t
      cout << endl;
    } // End if

    if(d4c_c.size()){
      cout << ">> D4Cs loaded by core indices <<" << endl;
      int count = 0;
      for(vector<SQterm>::iterator t = d4c_c.begin();t != d4c_c.end();++t,++count){
	cout << count << " : " << LT_d4c_c[count] << " += " << *t << endl;      
      } // End t
      cout << endl;
    } // End if

    if(d4c_o.size()){
      cout << ">> D4Cs loaded by active indices <<" << endl;
      int count = 0;
      for(vector<SQterm>::iterator t = d4c_o.begin();t != d4c_o.end();++t,++count){
	cout << count << " : " << LT_d4c_o[count] << " += " << *t << endl;      
      } // End t
      cout << endl;
    } // End if

    if(d4c_v.size()){
      cout << ">> D4Cs loaded by virtual indices <<" << endl;
      int count = 0;
      for(vector<SQterm>::iterator t = d4c_v.begin();t != d4c_v.end();++t,++count){
	cout << count << " : " << LT_d4c_v[count] << " += " << *t << endl;      
      } // End t
      cout << endl;
    } // End if

    if(noeri.size()+eri_c.size()+eri_o.size()+eri_v.size()+d4c_c.size()+d4c_o.size()+d4c_v.size() != inTerms_.size()){
      cout << "Some terms might be missing in classification steps ... " << endl;
      abort();
    } // End if

    // Compress to produce type_terms_
    type_terms_.insert(map<string, vector<SQterm> >::value_type("noeri", noeri));
    type_terms_.insert(map<string, vector<SQterm> >::value_type("eri_c", eri_c));
    type_terms_.insert(map<string, vector<SQterm> >::value_type("eri_o", eri_o));
    type_terms_.insert(map<string, vector<SQterm> >::value_type("eri_v", eri_v));
    type_terms_.insert(map<string, vector<SQterm> >::value_type("d4c_c", d4c_c));
    type_terms_.insert(map<string, vector<SQterm> >::value_type("d4c_o", d4c_o));
    type_terms_.insert(map<string, vector<SQterm> >::value_type("d4c_v", d4c_v));

    // Compress to produce type_LTensors_
    type_LTensors_.insert(map<string, vector<SQterm> >::value_type("noeri", LT_noeri));
    type_LTensors_.insert(map<string, vector<SQterm> >::value_type("eri_c", LT_eri_c));
    type_LTensors_.insert(map<string, vector<SQterm> >::value_type("eri_o", LT_eri_o));
    type_LTensors_.insert(map<string, vector<SQterm> >::value_type("eri_v", LT_eri_v));
    type_LTensors_.insert(map<string, vector<SQterm> >::value_type("d4c_c", LT_d4c_c));
    type_LTensors_.insert(map<string, vector<SQterm> >::value_type("d4c_o", LT_d4c_o));
    type_LTensors_.insert(map<string, vector<SQterm> >::value_type("d4c_v", LT_d4c_v));

  }                                               

}} //End Femto
