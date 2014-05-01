//
//  Set_Indices.cc
//  
//
//  Created by Masaaki Saitow on 13/10/25.
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
  void SQreaktor::set_and_regulate_indices(vector<vector<SQbinary> > &theBins)
  {
    if(LTensor_.get_indices().size() != 4 && isBareLHS_){
      cout << "Argument 0 cannnot be treated as a bareampack .... " << endl;
      abort();
    } // End if

    vector<SQindex*> Linds(LTensor_.get_indices());
    for(vector<SQindex*>::iterator i = Linds.begin();i != Linds.end();++i)
      if((*i)->get_isSummed()){
        cout << *i << " in 0th argument is a kind of dummy index." << endl;
        abort(); 
      } // End if

    if(is_sfGen(LTensor_.get_name())){
      cout << "Factorize: 1st argument cannot be a spin-free unitary group generator ... " << endl;
      abort();
    }

    if(!theBins.size()) return;

    // Convert term by term ....
    int numTerm(0); 
    for(auto thisBin = theBins.begin();thisBin != theBins.end();++thisBin, ++numTerm){

//--       // Achieve the best matching
//--       typedef pair<bool, SQtensor> conTen;
//--       conTen h2f, ampf;
//--       h2f.first  = false;
//--       ampf.first = false;
//-- 
//--       for(auto ts = thisBin->begin();ts != thisBin->end();++ts){
//-- 	for(size_t num_t = 0;num_t < ts->get_Rtensors().size();++num_t)
//-- 	  if      (ts->get_Rtensors()[num_t].get_name() == name_h2_){
//-- 	    h2f.first  = true;
//-- 	    h2f.second = (ts->get_Rtensors()[num_t]);
//-- 	  } // End if
//-- 	  else if (ts->get_Rtensors()[num_t].get_name() == name_amp_){
//-- 	    ampf.first  = true;
//-- 	    ampf.second = (ts->get_Rtensors()[num_t]);
//-- 	  } // End if
//--       } // End ts
//-- 
//--       // Match the external indices of the bareamp and eri, if possible.
//--       // If impossible, consider such the possibility of d4 and eri.
//--       bool SetT_V(false);
//--       bool SetL_V(false);
//--       bool SetL_T(false);
//--       bool SetV_v(false);
//-- 
//--       // LHS v.s. T2
//--       if(ampf.first && isBareLHS_){
//--         auto ampind(ampf.second.get_indices());
//--         auto LTind(thisBin->back().get_Ltensor().get_indices());
//--         for(vector<SQindex*>::iterator i = ampind.begin();i != ampind.end();++i){
//--           for(vector<SQindex*>::iterator j = LTind.begin();j != LTind.end();++j){
//--             if(**i == **j){
//--               IIvector ampperm(ampf.second.get_perms());
//--               IIvector LTperm(LTensors_[numTerm].get_tensors()[0].get_perms());
//--               typedef pair<bool, size_t> found;
//--               found amp, LT;
//--               amp.first = false;
//--               LT.first = false;
//--               amp.second = 0;
//--               LT.second = 0;
//--               for(size_t k = 0;k < ampperm.size();++k){
//--                 if(*ampind[ampperm[k][extamp_]] == **i) { amp.first = true; amp.second = k; break; } 
//--                 if(amp.first) break;
//-- 	      } //End k
//--               for(size_t k = 0;k < LTperm.size();++k){
//--                 if(*LTind[LTperm[k][extamp_]] == **j) { LT.first = true; LT.second = k; break; } 
//--                 if(LT.first) break;
//-- 	      } //End k
//-- 
//--               if(amp.first && LT.first) {
//-- 		vector<SQtensor*>::iterator amp_t;
//--                 vector<SQtensor*> ten_ptr;
//-- 		for(size_t num_t = 0;num_t < thisBin->size();++num_t)
//-- 		  for(size_t tt = 0;tt < thisBin->at(num_t).get_Rtensors_ptr().size();++tt) ten_ptr.push_back(thisBin->at(num_t).get_Rtensors_ptr()[tt]);
//-- 
//--                 for(auto t = ten_ptr.begin();t != ten_ptr.end();++t){
//--                   if     ((*t)->get_name() == name_amp_) 
//-- 		    { (*t)->rotateIndices(amp.second); amp_t = t; break; }
//-- 	        } // End t
//--                 //SQtensor* Lptr(LTensors_[numTerm].get_tensors()[0].rotateIndices(LT.second));
//--                 SQtensor* Lptr(thisBin->back().get_Ltensor_ptr());
//-- 		Lptr->rotateIndices(LT.second);
//--                 cout << ">> Indices of BareAmp are rotated to match with LHS.           -> " << numTerm << endl;
//-- 		cout << "   |--> ( " << *Lptr << "/ " << **amp_t << ")" << endl;
//--                 SetL_T = true;
//-- 	      } // End if
//-- 	      if(SetL_T) break;
//-- 	    } // End if
//-- 	  } // End j
//--           if(SetL_T) break;
//-- 	} // End i
//--       } // End if
//-- 
//--       // LHS v.s. V2
//--       if(h2f.first && isBareLHS_ && !SetL_T){
//--         vector<SQindex*> h2ind(h2f.second.get_indices());
//--         vector<SQindex*> LTind(thisBin->back().get_Ltensor().get_indices());
//-- 
//--         for(vector<SQindex*>::iterator i = h2ind.begin();i != h2ind.end();++i){
//--           for(vector<SQindex*>::iterator j = LTind.begin();j != LTind.end();++j){
//--             if(**i == **j){
//--               IIvector h2perm(h2f.second.get_perms());
//--               IIvector LTperm(LTensors_[numTerm].get_tensors()[0].get_perms());
//--               typedef pair<bool, size_t> found;
//--               found h2, LT;
//--               h2.first = false;
//--               LT.first = false;
//--               h2.second = 0;
//--               LT.second = 0;
//--               for(size_t k = 0;k < h2perm.size();++k){
//--                 if(*h2ind[h2perm[k][exth2_]] == **i) { h2.first = true; h2.second = k; break; } 
//--                 if(h2.first) break;
//-- 	      } //End k
//--               for(size_t k = 0;k < LTperm.size();++k){
//--                 if(*LTind[LTperm[k][extamp_]] == **j) { LT.first = true; LT.second = k; break; } 
//--                 if(LT.first) break;
//-- 	      } //End k
//-- 
//-- 	      SQtensor* ERI(NULL);
//--               if(h2.first && LT.first) {
//--                 vector<SQtensor*> ten_ptr(thisTerm->get_tensors_ptr()); 
//--                 for(vector<SQtensor*>::iterator t = ten_ptr.begin();t != ten_ptr.end();++t){
//--                   if     ((*t)->get_name() == name_h2_) 
//--                     { (*t)->rotateIndices(h2.second); ERI = (*t); }
//-- 	        } // End t
//--                 SQtensor* Lptr(LTensors_[numTerm].get_tensors_ptr()[0]);
//-- 		Lptr->rotateIndices(LT.second);
//--                 cout << ">> Indices of ERI are rotated to match with LHS.               -> " << numTerm << endl;
//-- 		cout << "   |--> ( " << *Lptr << "/ " << *ERI << ")" << endl;
//--                 SetL_V = true;
//-- 	      } // End if
//--               if(SetL_V) break;
//-- 	    } // End if
//-- 	  } // End j
//--           if(SetL_V) break;
//-- 	} // End i
//--       } // End if
//-- 
//--       // T2 v.s. V2
//--       if(h2f.first && ampf.first && !SetL_V || SetL_T){
//--         auto h2ind(h2f.second.get_indices());
//--         auto ampind(ampf.second.get_indices());
//-- 
//--         for(vector<SQindex*>::iterator i = h2ind.begin();i != h2ind.end();++i){
//--           for(vector<SQindex*>::iterator j = ampind.begin();j != ampind.end();++j){
//--             if(**i == **j){
//--               IIvector h2perm(h2f.second.get_perms());
//--               IIvector ampperm(ampf.second.get_perms());
//--               typedef pair<bool, size_t> found;
//--               found h2, amp;
//--               h2.first = false;
//--               amp.first = false;
//--               h2.second = 0;
//--               amp.second = 0;
//--               for(size_t k = 0;k < h2perm.size();++k){
//--                 if(*h2ind[h2perm[k][exth2_]] == **i) { h2.first = true; h2.second = k; break; } 
//--                 if(h2.first) break;
//-- 	      } //End k
//--               for(size_t k = 0;k < ampperm.size();++k){
//--                 if(*ampind[ampperm[k][extamp_]] == **j) { amp.first = true; amp.second = k; break; } 
//--                 if(amp.first) break;
//-- 	      } //End k
//-- 
//--               if(h2.first && amp.first) {
//--                 vector<SQtensor*> ten_ptr;
//-- 		for(size_t num_t = 0;num_t < thisBin->size();++num_t)
//-- 		  for(size_t tt = 0;tt < thisBin->at(num_t).get_Rtensors_ptr().size();++tt) ten_ptr.push_back(thisBin->at(num_t).get_Rtensors_ptr()[tt]);
//-- 
//--                 for(auto t = ten_ptr.begin();t != ten_ptr.end();++t){
//--                   if     ((*t)->get_name() == name_h2_) 
//--                     (*t)->rotateIndices(h2.second);
//--                   else if((*t)->get_name() == name_amp_) 
//--                     (*t)->rotateIndices(amp.second);
//-- 	        } // End t
//--                 cout << ">> Indices of ERI are rotated to match with Bareamp.           -> " << numTerm << endl;
//--                 SetT_V = true;
//-- 	      } // End if
//-- 	      if(SetT_V) break;
//-- 	    } // End if
//-- 	  } // End j
//--           if(SetT_V) break;
//-- 	} // End i
//--       } // End if
//-- 
//--       // Set loading index of ERI to become virtual if possible
//--       //if(!SetT_V && !SetL_V && /*!SetV_D &&*/ h2f.first){ // OLD setting (2013/01/08). This may be better
//--       if(!SetT_V && !SetL_V && h2f.first){ // <-- NEW 2013/01/08
//--         vector<SQindex*> h2ind(h2f.second.get_indices());
//--         for(size_t num_i = 0;num_i < h2f.second.get_indices().size();++num_i){
//--           SQindex* i = h2f.second.get_indices()[num_i];
//--           if(i->get_char() == Femto::virt){
//--             IIvector h2perm(h2f.second.get_perms());
//--             typedef pair<bool, size_t> found;
//--             found h2;
//--             h2.first = false;
//--             h2.second = 0;
//--             for(size_t k = 0;k < h2perm.size();++k){
//--               if(*h2ind[h2perm[k][exth2_]] == *i) { h2.first = true; h2.second = k; break; } 
//-- 	    } // End k
//-- 
//--             if(h2.first){
//-- 	      vector<SQtensor*> ten_ptr;
//-- 	      for(size_t num_t = 0;num_t < thisBin->size();++num_t)
//-- 		for(size_t tt = 0;tt < thisBin->at(num_t).get_Rtensors_ptr().size();++tt) ten_ptr.push_back(thisBin->at(num_t).get_Rtensors_ptr()[tt]);
//-- 	      
//--               for(vector<SQtensor*>::iterator t = ten_ptr.begin();t != ten_ptr.end();++t){
//--                 if((*t)->get_name() == name_h2_) (*t)->rotateIndices(h2.second);
//-- 	      } // End t
//--               cout << ">> The indices of ERI are rotated to became virtual.           -> " << numTerm << endl;
//--               SetV_v = true;              
//-- 	    } // End if
//-- 	  } // End if
//--           if(SetV_v) break;
//-- 	} // End num_i
//--       } // End if
    
      /////////////////////////////////////////////////
      // Set the loading indices as external
      /////////////////////////////////////////////////
      vector<SQindex*> inds_ptr;
      for(auto ts = thisBin->begin();ts != thisBin->end();++ts){
	vector<SQindex*> temp(ts->get_summedBody());
	inds_ptr.insert(inds_ptr.end(), temp.begin(), temp.end());
      } // End ts
      for(auto i = inds_ptr.begin();i != inds_ptr.end();++i) (*i)->switch_isExt(false);

      vector<SQindex> extInds;
      for(auto ts = thisBin->begin();ts != thisBin->end();++ts){
	for(size_t num_t = 0;num_t < ts->get_Rtensors().size();++num_t){
	  if     (ts->get_Rtensors()[num_t].get_name() == name_h2_)  extInds.push_back(*ts->get_Rtensors()[num_t].get_indices()[exth2_]);
	  else if(ts->get_Rtensors()[num_t].get_name() == name_amp_) extInds.push_back(*ts->get_Rtensors()[num_t].get_indices()[extamp_]);
	} // End size_t
	if(isBareLHS_ && ts->get_Ltensor().get_name() == LTensor_.get_name()) extInds.push_back(*(ts->get_Ltensor().get_indices()[extamp_]));
      } // End ts
 
      for(auto ts = thisBin->begin();ts != thisBin->end();++ts){
	auto tempInds(ts->get_summedBody());
	for(auto i = tempInds.begin();i != tempInds.end();++i) 
	  for(auto I = extInds.begin();I != extInds.end();++I)
	    if((**i) == (*I)) (*i)->switch_isExt(true);
      } // End ts

    } // End thisBins

    int count = 0;
    cout << endl;
    cout << "++ Indices are rotated to achive the best matching :: " << endl;
    for(auto t = theBins.begin();t != theBins.end();++t, ++count){
      SQcont<SQbinary> temp(*t);
      cout << temp << endl << endl;

//       cout << "[1]" << endl; 
//       temp[0].print_summedBody();
//       cout << "[2]" << endl; 
//       temp[1].print_summedBody();

    } // End t
    cout << endl;

    // Divide each term in inTerms_ into type_terms_
    vector<vector<SQbinary> > noeri; noeri.reserve(theBins.size());
    vector<vector<SQbinary> > eri_c; eri_c.reserve(theBins.size());
    vector<vector<SQbinary> > eri_o; eri_o.reserve(theBins.size());
    vector<vector<SQbinary> > eri_v; eri_v.reserve(theBins.size());
    vector<vector<SQbinary> > d4c_c; d4c_c.reserve(theBins.size());
    vector<vector<SQbinary> > d4c_o; d4c_o.reserve(theBins.size());
    vector<vector<SQbinary> > d4c_v; d4c_v.reserve(theBins.size());

    int num(0);
    for(auto t = theBins.begin();t != theBins.end();++t,++num){

      int flag(0);
      vector<SQtensor> myTens;
      for(size_t num = 0;num < t->size();++num){
	vector<SQtensor> temp(t->at(num).get_Rtensors());
	myTens.insert(myTens.end(), temp.begin(), temp.end());
      } // End num

      // 0 : No ERI
      // 1 : ERI core
      // 2 : ERI active
      // 3 : ERI virtual
      // 4 : D4C core
      // 5 : D4C active
      // 6 : D4C virtual
      for(auto myt = myTens.begin(); myt != myTens.end();++myt){
	if(myt->get_name() == name_h2_)
	  if     (myt->get_indices()[exth2_]->get_char() == core) flag = 1;
	  else if(myt->get_indices()[exth2_]->get_char() == act ) flag = 2;
	  else if(myt->get_indices()[exth2_]->get_char() == virt) flag = 3;
	if(is_D4C(myt->get_name()))
	  if     (myt->get_indices()[extd4c_]->get_char() == core) flag = 4;
	  else if(myt->get_indices()[extd4c_]->get_char() == act ) flag = 5;
	  else if(myt->get_indices()[extd4c_]->get_char() == virt) flag = 6;
      } // myt

      if     (flag == 0){
        noeri.push_back(*t);
      } // End if
      else if(flag == 1){
        eri_c.push_back(*t);
      } // End if
      else if(flag == 2){
        eri_o.push_back(*t);
      } // End if
      else if(flag == 3){
        eri_v.push_back(*t);
      } // End if
      else if(flag == 4){
        d4c_c.push_back(*t);
      } // End if
      else if(flag == 5){
        d4c_o.push_back(*t);
      } // End if
      else if(flag == 6){
        d4c_v.push_back(*t);
      } // End if
    } // End t

    cout << "++ Terms are classified into these groups :: " << endl;
    if(noeri.size()){
      cout << ">> No ERIs <<" << endl;
      int count = 0;
      for(auto t = noeri.begin();t != noeri.end();++t,++count){
	SQcont<SQbinary> temp(*t);
	cout << "<-- " << count << " -->" << endl;
	cout << temp << endl << endl;      
      } // End t
      cout << endl;
    } // End if

    if(eri_c.size()){
      cout << ">> ERIs loaded by core indices <<" << endl;
      int count = 0;
      for(auto t = eri_c.begin();t != eri_c.end();++t,++count){
	SQcont<SQbinary> temp(*t);
	cout << "<-- " << count << " -->" << endl;
	cout << temp << endl << endl;      
      } // End t
      cout << endl;
    } // End if

    if(eri_o.size()){
      cout << ">> ERIs loaded by active indices <<" << endl;
      int count = 0;
      for(auto t = eri_o.begin();t != eri_o.end();++t,++count){
	SQcont<SQbinary> temp(*t);
	cout << "<-- " << count << " -->" << endl;
	cout << temp << endl << endl;      
      } // End t
      cout << endl;
    } // End if

    if(eri_v.size()){
      cout << ">> ERIs loaded by virtual indices <<" << endl;
      int count = 0;
      for(auto t = eri_v.begin();t != eri_v.end();++t,++count){
	SQcont<SQbinary> temp(*t);
	cout << "<-- " << count << " -->" << endl;
	cout << temp << endl << endl;      
      } // End t
      cout << endl;
    } // End if

    if(d4c_c.size()){
      cout << ">> D4Cs loaded by core indices <<" << endl;
      int count = 0;
      for(auto t = d4c_c.begin();t != d4c_c.end();++t,++count){
	SQcont<SQbinary> temp(*t);
	cout << "<-- " << count << " -->" << endl;
	cout << temp << endl << endl;      
      } // End t
      cout << endl;
    } // End if

    if(d4c_o.size()){
      cout << ">> D4Cs loaded by active indices <<" << endl;
      int count = 0;
      for(auto t = d4c_o.begin();t != d4c_o.end();++t,++count){
	SQcont<SQbinary> temp(*t);
	cout << "<-- " << count << " -->" << endl;
	cout << temp << endl << endl;      
      } // End t
      cout << endl;
    } // End if

    if(d4c_v.size()){
      cout << ">> D4Cs loaded by virtual indices <<" << endl;
      int count = 0;
      for(auto t = d4c_v.begin();t != d4c_v.end();++t,++count){
	SQcont<SQbinary> temp(*t);
	cout << "<-- " << count << " -->" << endl;
	cout << temp << endl << endl;      
      } // End t
      cout << endl;
    } // End if

    if(noeri.size()+eri_c.size()+eri_o.size()+eri_v.size()+d4c_c.size()+d4c_o.size()+d4c_v.size() != theBins.size()){
      cout << "Some terms might be missing in classification steps ... " << endl;
      abort();
    } // End if

    // Compress to produce type_terms_
    FafBinaries_.insert(map<string, vector<vector<SQbinary> > >::value_type("noeri", noeri));
    FafBinaries_.insert(map<string, vector<vector<SQbinary> > >::value_type("eri_c", eri_c));
    FafBinaries_.insert(map<string, vector<vector<SQbinary> > >::value_type("eri_o", eri_o));
    FafBinaries_.insert(map<string, vector<vector<SQbinary> > >::value_type("eri_v", eri_v));
    FafBinaries_.insert(map<string, vector<vector<SQbinary> > >::value_type("d4c_c", d4c_c));
    FafBinaries_.insert(map<string, vector<vector<SQbinary> > >::value_type("d4c_o", d4c_o));
    FafBinaries_.insert(map<string, vector<vector<SQbinary> > >::value_type("d4c_v", d4c_v));

  }                                               

}} //End Femto
