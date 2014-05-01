//
//  Replace_Fock.cc
//  
//
//  Created by Masaaki Saitow on 12/10/24.
//  Copyright (c) 2012 Masaaki Saitow. All rights reserved.
//

#include <algorithm>
#include <fstream>
#include <cmath>
#include <string>
#include <boost/format.hpp>
#include <SQreaktor.hpp>

#define _FC_DEBUG
#define _FORCE_MODE

using namespace std;
using namespace Femto;
using namespace Femto::Core;

namespace Femto { namespace Reaktor { 

  // *********************************************************
  // 
  // *********************************************************
  void SQreaktor::replace_Fock()
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
      cout << "Replace_Fock: LTensors_ are not set. Do process_kDeltas first. if already done, something is wrong." << endl;
      abort();
    } // End if

    if(!inTerms_.size()) return;

    // If ERI is not of Mulliken form, transform it
    // But I don't care about RDMs (DIY)
    for(vector<SQterm>::iterator t = inTerms_.begin();t != inTerms_.end();++t){
      for(size_t num_t = 0;num_t < t->get_tensors().size();++num_t){
        if(t->get_tensors()[num_t].get_name() == name_h2_)
          if(t->get_tensors()[num_t].get_notation() == (notation)Dirac) 
	    t->get_tensors_ptr()[num_t]->convertD2M(); // Convert Dirac->Mulliken
      } // End num_t
    } // End t

    // Analyzing the index dependence to process the un-linked terms
    vector<string> Flags;
    vector<SQterm> Yeffs;
    int num_replaced = 0;
    for(vector<SQterm>::iterator t = inTerms_.begin();t != inTerms_.end();++t){
      vector<SQindex*>  indices(t->get_summedBody());
      vector<SQtensor>  tensors(t->get_tensors());
      vector<SQindex*>  summed;
      vector<int>       counters;
      vector<SQtensor*> holders;
      for(vector<SQindex*>::iterator i = indices.begin();i != indices.end();++i){
        if((*i)->get_isSummed()) summed.push_back(*i);
      } // End i
      // Count how many tensors share a dummy index i
      for(vector<SQindex*>::iterator i = summed.begin();i != summed.end();++i){
        int counter = 0;
        SQtensor* holder;
        for(vector<SQtensor>::iterator ten = tensors.begin();ten != tensors.end();++ten){
          vector<SQindex*> inds(ten->get_indices());
          if(find(inds.begin(), inds.end(), *i) != inds.end()) { ++counter; holder = &(*ten); }
	} // End ten
        counters.push_back(counter);
        holders.push_back(holder);
      } // End i
      for(size_t num_i = 0;num_i < summed.size();++num_i){
        // If only one tensor shares the index, summed[num_i], the tensor can be reduced possibly.
        if(counters[num_i] == 1){

          SQindex* i(summed[num_i]);
          vector<SQtensor> new_ten;
          for(vector<SQtensor>::iterator ten = tensors.begin();ten != tensors.end();++ten){
            if(!(*ten == *holders[num_i])) new_ten.push_back(*ten);  
	  } // End ten
          vector<SQindex*> new_inds(holders[num_i]->get_indices());
          for(vector<SQindex*>::iterator j = new_inds.begin();j != new_inds.end();){
            vector<SQindex*>::iterator j_ptr = find(summed.begin(), summed.end(), *j);
            
            if(j_ptr != summed.end()){
              if(counters[(size_t)(j_ptr-summed.begin())] == 1) 
                j = new_inds.erase(j);
              else ++j;
	    } // End if
            else ++j;
	  } // End j

          // >> Detremine the type of the reduced tensor from these possibilities below << 
          // [1] h1()    <-- h1(c1,c1)       :: one-body int of rank 0
          // [2] h2()    <-- V2(c1,c1,c2,c2) :: two-body J-type int of rank 0
          // [3] h3()    <-- V2(c1,c2,c1,c2) :: two-body K-type int of rank 0
          // [4] h4(*,*) <-- V2(c1,c1, *, *) :: two-body J-type int of rank 1 
          // [5] h5(*,*) <-- V2(c1, *,c1, *) :: two-body K-type int of rank 1
          // [6] h6()    <-- g1(c1,c1)       :: cas-fock matrix of rank 0

	  //cout << "mamamam1" << endl;
          string ten_name;
          // In case of [1]
          if(holders[num_i]->get_name() == name_h1_){
            if(!new_inds.size()) ten_name = "h1_int";
            else{
              cout << "Replace_Fock: I cannot handle this case" << endl;
              abort();
	    } // End else
	  } // End if
          // In case of [2],[3],[4],[5]
          else if(holders[num_i]->get_name() == name_h2_){
            // [2],[3]
            if(!new_inds.size()){
              vector<SQindex*> indices(holders[num_i]->get_indices());
              vector<SQindex*>::iterator i_ptr(find(indices.begin()+1, indices.end(), indices[0]));
              if     (i_ptr == indices.begin()+1) ten_name = "h2_J_int";
              else if(i_ptr == indices.begin()+2) ten_name = "h3_K_int";
              else{
                cout << "Something is wrong in algorithms ... " << endl;
                abort();
	      } // End else
	    } // End if
            // [4],[5]
            else if(new_inds.size() == 2){
              vector<SQindex*> indices(holders[num_i]->get_indices());
              vector<SQindex*>::iterator i1_ptr(find(indices.begin(), indices.end(), new_inds[0]));
              vector<SQindex*>::iterator i2_ptr(find(indices.begin(), indices.end(), new_inds[1]));
              if((i1_ptr+1) == (i2_ptr) || (i1_ptr-1) == (i2_ptr)) ten_name = "h4_J_int";
              // If i1_ptr == i2_ptr, whether the integral is J, or K type, is not determinable.
              // In such the case, determine by the position of the dummy core indices
              else if(i1_ptr == i2_ptr){
                size_t dum1_pos(-1);
                size_t dum2_pos(-1);
                size_t pos(0);
                for(vector<SQindex*>::iterator i = indices.begin();i != indices.end();++i,++pos)
                  if     (dum1_pos == -1 && dum2_pos == -1 && (*i)->get_isSummed() && (*i)->get_char() == core) dum1_pos = pos;
                  else if(dum1_pos != -1 && dum2_pos == -1 && (*i)->get_isSummed() && (*i)->get_char() == core) dum2_pos = pos;
                if(dum1_pos == -1 || dum2_pos == -1){
                  cout << "Replace_Fock: Intergarl type cannot be specified" << endl;
                  abort();
		} // End if
                if(dum2_pos-dum1_pos == 1) ten_name = "h4_J_int";
                else                       ten_name = "h5_K_int";
	      } // End else if
              else                         ten_name = "h5_K_int";
	    } // End if
            else{
              cout << "Replace_Fock: I don't know what to do with this, " << Yeffs.back().get_tensors()[0].get_name() << endl;
              abort();
	    } // End else
	  } // End if
          else if(holders[num_i]->get_name() == Fock_name()){
            if(!new_inds.size()) ten_name = "h6_int";
            else{
              cout << "Replace_Fock: I cannot handle this case" << endl;
              abort();
	    } // End else
	  } // End if
	  else {
              cout << "Replace_Fock: I don't know what to do with this, " << Yeffs.back().get_tensors()[0].get_name() << endl;
              abort();
	  } // End else          
	  //cout << "mamamam2" << endl; //*TEST*

          ostringstream stm;
          stm << num_replaced++;
          // Currently, the symmetry of the reduced integrals (J, K ints) are given as same to h1_symm
          if(new_inds.size()){
            Symmetry my_symm;
            //if     (new_inds.size() == 2) my_symm = u2_symm();
            if     (new_inds.size() == 2) my_symm = h1_symm();
            else if(new_inds.size() == 4) my_symm = u4_symm();
            new_ten.push_back(SQtensor(ten_name, new_inds, my_symm));
	  } // End is
          else{
            vector<string> my_coeff(t->get_Consts());
            my_coeff.push_back(ten_name);
            t->set_Consts(my_coeff);
	  } // End else
          t->set_tensors(new_ten);
          t->masquerade(); // Does it work correctly ?? 2012/10/22
          break;

	} // End if
      } // End num_i      
    } // End t

    cout << endl;
    cout << " >> " << num_replaced << " terms are replaced in the linking process .... " << endl << endl;;
    if(num_replaced){
      int count = 0;
      cout << " >> The linked formulas .... " << endl;
      for(vector<SQterm>::iterator t = inTerms_.begin();t != inTerms_.end();++t, ++count){
        cout << count << " : " << LTensors_[count] << " += "<< *t << endl;
      } // End t
      cout << endl;
      cout << "++ Contents of each effective integral :: " << endl;
      cout << "[1] h1_int()      <-- " + name_h1_     + "(c1,c1)"       << endl;
      cout << "[2] h2_J_int()    <-- " + name_h2_     + "(c1,c1,c2,c2)" << endl;
      cout << "[3] h3_K_int()    <-- " + name_h2_     + "(c1,c2,c1,c2)" << endl;
      cout << "[4] h4_J_int(P,Q) <-- " + name_h2_     + "(c1,c1, P, Q)" << endl;
      cout << "[5] h5_K_int(P,Q) <-- " + name_h2_     + "(c1, P,c1, Q)" << endl;
      cout << "[6] h6_int()      <-- " + Fock_name()  + "(c1,c1)"       << endl;
      cout << endl;
    } // End if

    // >> Definitions of the effective Fock matrices of rank 0 and 1
    // [1] Rank 0 (Fc0) : Fc      <-- 2h(c1,c1) + 2(c1,c1|c2,c2) - (c1,c2|c1,c2)
    // [2] Rank 1 (Fc1) : Fc(P,Q) <--  h( P, Q) + 2(c1,c1| P, Q) - (c1, P|c1, Q)

    // Now the unlinked integrals are replaced with the effective integrals, so let's form the core Fock matrix!
    int num_t1(0);
    for(vector<SQterm>::iterator t1 = inTerms_.begin();t1 != inTerms_.end();++t1,++num_t1){
      vector<string> Consts(t1->get_Consts());
      vector<string> ten_names;
      pair<bool, bool> found_JK; // Whether both J and K type terms are detected
      pair<vector<SQterm>::iterator , vector<SQterm>::iterator > points_JK; // Points where they are
      found_JK.first  = false;
      found_JK.second = false;

      for(size_t num_t = 0;num_t < t1->get_tensors().size();++num_t)
        ten_names.push_back(t1->get_tensors()[num_t].get_name());

      // If there's the h1_int, there should be J and K type ints anywhere in inTerms_ 
      if(find(Consts.begin(), Consts.end(), "h1_int") != Consts.end()){
        int num_t2(0);
	for(vector<SQterm>::iterator t2 = inTerms_.begin();t2 != inTerms_.end();++t2,++num_t2){
          vector<string> Consts2(t2->get_Consts());
          // Found J-term!
	  if     (isFactorizable(*t1, *t2) && t1 != t2
             && t1->get_numConst() == t2->get_numConst()
		  //&& find(Consts2.begin(), Consts2.end(), "h2_J_int") != Consts.end()
	     && find(Consts2.begin(), Consts2.end(), "h2_J_int") != Consts2.end()
             && LTensors_[num_t1] == LTensors_[num_t2]){
            if(found_JK.first){
              cout << "Something is wrong, maybe not combined yet?" << endl;
              abort();
	    } // End if
	    found_JK.first  = true;
            points_JK.first = t2;
	  } // End if
          // Found K-term!
	  else if(isFactorizable(*t1, *t2) && t1 != t2
	     && (-0.5) * t1->get_numConst() == t2->get_numConst()
		  //&& find(Consts2.begin(), Consts2.end(), "h3_K_int") != Consts.end()
	     && find(Consts2.begin(), Consts2.end(), "h3_K_int") != Consts2.end()
             && LTensors_[num_t1] == LTensors_[num_t2]){
            if(found_JK.second){
              cout << "Something is wrong, maybe not combined yet?" << endl;
              abort();
	    } // End if
	    found_JK.second  = true;
            points_JK.second = t2;
	  } // End if
	} // End t2

        if(found_JK.first && found_JK.second){
          vector<string>::iterator c_ptr(find(Consts.begin(), Consts.end(), "h1_int"));
          Consts.erase(c_ptr);
          Consts.push_back("Fc0");
          t1->set_numConst(t1->get_numConst()/2);
          t1->set_Consts(Consts); // Set trace of the core Fock matrix

          if     ((size_t)(points_JK.first-inTerms_.begin()) > (size_t)(points_JK.second-inTerms_.begin())){
	    inTerms_.erase(points_JK.first);  // Erase J
	    inTerms_.erase(points_JK.second); // Erase K
            size_t offset1(points_JK.first-inTerms_.begin()); 
            size_t offset2(points_JK.second-inTerms_.begin());
            vector<SQterm>::iterator ptr_J(LTensors_.begin()+offset1);
            vector<SQterm>::iterator ptr_K(LTensors_.begin()+offset2);
            LTensors_.erase(ptr_J); // Erase J
            LTensors_.erase(ptr_K); // Erase K 
	  } // End if
	  else if((size_t)(points_JK.first-inTerms_.begin()) < (size_t)(points_JK.second-inTerms_.begin())){
	    inTerms_.erase(points_JK.second); // Erase K
	    inTerms_.erase(points_JK.first);  // Erase J
            size_t offset1(points_JK.second-inTerms_.begin()); 
            size_t offset2(points_JK.first-inTerms_.begin());
            vector<SQterm>::iterator ptr_K(LTensors_.begin()+offset1);
            vector<SQterm>::iterator ptr_J(LTensors_.begin()+offset2);
            LTensors_.erase(ptr_K); // Erase K 
            LTensors_.erase(ptr_J); // Erase J
	  } // End if
          else{
            cout << "Something wrong in elimination step of terms on formation of Fc0" << endl;
            abort();
	  } // End else

	} // End if
#ifdef _FC_DEBUG
        else{
	  cout << "Failed to find either J, or K for " << *t1 << endl;
          cout << "Found J ? : " << (found_JK.first  ? "Yes" : "No") << endl;
          cout << "Found K ? : " << (found_JK.second ? "Yes" : "No") << endl;
	} // End else
#endif
      } // End if (Seek for h1_int)
      // If there's h1 int, there should be the J and K type reduced integrals
      else if(find(ten_names.begin(), ten_names.end(), name_h1_) != ten_names.end()){
        int num_t2(0);
	for(vector<SQterm>::iterator t2 = inTerms_.begin();t2 != inTerms_.end();++t2,++num_t2){
          vector<string> ten_name2;
          for(size_t num_t = 0;num_t < t2->get_tensors().size();++num_t) 
            ten_name2.push_back(t2->get_tensors()[num_t].get_name());
          // If t2 has possibilities to have J, or K type reduced integrals
          if((find(ten_name2.begin(), ten_name2.end(), "h4_J_int") != ten_name2.end()  || 
	      find(ten_name2.begin(), ten_name2.end(), "h5_K_int") != ten_name2.end()) &&
	      t1->get_Consts() == t2->get_Consts() && t1 != t2 && LTensors_[num_t1] == LTensors_[num_t2]){

	    vector<SQtensor> tens;
            for(size_t num_t = 0;num_t < t2->get_tensors().size();++num_t){ 
              if(t2->get_tensors()[num_t].get_name() != "h4_J_int" && t2->get_tensors()[num_t].get_name() != "h5_K_int") 
		tens.push_back(t2->get_tensors()[num_t]);
	      else{
                vector<SQindex*> inds(t2->get_tensors()[num_t].get_indices());
                tens.push_back(SQtensor(name_h1_, inds, h1_symm()));
	      } // End else
	    } // End num_t
            SQterm dummy(t2->get_numConst(), t2->get_Consts(), tens);
            // Find J!
            if     (isFactorizable(dummy, *t1) && (find(ten_name2.begin(), ten_name2.end(), "h4_J_int")!=ten_name2.end())
                    && t2->get_numConst() == t1->get_numConst()*2){
              if(found_JK.first){
                cout << "Something is wrong in searching J(P,Q) type integrals" << endl;
                cout << "[1] " << LTensors_[num_t1].get_tensors()[0] << " += " << *t1 << endl;
                cout << "[2] " << LTensors_[num_t2].get_tensors()[0] << " += " << *t2 << endl;
#ifndef _FORCE_MODE
                abort();
#endif
	      }	// End if
              found_JK.first  = true;
              points_JK.first = t2;    
	    } // End if
            // Find K!
            else if(isFactorizable(dummy, *t1) && (find(ten_name2.begin(), ten_name2.end(), "h5_K_int")!=ten_name2.end())
                    && t2->get_numConst() == t1->get_numConst()*(-1)){
              if(found_JK.second){
                cout << "Something is wrong in searching K(P,Q) type integrals" << endl;
                cout << "[1] " << LTensors_[num_t1].get_tensors()[0] << " += " << *t1 << endl;
                cout << "[2] " << LTensors_[num_t2].get_tensors()[0] << " += " << *t2 << endl;
#ifndef _FORCE_MODE
                abort();
#endif
	      }	// End if
              found_JK.second  = true;
              points_JK.second = t2;    
	    } // End else (possibility check2)
	  } // End if (possibility check1)

        } // End t2
        if(found_JK.first && found_JK.second){
	  vector<SQtensor> new_ten;
          for(size_t num_t = 0;num_t < t1->get_tensors().size();++num_t)
            if(t1->get_tensors()[num_t].get_name() != name_h1_) new_ten.push_back(t1->get_tensors()[num_t]);
            else new_ten.push_back(SQtensor("Fc1", t1->get_tensors()[num_t].get_indices(), h1_symm()));
          t1->set_tensors(new_ten);

          if     ((size_t)(points_JK.first-inTerms_.begin()) > (size_t)(points_JK.second-inTerms_.begin())){
	    inTerms_.erase(points_JK.first);  // Erase J
	    inTerms_.erase(points_JK.second); // Erase K
            size_t offset1(points_JK.first-inTerms_.begin()); 
            size_t offset2(points_JK.second-inTerms_.begin());
            vector<SQterm>::iterator ptr_J(LTensors_.begin()+offset1);
            vector<SQterm>::iterator ptr_K(LTensors_.begin()+offset2);
            LTensors_.erase(ptr_J); // Erase J
            LTensors_.erase(ptr_K); // Erase K 
	  } // End if
	  else if((size_t)(points_JK.first-inTerms_.begin()) < (size_t)(points_JK.second-inTerms_.begin())){
	    inTerms_.erase(points_JK.second); // Erase K
	    inTerms_.erase(points_JK.first);  // Erase J
            size_t offset1(points_JK.second-inTerms_.begin()); 
            size_t offset2(points_JK.first-inTerms_.begin());
            vector<SQterm>::iterator ptr_K(LTensors_.begin()+offset1);
            vector<SQterm>::iterator ptr_J(LTensors_.begin()+offset2);
            LTensors_.erase(ptr_K); // Erase K 
            LTensors_.erase(ptr_J); // Erase J
	  } // End if
          else{
            cout << "Something wrong in elimination step of terms on formation of Fc1" << endl;
            abort();
	  } // End else

	} // End if
#ifdef _FC_DEBUG
        else{
	  cout << "Failed to find either J, or K for " << *t1 << endl;
          cout << "Found J ? : " << (found_JK.first  ? "Yes" : "No") << endl;
          cout << "Found K ? : " << (found_JK.second ? "Yes" : "No") << endl;
	} // End else
#endif

      } // End if (Seek for h1)

    } // End t1

    int count = 0;
    cout << "++ Replaced terms by the Fock elements :: " << endl;
    for(vector<SQterm>::iterator t = inTerms_.begin();t != inTerms_.end();++t, ++count){
      cout << count << " : " << LTensors_[count] << " += "<< *t << endl;
    } // End t

    cout << endl;
    cout << "++ Contents of the core Fock matrices :: " << endl;
    cout << "[1] Fc0()    <-- 2 " + name_h1_ + "(c1,c1) + 2 " + name_h2_ + "(c1,c1,c2,c2) - " + name_h2_ + "(c1,c2,c1,c2)"<< endl;
    cout << "[2] Fc1(P,Q) <--   " + name_h1_ + "( P, Q) + 2 " + name_h2_ + "(c1,c1, P, Q) - " + name_h2_ + "(c1, P,c1, Q)"<< endl;
    cout << endl;

  }                                               

}} //End Femto
