//
//  SQfock.cc
//  
//
//  Created by Masaaki Saitow on 13/05/24.
//  Copyright (c) 2012 Masaaki Saitow. All rights reserved.
//

#include <algorithm>
#include <fstream>
#include <cmath>
#include <string>
#include <boost/format.hpp>
#include <SQterm.hpp>
#include <SQtensor.hpp>
#include <SQindex.hpp>
#include <SQints.hpp>

#define _FC_DEBUG
#define _FORCE_MODE

//#define _DEBUG // For debugging ..

using namespace std;

using namespace Femto::Core;

namespace Femto { namespace Core { 

  // *********************************************************
  // Function to make the core Fock matrix
  // *********************************************************
  void makeFock(vector<SQterm> &inTerms)
  {

    // If MO integrals are not of Mulliken form, transform them.
    for(auto t = inTerms.begin();t != inTerms.end();++t)
      for(size_t num_t = 0;num_t < t->get_tensors().size();++num_t)
        if(t->get_tensors()[num_t].get_name() == name_h2() || t->get_tensors()[num_t].get_name() == name_h1())
          if(t->get_tensors()[num_t].get_notation() == (notation)Dirac) 
	    t->get_tensors_ptr()[num_t]->convertD2M(); // Convert Dirac->Mulliken

    // Analyze the index dependence to process the un-linked terms
    vector<string> Flags;
    vector<SQterm> Yeffs;
    int num_replaced(0);
    for(auto t = inTerms.begin();t != inTerms.end();++t){
      vector<SQindex*>  indices(t->get_summedBody());
      vector<SQtensor>  tensors(t->get_tensors());
      vector<SQindex*>  summed;
      Ivector           counters;
      vector<SQtensor*> holders;

      for(auto i = indices.begin();i != indices.end();++i) if((*i)->get_isSummed()) summed.push_back(*i);

      for(auto i = summed.begin();i != summed.end();++i){
        int counter(0);
	SQtensor* holder;
	for(auto ten = tensors.begin();ten != tensors.end();++ten){
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
	  for(auto ten = tensors.begin();ten != tensors.end();++ten) 
	    if(!(*ten == *holders[num_i])) new_ten.push_back(*ten);
	  vector<SQindex*> new_inds(holders[num_i]->get_indices());
	  for(auto j = new_inds.begin();j != new_inds.end();){
            auto j_ptr = find(summed.begin(), summed.end(), *j);

	    if(j_ptr != summed.end()){
              if(counters[(size_t)(j_ptr-summed.begin())] == 1) j = new_inds.erase(j);
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

	  string ten_name;
          // In case of [1]
          if(holders[num_i]->get_name() == name_h1()){
            if(!new_inds.size()) ten_name = "h1_int";
            else{
              cout << " >>>>> SQfock: I cannot handle this case [0] <<<<< " << endl;
              abort();
	    } // End else
	  } // End if
          // In case of [2],[3],[4],[5]
          else if(holders[num_i]->get_name() == name_h2()){
            // [2],[3]
            if(!new_inds.size()){
              vector<SQindex*> indices(holders[num_i]->get_indices());
              auto i_ptr(find(indices.begin()+1, indices.end(), indices[0]));
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
              auto i1_ptr(find(indices.begin(), indices.end(), new_inds[0]));
              auto i2_ptr(find(indices.begin(), indices.end(), new_inds[1]));
              if((i1_ptr+1) == (i2_ptr) || (i1_ptr-1) == (i2_ptr)) ten_name = "h4_J_int";
              // If i1_ptr == i2_ptr, whether the integral is J, or K type, is not determinable.
              // In such the case, determine by the position of the dummy core indices
              else if(i1_ptr == i2_ptr){
                size_t dum1_pos(-1);
                size_t dum2_pos(-1);
                size_t pos(0);
                for(auto i = indices.begin();i != indices.end();++i,++pos)
                  if     (dum1_pos == -1 && dum2_pos == -1 && (*i)->get_isSummed() && (*i)->get_char() == core) dum1_pos = pos;
                  else if(dum1_pos != -1 && dum2_pos == -1 && (*i)->get_isSummed() && (*i)->get_char() == core) dum2_pos = pos;
                if(dum1_pos == -1 || dum2_pos == -1){
                  cout << " >>>>> SQfockock: Intergral type cannot be specified [1] <<<<< " << endl;
                  abort();
		} // End if
                if(dum2_pos-dum1_pos == 1) ten_name = "h4_J_int";
                else                       ten_name = "h5_K_int";
	      } // End else if
              else                         ten_name = "h5_K_int";
	    } // End if
            else{
              cout << " >>>>> SQfock: I don't know what to do with this, " << Yeffs.back().get_tensors()[0].get_name() << " <<<<< " << endl;
              abort();
	    } // End else
	  } // End if
          else if(holders[num_i]->get_name() == name_Fock()){
            if(!new_inds.size()) ten_name = "h6_int";
            else{
              cout << " >>>>> SQfock: I cannot handle this case <<<<< " << endl;
              abort();
	    } // End else
	  } // End if
	  else {
	    cout << " >>>>> SQfock: I don't know what to do with this, " << Yeffs.back().get_tensors()[0].get_name() << " <<<<< " << endl;
	    abort();
	  } // End else          

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
            auto my_coeff(t->get_Consts());
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
      int count(0);
      cout << " >> The linked formulas .... " << endl;
      for(vector<SQterm>::iterator t = inTerms.begin();t != inTerms.end();++t, ++count){
	cout << boost::format(" [%5d] ") % count << " " << *t << endl;
        //cout << count << " : " << LTensors_[count] << " += "<< *t << endl;
      } // End t
      cout << " >> Contents of each effective integral << " << endl;
      cout << "  [1] h1_int()      <-- " + name_h1()     + "(c1,c1)"       << endl;
      cout << "  [2] h2_J_int()    <-- " + name_h2()     + "(c1,c1,c2,c2)" << endl;
      cout << "  [3] h3_K_int()    <-- " + name_h2()     + "(c1,c2,c1,c2)" << endl;
      cout << "  [4] h4_J_int(P,Q) <-- " + name_h2()     + "(c1,c1, P, Q)" << endl;
      cout << "  [5] h5_K_int(P,Q) <-- " + name_h2()     + "(c1, P,c1, Q)" << endl;
      cout << "  [6] h6_int()      <-- " + name_Fock()   + "(c1,c1)"       << endl;
      cout << endl;
    } // End if

    // >> Definitions of the effective Fock matrices of rank 0 and 1
    // [1] Rank 0 (Fc0) : Fc      <-- 2h(c1,c1) + 2(c1,c1|c2,c2) - (c1,c2|c1,c2)
    // [2] Rank 1 (Fc1) : Fc(P,Q) <--  h( P, Q) + 2(c1,c1| P, Q) - (c1, P|c1, Q)

    // Now the unlinked integrals are replaced with the effective integrals, so let's form the core Fock matrix!
    for(auto t1 = inTerms.begin();t1 != inTerms.end();++t1){
      auto Consts(t1->get_Consts());
      vector<string> ten_names;
      pair<bool, bool> found_JK(false, false); // Whether both J and K type terms are detected
      pair<vector<SQterm>::iterator , vector<SQterm>::iterator > points_JK; // Points where they are
      //*// found_JK.first  = false;
      //*// found_JK.second = false;

      for(size_t num_t = 0;num_t < t1->get_tensors().size();++num_t) ten_names.push_back(t1->get_tensors()[num_t].get_name());

      // If there's the h1_int, there should be J and K type ints anywhere in inTerms 
      if(find(Consts.begin(), Consts.end(), "h1_int") != Consts.end()){
	for(auto t2 = inTerms.begin();t2 != inTerms.end();++t2){
	  if(t1 == t2) continue;
          auto Consts2(t2->get_Consts());
          // Found J-term!
	  if     (isFactorizable(*t1, *t2) && t1->get_numConst() == t2->get_numConst()
		  && find(Consts2.begin(), Consts2.end(), "h2_J_int") != Consts2.end()){
            if(found_JK.first){
              cout << "Something is wrong, maybe not combined yet?" << endl;
              abort();
	    } // End if
	    found_JK.first  = true;
            points_JK.first = t2;
	  } // End if
          // Found K-term!
	  else if(isFactorizable(*t1, *t2) && (-0.5) * t1->get_numConst() == t2->get_numConst()
		  && find(Consts2.begin(), Consts2.end(), "h3_K_int") != Consts2.end()){
            if(found_JK.second){
              cout << "Something is wrong, maybe not combined yet?" << endl;
              abort();
	    } // End if
	    found_JK.second  = true;
            points_JK.second = t2;
	  } // End if
	} // End t2

        if(found_JK.first && found_JK.second){
          auto c_ptr(find(Consts.begin(), Consts.end(), "h1_int"));
          Consts.erase(c_ptr);
          Consts.push_back(name_cFock0());
          t1->set_numConst(t1->get_numConst()/2);
          t1->set_Consts(Consts); // Set trace of the core Fock matrix

          if     ((size_t)(points_JK.first-inTerms.begin()) > (size_t)(points_JK.second-inTerms.begin())){
	    inTerms.erase(points_JK.first);  // Erase J
	    inTerms.erase(points_JK.second); // Erase K
	  } // End if
	  else if((size_t)(points_JK.first-inTerms.begin()) < (size_t)(points_JK.second-inTerms.begin())){
	    inTerms.erase(points_JK.second); // Erase K
	    inTerms.erase(points_JK.first);  // Erase J
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
      else if(find(ten_names.begin(), ten_names.end(), name_h1()) != ten_names.end()){
	for(auto t2 = inTerms.begin();t2 != inTerms.end();++t2){
	  if(t1 == t2) continue;
          vector<string> ten_name2;
          for(size_t num_t = 0;num_t < t2->get_tensors().size();++num_t) ten_name2.push_back(t2->get_tensors()[num_t].get_name());
          // If t2 has possibilities to have J, or K type reduced integrals
          if((find(ten_name2.begin(), ten_name2.end(), "h4_J_int") != ten_name2.end()  || 
	      find(ten_name2.begin(), ten_name2.end(), "h5_K_int") != ten_name2.end()) &&
	     t1->get_Consts() == t2->get_Consts()){
	  
	    vector<SQtensor> tens;
            for(size_t num_t = 0;num_t < t2->get_tensors().size();++num_t){ 
              if(t2->get_tensors()[num_t].get_name() != "h4_J_int" && t2->get_tensors()[num_t].get_name() != "h5_K_int") 
		tens.push_back(t2->get_tensors()[num_t]);
	      else{
                vector<SQindex*> inds(t2->get_tensors()[num_t].get_indices());
                tens.push_back(SQtensor(name_h1(), inds, h1_symm()));
                //tens.push_back(SQtensor("h", inds, h1_symm()));
	      } // End else
	    } // End num_t
            SQterm dummy(t2->get_numConst(), t2->get_Consts(), tens);
	    for(size_t num_h1 = 0;num_h1 < t1->get_tensors().size();++num_h1){
              if(t1->get_tensors()[num_h1].get_name() == name_h1()){
                if(t1->get_tensors()[num_h1].get_notation() == Femto::Dirac) break;
		else {
                  for(size_t num = 0;num < dummy.get_tensors().size();++num)
		    if(dummy.get_tensors()[num].get_name() == name_h1())
		      dummy.get_tensors_ptr()[num]->convertD2M();
		} // End else
	      } // End if
	    } // End num_h1

#ifdef _DEBUG
	    std::cout << " >> Are you ... " << *t2 << endl; //TEST	    
	    std::cout << " ====================================== " << endl;
	    std::cout << " [a] " << (isFactorizable(dummy, *t1) ? "Y" : "N") << endl; 
	    std::cout << " [b] " << ((find(ten_name2.begin(), ten_name2.end(), "h4_J_int")!=ten_name2.end()) ? "Y" : "N") << endl;
	    std::cout << " [c] " << (t2->get_numConst() == t1->get_numConst()*2 ? "Y" : "N") << endl;
	    std::cout << " [d] " << (tens == t1->get_tensors() ? "Y" : "N") << endl;
	    for(auto tt = tens.begin();tt != tens.end();++tt){
	      std::cout << "    >> " << *tt << (*tt == t1->get_tensors()[(size_t)(tt-tens.begin())] ? "Y" : "N") << " :: " << t1->get_tensors()[(size_t)(tt-tens.begin())] << endl;
	      if(!(*tt == t1->get_tensors()[(size_t)(tt-tens.begin())])) {
		std::cout << "       >> " << (tt->get_name() == t1->get_tensors()[(size_t)(tt-tens.begin())].get_name() ? "Y" : "N") << endl;
		std::cout << "       >> " << (*tt->get_indices()[0] == *t1->get_tensors()[(size_t)(tt-tens.begin())].get_indices()[0] ? "Y" : "N") << endl;
		std::cout << "       >> " << (*tt->get_indices()[1] == *t1->get_tensors()[(size_t)(tt-tens.begin())].get_indices()[1] ? "Y" : "N") << endl;
		std::cout << "       >> " << (tt->get_notation() == t1->get_tensors()[(size_t)(tt-tens.begin())].get_notation() ? "Y" : "N") << endl;
	      }
	    }
	    std::cout << " ====================================== " << endl;
	    std::cout << " >> Dummy ...  " << dummy << endl; //TEST
	    std::cout << " >> Target ... " << *t1   << endl; //TEST
	    cout << " --> ";
	    for(auto tt = tens.begin();tt != tens.end();++tt) cout << *tt << " "; //TEST
	    cout << endl << endl; //TEST
#endif

            // Find J!
            if     (isFactorizable(dummy, *t1) && (find(ten_name2.begin(), ten_name2.end(), "h4_J_int")!=ten_name2.end())
                    && t2->get_numConst() == t1->get_numConst()*2){
              if(found_JK.first){
                cout << "Something is wrong in searching J(P,Q) type integrals" << endl;
                cout << "[1] " << *t1 << endl;
                cout << "[2] " << *t2 << endl;
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
                cout << "[1] " << " += " << *t1 << endl;
                cout << "[2] " << " += " << *t2 << endl;
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
            if(t1->get_tensors()[num_t].get_name() != name_h1()) new_ten.push_back(t1->get_tensors()[num_t]);
            else new_ten.push_back(SQtensor(name_cFock1(), t1->get_tensors()[num_t].get_indices(), h1_symm()));
          t1->set_tensors(new_ten);

          if     ((size_t)(points_JK.first-inTerms.begin()) > (size_t)(points_JK.second-inTerms.begin())){
	    inTerms.erase(points_JK.first);  // Erase J
	    inTerms.erase(points_JK.second); // Erase K
	  } // End if
	  else if((size_t)(points_JK.first-inTerms.begin()) < (size_t)(points_JK.second-inTerms.begin())){
	    inTerms.erase(points_JK.second); // Erase K
	    inTerms.erase(points_JK.first);  // Erase J
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
    cout << " >> Replaced terms << " << endl;
    for(vector<SQterm>::iterator t = inTerms.begin();t != inTerms.end();++t, ++count){
      cout << boost::format(" [%5d] ") % count << " " << *t << endl;
    } // End t

    cout << endl;
    cout << " >> Contents of the core Fock matrices << " << endl;
    cout << "  [1] " + name_cFock0() + "()    <-- 2 " + name_h1() + "(c1,c1) + 2 " + name_h2() + "(c1,c1,c2,c2) - " + name_h2() + "(c1,c2,c1,c2)"<< endl;
    cout << "  [2] " + name_cFock1() + "(P,Q) <--   " + name_h1() + "( P, Q) + 2 " + name_h2() + "(c1,c1, P, Q) - " + name_h2() + "(c1, P,c1, Q)"<< endl;
    cout << endl;

  }

}} // Femto::core
