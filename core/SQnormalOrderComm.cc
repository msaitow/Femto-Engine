//
//  SQnormalOrderComm.cc
//  
//
//  Created by Masaaki Saitow on 12/10/31.
//  Copyright (c) 2012 Masaaki Saitow. All rights reserved.
//

#include <vector>
#include <algorithm>
#include <cstdlib>
#include <Femto.hpp>
#include <SQterm.hpp>
#include <SQtensor.hpp>
#include <SQindex.hpp>

#define _ZERO_CONT // Flag to avoid the operator contraction in terms of spin-free generators to give vanishing result 

using namespace std;

namespace Femto { namespace Core {

  // *****************************************************************
  // Code to achieve normal ordering in terms of multiple commutators
  // *****************************************************************
  // -----------------------------------------------------------------
  // *NOTE* :: Only spin-free generators in input are picked up first. 
  //           Then, they are normal ordered from right pairs, i.e.
  //
  //        inTerms = { ... En, En+1, En+2}  ----
  //                                            |
  //                                            /  Translated like this
  //           [ ... [En, [En+1, En+2]]] <------ 
  //                  ^      ^     ^
  //                  |      |     |
  //                  |      -------
  //                  |         | Contract first to prooduce new En+1
  //                  |         V
  //                  |         E'n+1
  //                  |         ^
  //                  |         |
  //                  -----------
  //                       | Contract next to produce new En
  //                       V
  //                       E'n
  //           and so forth.
  // -----------------------------------------------------------------
  void normalOrderComm(vector<SQterm> *inTerms, const bool OF_flag)
  {
    // Now, only one term is allowed as an input
    if(inTerms->size()!=1){
      cout << "Now, only one term is allowed as an input";
      abort();
    }

    // Now, only spin-free generator is habdlable
    for(size_t num_t = 0;num_t < inTerms->at(0).get_tensors().size();++num_t)
      if(inTerms->at(0).get_tensors()[num_t].get_name() == aCre_name() || inTerms->at(0).get_tensors()[num_t].get_name() == aDes_name()){
	cout << "Now, only spin-free generator can be handlable in SQnormalOrderComm" << endl;
	abort();
      }

    // If there's other indices than core, active and virtual, abort
    for(size_t num = 0;num < inTerms->at(0).get_summedBody().size();++num)
      if(inTerms->at(0).get_summedBody()[num]->get_char() != core && 
         inTerms->at(0).get_summedBody()[num]->get_char() != act  &&
         inTerms->at(0).get_summedBody()[num]->get_char() != virt){
        cout << "normalOrderComm: There's at least one index that doesn't belong to either of core, active or virtual orbital" << endl;
        abort();
      } // End if

    // Initialize number n, which is a number of the ramining sfGen
    if((*inTerms)[0].get_isInCanonical()) return;

    // If inTerm contains tensor in Mulliken notation, transform it to Dirac notation
    for(vector<SQterm>::iterator t = inTerms->begin();t != inTerms->end();++t){
      for(size_t num_t = 0;num_t < t->get_tensors().size();++num_t)
        if(t->get_tensors()[num_t].get_notation() == Mulliken){
          t->get_tensors()[num_t].convertM2D(); // Convert 2 Dirac notation
        } // End if
    } // End t

    // Make separate list of the sfGen and the others
    vector<SQtensor> sfGen_list;
    vector<SQtensor> other_list;

    vector<SQtensor> inTensors((*inTerms)[0].get_tensors());
    for(size_t i = 0;i < inTensors.size();++i){
      string t_name = inTensors[i].get_name();
      if(is_sfGen(t_name)) {
        sfGen_list.push_back(inTensors[i]);
      }
      else
        other_list.push_back(inTensors[i]);
    } // End i
    int n = (int)sfGen_list.size();

    // If members of sfGen_list give the zero contribution (in case that core and virtual indices appear different times on upper and lower cases),
    // no need to calculate such the contraction straightforwardly. Because it's just nothing.
#ifdef _ZERO_CONT
    typedef pair<int,int> cv_pair;
    cv_pair upp_cv; // <num core, num virt> for upper side
    cv_pair low_cv; // <num core, num virt> for lower side
    upp_cv.first = 0; upp_cv.second = 0;
    low_cv.first = 0; low_cv.second = 0;
    for(vector<SQtensor>::const_iterator e1 = sfGen_list.begin();e1 != sfGen_list.end();++e1){
      // Seek in upper side indices
      for(size_t num_i = 0;num_i < e1->get_indices().size()/2;++num_i)
	if     (e1->get_indices()[num_i]->get_char() == core) ++upp_cv.first;
	else if(e1->get_indices()[num_i]->get_char() == virt) ++upp_cv.second;
      // Seek in lower side indices
      for(size_t num_i = e1->get_indices().size()/2;num_i < e1->get_indices().size();++num_i)
	if     (e1->get_indices()[num_i]->get_char() == core) ++low_cv.first;
	else if(e1->get_indices()[num_i]->get_char() == virt) ++low_cv.second;
    } // End e1
    if(upp_cv.first != low_cv.first || upp_cv.second != low_cv.second) { 
      inTerms->erase(inTerms->begin());
      return;
    }
#endif    

    // Estimate the number of contractions arises from each step in the normal ordering
    int m = n;
    int totnum = 0;
    vector<int> orders;
    for(vector<SQtensor>::iterator gen = sfGen_list.begin();gen != sfGen_list.end();++gen)
      orders.push_back((gen->get_indices().size())/2);
    IIvector init_orders; init_orders.push_back(orders);
    IIvector out_orders;
    while(m > 1){
      out_orders.clear();
      for(vector<vector<int> >::iterator o = init_orders.begin();o != init_orders.end();++o){
        vector<int>::iterator o0(o->begin());
        vector<int>::iterator o2(o->end()); --o2;
        vector<int>::iterator o1(o2);       --o1;
        vector<int> except_last_two;
        except_last_two.insert(except_last_two.end(), o0, o1);
        // Count numbers of terms of order of newOrder
        for(int nc = 1;nc < min(*o1,*o2)+1;++nc){
          int newOrder = *o1 + *o2 - nc; 
          int nterm    = 1;
          if(nc!=0) {
            for(int i = 0;i < nc;++i) nterm = nterm * (*o1-i) * (*o2-i) / fact(i+1);
          } // End if
          nterm = nterm * 2; // Count - E2 E1 terms
          vector<int> o_list(except_last_two);
          o_list.push_back(newOrder); // + E1E1
          out_orders.insert(out_orders.end(), nterm, o_list);
        }        
      } // End o
      init_orders = out_orders;
      totnum += (int)init_orders.size();
      --m; 
    } // End while

    // Now we've already got number of contractions .... 
    vector<SQtensor> tensors(other_list);
    tensors.insert(tensors.end(), sfGen_list.begin(), sfGen_list.end());
    SQterm term1((*inTerms)[0].get_numConst(), (*inTerms)[0].get_Consts(), tensors);
    inTerms->reserve(init_orders.size());
    vector<SQterm> iter_terms; iter_terms.reserve(totnum);
    // Push back the conditioned term
    iter_terms.push_back(term1);
    // And erase the old one
    inTerms->erase(inTerms->begin());

   // Successively normal order the last two excitation operators until each term
   // has only one exitation operator left (at which point the term is normal ordered)

   // For each of this iteration's input terms, produce all terms resulting from normal ordering
   // the last two excitation operators
   size_t Count = 0;
   for(vector<SQterm>::iterator t = iter_terms.begin();t != iter_terms.end();++t){
     if(t->get_isInCanonical()) break;
     else ++Count;    
     
     //--* Enter the first product in the commutator, + E1 E2 *--//
     { 
     // Make a list of the term's tensors that excludes the last two excitation operators
     vector<SQtensor> tensors_except_last_two;
     vector<SQtensor> t_temp = t->get_tensors();
     vector<SQtensor>::iterator it    = t_temp.begin();
     vector<SQtensor>::iterator it_m1 = t_temp.end(); 
     vector<SQtensor>::iterator it_m2 = t_temp.end(); 
     --it_m1;
     --it_m2; --it_m2;
     tensors_except_last_two.insert(tensors_except_last_two.end(), it, it_m2);
     
     //  Give short names for the last two excitation operators and their orders
     SQtensor e1(*it_m2);
     SQtensor e2(*it_m1);
     int o1((int)(e1.get_indices().size()/2)); // cout << o1 << endl; //*TEST*
     int o2((int)(e2.get_indices().size()/2)); // cout << o2 << endl; //*TEST*

     // Loop over the number of contraction
     for(int nc = 1;nc < min(o1,o2)+1;++nc){

       int newOrder(o1 + o2 - nc);
       
       // Compute all nc-tuples of index numbers from e1 and e2, as well as all permutations of 
       // the order in which the tuples may be combined to form a contraction
       Ivector o1_list; o1_list.reserve(o1);
       Ivector o2_list; o2_list.reserve(o2);
       for(int i = 0;i < o1;++i) o1_list.push_back(i);
       for(int i = 0;i < o2;++i) o2_list.push_back(i);

       IIvector perms; 
       Ivector temp;
       temp.push_back(0);
       perms.push_back(temp);
       if(nc > 0)
         perms = makePermutations(nc);
       IIvector tups1(makeTuples2(nc, o1_list));
       IIvector tups2(makeTuples2(nc, o2_list));

       // For each contractions, compute the resulting term
       for(IIvector::iterator perm = perms.begin();perm != perms.end();++perm){
         for(IIvector::iterator tup1 = tups1.begin();tup1 != tups1.end();++tup1){
           for(IIvector::iterator tup2 = tups2.begin();tup2 != tups2.end();++tup2){
             // Initialize the term's tensor list
             vector<SQtensor> tensorList; 
             tensorList.insert(tensorList.end(), tensors_except_last_two.begin(),
                                                 tensors_except_last_two.end());

             bool isDead = false; // Whether kDelta vanises or not
             // Compute the pairs of indices to be contracted
             // Example: (conPairs[0][p], conPairs[1][p]) is the p-th contraction pair.
             IIvector conPairs;
             Ivector p; 
             conPairs.push_back(p);
             conPairs.push_back(p);
             for(int i = 0;i < nc;++i){
               conPairs[0].push_back((*tup1)[(*perm)[i]]);
               conPairs[1].push_back((*tup2)[i]);
             } // End i

             // Initialize the index list for the new excitation operator
             vector<SQindex*> indexList; indexList.reserve(2*newOrder);
             SQindex *q = NULL;
             for(int i = 0;i < 2*newOrder;++i) indexList.push_back(q);

             // Populate the index list for the new excitation operator
             // Also, create a kronecker delta function for each contraction pair
             vector<SQindex*> e1_indices(e1.get_indices());
             vector<SQindex*> e2_indices(e2.get_indices());

             for(int i = 0;i < o1;++i){
               indexList[i] = e1_indices[i];
	       Ivector::iterator i_index(find(conPairs[0].begin(), conPairs[0].end(), i));

               if(i_index != conPairs[0].end()){
                 int i1 = i + o1;
                 int i2 = conPairs[1][(size_t)(i_index-conPairs[0].begin())];
                 SQindex* ind1 = e1_indices[i1]; 
                 SQindex* ind2 = e2_indices[i2]; 

                 vector<SQindex*> ind_delta;
                 ind_delta.push_back(ind1);
                 ind_delta.push_back(ind2);
                 if(ind_delta[0]->get_char()!=ind_delta[1]->get_char()) isDead = true; 
                 tensorList.push_back(kDelta(ind_delta)); 
                 indexList[i+newOrder]  = e2_indices[o2+i2]; 

               } // End if
               else
                 indexList[i+newOrder]  = e1_indices[i+o1];

             } // End i
             int count = 0;
             for(int i = 0;i < o2;++i){

	       Ivector::iterator i_index(find(conPairs[1].begin(), conPairs[1].end(), i));

	       if(i_index==conPairs[1].end()){
                 indexList[o1+count]          = e2_indices[i];
                 indexList[o1+count+newOrder] = e2_indices[i+o2];
		 ++count;
               } // End if
             } // End i

             // Ensure that all slots in the index list have been filled
             // And check whether virtual exists or not
             for(size_t numInd = 0;numInd < indexList.size();++numInd){
               if(indexList[numInd] == NULL){
                 cout << "There is at least one unassigned index in the new spin-free unitary group generators" << endl;
                 abort();
               }
             } // End numInd
             // Add the new unitary generator into the tensorList
             tensorList.push_back(sfGen(indexList));

             // Add the resulting term to this iteration's list of output terms
             if(OF_flag){
               if(!isDead){
                 int Ecount=0;
                 vector<SQtensor>::iterator ten_E;
                 for(vector<SQtensor>::iterator ten = tensorList.begin();ten != tensorList.end();++ten){
                   if(is_sfGen(ten->get_name())) { ++Ecount; ten_E = ten; }
                 } // End ten
                 bool isVanished = false;
                 if(Ecount==0 || Ecount==1) {
                   pair<int, int> numCore; // Represents the number of core indices appear on the upper and lower case of the sfGen
                   numCore.first = 0; numCore.second = 0;
                   vector<SQindex*> ind(ten_E->get_indices());
                   for(vector<SQindex*>::iterator I = ind.begin();I != ind.end();++I){
                     if((*I)->get_char()==(char_state)2) isVanished = true;
                     else if((*I)->get_char()==(char_state)0 && (size_t)(I-ind.begin()) <  (ten_E->get_indices().size())/2) 
                       ++numCore.first;
                     else if((*I)->get_char()==(char_state)0 && (size_t)(I-ind.begin()) >= (ten_E->get_indices().size())/2) 
                       ++numCore.second;
                     if(isVanished) break;
		   } // End I
                   if(numCore.first != numCore.second) isVanished = true;
                 } // End if
                 if(!isVanished)
                   { iter_terms.push_back(SQterm(t->get_numConst(), term1.get_Consts(), tensorList)); }
	       } // End if
             } // End if
             else
               iter_terms.push_back(SQterm(t->get_numConst(), term1.get_Consts(), tensorList));
           } // End tup2
         } // End tup1
       } // End perm          
     } // End nc
     } // End first

     //--*  Enter the second product in the commutator, - E2 E1 *--//
     {
     // Make a list of the term's tensors that excludes the last two excitation operators
     vector<SQtensor> tensors_except_last_two;
     vector<SQtensor> t_temp = t->get_tensors();
     vector<SQtensor>::iterator it    = t_temp.begin();
     vector<SQtensor>::iterator it_m1 = t_temp.end(); 
     vector<SQtensor>::iterator it_m2 = t_temp.end(); 
     --it_m1;
     --it_m2; --it_m2;
     tensors_except_last_two.insert(tensors_except_last_two.end(), it, it_m2);
     
     //  Give short names for the last two excitation operators and their orders
     SQtensor e1(*it_m1);
     SQtensor e2(*it_m2);
     int o1((int)(e1.get_indices().size()/2)); // cout << o1 << endl; //*TEST*
     int o2((int)(e2.get_indices().size()/2)); // cout << o2 << endl; //*TEST*

     // Loop over the number of contraction
     for(int nc = 1;nc < min(o1,o2)+1;++nc){

       int newOrder(o1 + o2 - nc);
       
       // Compute all nc-tuples of index numbers from e1 and e2, as well as all permutations of 
       // the order in which the tuples may be combined to form a contraction
       Ivector o1_list; o1_list.reserve(o1);
       Ivector o2_list; o2_list.reserve(o2);
       for(int i = 0;i < o1;++i) o1_list.push_back(i);
       for(int i = 0;i < o2;++i) o2_list.push_back(i);

       IIvector perms; 
       Ivector temp;
       temp.push_back(0);
       perms.push_back(temp);
       if(nc > 0)
         perms = makePermutations(nc);
       IIvector tups1(makeTuples2(nc, o1_list));
       IIvector tups2(makeTuples2(nc, o2_list));

       // For each contractions, compute the resulting term
       for(IIvector::iterator perm = perms.begin();perm != perms.end();++perm){
         for(IIvector::iterator tup1 = tups1.begin();tup1 != tups1.end();++tup1){
           for(IIvector::iterator tup2 = tups2.begin();tup2 != tups2.end();++tup2){
             // Initialize the term's tensor list
             vector<SQtensor> tensorList; 
             tensorList.insert(tensorList.end(), tensors_except_last_two.begin(),
                                                 tensors_except_last_two.end());

             bool isDead = false; // Whether kDelta vanises or not
             // Compute the pairs of indices to be contracted
             // Example: (conPairs[0][p], conPairs[1][p]) is the p-th contraction pair.
             IIvector conPairs;
             Ivector p; 
             conPairs.push_back(p);
             conPairs.push_back(p);
             for(int i = 0;i < nc;++i){
               conPairs[0].push_back((*tup1)[(*perm)[i]]);
               conPairs[1].push_back((*tup2)[i]);
             } // End i

             // Initialize the index list for the new excitation operator
             vector<SQindex*> indexList; indexList.reserve(2*newOrder);
             SQindex *q(NULL);
             for(int i = 0;i < 2*newOrder;++i) indexList.push_back(q);

             // Populate the index list for the new excitation operator
             // Also, create a kronecker delta function for each contraction pair
             vector<SQindex*> e1_indices(e1.get_indices());
             vector<SQindex*> e2_indices(e2.get_indices());

             for(int i = 0;i < o1;++i){
               indexList[i] = e1_indices[i];
	       Ivector::iterator i_index(find(conPairs[0].begin(), conPairs[0].end(), i));

               if(i_index != conPairs[0].end()){
                 int i1 = i + o1;
                 int i2 = conPairs[1][(size_t)(i_index-conPairs[0].begin())];
                 SQindex* ind1(e1_indices[i1]); 
		 SQindex* ind2(e2_indices[i2]); 

                 vector<SQindex*> ind_delta;
                 ind_delta.push_back(ind1);
                 ind_delta.push_back(ind2);
                 if(ind_delta[0]->get_char()!=ind_delta[1]->get_char()) isDead = true; 
                 tensorList.push_back(kDelta(ind_delta)); 
                 indexList[i+newOrder]  = e2_indices[o2+i2]; 

               } // End if
               else
                 indexList[i+newOrder]  = e1_indices[i+o1];

             } // End i
             int count = 0;
             for(int i = 0;i < o2;++i){

	       Ivector::iterator i_index(find(conPairs[1].begin(), conPairs[1].end(), i));

	       if(i_index==conPairs[1].end()){
                 indexList[o1+count]          = e2_indices[i];
                 indexList[o1+count+newOrder] = e2_indices[i+o2];
		 ++count;
               } // End if
             } // End i

             // Ensure that all slots in the index list have been filled
             // And check whether virtual exists or not
             for(size_t numInd = 0;numInd < indexList.size();++numInd){
               if(indexList[numInd] == NULL){
                 cout << "There is at least one unassigned index in the new spin-free unitary group generators" << endl;
                 abort();
               }
             } // End numInd
             // Add the new unitary generator into the tensorList
             tensorList.push_back(sfGen(indexList));

             // Add the resulting term to this iteration's list of output terms
             if(OF_flag){
               if(!isDead){
                 int Ecount=0;
                 vector<SQtensor>::iterator ten_E;
                 for(vector<SQtensor>::iterator ten = tensorList.begin();ten != tensorList.end();++ten){
                   if(is_sfGen(ten->get_name())) { ++Ecount; ten_E = ten; }
                 } // End ten
                 bool isVanished = false;
                 if(Ecount==0 || Ecount==1) {
                   pair<int, int> numCore; // Represents the number of core indices appear on the upper and lower case of the sfGen
                   numCore.first = 0; numCore.second = 0;
                   vector<SQindex*> ind(ten_E->get_indices());
                   for(vector<SQindex*>::iterator I = ind.begin();I != ind.end();++I){
                     if((*I)->get_char()==(char_state)2) isVanished = true;
                     else if((*I)->get_char()==(char_state)0 && (size_t)(I-ind.begin()) <  (ten_E->get_indices().size())/2) 
                       ++numCore.first;
                     else if((*I)->get_char()==(char_state)0 && (size_t)(I-ind.begin()) >= (ten_E->get_indices().size())/2) 
                       ++numCore.second;
                     if(isVanished) break;
		   } // End I
                   if(numCore.first != numCore.second) isVanished = true;
                 } // End if
                 if(!isVanished)
                   { iter_terms.push_back(SQterm((-1.0) * t->get_numConst(), term1.get_Consts(), tensorList)); }
	       } // End if
             } // End if
             else
               iter_terms.push_back(SQterm((-1.0) * t->get_numConst(), term1.get_Consts(), tensorList));
           } // End tup2
         } // End tup1
       } // End perm          
     } // End nc
     } // End second

   } // End t
   
   // 
   inTerms->insert(inTerms->end(), iter_terms.begin()+Count, iter_terms.end());

   return;
  }

}} // Femto::core
