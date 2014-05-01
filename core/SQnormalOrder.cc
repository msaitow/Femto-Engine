//
//  SQnormalOrder.cc
//  
//
//  Created by Masaaki Saitow on 12/07/08.
//  Copyright (c) 2012 Masaaki Saitow. All rights reserved.
//

#include <vector>
#include <algorithm>
#include <cstdlib>
#include <Femto.hpp>
#include <SQterm.hpp>
#include <SQtensor.hpp>
#include <SQindex.hpp>

// * 2012/08/21 : If OF_flag is turned on, the normal-ordered sfGen that has diffent numbers of core indices on the upper and lower
//                side is not pushed back.

#define _ZERO_CONT // Flag to avoid the operator contraction in terms of spin-free generators to give vanishing result 

using namespace std;

namespace Femto { namespace Core { 

  // ****************************************************************
  // Bunch of code to achieve normal ordering. This is an inetrface
  // ****************************************************************
  void normalOrder(vector<SQterm> *inTerms, const bool OF_flag)
  {
    pair<bool, bool> so_sf; // <so-base, sf-base>
    so_sf.first  = false;
    so_sf.second = false;
    for(size_t num_t = 0;num_t < inTerms->at(0).get_tensors().size();++num_t)
      if(inTerms->at(0).get_tensors()[num_t].get_name() == aCre_name() || inTerms->at(0).get_tensors()[num_t].get_name() == aDes_name())
	so_sf.first = true;
      else if(is_sfGen(inTerms->at(0).get_tensors()[num_t].get_name()))
	so_sf.second = true;

    // In case only spin-free unitary group generator exists
    if     (!so_sf.first &&  so_sf.second) normalOrder_sf(inTerms, OF_flag);
//*NOT YET*     else if( so_sf.first && !so_sf.second) normalOrder_so(inTerms, OF_flag);
//*NOT YET*     else if( so_sf.first &&  so_sf.second) normalOrder_so_sf(inTerms, OF_flag);
    // In case only second quatization operator exists
    else if( so_sf.first && !so_sf.second) {cout << "Not implemented yet" << endl; abort();}
    // In case both second quatization operator and the generator exist
    else if( so_sf.first &&  so_sf.second) {cout << "Not implemented yet" << endl; abort();}

  }


  // ****************************************************************
  // Bunch of code to achieve normal ordering
  // ****************************************************************
  // ----------------------------------------------------------------
  // *NOTE* :: Normal ordering of the spin-free generator is more 
  //           straightforward procedure than that of the usual
  //           second quantization operators because the freedom of 
  //           permutation is already confined in some way. This code
  //           normal-orders a stream of the generator successively from
  //           right most pairs. In the contraction procedure, 
  //           if Kronecker's delta with indices of different orbital
  //           groups, or un-contractable core/virtual indices in the 
  //           generator are found, such the term is abondoned 
  //           at that moment to reduce manipulation and time-saving. 
  // ----------------------------------------------------------------
  void normalOrder_sf(vector<SQterm> *inTerms, const bool OF_flag)
  {
    // Now, only one term is allowed as an input
    if(inTerms->size()!=1){
      cout << "Now, only one term is allowed as an input";
      abort();
    }

    // If there's other indices than core, active and virtual, abort
    for(size_t num = 0;num < inTerms->at(0).get_summedBody().size();++num)
      if(inTerms->at(0).get_summedBody()[num]->get_char() != core && 
         inTerms->at(0).get_summedBody()[num]->get_char() != act  &&
         inTerms->at(0).get_summedBody()[num]->get_char() != virt){
        cout << "normalOrder: There's at least one index that doesn't belong to either of core, active or virtual orbital" << endl;
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
      //if(typeid(inTensors[i]) == typeid(sfGen)) // <---- Should be like this !!!!
      string t_name = inTensors[i].get_name(); //cout << inTensors[i] << endl;
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
        for(int nc = 0;nc < min(*o1,*o2)+1;++nc){
          int newOrder = *o1 + *o2 - nc; 
          int nterm = 1;
          if(nc!=0) {
            for(int i = 0;i < nc;++i) nterm = nterm * (*o1-i) * (*o2-i) / fact(i+1);
          } // End if
          vector<int> o_list(except_last_two);
          o_list.push_back(newOrder);
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
     SQtensor e1 = *it_m2;
     SQtensor e2 = *it_m1;
     int o1 = (int)(e1.get_indices().size()/2); // cout << o1 << endl; //*TEST*
     int o2 = (int)(e2.get_indices().size()/2); // cout << o2 << endl; //*TEST*

     // Loop over the number of contraction
     for(int nc = 0;nc < min(o1,o2)+1;++nc){

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

             bool isDead(false); // Whether kDelta vanises or not
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
	       Ivector::iterator i_index = find(conPairs[0].begin(), conPairs[0].end(), i);

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

	       Ivector::iterator i_index = find(conPairs[1].begin(), conPairs[1].end(), i);

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
                   { iter_terms.push_back(SQterm(term1.get_numConst(), term1.get_Consts(), tensorList)); }
	       } // End if
             } // End if
             else
               iter_terms.push_back(SQterm(term1.get_numConst(), term1.get_Consts(), tensorList));
           } // End tup2
         } // End tup1
       } // End perm          
     } // End nc
   } // End t
   inTerms->insert(inTerms->end(), iter_terms.begin()+Count, iter_terms.end());

   return;
  }


//*OLD*   // *********************************************************
//*OLD*   // Normal order the inTerm (mostly taken from sqaNormalOrder.py)
//*OLD*   // *********************************************************
//*OLD*   void normalOrder(vector<SQterm> *inTerms, const bool OF_flag)
//*OLD*   {
//*OLD*     // Now, only one term is allowed as an input
//*OLD*     if(inTerms->size()!=1){
//*OLD*       cout << "Now, only one term is allowed as an input";
//*OLD*       abort();
//*OLD*     }
//*OLD*     // Initialize number n, which is a number of the ramining sfGen
//*OLD*     if((*inTerms)[0].get_isInCanonical()) return;
//*OLD* 
//*OLD*     // Make separate list of the sfGen and the others
//*OLD*     vector<SQtensor> sfGen_list;
//*OLD*     vector<SQtensor> other_list;
//*OLD* 
//*OLD*     vector<SQtensor> inTensors((*inTerms)[0].get_tensors());
//*OLD*     for(size_t i = 0;i < inTensors.size();++i){
//*OLD*       //if(typeid(inTensors[i]) == typeid(sfGen)) // <---- Should be like this !!!!
//*OLD*       string t_name = inTensors[i].get_name(); //cout << inTensors[i] << endl;
//*OLD*       if(is_sfGen(t_name)) {
//*OLD*         sfGen_list.push_back(inTensors[i]);
//*OLD*       }
//*OLD*       else
//*OLD*         other_list.push_back(inTensors[i]);
//*OLD*     } // End i
//*OLD*     int n = (int)sfGen_list.size();
//*OLD* 
//*OLD*     // Estimate the number of contractions arise from each step in the normal ordering
//*OLD*     int m = n;
//*OLD*     int totnum = 0;
//*OLD*     vector<int> orders;
//*OLD*     for(vector<SQtensor>::iterator gen = sfGen_list.begin();gen != sfGen_list.end();++gen)
//*OLD*       orders.push_back((gen->get_indices().size())/2);
//*OLD*     IIvector init_orders; init_orders.push_back(orders);
//*OLD*     IIvector out_orders;
//*OLD*     while(m > 1){
//*OLD*       out_orders.clear();
//*OLD*       for(vector<vector<int> >::iterator o = init_orders.begin();o != init_orders.end();++o){
//*OLD*         vector<int>::iterator o0 = o->begin();
//*OLD*         vector<int>::iterator o2 = o->end(); --o2;
//*OLD*         vector<int>::iterator o1 = o2;       --o1;
//*OLD*         vector<int> except_last_two;
//*OLD*         except_last_two.insert(except_last_two.end(), o0, o1);
//*OLD*         // Count numbers of terms of order of newOrder
//*OLD*         for(int nc = 0;nc < min(*o1,*o2)+1;++nc){
//*OLD*           int newOrder = *o1 + *o2 - nc; 
//*OLD*           int nterm = 1;
//*OLD*           if(nc!=0) {
//*OLD*             for(int i = 0;i < nc;++i) nterm = nterm * (*o1-i) * (*o2-i) / fact(i+1);
//*OLD*           } // End if
//*OLD*           vector<int> o_list(except_last_two);
//*OLD*           o_list.push_back(newOrder);
//*OLD*           out_orders.insert(out_orders.end(), nterm, o_list);
//*OLD*         }        
//*OLD*       } // End o
//*OLD*       init_orders = out_orders;
//*OLD*       totnum += (int)init_orders.size();
//*OLD*       --m; 
//*OLD*     } // End while
//*OLD* 
//*OLD*     // Now we've already got number of contractions .... 
//*OLD*     vector<SQtensor> tensors(other_list);
//*OLD*     tensors.insert(tensors.end(), sfGen_list.begin(), sfGen_list.end());
//*OLD*     SQterm term1((*inTerms)[0].get_numConst(), (*inTerms)[0].get_Consts(), tensors);
//*OLD*     inTerms->reserve(init_orders.size());
//*OLD*     vector<SQterm> iter_terms; iter_terms.reserve(totnum);
//*OLD*     // Push back the conditioned term
//*OLD*     iter_terms.push_back(term1);
//*OLD*     // And erase the old one
//*OLD*     inTerms->erase(inTerms->begin());
//*OLD* 
//*OLD*    // Successively normal order the last two excitation operators until each term
//*OLD*    // has only one exitation operator left (at which point the term is normal ordered)
//*OLD* 
//*OLD*    // For each of this iteration's input terms, produce all terms resulting from normal ordering
//*OLD*    // the last two excitation operators
//*OLD*    size_t Count = 0;
//*OLD*    for(vector<SQterm>::iterator t = iter_terms.begin();t != iter_terms.end();++t){
//*OLD*      if(t->get_isInCanonical()) break;
//*OLD*      else ++Count;    
//*OLD*      
//*OLD*      // Make a list of the term's tensors that excludes the last two excitation operators
//*OLD*      vector<SQtensor> tensors_except_last_two;
//*OLD*      vector<SQtensor> t_temp = t->get_tensors();
//*OLD*      vector<SQtensor>::iterator it    = t_temp.begin();
//*OLD*      vector<SQtensor>::iterator it_m1 = t_temp.end(); 
//*OLD*      vector<SQtensor>::iterator it_m2 = t_temp.end(); 
//*OLD*      --it_m1;
//*OLD*      --it_m2; --it_m2;
//*OLD*      tensors_except_last_two.insert(tensors_except_last_two.end(), it, it_m2);
//*OLD*      
//*OLD*      //  Give short names for the last two excitation operators and their orders
//*OLD*      SQtensor e1 = *it_m2;
//*OLD*      SQtensor e2 = *it_m1;
//*OLD*      int o1 = (int)(e1.get_indices().size()/2); // cout << o1 << endl; //*TEST*
//*OLD*      int o2 = (int)(e2.get_indices().size()/2); // cout << o2 << endl; //*TEST*
//*OLD* 
//*OLD*      // Loop over the number of contraction
//*OLD*      for(int nc = 0;nc < min(o1,o2)+1;++nc){
//*OLD* 
//*OLD*        int newOrder = o1 + o2 - nc;
//*OLD*        
//*OLD*        // Compute all nc-tuples of index numbers from e1 and e2, as well as all permutations of 
//*OLD*        // the order in which the tuples may be combined to form a contraction
//*OLD*        Ivector o1_list; o1_list.reserve(o1);
//*OLD*           Ivector o2_list; o2_list.reserve(o2);
//*OLD*        for(int i = 0;i < o1;++i) o1_list.push_back(i);
//*OLD*        for(int i = 0;i < o2;++i) o2_list.push_back(i);
//*OLD* 
//*OLD*        IIvector perms; 
//*OLD*        Ivector temp;
//*OLD*        temp.push_back(0);
//*OLD*        perms.push_back(temp);
//*OLD*        if(nc > 0)
//*OLD*          perms = makePermutations(nc);
//*OLD*        IIvector tups1 = makeTuples2(nc, o1_list);
//*OLD*        IIvector tups2 = makeTuples2(nc, o2_list);
//*OLD* 
//*OLD*        // For each contractions, compute the resulting term
//*OLD*        for(IIvector::iterator perm = perms.begin();perm != perms.end();++perm){
//*OLD*          for(IIvector::iterator tup1 = tups1.begin();tup1 != tups1.end();++tup1){
//*OLD*            for(IIvector::iterator tup2 = tups2.begin();tup2 != tups2.end();++tup2){
//*OLD*              // Initialize the term's tensor list
//*OLD*              vector<SQtensor> tensorList; 
//*OLD*              tensorList.insert(tensorList.end(), tensors_except_last_two.begin(),
//*OLD*                                                  tensors_except_last_two.end());
//*OLD* 
//*OLD*              bool isDead = false; // Whether kDelta vanises or not
//*OLD*              // Compute the pairs of indices to be contracted
//*OLD*              // Example: (conPairs[0][p], conPairs[1][p]) is the p-th contraction pair.
//*OLD*              IIvector conPairs;
//*OLD*              Ivector p; 
//*OLD*              conPairs.push_back(p);
//*OLD*              conPairs.push_back(p);
//*OLD*              for(int i = 0;i < nc;++i){
//*OLD*                conPairs[0].push_back((*tup1)[(*perm)[i]]);
//*OLD*                conPairs[1].push_back((*tup2)[i]);
//*OLD*              } // End i
//*OLD* 
//*OLD*              // Initialize the index list for the new excitation operator
//*OLD*              vector<SQindex*> indexList; indexList.reserve(2*newOrder);
//*OLD*              SQindex *q = NULL;
//*OLD*              for(int i = 0;i < 2*newOrder;++i) indexList.push_back(q);
//*OLD* 
//*OLD*              // Populate the index list for the new excitation operator
//*OLD*              // Also, create a kronecker delta function for each contraction pair
//*OLD*              vector<SQindex*> e1_indices(e1.get_indices());
//*OLD*              vector<SQindex*> e2_indices(e2.get_indices());
//*OLD* 
//*OLD*              for(int i = 0;i < o1;++i){
//*OLD*                indexList[i] = e1_indices[i];
//*OLD* 
//*OLD*                bool i_flag = false;
//*OLD*                size_t i_index;
//*OLD*                for(size_t j = 0;j < conPairs[0].size();++j){
//*OLD*                  if(conPairs[0][j] == i) {
//*OLD*                    i_flag = true;
//*OLD*                    i_index = j;
//*OLD*                  }
//*OLD*                } // End j
//*OLD* 
//*OLD*                if(i_flag){
//*OLD*                  int i1 = i + o1;
//*OLD*                  int i2 = conPairs[1][i_index];
//*OLD*                  SQindex* ind1 = e1_indices[i1]; 
//*OLD*                  SQindex* ind2 = e2_indices[i2]; 
//*OLD* 
//*OLD*                  vector<SQindex*> ind_delta;
//*OLD*                  ind_delta.push_back(ind1);
//*OLD*                  ind_delta.push_back(ind2);
//*OLD*                  if(ind_delta[0]->get_char()!=ind_delta[1]->get_char()) isDead = true; 
//*OLD*                  tensorList.push_back(kDelta(ind_delta)); 
//*OLD*                  indexList[i+newOrder]  = e2_indices[o2+i2]; 
//*OLD* 
//*OLD*                } // End if
//*OLD*                else
//*OLD*                  indexList[i+newOrder]  = e1_indices[i+o1];
//*OLD* 
//*OLD*              } // End i
//*OLD*              int count = 0;
//*OLD*              for(int i = 0;i < o2;++i){
//*OLD*                bool i_flag = false;
//*OLD*                size_t i_index;
//*OLD*                for(size_t j = 0;j < conPairs[1].size();++j){
//*OLD*                  if(conPairs[1][j] == i) {
//*OLD*                    i_flag = true;
//*OLD*                    i_index = j;
//*OLD*                  }
//*OLD*                } // End j
//*OLD*                if(!i_flag){
//*OLD*                  indexList[o1+count]          = e2_indices[i];
//*OLD*                  indexList[o1+count+newOrder] = e2_indices[i+o2];
//*OLD*                     ++count;
//*OLD*                } // End if
//*OLD*              } // End i
//*OLD* 
//*OLD*              // Ensure that all slots in the index list have been filled
//*OLD*              // And check whether virtual exists or not
//*OLD*              for(size_t numInd = 0;numInd < indexList.size();++numInd){
//*OLD*                if(indexList[numInd] == NULL){
//*OLD*                  cout << "There is at least one unassigned index in the new spin-free unitary group generators" << endl;
//*OLD*                  abort();
//*OLD*                }
//*OLD*              } // End numInd
//*OLD*              // Add the new unitary generator into the tensorList
//*OLD*              tensorList.push_back(sfGen(indexList));
//*OLD* 
//*OLD*              // Add the resulting term to this iteration's list of output terms
//*OLD*              if(OF_flag){
//*OLD*                if(!isDead){
//*OLD*                  int Ecount=0;
//*OLD*                  vector<SQtensor>::iterator ten_E;
//*OLD*                  for(vector<SQtensor>::iterator ten = tensorList.begin();ten != tensorList.end();++ten){
//*OLD*                    if(is_sfGen(ten->get_name())) { ++Ecount; ten_E = ten; }
//*OLD*                  } // End ten
//*OLD*                  bool isVanished = false;
//*OLD*                  if(Ecount==0 || Ecount==1) {
//*OLD*                    vector<SQindex*> ind(ten_E->get_indices());
//*OLD*                    for(vector<SQindex*>::iterator I = ind.begin();I != ind.end();++I) 
//*OLD*                      if((*I)->get_char()==(char_state)2) isVanished = true;
//*OLD*                  } // End if
//*OLD*                  if(!isVanished)
//*OLD*                    { iter_terms.push_back(SQterm(term1.get_numConst(), term1.get_Consts(), tensorList)); }
//*OLD* 	       } // End if
//*OLD*              } // End if
//*OLD*              else
//*OLD*                iter_terms.push_back(SQterm(term1.get_numConst(), term1.get_Consts(), tensorList));
//*OLD*            } // End tup2
//*OLD*          } // End tup1
//*OLD*        } // End perm          
//*OLD*      } // End nc
//*OLD*    } // End t
//*OLD*    inTerms->insert(inTerms->end(), iter_terms.begin()+Count, iter_terms.end());
//*OLD* 
//*OLD*    return;
//*OLD*   }


//*OLD*   // *********************************************************
//*OLD*   // Normal order the inTerm (mostly taken from sqaNormalOrder.py)
//*OLD*   // *********************************************************
//*OLD*   vector<SQterm> normalOrder(const SQterm &inTerm, const bool OF_flag)
//*OLD*   {
//*OLD*     vector<SQterm> outTerms;
//*OLD* 
//*OLD*     // Make separate list of the sfGen and the others
//*OLD*     vector<SQtensor> sfGen_list;
//*OLD*     vector<SQtensor> other_list;
//*OLD* 
//*OLD*     vector<SQtensor> inTensors(inTerm.get_tensors());
//*OLD*     for(size_t i = 0;i < inTensors.size();++i){
//*OLD*       //if(typeid(inTensors[i]) == typeid(sfGen)) // <---- Should be like this !!!!
//*OLD*       string t_name = inTensors[i].get_name();
//*OLD*       if(is_sfGen(t_name)) {
//*OLD*         sfGen_list.push_back(inTensors[i]);
//*OLD*       }
//*OLD*       else
//*OLD*         other_list.push_back(inTensors[i]);
//*OLD*     } // End i
//*OLD* 
//*OLD*     // Initialize number n, which is a number of the ramining sfGen
//*OLD*     int n = (int)sfGen_list.size(); 
//*OLD*     if(inTerm.get_isInCanonical()){
//*OLD*       outTerms.push_back(inTerm);
//*OLD*       return outTerms;
//*OLD*     }
//*OLD* 
//*OLD*     // Estimate the number of contractions arise from each step in the normal ordering
//*OLD*     int m = n;
//*OLD*     vector<int> orders;
//*OLD*     for(vector<SQtensor>::iterator gen = sfGen_list.begin();gen != sfGen_list.end();++gen)
//*OLD*       orders.push_back((gen->get_indices().size())/2);
//*OLD*     IIvector init_orders; init_orders.push_back(orders);
//*OLD*     IIvector out_orders;
//*OLD*     while(m > 1){
//*OLD*       out_orders.clear();
//*OLD*       for(vector<vector<int> >::iterator o = init_orders.begin();o != init_orders.end();++o){
//*OLD*         vector<int>::iterator o0 = o->begin();
//*OLD*         vector<int>::iterator o2 = o->end(); --o2;
//*OLD*         vector<int>::iterator o1 = o2;       --o1;
//*OLD*         vector<int> except_last_two;
//*OLD*         except_last_two.insert(except_last_two.end(), o0, o1);
//*OLD*         // Count numbers of terms of order of newOrder
//*OLD*         for(int nc = 0;nc < min(*o1,*o2)+1;++nc){
//*OLD*           int newOrder = *o1 + *o2 - nc; 
//*OLD*           int nterm = 1;
//*OLD*           if(nc!=0) {
//*OLD*             for(int i = 0;i < nc;++i) nterm = nterm * (*o1-i) * (*o2-i) / fact(i+1);
//*OLD* 	  } // End if
//*OLD*           vector<int> o_list(except_last_two);
//*OLD*           o_list.push_back(newOrder);
//*OLD*           out_orders.insert(out_orders.end(), nterm, o_list);
//*OLD* 	}        
//*OLD*       } // End o
//*OLD*       init_orders = out_orders;
//*OLD*       --m; 
//*OLD*     } // End while
//*OLD* 
//*OLD*     // Now we've already got number of contractions .... 
//*OLD*     vector<SQterm> iter_input_terms;  iter_input_terms.reserve(init_orders.size()); // Replace with Nterms()
//*OLD*     vector<SQterm> iter_output_terms; iter_output_terms.reserve(init_orders.size()); // Replace with Nterms()
//*OLD* 
//*OLD*     vector<SQtensor> tensors(other_list);
//*OLD*     tensors.insert(tensors.end(), sfGen_list.begin(), sfGen_list.end());
//*OLD*     SQterm term1(inTerm.get_numConst(), inTerm.get_Consts(), tensors);
//*OLD*     iter_input_terms.push_back(term1);
//*OLD* 
//*OLD*     // Successively normal order the last two excitation operators until each term
//*OLD*     // has only one exitation operator left (at which point the term is normal ordered)
//*OLD*     while(n > 1){
//*OLD*  
//*OLD*       // Initialize the list to hold this iteration's output terms
//*OLD*       iter_output_terms.clear(); 
//*OLD* 
//*OLD*       // For each of this iteration's input terms, produce all terms resulting from normal ordering
//*OLD*       // the last two excitation operators
//*OLD*       vector<SQterm>::iterator t = iter_input_terms.begin();
//*OLD*       for(;t != iter_input_terms.end();++t){
//*OLD*         
//*OLD*         // Make a list of the term's tensors that excludes the last two excitation operators
//*OLD*         vector<SQtensor> tensors_except_last_two;
//*OLD*         vector<SQtensor> t_temp = t->get_tensors();
//*OLD*         vector<SQtensor>::iterator it    = t_temp.begin();
//*OLD*         vector<SQtensor>::iterator it_m1 = t_temp.end(); 
//*OLD*         vector<SQtensor>::iterator it_m2 = t_temp.end(); 
//*OLD*         --it_m1;
//*OLD* 	--it_m2; --it_m2;
//*OLD* 	tensors_except_last_two.insert(tensors_except_last_two.end(), it, it_m2);
//*OLD* 
//*OLD*         //  Give short names for the last two excitation operators and their orders
//*OLD*         SQtensor e1 = *it_m2;
//*OLD*         SQtensor e2 = *it_m1;
//*OLD*         int o1 = (int)(e1.get_indices().size()/2); // cout << o1 << endl; //*TEST*
//*OLD*         int o2 = (int)(e2.get_indices().size()/2); // cout << o2 << endl; //*TEST*
//*OLD* 
//*OLD*         // Loop over the number of contraction
//*OLD*         for(int nc = 0;nc < min(o1,o2)+1;++nc){
//*OLD* 
//*OLD*           int newOrder = o1 + o2 - nc;
//*OLD*           
//*OLD*           // Compute all nc-tuples of index numbers from e1 and e2, as well as all permutations of 
//*OLD*           // the order in which the tuples may be combined to form a contraction
//*OLD*           Ivector o1_list; o1_list.reserve(o1);
//*OLD* 	  Ivector o2_list; o2_list.reserve(o2);
//*OLD*           for(int i = 0;i < o1;++i) o1_list.push_back(i);
//*OLD*           for(int i = 0;i < o2;++i) o2_list.push_back(i);
//*OLD* 
//*OLD*           IIvector perms; 
//*OLD*           Ivector temp;
//*OLD*           temp.push_back(0);
//*OLD*           perms.push_back(temp);
//*OLD*           if(nc > 0)
//*OLD*             perms = makePermutations(nc);
//*OLD*           IIvector tups1 = makeTuples2(nc, o1_list);
//*OLD*           IIvector tups2 = makeTuples2(nc, o2_list);
//*OLD* 
//*OLD*           // For each contractions, compute the resulting term
//*OLD* 	  IIvector::iterator perm = perms.begin(); 
//*OLD* 	  IIvector::iterator tup1 = tups1.begin(); 
//*OLD* 	  IIvector::iterator tup2 = tups2.begin(); 
//*OLD*           for(IIvector::iterator perm = perms.begin();perm != perms.end();++perm){
//*OLD*             for(IIvector::iterator tup1 = tups1.begin();tup1 != tups1.end();++tup1){
//*OLD*               for(IIvector::iterator tup2 = tups2.begin();tup2 != tups2.end();++tup2){
//*OLD*                 // Initialize the term's tensor list
//*OLD*                 vector<SQtensor> tensorList; 
//*OLD*                 tensorList.insert(tensorList.end(), tensors_except_last_two.begin(),
//*OLD* 				                    tensors_except_last_two.end());
//*OLD* 
//*OLD*                 bool isDead = false; // Whether kDelta vanises or not
//*OLD*                 // Compute the pairs of indices to be contracted
//*OLD*                 // Example: (conPairs[0][p], conPairs[1][p]) is the p-th contraction pair.
//*OLD*                 IIvector conPairs;
//*OLD*                 Ivector p; 
//*OLD*                 conPairs.push_back(p);
//*OLD*                 conPairs.push_back(p);
//*OLD*                 for(int i = 0;i < nc;++i){
//*OLD*                   conPairs[0].push_back((*tup1)[(*perm)[i]]);
//*OLD*                   conPairs[1].push_back((*tup2)[i]);
//*OLD* 		} // End i
//*OLD* 
//*OLD*                 // Initialize the index list for the new excitation operator
//*OLD*                 vector<SQindex*> indexList; indexList.reserve(2*newOrder);
//*OLD*                 SQindex *q = NULL;
//*OLD*                 for(int i = 0;i < 2*newOrder;++i) indexList.push_back(q);
//*OLD* 
//*OLD*                 // Populate the index list for the new excitation operator
//*OLD*                 // Also, create a kronecker delta function for each contraction pair
//*OLD*                 vector<SQindex*> e1_indices(e1.get_indices());
//*OLD*                 vector<SQindex*> e2_indices(e2.get_indices());
//*OLD* 
//*OLD*                 for(int i = 0;i < o1;++i){
//*OLD*                   indexList[i] = e1_indices[i];
//*OLD* 
//*OLD*                   bool i_flag = false;
//*OLD*                   size_t i_index;
//*OLD*                   for(size_t j = 0;j < conPairs[0].size();++j){
//*OLD*                     if(conPairs[0][j] == i) {
//*OLD*                       i_flag = true;
//*OLD*                       i_index = j;
//*OLD* 		    }
//*OLD* 		  } // End j
//*OLD* 
//*OLD*                   if(i_flag){
//*OLD*                     int i1 = i + o1;
//*OLD*                     int i2 = conPairs[1][i_index];
//*OLD*                     SQindex* ind1 = e1_indices[i1]; 
//*OLD*                     SQindex* ind2 = e2_indices[i2]; 
//*OLD* 
//*OLD*                     vector<SQindex*> ind_delta;
//*OLD*                     ind_delta.push_back(ind1);
//*OLD*                     ind_delta.push_back(ind2);
//*OLD*                     if(ind_delta[0]->get_char()!=ind_delta[1]->get_char()) isDead = true;
//*OLD*                     tensorList.push_back(kDelta(ind_delta));
//*OLD* 
//*OLD*                     vector<SQtensor>::iterator it = tensorList.end();--it;
//*OLD* 
//*OLD*                     indexList[i+newOrder]  = e2_indices[o2+i2]; 
//*OLD* 
//*OLD* 		  } // End if
//*OLD*                   else
//*OLD*                     indexList[i+newOrder]  = e1_indices[i+o1];
//*OLD* 
//*OLD* 		} // End i
//*OLD*                 int count = 0;
//*OLD*                 for(int i = 0;i < o2;++i){
//*OLD*                   bool i_flag = false;
//*OLD*                   size_t i_index;
//*OLD*                   for(size_t j = 0;j < conPairs[1].size();++j){
//*OLD*                     if(conPairs[1][j] == i) {
//*OLD*                       i_flag = true;
//*OLD*                       i_index = j;
//*OLD* 		    }
//*OLD* 		  } // End j
//*OLD*                   if(!i_flag){
//*OLD*                     indexList[o1+count]          = e2_indices[i];
//*OLD*                     indexList[o1+count+newOrder] = e2_indices[i+o2];
//*OLD* 		    ++count;
//*OLD* 		  } // End if
//*OLD* 		} // End i
//*OLD* 
//*OLD*                 // Ensure that all slots in the index list have been filled
//*OLD*                 // And check whether virtual exists or not
//*OLD*                 for(size_t numInd = 0;numInd < indexList.size();++numInd){
//*OLD*                   if(indexList[numInd] == NULL){
//*OLD*                     cout << "There is at least one unassigned index in the new spin-free unitary generators" << endl;
//*OLD*                     abort();
//*OLD* 		  }
//*OLD* 		} // End numInd
//*OLD*                 // Add the new unitary generator into the tensorList
//*OLD*                 tensorList.push_back(sfGen(indexList));
//*OLD* 
//*OLD*                 // Add the resulting term to this iteration's list of output terms
//*OLD*                 if(OF_flag && !isDead)
//*OLD*                 iter_output_terms.push_back(SQterm(inTerm.get_numConst(), inTerm.get_Consts(), tensorList));
//*OLD* 
//*OLD* 	      } // End tup2
//*OLD* 	    } // End tup1
//*OLD* 	  } // End perm          
//*OLD* 
//*OLD* 	} // End nc
//*OLD* 
//*OLD*       } // End t
//*OLD* 
//*OLD*       // Set this iteration's list of output terms as the next iteration's input terms
//*OLD*       iter_input_terms = iter_output_terms;
//*OLD*       --n;
//*OLD* 
//*OLD*     } // End while    
//*OLD* 
//*OLD*     return iter_output_terms;
//*OLD*   }


//*OLD* //*OLD*   // *********************************************************
//*OLD* //*OLD*   // Normal order the inTerm (mostly taken from sqaNormalOrder.py)
//*OLD* //*OLD*   // *********************************************************
//*OLD* //*OLD*   vector<SQterm> normalOrder(const SQterm &inTerm, const bool OF_flag)
//*OLD* //*OLD*   {
//*OLD* //*OLD*     vector<SQterm> outTerms;
//*OLD* //*OLD* 
//*OLD* //*OLD*     // Make separate list of the sfGen and the others
//*OLD* //*OLD*     vector<SQtensor> sfGen_list;
//*OLD* //*OLD*     vector<SQtensor> other_list;
//*OLD* //*OLD* 
//*OLD* //*OLD* //     vector<SQindex*> temp_i;
//*OLD* //*OLD* //     SQindex tempInd1;
//*OLD* //*OLD* //     SQindex tempInd2;
//*OLD* //*OLD* //     temp_i.push_back(&tempInd1);
//*OLD* //*OLD* //     temp_i.push_back(&tempInd2);
//*OLD* //*OLD* //     sfGen Gen_ref(temp_i);
//*OLD* //*OLD* //     cout << typeid(Gen_ref).name() << endl << endl; //*TEST*
//*OLD* //*OLD* //     vector<SQtensor> test_v;
//*OLD* //*OLD* //     test_v.push_back(Gen_ref);
//*OLD* //*OLD* //     cout << typeid(test_v[0]).name() << endl << endl;    
//*OLD* //*OLD* 
//*OLD* //*OLD*     vector<SQtensor> inTensors(inTerm.get_tensors());
//*OLD* //*OLD*     for(size_t i = 0;i < inTensors.size();++i){
//*OLD* //*OLD*       //if(typeid(inTensors[i]) == typeid(sfGen)) // <---- Should be like this !!!!
//*OLD* //*OLD*       string t_name = inTensors[i].get_name();
//*OLD* //*OLD*       if(is_sfGen(t_name)) {
//*OLD* //*OLD*         sfGen_list.push_back(inTensors[i]);
//*OLD* //*OLD*       }
//*OLD* //*OLD*       else
//*OLD* //*OLD*         other_list.push_back(inTensors[i]);
//*OLD* //*OLD*     } // End i
//*OLD* //*OLD* 
//*OLD* //*OLD*     // Initialize number n, which is a number of the ramining sfGen
//*OLD* //*OLD*     int n = (int)sfGen_list.size(); 
//*OLD* //*OLD*     if(inTerm.get_isInCanonical()){
//*OLD* //*OLD*       outTerms.push_back(inTerm);
//*OLD* //*OLD*       return outTerms;
//*OLD* //*OLD*     }
//*OLD* //*OLD* 
//*OLD* //*OLD* //*TEST*    cout << n << endl; //*TEST*
//*OLD* //*OLD* //*TEST*    cout << other_list.size() << endl; //*TEST*
//*OLD* //*OLD* 
//*OLD* //*OLD*     vector<SQterm> iter_input_terms; iter_input_terms.reserve(Nterms());
//*OLD* //*OLD*     vector<SQtensor> tensors(other_list);
//*OLD* //*OLD*     tensors.insert(tensors.end(), sfGen_list.begin(), sfGen_list.end());
//*OLD* //*OLD*     SQterm term1(inTerm.get_numConst(), inTerm.get_Consts(), tensors);
//*OLD* //*OLD*     iter_input_terms.push_back(term1);
//*OLD* //*OLD* 
//*OLD* //*OLD*     vector<SQterm> iter_output_terms; iter_output_terms.reserve(Nterms());
//*OLD* //*OLD* 
//*OLD* //*OLD*     // Successively normal order the last two excitation operators until each term
//*OLD* //*OLD*     // has only one exitation operator left (at which point the term is normal ordered)
//*OLD* //*OLD*     while(n > 1){
//*OLD* //*OLD*  
//*OLD* //*OLD*       // Initialize the list to hold this iteration's output terms
//*OLD* //*OLD*       iter_output_terms.clear();
//*OLD* //*OLD* 
//*OLD* //*OLD*       // For each of this iteration's input terms, produce all terms resulting from normal ordering
//*OLD* //*OLD*       // the last two excitation operators
//*OLD* //*OLD*       vector<SQterm>::iterator t = iter_input_terms.begin();
//*OLD* //*OLD*       for(;t != iter_input_terms.end();++t){
//*OLD* //*OLD*         
//*OLD* //*OLD*         // Make a list of the term's tensors that excludes the last two excitation operators
//*OLD* //*OLD*         vector<SQtensor> tensors_except_last_two;
//*OLD* //*OLD*         vector<SQtensor> t_temp = t->get_tensors();
//*OLD* //*OLD*         vector<SQtensor>::iterator it    = t_temp.begin();
//*OLD* //*OLD*         vector<SQtensor>::iterator it_m1 = t_temp.end(); 
//*OLD* //*OLD*         vector<SQtensor>::iterator it_m2 = t_temp.end(); 
//*OLD* //*OLD*         --it_m1;
//*OLD* //*OLD* 	--it_m2;
//*OLD* //*OLD* 	--it_m2;
//*OLD* //*OLD* 	tensors_except_last_two.insert(tensors_except_last_two.end(), it, it_m2);
//*OLD* //*OLD* 
//*OLD* //*OLD* //*TEST*        cout << *it_m2 << ", " << *it_m1 << endl; //*TEST*
//*OLD* //*OLD* //*TEST*        cout << *t << endl;                       //*TEST*
//*OLD* //*OLD* 
//*OLD* //*OLD*         //  Give short names for the last two excitation operators and their orders
//*OLD* //*OLD*         SQtensor e1 = *it_m2;
//*OLD* //*OLD*         SQtensor e2 = *it_m1;
//*OLD* //*OLD*         int o1 = (int)(e1.get_indices().size()/2); // cout << o1 << endl; //*TEST*
//*OLD* //*OLD*         int o2 = (int)(e2.get_indices().size()/2); // cout << o2 << endl; //*TEST*
//*OLD* //*OLD* 
//*OLD* //*OLD*         // Loop over the number of contraction
//*OLD* //*OLD*         for(int nc = 0;nc < min(o1,o2)+1;++nc){
//*OLD* //*OLD* 
//*OLD* //*OLD*           int newOrder = o1 + o2 - nc;// cout << "made " <<  newOrder << endl;
//*OLD* //*OLD*           
//*OLD* //*OLD*           // Compute all nc-tuples of index numbers from e1 and e2, as well as all permutations of 
//*OLD* //*OLD*           // the order in which the tuples may be combined to form a contraction
//*OLD* //*OLD*           Ivector o1_list;
//*OLD* //*OLD* 	  Ivector o2_list;
//*OLD* //*OLD*           for(int i = 0;i < o1;++i) o1_list.push_back(i);
//*OLD* //*OLD*           for(int i = 0;i < o2;++i) o2_list.push_back(i);
//*OLD* //*OLD* 
//*OLD* //*OLD*           IIvector perms; //perms.push_back(0);
//*OLD* //*OLD*           Ivector temp;
//*OLD* //*OLD*           temp.push_back(0);
//*OLD* //*OLD*           perms.push_back(temp);
//*OLD* //*OLD*           if(nc > 0)
//*OLD* //*OLD*             perms = makePermutations(nc);
//*OLD* //*OLD*           IIvector tups1 = makeTuples2(nc, o1_list);
//*OLD* //*OLD*           IIvector tups2 = makeTuples2(nc, o2_list);
//*OLD* //*OLD* 
//*OLD* //*OLD*           // For each contractions, compute the resulting term
//*OLD* //*OLD* 	  IIvector::iterator perm = perms.begin(); //cout << newOrder << " " << perms.size() << endl;
//*OLD* //*OLD* 	  IIvector::iterator tup1 = tups1.begin(); //cout << newOrder << " " << tups1.size() << endl;
//*OLD* //*OLD* 	  IIvector::iterator tup2 = tups2.begin(); //cout << newOrder << " " << tups2.size() << endl;
//*OLD* //*OLD*           for(IIvector::iterator perm = perms.begin();perm != perms.end();++perm){
//*OLD* //*OLD*             //if(nc==0) cout << "nccccc" << endl;
//*OLD* //*OLD*             for(IIvector::iterator tup1 = tups1.begin();tup1 != tups1.end();++tup1){
//*OLD* //*OLD*               for(IIvector::iterator tup2 = tups2.begin();tup2 != tups2.end();++tup2){
//*OLD* //*OLD* 		//		cout << newOrder << " mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm " << endl;
//*OLD* //*OLD*                 // Initialize the term's tensor list
//*OLD* //*OLD*                 vector<SQtensor> tensorList;
//*OLD* //*OLD*                 tensorList.insert(tensorList.end(), tensors_except_last_two.begin(),
//*OLD* //*OLD* 				                    tensors_except_last_two.end());
//*OLD* //*OLD* 
//*OLD* //*OLD*                 // Compute the pairs of indices to be contracted
//*OLD* //*OLD*                 // Example: (conPairs[0][p], conPairs[1][p]) is the p-th contraction pair.
//*OLD* //*OLD*                 IIvector conPairs;
//*OLD* //*OLD*                 Ivector p;
//*OLD* //*OLD*                 conPairs.push_back(p);
//*OLD* //*OLD*                 conPairs.push_back(p);
//*OLD* //*OLD*                 for(int i = 0;i < nc;++i){
//*OLD* //*OLD*                   conPairs[0].push_back((*tup1)[(*perm)[i]]);
//*OLD* //*OLD*                   conPairs[1].push_back((*tup2)[i]);
//*OLD* //*OLD* 		} // End i
//*OLD* //*OLD* 
//*OLD* //*OLD* 
//*OLD* //*OLD* //*TEST*                cout << "*TEST*TEST*TES* ......... " << endl;
//*OLD* //*OLD* //*TEST*                if(nc==0){
//*OLD* //*OLD* //*TEST*                  for(size_t i = 0;i < nc;++i) cout << conPairs[0][i] << " ";
//*OLD* //*OLD* //*TEST*                  cout << endl;
//*OLD* //*OLD* //*TEST*		}
//*OLD* //*OLD* 
//*OLD* //*OLD*                 // Initialize the index list for the new excitation operator
//*OLD* //*OLD*                 vector<SQindex*> indexList;
//*OLD* //*OLD*                 SQindex *q = NULL;
//*OLD* //*OLD*                 for(int i = 0;i < 2*newOrder;++i) indexList.push_back(q);
//*OLD* //*OLD* 
//*OLD* //*OLD*                 // Populate the index list for the new excitation operator
//*OLD* //*OLD*                 // Also, create a kronecker delta function for each contraction pair
//*OLD* //*OLD*                 vector<SQindex*> e1_indices(e1.get_indices());
//*OLD* //*OLD*                 vector<SQindex*> e2_indices(e2.get_indices());
//*OLD* //*OLD* 
//*OLD* //*OLD* //*TEST*                for(size_t i = 0;i < e1_indices.size();++i){
//*OLD* //*OLD* //*TEST*                  cout << *e1_indices[i] << endl;
//*OLD* //*OLD* //*TEST*		}
//*OLD* //*OLD* //*TEST*                for(size_t i = 0;i < e2_indices.size();++i){
//*OLD* //*OLD* //*TEST*                  cout << *e2_indices[i] << endl;
//*OLD* //*OLD* //*TEST*		}
//*OLD* //*OLD* 
//*OLD* //*OLD*                 for(int i = 0;i < o1;++i){
//*OLD* //*OLD*                   indexList[i] = e1_indices[i]; //cout << "e1 " << *e1_indices[i] << endl;
//*OLD* //*OLD* 
//*OLD* //*OLD*                   bool i_flag = false;
//*OLD* //*OLD*                   size_t i_index;
//*OLD* //*OLD*                   for(size_t j = 0;j < conPairs[0].size();++j){
//*OLD* //*OLD*                     if(conPairs[0][j] == i) {
//*OLD* //*OLD*                       i_flag = true;
//*OLD* //*OLD*                       i_index = j;
//*OLD* //*OLD* 		    }
//*OLD* //*OLD* 		  } // End j
//*OLD* //*OLD* 
//*OLD* //*OLD*                   if(i_flag){
//*OLD* //*OLD*                     int i1 = i + o1;
//*OLD* //*OLD*                     //cout << "YES 1  " << i_index << " " << o1 << endl;
//*OLD* //*OLD*                     int i2 = conPairs[1][i_index];
//*OLD* //*OLD*                     SQindex* ind1 = e1_indices[i1]; //cout << *ind1 << " ";
//*OLD* //*OLD*                     SQindex* ind2 = e2_indices[i2]; //cout << *ind2 << " ";
//*OLD* //*OLD* 
//*OLD* //*OLD*                     vector<SQindex*> ind_delta;
//*OLD* //*OLD*                     ind_delta.push_back(ind1);
//*OLD* //*OLD*                     ind_delta.push_back(ind2);
//*OLD* //*OLD* 		    //                    tensorList.insert(tensorList.begin(), kDelta(ind_delta));
//*OLD* //*OLD*                     tensorList.push_back(kDelta(ind_delta));
//*OLD* //*OLD* 
//*OLD* //*OLD*                     vector<SQtensor>::iterator it = tensorList.end();--it;
//*OLD* //*OLD*                     //cout << "pushed ... " << *it << endl;
//*OLD* //*OLD* 
//*OLD* //*OLD*                     indexList[i+newOrder]  = e2_indices[o2+i2]; //cout << "e2" << *e2_indices[o2+i2] << "," << newOrder<< endl;
//*OLD* //*OLD* 
//*OLD* //*OLD* //*TEST*                    for(size_t l = 0;l < indexList.first.size();++l) cout << *indexList.first[i] << " ";
//*OLD* //*OLD* //*TEST*                    cout << endl;
//*OLD* //*OLD* 
//*OLD* //*OLD* 		  } // End if
//*OLD* //*OLD*                   else{
//*OLD* //*OLD*                     indexList[i+newOrder]  = e1_indices[i+o1];
//*OLD* //*OLD* 		    //                    indexList.second[i+newOrder] = true;
//*OLD* //*OLD* 		  }
//*OLD* //*OLD* 		} // End i
//*OLD* //*OLD*                 int count = 0;
//*OLD* //*OLD*                 for(int i = 0;i < o2;++i){
//*OLD* //*OLD*                   bool i_flag = false;
//*OLD* //*OLD*                   size_t i_index;
//*OLD* //*OLD*                   for(size_t j = 0;j < conPairs[1].size();++j){
//*OLD* //*OLD*                     if(conPairs[1][j] == i) {
//*OLD* //*OLD*                       i_flag = true;
//*OLD* //*OLD*                       i_index = j;
//*OLD* //*OLD* 		    }
//*OLD* //*OLD* 		  } // End j
//*OLD* //*OLD*                   if(!i_flag){
//*OLD* //*OLD* 		    //                    cout << "YES 1" << endl;
//*OLD* //*OLD*                     indexList[o1+count]  = e2_indices[i];
//*OLD* //*OLD*                     indexList[o1+count+newOrder]  = e2_indices[i+o2];
//*OLD* //*OLD* 		    ++count;
//*OLD* //*OLD* 		  } // End if
//*OLD* //*OLD* 		} // End i
//*OLD* //*OLD* 
//*OLD* //*OLD* //*TEST*                 cout << "XXXXXXX " << endl;
//*OLD* //*OLD* //*TEST*                 for(size_t i = 0;i < tensorList.size();++i) cout << tensorList[i] << endl;
//*OLD* //*OLD* //*TEST*                 cout << "YYYYYYY " << endl;
//*OLD* //*OLD* 
//*OLD* //*OLD* //*TEST*                cout << "newOrder1 " << indexList.first.size() << endl;
//*OLD* //*OLD* //*TEST*                cout << "newOrder2 " << indexList.second.size() << endl;
//*OLD* //*OLD* //*TEST*                vector<SQindex*>::iterator it = indexList.first.begin();
//*OLD* //*OLD* //*TEST*                for(;it != indexList.first.end();++it) cout << "ind .. " << *it << endl;
//*OLD* //*OLD* 
//*OLD* //*OLD*                 // Ensure that all slots in the index list have been filled
//*OLD* //*OLD*                 for(size_t numInd = 0;numInd < indexList.size();++numInd){
//*OLD* //*OLD*                   //if(nc==0) cout << indexList.at(numInd) << ". ";
//*OLD* //*OLD*  
//*OLD* //*OLD*                   if(indexList[numInd] == NULL){
//*OLD* //*OLD*                     cout << "There is at least one unassigned index in the new spin-free unitary generators" << endl;
//*OLD* //*OLD*                     abort();
//*OLD* //*OLD* 		  }
//*OLD* //*OLD* 		} // End numInd
//*OLD* //*OLD* 
//*OLD* //*OLD*                 // Add the new unitary generator into the tensorList
//*OLD* //*OLD*                 tensorList.push_back(sfGen(indexList));
//*OLD* //*OLD* 
//*OLD* //*OLD* //*TEST*                 cout << "XXXXXXX2 " << endl;
//*OLD* //*OLD* //*TEST*                 for(size_t i = 0;i < tensorList.size();++i) cout << tensorList[i] << endl;
//*OLD* //*OLD* //*TEST*                 cout << "YYYYYYY2 " << endl;
//*OLD* //*OLD* 
//*OLD* //*OLD* //*TEST*                vector<SQtensor>::iterator it = tensorList.end();--it;
//*OLD* //*OLD* //*TEST*                cout << *it << endl;
//*OLD* //*OLD* 
//*OLD* //*OLD*                 // Add the resulting term to this iteration's list of output terms
//*OLD* //*OLD*                 iter_output_terms.push_back(SQterm(inTerm.get_numConst(), inTerm.get_Consts(), tensorList));
//*OLD* //*OLD* 
//*OLD* //*OLD*                 if(OF_flag){
//*OLD* //*OLD* 		  vector<SQterm>::iterator iter = iter_output_terms.begin();
//*OLD* //*OLD*                   for(;iter != iter_output_terms.begin();++iter) {iter->contractkDeltas(); iter->decomposeRDM();}
//*OLD* //*OLD*                   iter_output_terms = termChop(iter_output_terms);
//*OLD* //*OLD* 		} // End if
//*OLD* //*OLD* 
//*OLD* //*OLD* //*TEST*                 cout << "ZZZZZZZZZZZZZZZZZ " << endl;
//*OLD* //*OLD* //*TEST*                 for(size_t i = 0;i < iter_output_terms.size();++i) cout << iter_output_terms[i] << endl;
//*OLD* //*OLD* //*TEST*                 cout << "ZZZZZZZZZZZZZZZZZ " << endl;
//*OLD* //*OLD* 
//*OLD* //*OLD* 	      } // End tup2
//*OLD* //*OLD* 	    } // End tup1
//*OLD* //*OLD* 	  } // End perm          
//*OLD* //*OLD* 
//*OLD* //*OLD* 	} // End nc
//*OLD* //*OLD* 
//*OLD* //*OLD*       } // End t
//*OLD* //*OLD* 
//*OLD* //*OLD*       // Set this iteration's list of output terms as the next iteration's input terms
//*OLD* //*OLD*       iter_input_terms = iter_output_terms;
//*OLD* //*OLD*       --n;
//*OLD* //*OLD* 
//*OLD* //*OLD*     } // End while    
//*OLD* //*OLD* 
//*OLD* //*OLD*     // Set the return value as the final iteration's output terms
//*OLD* //*OLD*     outTerms = iter_output_terms;
//*OLD* //*OLD* 
//*OLD* //*OLD*     return outTerms;
//*OLD* //*OLD*   }

}} // Femto::core
