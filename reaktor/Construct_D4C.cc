//
//  Constuct_D4C.cc
//  
//
//  Created by Masaaki Saitow on 12/12/11.
//  Copyright (c) 2012 Masaaki Saitow. All rights reserved.
//

#include <algorithm>
#include <fstream>
#include <cmath>
#include <string>
#include <boost/format.hpp>
#include <SQreaktor.hpp>

#define _IND_DEBUG
//#define _FORCE_MODE // *UNRECOMMENDABLE*
#define _REPLACE_D4C // Turn on to replace ERI and 4-RDM by D4C tensor
#define _REPLACE_C4  // Turn on to replace ERI and 4-RDM by C4 tensor
#define _REPLACE_C6  // Turn on to replace ERI and 4-RDM by C6 tensor

using namespace std;
using namespace Femto;
using namespace Femto::Core;

namespace Femto { namespace Reaktor { 

  // *********************************************************
  // Returns symmetry for D4C tensor in Mulliken notation
  // *********************************************************
  Symmetry SQreaktor::d4c_symm()
  {
    Symmetry D4C_symm;

    Femto::Ivector S0; // {0,1,2,3,4,5}
    S0.push_back(0);
    S0.push_back(1); 
    S0.push_back(2); 
    S0.push_back(3);
    S0.push_back(4); 
    S0.push_back(5);   

    Femto::Ivector S1; // {2,3,0,1,4,5}
    S1.push_back(2); 
    S1.push_back(3);
    S1.push_back(0);
    S1.push_back(1); 
    S1.push_back(4); 
    S1.push_back(5);   

    D4C_symm.first.push_back(S0);
    D4C_symm.first.push_back(S1);

    D4C_symm.second.push_back(1);
    D4C_symm.second.push_back(1);

    return D4C_symm;
  }

  // *********************************************************
  // Returns symmetry for C6 tensor in Mulliken notation
  // *********************************************************
  Symmetry SQreaktor::c6_symmM()
  {
    Symmetry C6_symm;
    SQindex i("i", act);
    vector<SQindex*> inds;
    inds.push_back(&i);
    inds.push_back(&i);
    inds.push_back(&i);
    inds.push_back(&i);
    inds.push_back(&i);
    inds.push_back(&i);
    RDM D3(inds);
    D3.convertD2M(); // Convert Dirac 2 Mulliken 

    C6_symm.first = D3.get_perms(); // Not cool, but it's OK!

    for(size_t num = 0;num < C6_symm.first.size();++num) C6_symm.second.push_back(1);

    return C6_symm;
  }

  // *********************************************************
  // 
  // *********************************************************
  void SQreaktor::construct_D4C()
  {
    cout << ">> Contract_D4C is called << " << endl;
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
      cout << "Construct_D4C: LTensors_ are not set. Do process_kDeltas first. if already done, something is wrong." << endl;
      abort();
    } // End if

    if(!inTerms_.size()) return;

    // If ERI or RDM is not of Mulliken form, transform it
    for(vector<SQterm>::iterator t = inTerms_.begin();t != inTerms_.end();++t){
      for(size_t num_t = 0;num_t < t->get_tensors().size();++num_t){
        if(t->get_tensors()[num_t].get_name() == name_h2_ || is_RDM(t->get_tensors()[num_t].get_name()))
          if(t->get_tensors()[num_t].get_notation() == (notation)Dirac) 
	    t->get_tensors_ptr()[num_t]->convertD2M(); // Convert Dirac->Mulliken
      } // End num_t
    } // End t

    pair<int, int> num_cont_d4cs;
    num_cont_d4cs.first  = 0; // Number of D4C tensor
    num_cont_d4cs.second = 0; // Number of C4 tensor
    int num_c6s(0);           // Number of C6 tensor
    for(vector<SQterm>::iterator t = inTerms_.begin();t != inTerms_.end();++t){
      int num_commact(0);
      size_t V2_num(-1);
      size_t Fk_num(-1);
      size_t D4_num(-1);
      for(size_t num_t = 0;num_t < t->get_tensors().size();++num_t){
        if     (t->get_tensors()[num_t].get_name() == name_h2_)    V2_num = num_t;
        else if(t->get_tensors()[num_t].get_name() == name_d4_)    D4_num = num_t;
        else if(t->get_tensors()[num_t].get_name() == Fock_name()) Fk_num = num_t;
      } // End num_t

      // In case of V2 and D4 have several indices in common
      if(V2_num != -1 && D4_num != -1){ 
        vector<SQindex*> D4inds; D4inds.reserve(t->get_tensors()[D4_num].get_indices().size());
        vector<SQindex*> V2inds; V2inds.reserve(t->get_tensors()[V2_num].get_indices().size());

	for(size_t num_i = 0;num_i < t->get_tensors()[D4_num].get_indices().size();++num_i) 
	  if(find(D4inds.begin(), D4inds.end(), t->get_tensors()[D4_num].get_indices()[num_i]) == D4inds.end()) 
	    D4inds.push_back(t->get_tensors()[D4_num].get_indices()[num_i]);
	for(size_t num_i = 0;num_i < t->get_tensors()[V2_num].get_indices().size();++num_i) 
	  if(find(V2inds.begin(), V2inds.end(), t->get_tensors()[V2_num].get_indices()[num_i]) == V2inds.end()) 
	    V2inds.push_back(t->get_tensors()[V2_num].get_indices()[num_i]);

        // Count how many active indices are shared by D4 and V2
        for(vector<SQindex*>::iterator i = V2inds.begin();i != V2inds.end();++i) 
          if((*i)->get_char() == act && find(D4inds.begin(), D4inds.end(), *i) != D4inds.end()) ++num_commact;
      } // End if
      // In case of common Fock and D4 have several indices in common
      else if(Fk_num != -1 && D4_num != -1){
        vector<SQindex*> D4inds; D4inds.reserve(t->get_tensors()[D4_num].get_indices().size());
        vector<SQindex*> Fkinds; Fkinds.reserve(t->get_tensors()[Fk_num].get_indices().size());

	for(size_t num_i = 0;num_i < t->get_tensors()[D4_num].get_indices().size();++num_i) 
	  if(find(D4inds.begin(), D4inds.end(), t->get_tensors()[D4_num].get_indices()[num_i]) == D4inds.end()) 
	    D4inds.push_back(t->get_tensors()[D4_num].get_indices()[num_i]);
	for(size_t num_i = 0;num_i < t->get_tensors()[Fk_num].get_indices().size();++num_i) 
	  if(find(Fkinds.begin(), Fkinds.end(), t->get_tensors()[Fk_num].get_indices()[num_i]) == Fkinds.end()) 
	    Fkinds.push_back(t->get_tensors()[Fk_num].get_indices()[num_i]);

        // Count how many active indices are shared by D4 and Fk
        for(vector<SQindex*>::iterator i = Fkinds.begin();i != Fkinds.end();++i) 
          if((*i)->get_char() == act && find(D4inds.begin(), D4inds.end(), *i) != D4inds.end()) ++num_commact;        
      } // End if

      // If *t has both V2 and D4 with three active indices in common, they may possibly be contractable to give D4C(o,o,o,o,o,*) type contracted D4 tensor
      // e.g.)
      //
      //  D4C(o,o,o,o,o,*) <--- D4(o,o,o,o,o,o1,o2,o3) V(*,o1,o2,o3)
      // 
      // Farewell to costly and bitchy, cumulant !
      /*- In case of D4C tensor -*/ 
      if(num_commact == 3 && V2_num != -1 && D4_num != -1 && Fk_num == -1){

        vector<SQtensor>::iterator V2_ptr; // Points where V2 is in tens
        vector<SQtensor>::iterator D4_ptr; // Points where D4 is in tens
        vector<SQtensor> tens(t->get_tensors());
        for(vector<SQtensor>::iterator tt = tens.begin();tt != tens.end();++tt){
          if     (tt->get_name() == name_h2_) V2_ptr = tt;
          else if(tt->get_name() == name_d4_) D4_ptr = tt;
	} // End tt

        // First, extract three active indices in common with D4 and V2
        vector<SQindex*> V2inds(V2_ptr->get_indices());
	vector<SQindex*> D4inds(D4_ptr->get_indices());
        vector<SQindex*> comminds;
        SQindex* leftind(NULL);
        for(vector<SQindex*>::iterator i = V2inds.begin();i != V2inds.end();++i){
          if(find(D4inds.begin(), D4inds.end(), *i) != D4inds.end()) comminds.push_back(*i);
	  else leftind = *i;
	} // End i
        sort(comminds.begin(), comminds.end());

#ifdef _IND_DEBUG
        cout << "++ size D4inds : " << D4inds.size() << endl; //*TEST* 
        for(vector<SQindex*>::iterator i = D4inds.begin();i != D4inds.end();++i)
          cout << " " << **i << ", " << *i << endl; 
	cout << endl; 
        cout << "++ size V2inds : " << V2inds.size() << endl; //*TEST* 
        for(vector<SQindex*>::iterator i = V2inds.begin();i != V2inds.end();++i)
          cout << " " << **i << ", " << *i << endl; 
	cout << endl;
        cout << "++ size commin : " << comminds.size() << endl; //*TEST* 
        for(vector<SQindex*>::iterator i = comminds.begin();i != comminds.end();++i)
          cout << " " << **i << ", " << *i << endl;         
#endif

	//:::::::::::::::::::::::::::: NEW ::::::::::::::::::::::::::::::::::
	bool set_OK(false);
        IIvector D4perms(D4_ptr->get_perms());
        IIvector V2perms(V2_ptr->get_perms());
        for(size_t k = 0;k < D4perms.size();++k){
          vector<SQindex*> D4_inds_part;
          D4_inds_part.push_back(D4inds[D4perms[k][5]]);
          D4_inds_part.push_back(D4inds[D4perms[k][6]]);
          D4_inds_part.push_back(D4inds[D4perms[k][7]]);
          sort(D4_inds_part.begin(), D4_inds_part.end());
          // Need D4(o,o,o,o,o,o1,o2,o3) V2(g,o1,o2,o3) like conf.
          if(D4_inds_part == comminds){
	    for(size_t l = 0;l < V2perms.size();++l){

	      if((D4inds[D4perms[k][5]] == V2inds[V2perms[l][1]] &&
		  D4inds[D4perms[k][6]] == V2inds[V2perms[l][2]] &&
		  D4inds[D4perms[k][7]] == V2inds[V2perms[l][3]])){

		D4_ptr->rotateIndices(k); // -> V2(g,o1,o2,o3)
		V2_ptr->rotateIndices(l); // -> D4(o,o,o,o,o,o1,o2,o3)

#ifdef _IND_DEBUG
            cout << " -- Indices of ERI and D4 are rotated to form D4C :: " << endl;
	    cout << *D4_ptr << endl;
	    cout << *V2_ptr << endl;
#endif
	        set_OK = true;
		break;

              } // End if
	    } // End l
            if(set_OK) break;
	  } // End if

	} // End k
	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        if(!set_OK){
	  cout << endl;
          cout << "Construct_D4C: Something is wrong in reordering indices of ERI and D4 << " << *t << ">> for D4C contracted tensor" << endl;
#ifndef _FORCE_MODE
          abort();
#endif
	} // End if

        // If everything is OK, form D4C
        if(set_OK){
	  vector<SQindex*> D4inds_new(D4_ptr->get_indices());
	  vector<SQindex*> V2inds_new(V2_ptr->get_indices());
	  D4inds_new.erase(D4inds_new.begin()+5, D4inds_new.end());
	  D4inds_new.push_back(V2inds_new.at(0)); // This should be the new indices of D4C !

	  vector<SQtensor> new_ten; new_ten.reserve(tens.size());
#ifdef _REPLACE_D4C
	  for(vector<SQtensor>::iterator tt = tens.begin();tt != tens.end();++tt){
	    if (tt != V2_ptr && tt != D4_ptr) new_ten.push_back(*tt);
	  } // End tt
	  new_ten.push_back(SQtensor(D4C_name(), D4inds_new, d4c_symm(), true, Mulliken));
#else
	  for(vector<SQtensor>::iterator tt = tens.begin();tt != tens.end();++tt) new_ten.push_back(*tt);
#endif

#ifdef _REPLACE_D4C	  
	  cout << " -- One 4-RDM and ERI are contracted to form D4C :: " << endl;
#endif
	  cout << "    + Before: " << *t << endl;
	  t->set_tensors(new_ten);
	  cout << "    + After : " << *t << endl << endl;
          num_cont_d4cs.first++;                         

	} // End if

      } // End if

      // If *t has both V2 and D4 with four active indices in common, they may possibly be contractable to give C4(o,o,o,o) type contracted D4 tensor
      // e.g.)
      //
      //  C4(o,o,o,o) <--- D4(o,o,o,o,o1,o2,o3,o4) V(o1,o2,o3,o4)
      // 
      // Farewell to costly and bitchy, cumulant !     
      /*- In case of C4 tensor  -*/
      else if(num_commact == 4 && V2_num != -1 && D4_num != -1 && Fk_num == -1){
        vector<SQtensor>::iterator V2_ptr; // Points where V2 is in tens
        vector<SQtensor>::iterator D4_ptr; // Points where D4 is in tens
        vector<SQtensor> tens(t->get_tensors());
        for(vector<SQtensor>::iterator tt = tens.begin();tt != tens.end();++tt){
          if     (tt->get_name() == name_h2_) V2_ptr = tt;
          else if(tt->get_name() == name_d4_) D4_ptr = tt;
	} // End tt
        
        // First, take all the indices of V2 (all of which are assumed to be active and commonly shared by D4)
        vector<SQindex*> V2inds(V2_ptr->get_indices()); // <-- [o1,o2,o3,o4]
	vector<SQindex*> D4inds(D4_ptr->get_indices());
        vector<SQindex*> comminds;
        for(vector<SQindex*>::iterator i = V2inds.begin();i != V2inds.end();++i){
          if((*i)->get_char() != act) {
            cout << "Construct_D4C: Non-active index detected in ERI in forming C4 contracted tensor" << endl;
            abort();
	  } // End if
	  else comminds.push_back(*i);
	} // End i
        sort(comminds.begin(), comminds.end());

	//:::::::::::::::::::::::::::: NEW ::::::::::::::::::::::::::::::::::
	bool set_OK(false);
        IIvector D4perms(D4_ptr->get_perms());
        IIvector V2perms(V2_ptr->get_perms());
        for(size_t k = 0;k < D4perms.size();++k){
          vector<SQindex*> D4_inds_part;
          D4_inds_part.push_back(D4inds[D4perms[k][4]]);
          D4_inds_part.push_back(D4inds[D4perms[k][5]]);
          D4_inds_part.push_back(D4inds[D4perms[k][6]]);
          D4_inds_part.push_back(D4inds[D4perms[k][7]]);
          sort(D4_inds_part.begin(), D4_inds_part.end());
          // Need D4(o,o,o,o,o1,o2,o3,o4) V2(o1,o2,o3,o4) like conf.
          if(D4_inds_part == comminds){
	    for(size_t l = 0;l < V2perms.size();++l){

	      if((D4inds[D4perms[k][4]] == V2inds[V2perms[l][0]] &&
		  D4inds[D4perms[k][5]] == V2inds[V2perms[l][1]] &&
		  D4inds[D4perms[k][6]] == V2inds[V2perms[l][2]] &&
		  D4inds[D4perms[k][7]] == V2inds[V2perms[l][3]])){

		D4_ptr->rotateIndices(k); // -> D4(o,o,o,o,o1,o2,o3,o4) 
		V2_ptr->rotateIndices(l); // -> V2(o1,o2,o3,o4)

#ifdef _IND_DEBUG
            cout << " -- Indices of ERI and D4 are rotated to form C4 :: " << endl;
	    cout << *D4_ptr << endl;
	    cout << *V2_ptr << endl;
#endif
	        set_OK = true;
		break;

              } // End if
	    } // End l
            if(set_OK) break;
	  } // End if

	} // End k
	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        if(!set_OK){
	  cout << endl;
          cout << "Construct_D4C: Something is wrong in reordering indices of ERI and D4 << " << *t << ">> for C4 contracted tensor" << endl;
#ifndef _FORCE_MODE
          abort();
#endif
	} // End if

        // At last, orm the C4 tensor to update the inTerms_ (if everything is OK)
        if(set_OK){
	  vector<SQindex*> D4inds_new(D4_ptr->get_indices());
	  D4inds_new.erase(D4inds_new.begin()+4, D4inds_new.end());

	  vector<SQtensor> new_ten; new_ten.reserve(tens.size());
#ifdef _REPLACE_C4
	  for(vector<SQtensor>::iterator tt = tens.begin();tt != tens.end();++tt){
	    if (tt != V2_ptr && tt != D4_ptr) new_ten.push_back(*tt);
	  } // End tt
	  new_ten.push_back(SQtensor(C4_name(), D4inds_new, Femto::h2_symmM(), true, Mulliken));
#else
	  for(vector<SQtensor>::iterator tt = tens.begin();tt != tens.end();++tt) new_ten.push_back(*tt);
#endif

#ifdef _REPLACE_C4	  
	  cout << " -- One 4-RDM and V2 are contracted to form C4 :: " << endl;
#endif
	  cout << "    + Before: " << *t << endl;
	  t->set_tensors(new_ten);
	  cout << "    + After : " << *t << endl << endl;        
          num_cont_d4cs.second++;                                          

	} // End if

      } // End if

      // If *t has both Fock and D4 with four active indices in common, they may possibly be contractable to give C6(o,o,o,o,o,o) type contracted D4 tensor
      // e.g.)
      //
      //  C6(o,o,o,o,o,o) <--- D4(o,o,o,o,o,o,o1,o2) F(o1,o2)
      // 
      // Farewell to costly and bitchy, cumulant !     
      /*- In case of C6 tensor  -*/
      else if(num_commact == 2 && V2_num == -1 && D4_num != -1 && Fk_num != -1){
        vector<SQtensor>::iterator Fk_ptr; // Points where Fock is in tens
        vector<SQtensor>::iterator D4_ptr; // Points where D4   is in tens
        vector<SQtensor> tens(t->get_tensors());
        for(vector<SQtensor>::iterator tt = tens.begin();tt != tens.end();++tt){
          if     (tt->get_name() == Fock_name()) Fk_ptr = tt;
          else if(tt->get_name() == name_d4_)    D4_ptr = tt;
	} // End tt
        
        // First, take all the indices of Fock (all of which are assumed to be active and commonly shared by D4)
        vector<SQindex*> Fkinds(Fk_ptr->get_indices()); // <-- [o1,o2]
	vector<SQindex*> D4inds(D4_ptr->get_indices());
        vector<SQindex*> comminds;
        for(vector<SQindex*>::iterator i = Fkinds.begin();i != Fkinds.end();++i){
          if((*i)->get_char() != act) {
            cout << "Construct_D4C: Non-active index detected in ERI in forming C6 contracted tensor" << endl;
            abort();
	  } // End if
	  else comminds.push_back(*i);
	} // End i
        sort(comminds.begin(), comminds.end());

	//:::::::::::::::::::::::::::: NEW ::::::::::::::::::::::::::::::::::
	bool set_OK(false);
        IIvector D4perms(D4_ptr->get_perms());
        IIvector Fkperms(Fk_ptr->get_perms());
        for(size_t k = 0;k < D4perms.size();++k){
          vector<SQindex*> D4_inds_part;
          D4_inds_part.push_back(D4inds[D4perms[k][6]]);
          D4_inds_part.push_back(D4inds[D4perms[k][7]]);
          sort(D4_inds_part.begin(), D4_inds_part.end());
          // Need D4(o,o,o,o,o,o,o1,o2) Fock(o1,o2) like conf.
          if(D4_inds_part == comminds){
	    for(size_t l = 0;l < Fkperms.size();++l){

	      if((D4inds[D4perms[k][6]] == Fkinds[Fkperms[l][0]] &&
		  D4inds[D4perms[k][7]] == Fkinds[Fkperms[l][1]])){

		D4_ptr->rotateIndices(k); // -> D4(o,o,o,o,o,o,o1,o2)
		Fk_ptr->rotateIndices(l); // -> Fk(o1,o2) 

#ifdef _IND_DEBUG
            cout << " -- Indices of ERI and D4 are rotated to form C6 :: " << endl;
	    cout << *D4_ptr << endl;
	    cout << *Fk_ptr << endl;
#endif
	        set_OK = true;
		break;

              } // End if
	    } // End l
            if(set_OK) break;
	  } // End if

	} // End k
	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        if(!set_OK){
	  cout << endl;
          cout << "Construct_D4C: Something is wrong in reordering indices of Fock matrix and D4 << " << *t << ">> for C6 contracted tensor" << endl;
#ifndef _FORCE_MODE
          abort();
#endif
	} // End if

        // At last, orm the C6 tensor to update the inTerms_ (if everything is OK)
        if(set_OK){
	  vector<SQindex*> D4inds_new(D4_ptr->get_indices());
	  D4inds_new.erase(D4inds_new.begin()+6, D4inds_new.end());

	  vector<SQtensor> new_ten; new_ten.reserve(tens.size());
#ifdef _REPLACE_C4
	  for(vector<SQtensor>::iterator tt = tens.begin();tt != tens.end();++tt){
	    if (tt != Fk_ptr && tt != D4_ptr) new_ten.push_back(*tt);
	  } // End tt
	  new_ten.push_back(SQtensor(C6_name(), D4inds_new, c6_symmM(), true, Mulliken));
#else
	  for(vector<SQtensor>::iterator tt = tens.begin();tt != tens.end();++tt) new_ten.push_back(*tt);
#endif

#ifdef _REPLACE_C6	  
	  cout << " -- One 4-RDM and V2 are contracted to form C6 :: " << endl;
#endif
	  cout << "    + Before: " << *t << endl;
	  t->set_tensors(new_ten);
	  cout << "    + After : " << *t << endl << endl;
          num_c6s++;

	} // End if        

      } // End if

    } // End t

    cout << boost::format("  -- %d terms are replaced to form d4c tensor ") % num_cont_d4cs.first  << endl;
    cout << boost::format("  -- %d terms are replaced to form c4 tensor  ") % num_cont_d4cs.second << endl;
    cout << boost::format("  -- %d terms are replaced to form c6 tensor  ") % num_c6s              << endl << endl;

  }                                               

}} //End Femto
