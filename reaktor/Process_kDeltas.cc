//
//  Process_kDeltas.cc
//  
//
//  Created by Masaaki Saitow on 12/11/08.
//  Copyright (c) 2012 Masaaki Saitow. All rights reserved.
//

#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <SQreaktor.hpp>
#include <boost/format.hpp>

#define _KD_DEBUG

using namespace std;
using namespace Femto;
using namespace Femto::Core;

namespace Femto { namespace Reaktor { 

  // *********************************************************
  // 
  // *********************************************************
  void SQreaktor::process_kDeltas()
  {

    cout << ">> Check there is any unprocessed Kronecker's delta .... <<" << endl << endl;

    // This assumes threre are no Kronecker deltas with indices of different orbital type,
    // and those of dummy indices
    LTensors_.reserve(inTerms_.size());
    int numTerm(0);
    for(vector<SQterm>::iterator t = inTerms_.begin();t != inTerms_.end();++t,++numTerm){
      vector<SQtensor> tens(t->get_tensors());

      // Preserve LTensor_ for this term (In fact, it is suficient that only indices of LTensor_ are preserved)
      vector<SQtensor> temp;
      temp.push_back(LTensor_);
      vector<string> coeff;
      LTensors_.push_back(SQterm(1.0, coeff, temp));

#ifdef _FOR_PUBLICATION
      cout << "THIS EXECUTION IS JUST FOR PREPARATION OF THE SUPPLIMENTAL MATERIAL" << endl;
      //return;
#else
      for(vector<SQtensor>::iterator tt = tens.begin();tt != tens.end();){
        if(tt->get_name() == kDelta_name()){
          if(tt->get_indices()[0]->get_isSummed() || tt->get_indices()[1]->get_isSummed()){
            cout << "Process_kDeltas: There's kDelta with dummy indices in " << *t << endl;
            abort();
	  } // End if
          else if(tt->get_indices()[0]->get_char() != tt->get_indices()[1]->get_char()){
            cout << "Process_kDeltas: There's kDelta which is definitely zero" << *t << endl;
            abort();
	  } // End if
          // For LTensor (this should be the first)
	  {
            vector<SQindex*> Linds(LTensors_.back().get_summedBody());
            vector<SQindex*> Dinds(tt->get_indices());
	    SQindex* killer(NULL); // Points index to kill
	    SQindex* killed(NULL); // Points index to be killed
	    for(vector<SQindex*>::iterator i = Dinds.begin();i != Dinds.end();++i)
	      for(vector<SQindex*>::iterator j = Linds.begin();j != Linds.end();++j)
		if     (**i == **j && *i == Dinds[0] && *Dinds[0] < *Dinds[1]) killer = *j;
		else if(**i == **j && *i == Dinds[0] && *Dinds[0] > *Dinds[1]) killed = *j;
		else if(**i == **j && *i == Dinds[1] && *Dinds[0] < *Dinds[1]) killed = *j;
		else if(**i == **j && *i == Dinds[1] && *Dinds[0] > *Dinds[1]) killer = *j;
	    if(killer == NULL || killed == NULL){
	      cout << "Binary_Decomposition: Algorithmic error occured for .... " << endl;
	      cout << "  " << LTensors_.back() << " += " << *t << endl;
              if(killer == NULL) cout << "killer not found" << endl;
              else               cout << "killer : " << *killer << endl;
              if(killed == NULL) cout << "killed not found" << endl;
              else               cout << "killed : " << *killed << endl; 
	      abort();
	    } // End if

#ifdef _KD_DEBUG
            cout << " >> Debug mode of Process_kDeltas << " << endl;
	    cout << LTensors_.back() << " += " << *t << endl;
            cout << "killed : " << *killed << ", killer : " << killer << endl;
#endif
	    *killed = *killer;
	  } // End scope
          // For this term
	  {
	    SQindex* killer((*(tt->get_indices()[0]) > *(tt->get_indices()[1]) ? (tt->get_indices()[1]) : (tt->get_indices()[0])));
	    SQindex* killed((*(tt->get_indices()[0]) > *(tt->get_indices()[1]) ? (tt->get_indices()[0]) : (tt->get_indices()[1])));
	    *killed = *killer;   // Kill index to be killed
	  } // End scope

#ifdef _KD_DEBUG
          cout << "< After KD is killed > " << endl;
          cout << LTensors_.back() << endl; //*TEST*
#endif
 
          tt = tens.erase(tt); // Then, kill the kDelta
	} // End if
        else ++tt;
      } // End tt
      if(tens.size() != t->get_tensors().size()) t->set_tensors(tens);
#endif
    } // End t

#ifdef _KD_DEBUG
    if(LTensors_.size() != inTerms_.size()){
      cout << "Process_kDeltas: Something is wrong" << endl;
      cout << "length of LTensors_ : " << LTensors_.size() << endl;
      cout << "length of inTerms_  : " << inTerms_.size() << endl;
      abort();
    } // End if
    int count(0);
    for(vector<SQterm>::iterator t = inTerms_.begin();t != inTerms_.end();++t,++count){
      cout << count << " : "<< LTensors_[count].get_tensors()[0] << " += " << *t << endl;;
    } // End t
#endif
  }

}} // Femto::
