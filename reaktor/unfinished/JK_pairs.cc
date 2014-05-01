//
//  JK_pairs.cc
//  
//
//  Created by Masaaki Saitow on 12/09/05.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <vector>
#include <cmath>
#include <femto.hpp>
#include <SQterm.hpp>
#include <SQtensor.hpp>
#include <SQindex.hpp>
#include <SQreaktor.hpp>

namespace femto {

  // This is applicable only if each term in inTerms_ contains only one integrals
  void SQreaktor::find_JK_pairs(vector<SQterm> &outTerms)
  {
    size_t init_size = inTerms_.size();
    outTerms.reserve(init_size);
    for(vector<SQterm>::iterator t1 = inTerms_.begin();t1 != inTerms_.end();){
      outTerms.push_back(*t1);
      t1 = inTerms_.erase(t1);
      vector<SQterm>::iterator this_t = outTerms.end(); --this_t;
      vector<SQtensor> inTensor = this_t->get_tensors();
      bool Intflag = false;
      for(vector<SQtensor>::iterator t = inTensor.begin();t != inTensor.end();++t)
        if(t->get_name() == name_h2_ || t->get_name() == name_h1_) Intflag = true;

      // t1 may possibly be J, so find K!
      if((int)fabs(this_t->get_numConst())==2 && Intflag){
        for(vector<SQterm>::iterator t2 = inTerms_.begin();t2 != inTerms_.end();){
          //if(*this_t==*t2) { ++t2; continue; }
          if(this_t->get_tensors().size() != t2->get_tensors().size()) { ++t2; continue; }
          bool allOK = true;
          vector<SQtensor> tb = t2->get_tensors();
          vector<SQtensor>::iterator ta = inTensor.begin(); 
          for(vector<SQtensor>::iterator t = tb.begin();t != tb.end();++t,++ta){
            if(t->get_name() == name_h1_ || t->get_name() == name_h2_) continue;
            if(!(*ta==*t)) allOK = false;
	  } // End t
          if(allOK){
            outTerms.push_back(*t2);
            t2 = inTerms_.erase(t2);
	  } // End if
          else
            ++t2;
	} // End t2
      }// End if
      
    } // End t1
    if(outTerms.size() != init_size) {
      cout << "Something is wrong in find_JK_pairs .... " << endl;
      abort();
    }
  }

} // femto::
