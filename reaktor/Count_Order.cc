//
//  Count_Order.cc
//  
//
//  Created by Masaaki Saitow on 12/10/26.
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
  void SQreaktor::count_order(vector<SQbinary> &bin)
  {

    if(bin.size() >= 3){
      cout << "Size of the input vector is expected to be less than 3" << endl;
      abort();
    } // End if

    cout << endl;
    cout << "1: " << bin[0] << endl; 
    if(bin.size() == 2) 
    cout << "2: " << bin[1]  << endl;
    cout << endl;  
    
    vector<SQindex*> cont1(bin[0].get_summedBody());
    vector<SQindex*> cont2;
    if(bin.size() == 2) cont2 = bin[1].get_summedBody();
    vector<SQindex*> conti(bin[0].get_Ltensor().get_indices());
    vector<int> o1(3); 
    vector<int> o2(3); 
    vector<int> oi(3); 

    for(vector<SQindex*>::iterator i1 = cont1.begin();i1 != cont1.end();++i1)
      if     ((*i1)->get_char() == core) ++o1[0];
      else if((*i1)->get_char() == act ) ++o1[1];
      else if((*i1)->get_char() == virt) ++o1[2];

    for(vector<SQindex*>::iterator i2 = cont2.begin();i2 != cont2.end();++i2)
      if     ((*i2)->get_char() == core) ++o2[0];
      else if((*i2)->get_char() == act ) ++o2[1];
      else if((*i2)->get_char() == virt) ++o2[2];

    for(vector<SQindex*>::iterator ii = conti.begin();ii != conti.end();++ii)
      if     ((*ii)->get_char() == core) ++oi[0];
      else if((*ii)->get_char() == act ) ++oi[1];
      else if((*ii)->get_char() == virt) ++oi[2];

    vector<int> max_o(3);
    if(pow(numCore_, o1[0])*pow(numOcc_, o1[1])*pow(numVirt_, o1[2]) > pow(numCore_, o2[0])*pow(numOcc_, o2[1])*pow(numVirt_, o2[2])){
      max_o = o1;
    } // End if
    else max_o = o2;

    // Print all the information on the contraction ....
    string scaling(  "! Scaling       : O("); 
    string size_mem("! Max size of X : ");
    if(max_o[0]){
      ostringstream stm;
      stm << max_o[0];
      scaling += "c^" + stm.str();
    }
    if(max_o[1]){
      ostringstream stm;
      stm << max_o[1];
      scaling += "o^" + stm.str();
    }
    if(max_o[2]){
      ostringstream stm;
      stm << max_o[2];
      scaling += "v^" + stm.str();
    }
    scaling += ")";
    if(oi[0]){
      ostringstream stm;
      stm << oi[0];
      size_mem += "c^" + stm.str();
    }
    if(oi[1]){
      ostringstream stm;
      stm << oi[1];
      size_mem += "o^" + stm.str();
    }
    if(oi[2]){
      ostringstream stm;
      stm << oi[2];
      size_mem += "v^" + stm.str();
    }
    cout << scaling << endl;
    if(bin.size() == 2) cout << size_mem << endl;
    cout << endl;
    
  }                                               

}} //End Femto
