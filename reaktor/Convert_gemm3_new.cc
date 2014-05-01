//
//  Convert_gemm3_new.cc
//  
//
//  Created by Masaaki Saitow on 13/11/1.
//  Copyright (c) 2012 Masaaki Saitow. All rights reserved.
//

#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <SQreaktor.hpp>

#include <boost/format.hpp>
#include <boost/tuple/tuple.hpp>

#define _DEBUG_BLAS // Debug option for use of blas
#define _DEBUG_OPT  // Debug option for optimization 

#define _DGEMM_OPT // Whether to do optimization or not
#define _OPT_RANK 3 // Rank of RDM for which the loops are optimized

using namespace std;
using namespace Femto;
using namespace Femto::Core;

namespace Femto { namespace Reaktor {

  // *********************************************************
  // 
  // *********************************************************
  void SQreaktor::binary_contract3_new(SQbinary &bin, string title, ofstream &f, string s, contDecl c, bool isBareLHS)
  {

    if(bin.get_Rtensors().size() != 2){
      cout << "Number of terms pushed into binary_contract should be 2." << endl;
      abort();
    } 

    auto ExtInd (c["DecInd"]   );
    auto Consts (c["DecConst"] );
    auto NameTen(c["DecTensor"]);

    f <<  "\n"                                                                              ;     
    f <<  "!                  >> binary_contract3_new subroutine is called <<            \n";
    f <<  "! ---------------------------- Parameters used -------------------------------\n" ;
    f <<  "!                                                                             \n" ;    
    f <<  "! Whether the LHS is a BareAmpPack ....... " << (isBareLHS ? "Yes" : "No") <<"\n" ;
    f <<  "! Name of ERI ............................ " + name_h2_                    +  "\n" ;
    f <<  "! Name of BareAmpPack appearing in RHS.... " + name_amp_                   +  "\n" ;
    f <<  "!                                                                             \n" ;
    f <<  "!-----------------------------------------------------------------------------\n" ;  
    
    string printName("");
    printName  = "! subroutine " + title + "(";
    string printName2("subroutine " + title + " &\n  (");
    
    string printEnd("end subroutine " + title + "\n\n");

    string DecComm8;
    DecComm8  = "! FEMTO BEGIN  **************************************************************\n";
    DecComm8 += "use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6\n";
    DecComm8 += "\n";
    if(use_omp_)
    DecComm8 += "use omp_lib\n";
    DecComm8 += "implicit none\n\n";
    
    string DecExtInd("");
    auto tensors(bin.get_Rtensors());
    auto summedBody(bin.get_summedBody());
    auto LTinds(bin.get_Ltensor().get_indices());
    for(auto i = LTinds.begin();i != LTinds.end();++i){
      bool i_flag(true);
      for(auto j = summedBody.begin();j != summedBody.end();++j)
        if(**i == **j) i_flag = false;
      if(i_flag) summedBody.push_back(*i);
    } // End i

    for(auto i = summedBody.begin();i != summedBody.end();++i){
      if((*i)->get_isExt()) {
        DecExtInd += "integer, intent(in) :: i_" +  (*i)->get_index() + ", s_" + (*i)->get_index() + "\n";
        printName += "i_" + (*i)->get_index() + ", s_" + (*i)->get_index() + ", ";
      } // End if
      //f << "! OReORE : " << *i << endl; //TEST
    } // End i

    string DecSymm("");
    DecSymm  = "! Information of the Irreps ....\n" ;
    DecSymm += "integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)\n";  
    DecSymm += "! Flop count\n";
    DecSymm += "real*8, intent(inout) :: flops\n";

    // Print all the names of parameters 
    string DecConsts("");
    if(Consts.size())
      DecConsts = "! Declaration of numerical constants .... \n";
    for(auto s = Consts.begin();s != Consts.end();++s){
      DecConsts += "real(kind=8), intent(inout) :: " + *s;
      printName += *s + ", ";
    } // End s

    // Set-up the X_indices, which is true associated indices to intermediate array.
    vector<SQindex*>         X_indicesL; // X_indices for the left-hand side tensor
    vector<vector<SQindex*> >X_indicesR; // X_indices for the first-right hand side tensor
    X_indicesR.resize(2);
    X_indicesL = bin.get_LInnerIndices();
    for(size_t i = 0;i < bin.get_Rtensors().size();++i) X_indicesR[i] = bin.get_RInnerIndices(i);

    // Print all the names of tensors
    map<string, int> mapTen; // Correspondence between tensor name and number of the associated internal indices
    if(isBareLHS){
      mapTen.insert( map<string, int>::value_type(bin.get_Ltensor().get_name(), 3) );
    } // End if
    else if(is_Interm(bin.get_Ltensor().get_name())){
      mapTen.insert( map<string, int>::value_type(bin.get_Ltensor().get_name(), X_indicesL.size()) );
    } // End if
    else if(bin.get_Ltensor().get_name() == D4C_nameL()){
      mapTen.insert( map<string, int>::value_type(bin.get_Ltensor().get_name(), 5) );
    } // End if
    else{
      mapTen.insert( map<string, int>::value_type(bin.get_Ltensor().get_name(), bin.get_Ltensor().get_indices().size()) );
    } // End else
    int myNum(0);
    for(auto t = tensors.begin();t != tensors.end();++t,++myNum){
      auto key(mapTen.find(t->get_name()));
      if(key == mapTen.end()){
        if(t->get_name() == name_amp_ || t->get_name() == name_h2_){
          mapTen.insert( map<string, int>::value_type(t->get_name(), 3) );
	} // End if
        else if(t->get_name() == RDM_name() + "4"){
          mapTen.insert( map<string, int>::value_type(t->get_name(), 6) );
	} // End if
        else if(is_Interm(t->get_name())){
	  mapTen.insert( map<string, int>::value_type(t->get_name(), X_indicesR[myNum].size()) );
	} // End if
        else if(is_D4C(t->get_name())){
          mapTen.insert( map<string, int>::value_type(t->get_name(), 5) );
	} // End if
        else{
          mapTen.insert( map<string, int>::value_type(t->get_name(), t->get_indices().size()) );
	}
      } // End if
    } // End t

    // Declration of the external variables
    for(vector<string>::iterator i = ExtInd.begin();i != ExtInd.end();++i)
      printName2 += "s_" + *i + ", i_" + *i + ", ";
    for(vector<string>::iterator C = Consts.begin();C != Consts.end();++C)
      printName2 +=  *C + ", ";    
    for(vector<string>::iterator t = NameTen.begin();t != NameTen.end();++t)
      printName2 += *t + "_, ";

    string DecTensor("\n! Declare tensors used ...\n");
    string fraction;
    // Guys in module (such as RDMs) are declared here
    // This bunch of codes should be same to Convert2.cc.
    for(auto t = mapTen.begin();t != mapTen.end();++t){
      if(is_RDM(t->first) || t->first == "Fc1" || is_C4(t->first) || is_C6(t->first) || t->first == C2_name()) printName2 += t->first + "_, ";
      size_t dimTen(t->second);
      printName += t->first + "_, ";
      if(!dimTen){
        DecTensor += "real(kind=8)                   :: " + t->first + "_\n";
        continue;
      } // End if
      else{
        ostringstream stm;
        stm << dimTen;
        fraction = "type(symblock" + stm.str() + "), intent(inout) :: " + t->first + "_(";
      } // End else
      for(size_t i = 0;i < dimTen;++i){
        fraction += "0:nir-1";
        if(i == dimTen-1) { fraction += ")"; DecTensor += fraction + "\n"; }
        else              fraction += ", ";
      } // End i
    } // End t
    DecTensor += "\n";
    f << endl;

    printName  += "nir, nsym, psym, flops)\n\n";
    printName2 += "nir, nsym, psym, flops)\n\n";

    // Determine the order of the tensors in dgemm
    vector<vector<SQindex*> > X_indicesR2;
    X_indicesR2.resize(2);    
    vector<SQtensor> the_tensors;    
    vector<SQindex*> L_nd_ind;
    vector<SQindex*> LTind(bin.get_Ltensor().get_indices());
    for(vector<SQindex*>::reverse_iterator i = LTind.rbegin();i != LTind.rend();++i){
      if(find(L_nd_ind.begin(), L_nd_ind.end(), *i) == L_nd_ind.end() &&!(*i)->get_isExt()) 
        L_nd_ind.push_back(*i);
    } // End n_ind
    if(!LTind.size() || !L_nd_ind.size()){
      the_tensors.push_back(bin.get_Rtensors()[0]);
      the_tensors.push_back(bin.get_Rtensors()[1]);
      X_indicesR2[0] = X_indicesR[0];
      X_indicesR2[1] = X_indicesR[1];
    } // End if
    else{
      vector<SQindex*> temp_i(bin.get_Rtensors()[0].get_indices());
      if(find(temp_i.begin(), temp_i.end(), L_nd_ind[0]) != temp_i.end()){
        the_tensors.push_back(bin.get_Rtensors()[0]);
        the_tensors.push_back(bin.get_Rtensors()[1]);
	X_indicesR2[0] = X_indicesR[0];
	X_indicesR2[1] = X_indicesR[1];
      } // End if
      else{
        the_tensors.push_back(bin.get_Rtensors()[1]);
        the_tensors.push_back(bin.get_Rtensors()[0]);
	X_indicesR2[1] = X_indicesR[0];
	X_indicesR2[0] = X_indicesR[1];
      } // End if
    } // End else

    // If Z1 is composed of only external indices, swap Z1 and Z2
    bool external_only = true;
    for(size_t num_i = 0;num_i < the_tensors[0].get_indices().size();++num_i)
      if(!the_tensors[0].get_indices()[num_i]->get_isExt()) external_only = false;
    if(external_only){
      SQtensor temp(the_tensors[1]);
      the_tensors[1] = the_tensors[0];
      the_tensors[0] = temp;
    } // End if

    // Dummy indices to be gemmed!
    vector<SQindex*> rowInd;    // m in dgemm
    vector<SQindex*> colInd;    // n in dgemm
    vector<SQindex*> summedInd; // k in dgemm
    vector<SQindex*> extInd;    // non-dummy indices shared by several tensors
    vector<SQindex*> I1(the_tensors[0].get_indices());
    vector<SQindex*> I2(the_tensors[1].get_indices());
    // Define group of commonly shared, dummy indices which doesn't include external ones, as summedInd
    for(vector<SQindex*>::iterator i1 = I1.begin();i1 != I1.end();++i1){
      if((find(I2.begin(), I2.end(), *i1) != I2.end()/* || find(I3.begin(), I3.end(), *i1) != I3.end()*/) && 
	  find(extInd.begin(), extInd.end(), *i1) == extInd.end() && !(*i1)->get_isExt() && !(*i1)->get_isSummed())
        extInd.push_back(*i1);
      else if(find(I2.begin(), I2.end(), *i1) != I2.end() && 
              find(summedInd.begin(), summedInd.end(), *i1) == summedInd.end() && !(*i1)->get_isExt() && (*i1)->get_isSummed())
        summedInd.push_back(*i1);
      else if(find(rowInd.begin(), rowInd.end(), *i1) == rowInd.end() && !(*i1)->get_isExt())
        rowInd.push_back(*i1);
    } // End i1
    for(vector<SQindex*>::iterator i2 = I2.begin();i2 != I2.end();++i2){
      if(find(summedInd.begin(), summedInd.end(), *i2) == summedInd.end() && 
         find(colInd.begin(), colInd.end(),       *i2) == colInd.end()    && !(*i2)->get_isExt() &&
         find(extInd.begin(), extInd.end(),       *i2) == extInd.end()    &&
         find(rowInd.begin(), rowInd.end(),       *i2) == rowInd.end())
        colInd.push_back(*i2);
    } // End i2

#ifndef _DEBUG_BLAS
    cout << " - extInd :: " << endl << "    ";
    for(vector<SQindex*>::iterator i = extInd.begin();i != extInd.end();++i)
      cout << **i << ", ";
    cout << endl;
    cout << " - rowInd :: " << endl << "    ";
    for(vector<SQindex*>::iterator i = rowInd.begin();i != rowInd.end();++i)
      cout << **i << ", ";
    cout << endl;
    cout << " - colInd :: " << endl << "    ";
    for(vector<SQindex*>::iterator i = colInd.begin();i != colInd.end();++i)
      cout << **i << ", ";
    cout << endl;
    cout << " - summedInd :: " << endl << "    ";
    for(vector<SQindex*>::iterator i = summedInd.begin();i != summedInd.end();++i)
      cout << **i << ", ";
    cout << endl << endl;
#endif

    // If L_nd_ind has indices of extInd, remove these
    for(vector<SQindex*>::iterator i = L_nd_ind.begin();i != L_nd_ind.end();)
      if(find(extInd.begin(), extInd.end(), *i) != extInd.end()) i = L_nd_ind.erase(i);
      else ++i;

    pair<bool, bool> bareRDM; // <the_tensor[0], the_tensor[1]> 
    bareRDM.first  = false;
    bareRDM.second = false;
#ifdef _DGEMM_OPT
    /////////////////////////////////////////////////////////////////////////////////////////////
    { // Optimize the ordering of the indices and rotate the indices of the RDM if possible
      int precedence(-1);
      if(is_RDM(the_tensors[0].get_name()) && is_RDM(the_tensors[1].get_name())){
        if(the_tensors[0].get_indices().size() < the_tensors[1].get_indices().size())
	     precedence = 1;
	else precedence = 0;
      } // End if

#ifdef _OPT_RANK
      const int opt_rank(2*_OPT_RANK);
#else
      const int opt_rank(0);
#endif
      pair<bool, bool> is4RDM;
      is4RDM.first  = false;
      is4RDM.second = false;
      if(is_RDM(the_tensors[0].get_name()) && the_tensors[0].get_indices().size() == 8 || 
	 is_RDM(the_tensors[0].get_name()) && the_tensors[0].get_indices().size()  < opt_rank) is4RDM.first  = true;
      if(is_RDM(the_tensors[1].get_name()) && the_tensors[1].get_indices().size() == 8 || 
	 is_RDM(the_tensors[1].get_name()) && the_tensors[1].get_indices().size()  < opt_rank) is4RDM.second = true;
      
      pair<bool,bool> allRDM;
      allRDM.first  = false;
      allRDM.second = false;
      if(is_RDM(the_tensors[0].get_name()) && rowInd.size()+summedInd.size() == the_tensors[0].get_indices().size()) allRDM.first  = true;
      if(is_RDM(the_tensors[1].get_name()) && colInd.size()+summedInd.size() == the_tensors[1].get_indices().size()) allRDM.second = true;

      long unsigned int maxScore(0); // Mas score for the best matching

      if(((is_RDM(the_tensors[0].get_name()) && precedence == -1) || precedence == 0) && !is4RDM.first && allRDM.first){
	// Compare,
	//        [[ rowInd ], [ summedInd ]] and 
	//        [ Reversed indices of t0  ]
	// to determine the best matching
	IIvector rowPerm(makePermutations(rowInd.size()));
	IIvector summedPerm(makePermutations(summedInd.size()));
	IIvector tPerm(the_tensors[0].get_perms());
	size_t rI = (rowPerm.size()    ? rowPerm.size()    : 1);
	size_t sI = (summedPerm.size() ? summedPerm.size() : 1);
	size_t tI = (tPerm.size()      ? tPerm.size()      : 1);
	
	typedef tuple<Ivector, Ivector, size_t, long unsigned int> combis; // <rowwPerm, summedPerm, tPerm, (ULI)score> >
	vector<combis> combination; 
	
	vector<SQindex*> tinds(the_tensors[0].get_indices());
	for(size_t num = 0;num < tinds.size();++num) maxScore += pow(10,num);

	for(size_t ir = 0;ir < rI;++ir){
	  for(size_t is = 0;is < sI;++is){
	    for(size_t it = 0;it < tI;++it){

	      long unsigned int score(0);
	      for(size_t num = 0; num < tinds.size();++num){
		size_t num_ti(tPerm[it][num]);                  // Usual order
		size_t num_ri(tPerm[it][tinds.size()-(num+1)]); // Reversed order
                if     (num < rowInd.size() && rowPerm.size()){ // if the index belongs to rowInd
                  if(*(rowInd[rowPerm[ir][num]]) == *(tinds[num_ri])) score += pow(10, num_ti);
		} // End if
		else if(summedPerm.size()){ // if the index belongs to the colInd
		  if(*(summedInd[summedPerm[is][num-rowInd.size()]]) == *(tinds[num_ri])) score += pow(10, num_ti);
		} // End else
	      } // End num_ti

#ifdef _DEBUG_OPT
	      cout << " -- title : " << title << endl;
	      cout << " -- Score " << score << " " << sI << " " << rI << " " << tI << endl;
	      cout << " -- rowInd :: " << endl << "    ";
	      for(size_t num_i = 0;num_i < rowInd.size();++num_i)
		cout << *rowInd[rowPerm[ir][num_i]] << ", ";
	      cout << endl;
	      cout << " -- summedInd :: " << endl << "    ";
	      for(size_t num_i = 0;num_i < summedInd.size();++num_i)
		cout << *summedInd[summedPerm[is][num_i]] << ", ";
	      cout << endl;
	      cout << " -- tenInd :: " << endl << "    ";
	      for(size_t num_ti = 0; num_ti < tinds.size();++num_ti)
		cout << *tinds[tinds.size()-(tPerm[it][num_ti]+1)] << ", ";
	      cout << endl << endl;
#endif
	      Ivector dummy;
	      if(!rowPerm.size())    rowPerm.push_back(dummy);
	      if(!summedPerm.size()) summedPerm.push_back(dummy);
	      combination.push_back(make_tuple(rowPerm[ir], summedPerm[is], it, score));
	      
	    } // End it
	  } // End is
	} // End ir
	
	cout << bin << ", " << "madekita" << endl; //abort();
	
	// Find the one with the largest score
	pair<vector<combis>::iterator, long unsigned int> the_best; // <the_best_one, the_largest_score>
	the_best.second = 0;
	for(auto c = combination.begin();c != combination.end();++c){
	  if(get<3>(*c) > the_best.second) { the_best.first = c; the_best.second = get<3>(*c); }
	} // End c
	
	if(the_best.second){
	  cout << " >> The best combinations :: " << the_best.second << endl;
	  f << "! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
	  f << "! Information for debugging " << endl;
	  f << "! >> Score :: " << the_best.second << endl;
	  if(the_best.second == maxScore) 
	    { f << "! The best shot!" << endl; bareRDM.first = true; /* Turn this on! */ }
	  
	  cout << "madekita - 2" << endl;
	  
	  // Re-set all the indices
	  auto rbest(get<0>(*the_best.first));
	  auto sbest(get<1>(*the_best.first));
	  auto tbest(get<2>(*the_best.first));
	  
	  cout << "madekita - 2'" << endl;
	  
	  f << "! RDM is rotated :: " << the_tensors[0] << " >> ";
	  the_tensors[0].rotateIndices(tbest);
	  f << the_tensors[0] << endl;
	  
	  cout << "madekita - 3" << endl;
	  
	  vector<SQindex*> rtemp;
	  vector<SQindex*> stemp;
	  for(size_t num_i = 0;num_i < rowInd.size();   ++num_i) rtemp.push_back(rowInd[rbest[num_i]]   );

	  f << "! rowInd : ";
	  for(size_t num_i = 0;num_i < rowInd.size();   ++num_i) f << *rtemp[num_i] << " ";
	  f<< endl;
	  
	  cout << "madekita - 4" << endl;
	  
	  f << "! summedInd : ";
	  for(size_t num_i = 0;num_i < summedInd.size();++num_i) stemp.push_back(summedInd[sbest[num_i]]);
	  for(size_t num_i = 0;num_i < summedInd.size();++num_i) f << *stemp[num_i] << " ";
	  f << endl;
	  
	  cout << "madekita - 5" << endl;
	  
	  rowInd    = rtemp;
	  summedInd = stemp;
	} // End if
	
      } // End if
      else { 
	f << "! -- Check1 is skipped " << endl;
	f << "!    -- is4RDM.first " << (is4RDM.first ? "true" : "false") << endl; 
	f << "!    -- allRDM.first " << (allRDM.first ? "true" : "false") << endl; 
	f << "!    -- precedence   " <<  precedence                        << endl; 
      } // End else

      if(((is_RDM(the_tensors[1].get_name()) && precedence == -1) || precedence == 1) && !is4RDM.second && allRDM.second){
	//f << "! NOT YET IMPLEMENTED" << endl;
	f << "! IS IT OK??" << endl;
	// Compare,
	//        [[ summedInd ], [ colInd ]] and 
	//        [ Reversed indices of t1  ]
	// to determine the best matching
	IIvector summedPerm(makePermutations(summedInd.size()));
	IIvector colPerm(makePermutations(colInd.size()));
	IIvector tPerm(the_tensors[1].get_perms());
	size_t sI = (summedPerm.size() ? summedPerm.size() : 1);
	size_t cI = (colPerm.size()    ? colPerm.size()    : 1);
	size_t tI = (tPerm.size()      ? tPerm.size()      : 1);
	
	typedef tuple<Ivector, Ivector, size_t, long unsigned int> combis; // <summedPerm, colPerm, tPerm, (ULI)score> >
	vector<combis> combination; 
	
	vector<SQindex*> tinds(the_tensors[1].get_indices());
	for(size_t num = 0;num < tinds.size();++num) maxScore += pow(10,num);

	for(size_t is = 0;is < sI;++is){
	  for(size_t ic = 0;ic < cI;++ic){
	    for(size_t it = 0;it < tI;++it){

	      long unsigned int score(0);
	      for(size_t num = 0; num < tinds.size();++num){
		size_t num_ti(tPerm[it][num]);                  // Usual order
		size_t num_ri(tPerm[it][tinds.size()-(num+1)]); // Reversed order
                if     (num < summedInd.size() && summedPerm.size()){ // if the index belongs to summedInd
                  if(*(summedInd[summedPerm[is][num]]) == *(tinds[num_ri])) score += pow(10, num_ti);
		} // End if
		else if(colPerm.size()){ // if the index belongs to the colInd
		  if(*(colInd[colPerm[ic][num-summedInd.size()]]) == *(tinds[num_ri])) score += pow(10, num_ti);
		} // End else
	      } // End num_ti

#ifdef _DEBUG_OPT
	      cout << " -- title : " << title << endl;
	      cout << " -- Score " << score << " " << sI << " " << cI << " " << tI << endl;
	      cout << " -- summedInd :: " << endl << "    ";
	      for(size_t num_i = 0;num_i < summedInd.size();++num_i)
		cout << *summedInd[summedPerm[is][num_i]] << ", ";
	      cout << endl;
	      cout << " -- colInd :: " << endl << "    ";
	      for(size_t num_i = 0;num_i < colInd.size();++num_i)
		cout << *colInd[colPerm[ic][num_i]] << ", ";
	      cout << endl;
	      cout << " -- tenInd :: " << endl << "    ";
	      for(size_t num_ti = 0; num_ti < tinds.size();++num_ti)
		cout << *tinds[tinds.size()-(tPerm[it][num_ti]+1)] << ", ";
	      cout << endl << endl;
#endif
	      Ivector dummy;
	      if(!summedPerm.size()) summedPerm.push_back(dummy);
	      if(!colPerm.size())    colPerm.push_back(dummy);
	      combination.push_back(make_tuple(summedPerm[is], colPerm[ic], it, score));
	      
	    } // End it
	  } // End ic
	} // End is
	
	cout << bin << ", " << "madekita" << endl; //abort();
	
	// Find the one with the largest score
	pair<vector<combis>::iterator, long unsigned int> the_best; // <the_best_one, the_largest_score>
	the_best.second = 0;
	for(auto c = combination.begin();c != combination.end();++c){
	  if(get<3>(*c) > the_best.second) { the_best.first = c; the_best.second = get<3>(*c); }
	} // End c
	
	if(the_best.second){
	  cout << " >> The best combinations :: " << the_best.second << endl;
	  f << "! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
	  f << "! Information for debugging " << endl;
	  f << "! >> Score :: " << the_best.second << endl;
	  if(the_best.second == maxScore) 
	    { f << "! The best shot!" << endl; bareRDM.second = true; /* Turn this on! */ }
	  
	  cout << "madekita - 2" << endl;
	  
	  // Re-set all the indices
	  auto sbest(get<0>(*the_best.first));
	  auto cbest(get<1>(*the_best.first));
	  auto tbest(get<2>(*the_best.first));
	  
	  cout << "madekita - 2'" << endl;
	  
	  f << "! RDM is rotated :: " << the_tensors[1] << " >> ";
	  the_tensors[1].rotateIndices(tbest);
	  f << the_tensors[1] << endl;
	  
	  cout << "madekita - 3" << endl;
	  
	  vector<SQindex*> stemp;
	  vector<SQindex*> ctemp;
	  
	  f << "! summedInd : ";
	  for(size_t num_i = 0;num_i < summedInd.size();++num_i) stemp.push_back(summedInd[sbest[num_i]]);
	  for(size_t num_i = 0;num_i < summedInd.size();++num_i) f << *stemp[num_i] << " ";
	  f << endl;

	  cout << "madekita - 4" << endl;

	  f << "! colInd : ";
	  for(size_t num_i = 0;num_i < colInd.size();   ++num_i) ctemp.push_back(colInd[cbest[num_i]]   );
	  for(size_t num_i = 0;num_i < colInd.size();   ++num_i) f << *ctemp[num_i] << " ";
	  f<< endl;
	  
	  cout << "madekita - 5" << endl;
	  
	  summedInd = stemp;
	  colInd    = ctemp;
	} // End if

      } // End if
      else { 
	f << "! -- Check2 is skipped " << endl;
	f << "!    -- is4RDM.second " << (is4RDM.second ? "true" : "false") << endl; 
	f << "!    -- allRDM.second " << (allRDM.second ? "true" : "false") << endl; 
	f << "!    -- precedence    " <<  precedence                        << endl; 
      } // End else
      
    } // End scope

    vector<SQindex*> Z1inds(rowInd);    // Indices of Z1 tensor
    Z1inds.insert(Z1inds.end(), summedInd.begin(), summedInd.end());
    vector<SQindex*> Z2inds(summedInd); // Indices of Z2 tensor
    Z2inds.insert(Z2inds.end(), colInd.begin(), colInd.end());
    vector<SQindex*> Z3inds(rowInd);    // Indices of Z3 tensor
    Z3inds.insert(Z3inds.end(), colInd.begin(), colInd.end());

#ifndef _DEBUG_OPT
    cout << " - rowInd :: " << endl << "    ";
    for(vector<SQindex*>::iterator i = rowInd.begin();i != rowInd.end();++i)
      cout << **i << ", ";
    cout << endl;
    cout << " - colInd :: " << endl << "    ";
    for(vector<SQindex*>::iterator i = colInd.begin();i != colInd.end();++i)
      cout << **i << ", ";
    cout << endl;
    cout << " - summedInd :: " << endl << "    ";
    for(vector<SQindex*>::iterator i = summedInd.begin();i != summedInd.end();++i)
      cout << **i << ", ";
    cout << endl << endl;
#endif

    /////////////////////////////////////////////////////////////////////////////////////////////
#endif

    // Variables that define which intermediate is allocated 
    vector<bool> is_allocated(3);

    int count = 1;
    string DecInterm("! Intermediate arrays                \n");
    for(vector<SQtensor>::iterator t = the_tensors.begin();t != the_tensors.end();++t, ++count){

      if(count == 1 && bareRDM.first ) continue; 
      if(count == 2 && bareRDM.second) continue;

      ostringstream stm;
      stm << count;
      size_t id_count = (count == 1 ? rowInd.size() : colInd.size()) + summedInd.size();
      DecInterm += "real*8";
      if(id_count) { DecInterm += ", allocatable"; is_allocated[count-1] = true; }
      else is_allocated[count-1] = false;
      DecInterm += " :: Z" + stm.str() + "_" + (id_count ? "(" : "\n");
      //*TEST* DecInterm += "real*8, allocatable :: Z" + stm.str() + "_(";
      for(size_t i = 0;i < id_count;++i)
        if(i == id_count-1) DecInterm += ":)\n";
        else                DecInterm += ":,";
    } // End t
    if(rowInd.size() || colInd.size()){
      size_t id_count = rowInd.size() + colInd.size();
      DecInterm += "real*8, allocatable :: Z3_(";
      for(size_t i = 0;i < id_count;++i)
        if(i == id_count-1) DecInterm += ":)\n";
        else                DecInterm += ":,";
      is_allocated[2] = true;
    } // End if
    else{
      DecInterm += "real*8 :: Z3_\n";
      is_allocated[2] = false;
    } // End if

    //f << printData;
    //f << printName;
    f << printName2;
    f << DecComm8;
    f << DecExtInd;
    f << DecSymm;
    f << DecConsts;
    f << DecTensor;
    f << DecInterm;

    // Entering the body of the tensorial contraction part
    int LoopCount = 0;
    vector<SQindex*> NDindices;
    for(vector<SQindex*>::iterator i = summedBody.begin();i != summedBody.end();++i){
      if(!(*i)->get_isExt()) NDindices.push_back(*i);
    }// End i
    
    f << "! Indices used in the contractions as dummy ... \n";
    
    for(size_t num_i = 0;num_i < NDindices.size();++num_i){
      SQindex* i = NDindices[num_i];
      if(num_i % 5 != 0) f << ", s_" + i->get_index() + ", i_" + i->get_index();
      else{
        f << endl;
        f << "integer :: ";
        f << "s_" + i->get_index() + ", i_" + i->get_index();
      } // End else
    } // End num_i
    f << endl;
    f << "! " << bin << endl;

    for(vector<SQindex*>::iterator i = NDindices.begin();i != NDindices.end();++i){
      f << "do s_" + (*i)->get_index() + " = 0, nir-1" + "\n";
      ++LoopCount;
    } // End i

    // Symmetry constraints
    f << "if( &" << endl;
    if     (bin.get_Ltensor().get_indices().size() == 1){ // In case of unity
      vector<SQindex*> inds(bin.get_Ltensor().get_indices());
      f << "s_" + inds[0]->get_index() + " == 0";
    } // End if
    else if(bin.get_Ltensor().get_indices().size() == 2){ // In case of two
      vector<SQindex*> inds(bin.get_Ltensor().get_indices());
      f << "IEOR(s_" + inds[0]->get_index() + ",s_" + inds[1]->get_index() + ") == 0";
    } // End if
    else if(bin.get_Ltensor().get_indices().size() == 3){ // In case of three
      vector<SQindex*> inds(bin.get_Ltensor().get_indices());
      string cstr;
      cstr  = "IEOR(IEOR(s_" + inds[0]->get_index() + ",s_" + inds[1]->get_index() + "),";
      cstr += "s_" + inds[2]->get_index() + ") == 0";
      f << cstr;
    } // End if
    else if(bin.get_Ltensor().get_indices().size() == 4){ // In case of four
      vector<SQindex*> inds(bin.get_Ltensor().get_indices());
      string cstr;
      cstr  = "IEOR(s_" + inds[0]->get_index() + ",s_" + inds[1]->get_index() + ") == ";
      cstr += "IEOR(s_" + inds[2]->get_index() + ",s_" + inds[3]->get_index() + ")";
      f << cstr;
    } // End if
    else if(bin.get_Ltensor().get_indices().size() == 5){ // In case of five
      vector<SQindex*> inds(bin.get_Ltensor().get_indices());
      string cstr;
      cstr  = "IEOR(IEOR(s_" + inds[0]->get_index() + ",s_" + inds[1]->get_index() + "),s_" + inds[2]->get_index() + ") == ";
      cstr += "IEOR(s_" + inds[3]->get_index() + ",s_" + inds[4]->get_index() + ")";
      f << cstr;
    } // End if
    else if(bin.get_Ltensor().get_indices().size() == 6){ // In case of six
      vector<SQindex*> inds(bin.get_Ltensor().get_indices());
      string cstr;
      cstr  = "IEOR(IEOR(s_" + inds[0]->get_index() + ",s_" + inds[1]->get_index() + "),s_" + inds[2]->get_index() + ") == ";
      cstr += "IEOR(IEOR(s_" + inds[3]->get_index() + ",s_" + inds[4]->get_index() + "),s_" + inds[5]->get_index() + ")";
      f << cstr;
    } // End if
    else if(bin.get_Ltensor().get_indices().size() == 7){ // In case of seven
      vector<SQindex*> inds(bin.get_Ltensor().get_indices());
      string cstr;
      cstr  = "IEOR(IEOR(s_" + inds[0]->get_index() + ",s_" + inds[1]->get_index() + "),IEOR(s_" + inds[2]->get_index() + ",s_" + inds[3]->get_index() + ")) == ";
      cstr += "IEOR(IEOR(s_" + inds[4]->get_index() + ",s_" + inds[5]->get_index() + "),s_" + inds[6]->get_index()+ ")";
      f << cstr;
    } // End if
    else if(bin.get_Ltensor().get_indices().size() == 8){ // In case of eight
      vector<SQindex*> inds(bin.get_Ltensor().get_indices());
      string cstr;
      cstr  = "IEOR(IEOR(s_" + inds[0]->get_index() + ",s_" + inds[1]->get_index() + "),IEOR(s_" + inds[2]->get_index() + ",s_" + inds[3]->get_index() + ")) == ";
      cstr += "IEOR(IEOR(s_" + inds[4]->get_index() + ",s_" + inds[5]->get_index() + "),IEOR(s_" + inds[6]->get_index() + ",s_" + inds[7]->get_index() + "))";
      f << cstr;
    } // End if
    else if(bin.get_Ltensor().get_indices().size() != 0){
      cout << "Factorize: Number of indices in tensors must be 0 - 8" << endl;
      abort();
    } // End if
    if (bin.get_Ltensor().get_indices().size() && tensors.size())
      f << " .and. & " << endl;

    vector<SQtensor>::iterator t_end = tensors.end(); --t_end;
    for(vector<SQtensor>::iterator t = tensors.begin();t != tensors.end();++t){
      if     (t->get_indices().size() == 0){ // In case of scalar
        vector<SQindex*> inds(t->get_indices());
        f << ".True.";
      } // End if
      else if(t->get_indices().size() == 1){ // In case of unity
        vector<SQindex*> inds(t->get_indices());
        f << "s_" + inds[0]->get_index() + " == 0";
      } // End if
      else if(t->get_indices().size() == 2){ // In case of two
        vector<SQindex*> inds(t->get_indices());
        f << "IEOR(s_" + inds[0]->get_index() + ",s_" + inds[1]->get_index() + ") == 0";
      } // End if
      else if(t->get_indices().size() == 3){ // In case of three
        vector<SQindex*> inds(t->get_indices());
        string cstr;
        cstr  = "IEOR(IEOR(s_" + inds[0]->get_index() + ",s_" + inds[1]->get_index() + "),";
        cstr += "s_" + inds[2]->get_index() + ") == 0";
        f << cstr;
      } // End if
      else if(t->get_indices().size() == 4){ // In case of four
        vector<SQindex*> inds(t->get_indices());
        string cstr;
        cstr  = "IEOR(s_" + inds[0]->get_index() + ",s_" + inds[1]->get_index() + ") == ";
        cstr += "IEOR(s_" + inds[2]->get_index() + ",s_" + inds[3]->get_index() + ")";
        f << cstr;
      } // End if
      else if(t->get_indices().size() == 5){ // In case of five
        vector<SQindex*> inds(t->get_indices());
        string cstr;
        cstr  = "IEOR(IEOR(s_" + inds[0]->get_index() + ",s_" + inds[1]->get_index() + "),s_" + inds[2]->get_index() + ") == ";
        cstr += "IEOR(s_" + inds[3]->get_index() + ",s_" + inds[4]->get_index() + ")";
        f << cstr;
      } // End if
      else if(t->get_indices().size() == 6){ // In case of six
        vector<SQindex*> inds(t->get_indices());
        string cstr;
        cstr  = "IEOR(IEOR(s_" + inds[0]->get_index() + ",s_" + inds[1]->get_index() + "),s_" + inds[2]->get_index() + ") == ";
        cstr += "IEOR(IEOR(s_" + inds[3]->get_index() + ",s_" + inds[4]->get_index() + "),s_" + inds[5]->get_index() + ")";
        f << cstr;
      } // End if
      else if(t->get_indices().size() == 7){ // In case of seven
        vector<SQindex*> inds(t->get_indices());
        string cstr;
        cstr  = "IEOR(IEOR(s_" + inds[0]->get_index() + ",s_" + inds[1]->get_index() + "),IEOR(s_" + inds[2]->get_index() + ",s_" + inds[3]->get_index() + ")) == ";
        cstr += "IEOR(IEOR(s_" + inds[4]->get_index() + ",s_" + inds[5]->get_index() + "),s_" + inds[6]->get_index()+ ")";
        f << cstr;
      } // End if
      else if(t->get_indices().size() == 8){ // In case of eight
        vector<SQindex*> inds(t->get_indices());
        string cstr;
        cstr  = "IEOR(IEOR(s_" + inds[0]->get_index() + ",s_" + inds[1]->get_index() + "),IEOR(s_" + inds[2]->get_index() + ",s_" + inds[3]->get_index() + ")) == ";
        cstr += "IEOR(IEOR(s_" + inds[4]->get_index() + ",s_" + inds[5]->get_index() + "),IEOR(s_" + inds[6]->get_index() + ",s_" + inds[7]->get_index() + "))";
        f << cstr;
      } // End if
      else if(t->get_indices().size() != 0){
        cout << "Factorize: Number of indices in tensors must be 0 - 8" << endl;
        abort();
      } // End if
      if(t != t_end)
        f << " .and. &" << endl;
    } // End t
    f << ") then" << endl;

    // Dim constraints (these are necessary for MKL)
    if(rowInd.size() || colInd.size() || summedInd.size()) f << endl;
    bool dimCount = false;
    string m_label("");
    for(vector<SQindex*>::iterator i = rowInd.begin();i != rowInd.end();++i){
      string sname("s_" + (*i)->get_index());
      string otype;
      if     ((*i)->get_char() == (char_state)0) otype = "I_C";
      else if((*i)->get_char() == (char_state)1) otype = "I_O";
      else if((*i)->get_char() == (char_state)2) otype = "I_V";
      
      m_label += "psym(I_LENGTH," + otype + ", " + sname + ")";

      if(*i != rowInd.back()) m_label += "*"; 
    } // End i
    if(rowInd.size() || colInd.size() || summedInd.size()){
      f << "if(";
      if(rowInd.size()){
        f << m_label + " > 0";
        if(colInd.size() || summedInd.size()) f << " .and. &\n";
        else                                  f << ") then" << endl;
      }// End if
      dimCount = true;
    } // End if

    string n_label("");
    for(vector<SQindex*>::iterator i = colInd.begin();i != colInd.end();++i){
      string sname("s_" + (*i)->get_index());
      string otype;
      if     ((*i)->get_char() == (char_state)0) otype = "I_C";
      else if((*i)->get_char() == (char_state)1) otype = "I_O";
      else if((*i)->get_char() == (char_state)2) otype = "I_V";
      
      n_label += "psym(I_LENGTH," + otype + ", " + sname + ")";

      if(*i != colInd.back()) n_label += "*"; 
    } // End i
    if(colInd.size() || summedInd.size()){
      if(colInd.size()){
        f << (rowInd.size() ? "   " : "" ) + n_label + " > 0";
        if(summedInd.size()) f << " .and. &\n";
        else                 f << ") then" << endl;
      } // End if
    } // End if

    string k_label("");
    for(vector<SQindex*>::iterator i = summedInd.begin();i != summedInd.end();++i){
      string sname("s_" + (*i)->get_index());
      string otype;
      if     ((*i)->get_char() == (char_state)0) otype = "I_C";
      else if((*i)->get_char() == (char_state)1) otype = "I_O";
      else if((*i)->get_char() == (char_state)2) otype = "I_V";
      
      k_label += "psym(I_LENGTH," + otype + ", " + sname + ")";

      if(*i != summedInd.back()) k_label += "*";
    } // End i
    if(summedInd.size()){
      f << (rowInd.size() || colInd.size() ? "   " : "") + k_label + " > 0) then" << endl;
    } // End if
    if(rowInd.size() || colInd.size() || summedInd.size()) f << endl;

    // Declare loops over external indices
    if(extInd.size()) f << endl;
    for(vector<SQindex*>::iterator i = extInd.begin();i != extInd.end();++i){
      string iname = "i_" + (*i)->get_index();
      string sname = "s_" + (*i)->get_index();
      string otype;
      if     ((*i)->get_char() == (char_state)0) otype = "I_C";
      else if((*i)->get_char() == (char_state)1) otype = "I_O";
      else if((*i)->get_char() == (char_state)2) otype = "I_V";
      f << "do " + iname + " = psym(I_BEGIN, " + otype + ", " + sname + "), psym(I_END, " + otype + ", " + sname+ ")" << endl;
    } // End i
    if(extInd.size()) f << endl;

    // Allocate the first tensor (if necessary)
    if(!bareRDM.first){
      f << "! Z1 <-- " << the_tensors[0] << endl;
      if(is_allocated[0]) f << "allocate(Z1_(";
      for(vector<SQindex*>::iterator i = Z1inds.begin();i != Z1inds.end();++i){
	string sname("s_" + (*i)->get_index());
	string otype;
	if     ((*i)->get_char() == (char_state)0) otype = "I_C";
	else if((*i)->get_char() == (char_state)1) otype = "I_O";
	else if((*i)->get_char() == (char_state)2) otype = "I_V";
	
	if(i == Z1inds.begin())
	  f << "psym(I_BEGIN," + otype + ", " + sname + "):psym(I_END," + otype + ", " + sname + ")";
	else
	  f << "             psym(I_BEGIN," + otype + ", " + sname + "):psym(I_END," + otype + ", " + sname + ")";
	
	if(*i == Z1inds.back()) f << "))\n";
	else                   f << ", &\n";
      } // End i
      
      // Sort the first tensor
      int count1 = 0;
      if(use_omp_)
	f << endl << "!$omp parallel do" << endl;
      //for(vector<SQindex*>::iterator i = Z1inds.begin();i != Z1inds.end();++i, ++count1){
      for(vector<SQindex*>::reverse_iterator i = Z1inds.rbegin();i != Z1inds.rend();++i, ++count1){ // Modified: 2013/05/22
	string iname = "i_" + (*i)->get_index();
	string sname = "s_" + (*i)->get_index();
	string otype;
	if     ((*i)->get_char() == (char_state)0) otype = "I_C";
	else if((*i)->get_char() == (char_state)1) otype = "I_O";
	else if((*i)->get_char() == (char_state)2) otype = "I_V";
	f << "do " + iname + " = psym(I_BEGIN, " + otype + ", " + sname + "), psym(I_END, " + otype + ", " + sname+ ")" << endl;
      } // End i
      
      { // Begin scope 1
	string nterm("Z1_");
	if(is_allocated[0]) nterm += "(";
	for(size_t num_i = 0;num_i < Z1inds.size();++num_i){
	  nterm += "i_" + Z1inds[num_i]->get_index();
	  if(num_i != Z1inds.size()-1) nterm += ", ";
	  else nterm += ")";
	} // End num_i
	f << nterm + " =  &\n";

	vector<SQindex*> Rindlist;
	if(the_tensors[0].get_name() == name_h2_){
	  for(size_t num_i = 0;num_i < the_tensors[0].get_indices().size();++num_i){
	    SQindex* i(the_tensors[0].get_indices()[num_i]);
	    if(num_i != exth2_) Rindlist.push_back(i);
	  } // End num_i
	} // End if
	else if(the_tensors[0].get_name() == name_amp_){
	  for(size_t num_i = 0;num_i < the_tensors[0].get_indices().size();++num_i){
	    SQindex* i(the_tensors[0].get_indices()[num_i]);
	    if(num_i != extamp_) Rindlist.push_back(i);
	  } // End num_i
	} // End if
	else if(the_tensors[0].get_name() == RDM_name() + "4"){
	  int external_count = 0;
	  for(size_t num_i = 0;num_i < the_tensors[0].get_indices().size();++num_i){
	    SQindex* i(the_tensors[0].get_indices()[num_i]);
	    if(num_i != 0 && num_i != 1) Rindlist.push_back(i);
	    else ++external_count;
	  } // End num_i
	  if(external_count != 2) {
	    cout << "Factorize: Something is wrong in treatment of 4-RDM .... " << endl;
	    cout << "external_count: " << external_count << endl; 
	    abort(); 
	  } // End if
	} // End if
	//else if(the_tensors[0].get_name() == "X" && !X_indices.size()){
	else if(is_Interm(the_tensors[0].get_name()) && !X_indicesR2[0].size()){
	  for(size_t num_i = 0;num_i < the_tensors[0].get_indices().size();++num_i){
	    SQindex* i(the_tensors[0].get_indices()[num_i]);
	    if(!i->get_isExt()) Rindlist.push_back(i);
	  } // End num_i
	} // End if
	//else if(the_tensors[0].get_name() == "X" &&  X_indices.size()){
	else if(is_Interm(the_tensors[0].get_name()) &&  X_indicesR2[0].size()){
	  Rindlist.insert(Rindlist.end(), X_indicesR2[0].begin(), X_indicesR2[0].end());
	} // End if
	else if(is_D4C(the_tensors[0].get_name())){
	  for(size_t num_i = 0;num_i < the_tensors[0].get_indices().size();++num_i){
	    SQindex* i(the_tensors[0].get_indices()[num_i]);
	    if(num_i != extd4c_) Rindlist.push_back(i);
	  } // End num_i
	} // End else
	else{
	  Rindlist = the_tensors[0].get_indices();
	} // End else
	
	string nterm2("  " + the_tensors[0].get_name() + "_");
	nterm2 += "(";
	for(size_t num_i = 0;num_i < Rindlist.size();++num_i){
	  nterm2 += "s_" + Rindlist[Rindlist.size()-(num_i+1)]->get_index();
	  if(num_i != Rindlist.size()-1) nterm2 += ", ";
	  else nterm2 += ")";
	} // End num_i
	
	nterm2 += "%array("; 
	for(size_t num_i = 0;num_i < Rindlist.size();++num_i){
	  nterm2 += "i_" + Rindlist[Rindlist.size()-(num_i+1)]->get_index();
	  if(num_i != Rindlist.size()-1) nterm2 += ", ";
	  else nterm2 += ")";
	} // End num_i
	f << nterm2 << endl;
      } // End scope 1    
      for(int n = 0;n < count1;++n) f << "end do"   << endl;
    } // End bareRDM
    else is_allocated[0] = false;

    if(!bareRDM.second){    
      // Allocate the second tensor
      f << "! Z2 <-- " << the_tensors[1] << endl;
      if(is_allocated[1]) f << "allocate(Z2_(";
      for(vector<SQindex*>::iterator i = Z2inds.begin();i != Z2inds.end();++i){
	string sname("s_" + (*i)->get_index());
	string otype;
	if     ((*i)->get_char() == Femto::core) otype = "I_C";
	else if((*i)->get_char() == Femto::act ) otype = "I_O";
	else if((*i)->get_char() == Femto::virt) otype = "I_V";
	
	if(i == Z2inds.begin())      
	  f << "psym(I_BEGIN," + otype + ", " + sname + "):psym(I_END," + otype + ", " + sname + ")";
	else
	  f << "             psym(I_BEGIN," + otype + ", " + sname + "):psym(I_END," + otype + ", " + sname + ")";
	
	if(*i == Z2inds.back()) f << "))\n";
	else                    f << ", &\n";
      } // End i
      
      // Sort the second tensor
      int count2 = 0;
      if(use_omp_)
	f << endl << "!$omp parallel do" << endl;
      //for(vector<SQindex*>::iterator i = Z2inds.begin();i != Z2inds.end();++i, ++count2){
      for(vector<SQindex*>::reverse_iterator i = Z2inds.rbegin();i != Z2inds.rend();++i, ++count2){ // Modified: 2013/05/22
	string iname = "i_" + (*i)->get_index();
	string sname = "s_" + (*i)->get_index();
	string otype;
	if     ((*i)->get_char() == Femto::core) otype = "I_C";
	else if((*i)->get_char() == Femto::act ) otype = "I_O";
	else if((*i)->get_char() == Femto::virt) otype = "I_V";
	f << "do " + iname + " = psym(I_BEGIN, " + otype + ", " + sname + "), psym(I_END, " + otype + ", " + sname+ ")" << endl;
      } // End i
      
      { // Begin scope 2
	string nterm("Z2_");
	if(is_allocated[1]) nterm += "(";
	for(size_t num_i = 0;num_i < Z2inds.size();++num_i){
	  nterm += "i_" + Z2inds[num_i]->get_index();
	  if(num_i != Z2inds.size()-1) nterm += ", ";
	  else nterm += ")";
	} // End num_i
	f << nterm + " =  &\n";
	
	vector<SQindex*> Rindlist;
	if(the_tensors[1].get_name() == name_h2_){
	  for(size_t num_i = 0;num_i < the_tensors[1].get_indices().size();++num_i){
	    SQindex* i(the_tensors[1].get_indices()[num_i]);
	    if(num_i != exth2_) Rindlist.push_back(i);
	  } // End num_i
	} // End if
	else if(the_tensors[1].get_name() == name_amp_){
	  for(size_t num_i = 0;num_i < the_tensors[1].get_indices().size();++num_i){
	    SQindex* i(the_tensors[1].get_indices()[num_i]);
	    if(num_i != extamp_) Rindlist.push_back(i);
	  } // End num_i
	} // End if
	else if(the_tensors[1].get_name() == RDM_name() + "4"){
	  int external_count = 0;
	  for(size_t num_i = 0;num_i < the_tensors[1].get_indices().size();++num_i){
	    SQindex* i(the_tensors[1].get_indices()[num_i]);
	    if(num_i != 0 && num_i != 1) Rindlist.push_back(i);
	    else ++external_count;
	  } // End num_i
	  if(external_count != 2) {
	    cout << "Factorize: Something is wrong in treatment of 4-RDM .... " << endl;
	    cout << "external_count: " << external_count << endl; 
	    abort(); 
	  } // End if
	} // End if
	//else if(the_tensors[1].get_name() == "X" && !X_indices.size()){
	else if(is_Interm(the_tensors[1].get_name()) && !X_indicesR2[1].size()){
	  for(size_t num_i = 0;num_i < the_tensors[1].get_indices().size();++num_i){
	    SQindex* i(the_tensors[1].get_indices()[num_i]);
	    if(!i->get_isExt()) Rindlist.push_back(i);
	  } // End num_i
	} // End if
	//else if(the_tensors[1].get_name() == "X" &&  X_indices.size()){
	else if(is_Interm(the_tensors[1].get_name()) &&  X_indicesR2[1].size()){
	  Rindlist.insert(Rindlist.end(), X_indicesR2[1].begin(), X_indicesR2[1].end());
	} // End if
	else if(is_D4C(the_tensors[1].get_name())){
	  for(size_t num_i = 0;num_i < the_tensors[1].get_indices().size();++num_i){
	    SQindex* i(the_tensors[1].get_indices()[num_i]);
	    if(num_i != extd4c_) Rindlist.push_back(i);
	  } // End num_i
	} // End else
	else{
	  Rindlist = the_tensors[1].get_indices();
	} // End else
	
	string nterm2("  " + the_tensors[1].get_name() + "_");
	if(Rindlist.size()){
	  nterm2 += "(";
	  for(size_t num_i = 0;num_i < Rindlist.size();++num_i){
	    nterm2 += "s_" + Rindlist[Rindlist.size()-(num_i+1)]->get_index();
	    if(num_i != Rindlist.size()-1) nterm2 += ", ";
	    else nterm2 += ")";
	  } // End num_i
	  
	  nterm2 += "%array("; 
	  for(size_t num_i = 0;num_i < Rindlist.size();++num_i){
	    nterm2 += "i_" + Rindlist[Rindlist.size()-(num_i+1)]->get_index();
	    if(num_i != Rindlist.size()-1) nterm2 += ", ";
	    else nterm2 += ")";
	  } // End num_i
	} // End if
	f << nterm2 << endl;
      } // End scope 2
      for(int n = 0;n < count2;++n) f << "end do"   << endl;
      f << endl;
    } // End if
    else is_allocated[1] =false;

    if(Z3inds.size()){
    // Allocate the resultant tensor
    f << "! Z3 <-- " << bin.get_Ltensor() << endl;
    f << "allocate(Z3_(";
    for(vector<SQindex*>::iterator i = Z3inds.begin();i != Z3inds.end();++i){
      string sname("s_" + (*i)->get_index());
      string otype;
      if     ((*i)->get_char() == Femto::core) otype = "I_C";
      else if((*i)->get_char() == Femto::act ) otype = "I_O";
      else if((*i)->get_char() == Femto::virt) otype = "I_V";
      
      if(i == Z3inds.begin())
      f << "psym(I_BEGIN," + otype + ", " + sname + "):psym(I_END," + otype + ", " + sname + ")";
      else
      f << "             psym(I_BEGIN," + otype + ", " + sname + "):psym(I_END," + otype + ", " + sname + ")";

      if(*i == Z3inds.back()) f << "))\n";
      else                    f << ", &\n";
    } // End i
    // Zero out Z3_, if neccessary
    f << endl;
    } // End if

    // Flop count
    string flopcount("flops = flops + ");

    // Call dgemm!
    if(bin.get_Ltensor().get_indices().size())
    f << "! Gemm Z1 * Z2 to form Z3\n";
    else
    f << "! Gemm Z1 * Z2 to form " << bin.get_Ltensor() << "\n";
    f << "call dgemm('n', 'n', ";

    // Declare m
    if(!rowInd.size()) m_label += "1";
    f << m_label << ",&" <<endl;
    // Declare n
    if(!colInd.size()) n_label += "1";
    f << "                     " << n_label << ",&" << endl;
    // Declare k
    if(!summedInd.size()) k_label += "1";
    f << "                     " << k_label << ",&" << endl;

    // Form flop countter
    flopcount += m_label + " * &\n";
    flopcount += "                " + n_label + " * &\n";   
    flopcount += "                " + k_label + " * 2.0d+00\n";

    string alpha("                     ");    
    // Generate the numerical coefficients
    if(bin.get_numConst() > 0) {
      ostringstream stm;
      stm << (boost::format("%10.8f") % bin.get_numConst());
      alpha += stm.str() + "d+00";
    } // end if
    else {
      ostringstream stm;
      stm << (boost::format("%10.8f") % fabs(bin.get_numConst()));
      alpha += "- " + stm.str() + "d+00";
    } // end if

    // Print all the Coeffient
    for(size_t i = 0;i < bin.get_Consts().size();++i){
      if(bin.get_Consts()[i] != "") alpha += "*" + bin.get_Consts()[i];
    } // End if
    alpha += ", &";
    f << alpha << endl;
    ////////////////////////////////////////////////////////////////////////////////////
    // ! Gemm !
    ////////////////////////////////////////////////////////////////////////////////////
    if(!bareRDM.first)
      f << "                     Z1_,&"               << endl;
    else{
      vector<SQindex*> Rindlist(the_tensors[0].get_indices());
	string nterm2("                     " + the_tensors[0].get_name() + "_");
	nterm2 += "(";
	for(size_t num_i = 0;num_i < Rindlist.size();++num_i){
	  nterm2 += "s_" + Rindlist[Rindlist.size()-(num_i+1)]->get_index();
	  if(num_i != Rindlist.size()-1) nterm2 += ", ";
	  else nterm2 += ")";
	} // End num_i
	nterm2 += "%array,&";
	f << nterm2 << endl;
    } // End else
    f << "                     " << m_label << ",&" << endl;
    if(!bareRDM.second)
    f << "                     Z2_,&"               << endl;
    else{
      vector<SQindex*> Rindlist(the_tensors[1].get_indices());
      string nterm2("                     " + the_tensors[1].get_name() + "_");
      nterm2 += "(";
      for(size_t num_i = 0;num_i < Rindlist.size();++num_i){
	nterm2 += "s_" + Rindlist[Rindlist.size()-(num_i+1)]->get_index();
	if(num_i != Rindlist.size()-1) nterm2 += ", ";
	else nterm2 += ")";
      } // End num_i
      nterm2 += "%array,&";
      f << nterm2 << endl;	
    } // End else
    f << "                     " << k_label << ",&" << endl;
    if(Z3inds.size() || bin.get_Ltensor().get_indices().size()){        
    f << "                     0.0d+00,&"           << endl;        
    f << "                     Z3_,&"               << endl;
    }
    else{
    f << "                     1.0d+00,&"           << endl;        
    f << "                     " << bin.get_Ltensor().get_name() << "_,&" << endl;
    }
    f << "                     " << m_label << ")" << endl << endl;        
    ////////////////////////////////////////////////////////////////////////////////////

    if(bin.get_Ltensor().get_indices().size()){
    // Sort Z3 to form bin.get_Ltensor()
    f << "! " << bin.get_Ltensor() << " <-- Z3" << endl;
    int count3 = 0;
    if(use_omp_)
      f << endl << "!$omp parallel do" << endl;
    //for(vector<SQindex*>::iterator i = L_nd_ind.begin();i != L_nd_ind.end();++i, ++count3){
    for(vector<SQindex*>::reverse_iterator i = L_nd_ind.rbegin();i != L_nd_ind.rend();++i, ++count3){ // Modified: 2013/05/22
      string iname = "i_" + (*i)->get_index();
      string sname = "s_" + (*i)->get_index();
      string otype;
      if     ((*i)->get_char() == (char_state)0) otype = "I_C";
      else if((*i)->get_char() == (char_state)1) otype = "I_O";
      else if((*i)->get_char() == (char_state)2) otype = "I_V";
      f << "do " + iname + " = psym(I_BEGIN, " + otype + ", " + sname + "), psym(I_END, " + otype + ", " + sname+ ")" << endl;
    } // End i

    // Generate the left-hand side tensor
    vector<SQindex*> Lindlist;
    //if     (!isBareLHS && bin.get_Ltensor().get_name() != "X"){
    if     (!isBareLHS && !is_Interm(bin.get_Ltensor().get_name()) && bin.get_Ltensor().get_name() != D4C_nameL()){
      Lindlist = bin.get_Ltensor().get_indices();
    } // End if
    else if(!isBareLHS && is_Interm(bin.get_Ltensor().get_name())){
      if(!X_indicesL.size()){
        for(size_t num_i = 0;num_i < bin.get_Ltensor().get_indices().size();++num_i){
          if(!bin.get_Ltensor().get_indices()[num_i]->get_isExt()) 
            Lindlist.push_back(bin.get_Ltensor().get_indices()[num_i]);
        } // End num_i
      } // End if
      else{
        Lindlist.insert(Lindlist.end(), X_indicesL.begin(), X_indicesL.end());
      } // End else
    } // End if
    else if(bin.get_Ltensor().get_name() == D4C_nameL()){
      for(size_t num_i = 0;num_i < bin.get_Ltensor().get_indices().size();++num_i){
	if(!bin.get_Ltensor().get_indices()[num_i]->get_isExt()) 
	  Lindlist.push_back(bin.get_Ltensor().get_indices()[num_i]);
      } // End num_i      
    } // End if
    else{
      for(size_t num_i = 0;num_i < bin.get_Ltensor().get_indices().size();++num_i){
        if(num_i != extamp_) 
          Lindlist.push_back(bin.get_Ltensor().get_indices()[num_i]);
      } // End num_i      
    } // End 

    string nterm(bin.get_Ltensor().get_name() + "_");
    if(Lindlist.size()){
      nterm += "(";
      for(size_t num_i = 0;num_i < Lindlist.size();++num_i){
        nterm += "s_" + Lindlist[Lindlist.size()-(num_i+1)]->get_index();
        if(num_i != Lindlist.size()-1) nterm += ", ";
        else nterm += ")";
      } // End num_i
      nterm += "%array(";
      for(size_t num_i = 0;num_i < Lindlist.size();++num_i){
        nterm += "i_" + Lindlist[Lindlist.size()-(num_i+1)]->get_index();
        if(num_i != Lindlist.size()-1) nterm += ", ";
        else nterm += ")";
      } // End num_i
    } // End if

    f << nterm + " = &\n    " + nterm + " &" << endl;
    f << "  + Z3_";
    if(Z3inds.size()) f << "(";
    else              f << endl;
    for(vector<SQindex*>::iterator i = Z3inds.begin();i != Z3inds.end();++i){
      string sname("i_" + (*i)->get_index());
      f << sname;
      if(*i == Z3inds.back()) f << ")\n";
      else                    f << ", ";
    } // End i
    for(int n = 0;n < count3;++n) f << "end do"   << endl;
    f << endl;
    } // End if

    // Create flop counter
    f << "! Flop count" << endl;
    f << flopcount << endl;

    // Deallocate the Z intermediates 
    vector<string> Z_names; Z_names.reserve(3);
    if(is_allocated[0]) Z_names.push_back("Z1_");
    if(is_allocated[1]) Z_names.push_back("Z2_");
    if(is_allocated[2]) Z_names.push_back("Z3_");
    if(Z_names.size()) f << "deallocate(";
    for(vector<string>::iterator s = Z_names.begin();s != Z_names.end();++s)
      if(s != Z_names.end()-1) f << *s << ", ";
      else                     f << *s << ")\n";

    if(extInd.size())                f << "end do ! Orbital Loop" << endl << endl;
    if(dimCount)                     f << "end if ! Dim Const"    << endl;
    f << endl;
                                     f << "end if ! Irrep Cond"   << endl;
    for(int n = 0;n < LoopCount;++n) f << "end do ! Irrep Loop"   << endl;
    f << "! FEMTO END  ****************************************************************\n" << endl;
    f << printEnd;
    f.flush(); 
    //X_indices.clear();
    //abort(); //*TEST*
  }


}} // Femto::

