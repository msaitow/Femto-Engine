//
//  Generate_Contract_new.cc
//  
//
//  Created by Masaaki Saitow on 13/10/29.
//  Copyright (c) 2012 Masaaki Saitow. All rights reserved.
//

#include <algorithm>
#include <fstream>
#include <cmath>
#include <string>
#include <boost/format.hpp>
#include <SQreaktor.hpp>

#define _DEBUG_LVL1

using namespace std;
using namespace Femto;
using namespace Femto::Core;

namespace Femto { namespace Reaktor {

  // *********************************************************
  // 
  // *********************************************************
  void SQreaktor::generate_contract_new(ofstream &CPfile, ofstream &CHfile, ofstream &F90file, vector<vector<SQbinary> > &theBins, string TheLabel)
  {
    if(LTensor_.get_indices().size() != 4 && isBareLHS_){
      cout << "Argument cannnot be treated as a bareampack .... " << endl;
      abort();
    } // End if

    vector<SQindex*> Linds(LTensor_.get_indices());
    for(vector<SQindex*>::iterator i = Linds.begin();i != Linds.end();++i)
      if((*i)->get_isSummed()){
        cout << *i << " in argument is a kind of dummy index." << endl;
        abort(); 
      } // End if

    if(is_sfGen(LTensor_.get_name())){
      cout << "Generate_Contract_new: Argument cannot be a spin-free unitary group generator ... " << endl;
      abort();
    }   

    bool ERIflag(false);
    bool D4Cflag(false);
    pair<char_state, string> nature;
    for(auto t = theBins.begin();t != theBins.end();++t){
      for(auto tt = t->begin();tt != t->end();++tt){
	for(size_t num_t = 0;num_t < tt->get_Rtensors().size();++num_t){
	  if(is_sfGen(tt->get_Rtensors()[num_t].get_name())){
	    cout << "Generate_Contract: A spin-free unitary group generator is detected in the argument .... " << endl;
	    abort();
	  } // End if
          else if(tt->get_Rtensors()[num_t].get_name() == name_h2_){
	    ERIflag = true;
            if     (tt->get_Rtensors()[num_t].get_indices()[exth2_]->get_char() == core) { nature.first = core; nature.second = "core"; }
            else if(tt->get_Rtensors()[num_t].get_indices()[exth2_]->get_char() == act ) { nature.first = act ; nature.second = "active"; }
	    else if(tt->get_Rtensors()[num_t].get_indices()[exth2_]->get_char() == virt) { nature.first = virt; nature.second = "virtual"; }
	  } // End else
          else if(is_D4C(tt->get_Rtensors()[num_t].get_name())){
	    D4Cflag = true;
            if     (tt->get_Rtensors()[num_t].get_indices()[extd4c_]->get_char() == core) { nature.first = core; nature.second = "core"; }
            else if(tt->get_Rtensors()[num_t].get_indices()[extd4c_]->get_char() == act ) { nature.first = act ; nature.second = "active"; }
	    else if(tt->get_Rtensors()[num_t].get_indices()[extd4c_]->get_char() == virt) { nature.first = virt; nature.second = "virtual"; }
	  } // End else
	} // End num_t
      } // End tt
    } // End t

    if(ERIflag)
      cout << ">> ERI is detected. So, I assume all the terms contain ERI with same fixed loading index[" << nature.second << "]." << endl;
    else if(D4Cflag)
      cout << ">> D4C is detected. So, I assume all the terms contain D4C with same fixed loading index[" << nature.second << "]." << endl;
    else 
      cout << ">> There is at least one term that doesn't contain neither ERIs, nor D4C. So, I assume all the terms are like that." << endl;

    if(!theBins.size()) return;

    vector<vector<SQbinary> > workObjs0; // :: Actual contraction groups (calculated outside ERI loop and prior to the workObjs1)
    vector<vector<SQbinary> > workObjs1; // :: Actual contraction groups (calculated within ERI loop)
    vector<vector<SQbinary> > workObjs2; // :: Actual contraction groups (calculated outside of ERI loop)
    vector<SQindex*> eriIndices;         // :: Loading indices of the ERI
    vector<int>      eriPos;             // :: Position of the contraction with the loading index of ERI
    SQcont<SQtensor> extTensor;          // :: Tensors that should be called outside of the grand loop
    if(!ERIflag && !D4Cflag) workObjs0 = theBins; // In case of no-eri and no-d4c terms, it's easy
    else{
      ///////////////////////////////////////////////////////
      // Detect which index is the loading index of the ERI
      ///////////////////////////////////////////////////////
      for(auto t = theBins.begin();t != theBins.end();++t){
	int myERI(0);
	for(auto tt = t->begin();tt != t->end();++tt, ++myERI){
	  vector<SQtensor*> ten_ptr(tt->get_Rtensors_ptr());
	  for(auto myt = ten_ptr.begin();myt != ten_ptr.end();++myt)
	    if     ((*myt)->get_name() == name_h2_){ 
	      eriPos.push_back(myERI);
	      eriIndices.push_back((*myt)->get_indices()[exth2_ ]);
	    } // End if
	    else if(is_D4C((*myt)->get_name())){ 
	      eriPos.push_back(myERI);
	      eriIndices.push_back((*myt)->get_indices()[extd4c_]);
	    } // End if
	} // End tt
      } // End t

      if(eriIndices.size() != theBins.size()) { 
	cout << "generate_contract_new: Is there several ERIs?" << endl;
	cout << "Debug info: -->> Contents if theBins <<--" << endl;
	int count(0);
	for(auto t = theBins.begin();t != theBins.end();++t,++count){
	  cout << boost::format(" <<-- [%7d] -->> \n") % count << endl;
	  SQcont<SQbinary> temp(*t);
	  for(auto tt = temp.begin();tt != temp.end();++tt) cout << "     >> " << *tt << endl;
	  cout << endl;
	} // End t
	abort(); 
      } // End if

      for(auto t = theBins.begin();t != theBins.end();++t){

	if(t->size() == 1) { workObjs1.push_back(*t); continue; }
	vector<vector<SQindex*> > ind_list;
	for(auto tt = t->begin();tt != t->end();++tt) ind_list.push_back(tt->get_summedBody());

	size_t contraNum(t-theBins.begin());
	vector<bool> hasERIind; // Whether the contraction has the loading index of the ERI
	for(auto inds = ind_list.begin();inds != ind_list.end();++inds){
	  bool flag(false);
	  for(auto i = inds->begin();i != inds->end();++i) 
	    if(**i == *eriIndices[contraNum]) flag = true;
	  if(flag) hasERIind.push_back(true);
	  else     hasERIind.push_back(false);
	} // End inds

	int myNum(0);
	vector<SQbinary> myObjs0;
	vector<SQbinary> myObjs1;
	vector<SQbinary> myObjs2;
	for(auto tt = t->begin();tt != t->end();++tt,++myNum){
	  if     (hasERIind[myNum])           myObjs1.push_back(*tt);
	  else if(myNum <  eriPos[contraNum]) myObjs0.push_back(*tt);
	  else                                myObjs2.push_back(*tt);
	} // End tt	
	if(myObjs0.size()) workObjs0.push_back(myObjs0);
	if(myObjs1.size()) workObjs1.push_back(myObjs1);
	if(myObjs2.size()) workObjs2.push_back(myObjs2);

      } // End t

    } // End else

    //////////////////////////////////////////////////////////////////////////////////
    // Cut the linkage between the contractions classified into the different groups
    //////////////////////////////////////////////////////////////////////////////////
    // -- [0] In case of the type0 contractions
    for(auto objs = workObjs0.begin();objs != workObjs0.end();++objs){
      if(objs->size()){
	auto myEnd(objs->end()); --myEnd;
	if(myEnd->get_Ltensor().get_name() != LTensor_.get_name()){
	  extTensor <= myEnd->get_Ltensor();
	  myEnd->set_Lindices(true);
	} // End if

	//------------------------------------------------
	// Set all the indices as internal once
	for(auto ts = objs->begin();ts != objs->end();++ts){
	  auto summedInds(ts->get_summedBody());
	  for(auto i = summedInds.begin();i != summedInds.end();++i)
	    (*i)->switch_isExt(false);
	} // End ts
	// Re-set the indices as external if it is necessary
	SQcont<SQindex> tempI;
	for(auto ts = objs->begin();ts != objs->end();++ts){
	  auto temp2(ts->get_Rtensors_ptr());
	  for(auto t = temp2.begin();t != temp2.end();++t)
	    if     ((*t)->get_name() == name_h2_ && tempI.count(*(*t)->get_indices()[exth2_])  == 0) tempI <= *(*t)->get_indices()[exth2_];
	    else if((*t)->get_name() == name_amp_&& tempI.count(*(*t)->get_indices()[extamp_]) == 0) tempI <= *(*t)->get_indices()[extamp_];
	  if(ts->get_Ltensor().get_name() == LTensor_.get_name() && isBareLHS_ && tempI.count(*ts->get_Ltensor().get_indices()[extamp_]) == 0) 
	    tempI <= *ts->get_Ltensor().get_indices()[extamp_]; 
	} // End ts
	for(auto ts = objs->begin();ts != objs->end();++ts){
	  SQcont<SQindex*> summedI(ts->get_summedBody());
	  for(auto i = summedI.begin();i != summedI.end();++i) if(tempI.count(**i)) (*i)->switch_isExt(true);
	} // End ts
	//------------------------------------------------

      } // End if
    } // End objs

    // -- [1] In case of the type1 contractions
    for(auto objs = workObjs1.begin();objs != workObjs1.end();++objs){
      if(objs->size()){
	auto myBegin(objs->begin());
	auto myEnd(objs->end()); --myEnd;
	auto bRtensors(myBegin->get_Rtensors_ptr());
	auto eLtensor(myEnd->get_Ltensor_ptr());
	// If the first contraction contains the intermediate on the right-hand side, it should be treated as type2 intermediate
	for(auto t = bRtensors.begin();t != bRtensors.end();++t){
	  if(is_Interm((*t)->get_name())) {
	    myBegin->set_Rindices((size_t)(t-bRtensors.begin()), true);
	    if(extTensor.count(**t) == 0) extTensor <= **t;
	  } // End if
	} // End t
	// If the last contraction contains the intermediate on the keft-hand side, it should be treated as te type2 intermediate
	if(is_Interm(eLtensor->get_name())) {
	  myEnd->set_Lindices(true);
	  extTensor <= *eLtensor;
	} // End if

	//------------------------------------------------
	// Set all the indices as internal once
	for(auto ts = objs->begin();ts != objs->end();++ts){
	  auto summedInds(ts->get_summedBody());
	  for(auto i = summedInds.begin();i != summedInds.end();++i)
	    (*i)->switch_isExt(false);
	} // End ts
	// Re-set the indices as external if it is necessary
	SQcont<SQindex> tempI;
	for(auto ts = objs->begin();ts != objs->end();++ts){
	  auto temp2(ts->get_Rtensors_ptr());
	  for(auto t = temp2.begin();t != temp2.end();++t)
	    if     ((*t)->get_name() == name_h2_ && tempI.count(*(*t)->get_indices()[exth2_])  == 0) tempI <= *(*t)->get_indices()[exth2_];
	    else if((*t)->get_name() == name_amp_&& tempI.count(*(*t)->get_indices()[extamp_]) == 0) tempI <= *(*t)->get_indices()[extamp_];
	    else if(is_D4C((*t)->get_name())     && tempI.count(*(*t)->get_indices()[extd4c_]) == 0) tempI <= *(*t)->get_indices()[extd4c_];
	  if(ts->get_Ltensor().get_name() == LTensor_.get_name() && isBareLHS_ && tempI.count(*ts->get_Ltensor().get_indices()[extamp_]) == 0) 
	    tempI <= *ts->get_Ltensor().get_indices()[extamp_]; 
	} // End ts
	for(auto ts = objs->begin();ts != objs->end();++ts){
	  SQcont<SQindex*> summedI(ts->get_summedBody());
	  for(auto i = summedI.begin();i != summedI.end();++i) if(tempI.count(**i)) (*i)->switch_isExt(true);
	} // End ts
	//------------------------------------------------

      } // End if
    } // End objs

    // -- [2] In case of the type2 contractions
    for(auto objs = workObjs2.begin();objs != workObjs2.end();++objs){
      if(objs->size()){
	auto myBegin(objs->begin());
	auto bRtensors(myBegin->get_Rtensors_ptr());
	// If the first contraction contains the intermediate on the right-hand side, it should be treated as type2 intermediate
	for(auto t = bRtensors.begin();t != bRtensors.end();++t){
	  if(is_Interm((*t)->get_name())) myBegin->set_Rindices((size_t)(t-bRtensors.begin()), true);
	} // End t

	//------------------------------------------------
	// Set all the indices as internal once
	for(auto ts = objs->begin();ts != objs->end();++ts){
	  auto summedInds(ts->get_summedBody());
	  for(auto i = summedInds.begin();i != summedInds.end();++i)
	    (*i)->switch_isExt(false);
	} // End ts
	// Re-set the indices as external if it is necessary
	SQcont<SQindex> tempI;
	for(auto ts = objs->begin();ts != objs->end();++ts){
	  auto temp2(ts->get_Rtensors_ptr());
	  for(auto t = temp2.begin();t != temp2.end();++t)
	    if     ((*t)->get_name() == name_h2_ && tempI.count(*(*t)->get_indices()[exth2_])  == 0) tempI <= *(*t)->get_indices()[exth2_];
	    else if((*t)->get_name() == name_amp_&& tempI.count(*(*t)->get_indices()[extamp_]) == 0) tempI <= *(*t)->get_indices()[extamp_];
	  if(ts->get_Ltensor().get_name() == LTensor_.get_name() && isBareLHS_ && tempI.count(*ts->get_Ltensor().get_indices()[extamp_]) == 0) 
	    tempI <= *ts->get_Ltensor().get_indices()[extamp_]; 
	} // End ts
	for(auto ts = objs->begin();ts != objs->end();++ts){
	  SQcont<SQindex*> summedI(ts->get_summedBody());
	  for(auto i = summedI.begin();i != summedI.end();++i) if(tempI.count(**i)) (*i)->switch_isExt(true);
	} // End ts
	//------------------------------------------------

      } // End if
    } // End objs
    //////////////////////////////////////////////////////////////////////////////////

#ifdef _DEBUG_LVL1
    {
      int count(0);
      cout << "++ Contents of workObjs0 :: " << endl;
      for(auto t = workObjs0.begin();t != workObjs0.end();++t,++count){
	SQcont<SQbinary> temp(*t);
	if(temp.size()) cout << temp << endl;
      } // End t
    }
    //cout << "korekore" << endl; abort();
    {
      int count(0);
      cout << "++ Contents of workObjs1 :: " << endl;
      for(vector<vector<SQbinary> >::iterator t = workObjs1.begin();t != workObjs1.end();++t,++count){
	SQcont<SQbinary> temp(*t);
	if(temp.size()) cout << temp << endl;
      } // End t
    }
    {
      int count(0);
      cout << "++ Contents of workObjs2 :: " << endl;
      for(vector<vector<SQbinary> >::iterator t = workObjs2.begin();t != workObjs2.end();++t,++count){
	SQcont<SQbinary> temp(*t);
	if(temp.size()) cout << temp << endl;
      } // End t
    }
#endif

    // Print the flag
    if(ERIflag || D4Cflag){
      SQindex i_eri("_eri", nature.first);
      string sign;
      if     (nature.first == Femto::core) sign = "c";
      else if(nature.first == Femto::act ) sign = "a";
      else if(nature.first == Femto::virt) sign = "v";

      if      (ERIflag) CPfile << "//-@loadERI(" << sign << ",begin)" << endl;
      else if (D4Cflag)	CPfile << "//-@loadD4C(" << sign << ",begin)" << endl;
    } // End if

    // Let's begin to generate!!!
    CPfile << "  //*-- FEMTO begins --//*" << endl;
    CPfile << "  // Label : " << TheLabel << endl;
    CPfile << "  {" << endl << endl;

    size_t LoopCount(0); // :: Represents the depth of the loop of the current position
    string indent("");
    SQcont<string> Indents; Indents.reserve(100);
    for(int i = 0;i < 100;++i){
      Indents.push_back(indent);
      indent += "  ";
    } // End i 

    if(ERIflag || D4Cflag){
      if(extTensor.size()) {
	CPfile << "//-@type(2).declaration(begin)" << endl;
	CPfile << "  // --  Title : " << title_ << endl;            
	CPfile << "  //  >> Intermediates for the external contractions are defined here << " << endl;
      } // End if

      for(auto et = extTensor.begin();et != extTensor.end();++et)
	declareInterm(LoopCount, CPfile, *et, SQcont<SQindex*>(et->get_indices()), Femto::Reaktor::External);

      if(extTensor.size()) { // In case of the external contraction or the one-body contractions
	CPfile << "//-@type(2).declaration(end)" << endl;

	// Make the type 0 contractions
	if(workObjs0.size()){ 
	  CPfile << endl;
	  int numTerm0(0);
	  CPfile << "  //*-- Entering to take the type 0 contractions --*//" << endl;
	  CPfile << "//-@type(0).contraction(begin)" << endl;
	  for(auto b = workObjs0.begin();b != workObjs0.end();++b,++numTerm0){
	    CPfile << endl;
	    ostringstream stm;
	    stm << (int)(b-workObjs0.begin());
	    //CPfile << "  if(myrank == 0)" << endl;
	    CPfile << "  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< " << endl;
	    CPfile << boost::format("  // -- No. [%5d] -- ") % numTerm0 << endl;
	    for(auto tt = b->begin();tt != b->end();++tt){
	      CPfile << boost::format("  // |-- [%5d] --| ") % ((int)(tt-b->begin())) << *tt << endl;
	    } // End tt
	    CPfile << "  double flops(0); // Flop count  " << endl;
	    makeContractions(LoopCount, CPfile, CHfile, F90file, stm.str(), "_type0_"+TheLabel, *b);
	    CPfile << "  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> " << endl;
	    CPfile << endl;
	  } // End b
	  CPfile << "//-@type(0).contraction(end)" << endl; 
	} // End if

      } // End if

      CPfile << endl;
      SQindex i_eri("_eri", nature.first);
      CLoop(Indents[LoopCount], i_eri, CPfile);

      if(ERIflag){
	vector<SQindex*> inds(4, &i_eri);
	SQtensor ERI(name_h2_, inds, h2_symm());
	ReadERI(Indents[LoopCount], ERI, CPfile);
      }
      else if(D4Cflag){
	vector<SQindex*> inds(6, &i_eri);
        SQtensor D4C(D4C_name(), inds, d4c_symm());       
        MPIiproc(Indents[LoopCount], i_eri, CPfile);
	ReadD4C(Indents[LoopCount], D4C, CPfile);
      }
      CPfile << endl;
    } // End if

    // Contract workObjs1 here
    if(ERIflag || D4Cflag){
      CPfile << "  //*-- Entering to take the type 1 contractions --*//" << endl;
      CPfile << "//-@type(1).contraction(begin)" << endl;
      CPfile << "  // -- Title : " << title_ << endl;      
    }    
    else{
      // In case of the one-body contractions
      if(workObjs0.size()){
	CPfile << endl;
	CPfile << "  //*-- Entering to take the type 0 contractions --*//" << endl;
	//CPfile << "//-@type(0).contraction(begin)" << endl;
	int numTerm0(0);
	for(auto b = workObjs0.begin();b != workObjs0.end();++b,++numTerm0){
	  CPfile << endl;
	  ostringstream stm;
	  stm << (int)(b-workObjs0.begin());
	  CPfile << "  if(myrank == 0)" << endl;
	  CPfile << "  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< " << endl;
	  CPfile << boost::format("  // -- No. [%5d] -- ") % numTerm0 << endl;
	  for(auto tt = b->begin();tt != b->end();++tt){
	    CPfile << boost::format("  // |-- [%5d] --| ") % ((int)(tt-b->begin())) << *tt << endl;
	  } // End tt
	  CPfile << "  double flops(0); // Flop count  " << endl;
	  makeContractions(LoopCount, CPfile, CHfile, F90file, stm.str(), "_type0_"+TheLabel, *b);
	  CPfile << "  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> " << endl;
	  CPfile << endl;
	} // End b
	//CPfile << "//-@type(0).contraction(end)" << endl; 
      } // End if

    } // End else
    
    int numTerm1(0);
    for(auto t = workObjs1.begin();t != workObjs1.end();++t,++numTerm1){

      CPfile << endl;
      ostringstream stm;
      stm << (int)(t-workObjs1.begin());
      cout << "STM1 << " << stm.str() << endl;
      cout << SQcont<SQbinary>(*t) << endl;
      CPfile << "  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< " << endl;
      CPfile << boost::format("  // -- No. [%5d] -- ") % numTerm1 << endl;
      for(auto tt = t->begin();tt != t->end();++tt){
	CPfile << boost::format("  // |-- [%5d] --| ") % ((int)(tt-t->begin())) << *tt << endl;
      } // End tt
      SQindex i_eri(ERIflag ? findERIIndex(*t) : findD4CIndex(*t));
      CPfile << "  int s" + i_eri.get_index() + "(s_eri);" << endl;
      CPfile << "  int i" + i_eri.get_index() + "(i_eri);" << endl;
      CPfile << "  double flops(0); // Flop count  " << endl;
      makeContractions(LoopCount, CPfile, CHfile, F90file, stm.str(), "_type1_"+TheLabel, *t);
      if(t->back().get_Ltensor().get_name() == LTensor_.get_name() && isBareLHS_)
	if(*(t->back().get_Ltensor().get_indices()[extamp_]) == i_eri)
	  AccAmp(Indents[LoopCount], t->back().get_Ltensor(), CPfile);
      CPfile << "  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> " << endl;
      cout << "STM2 >> " << stm.str() << endl;
      
#ifdef _VERBOSE_MODE
      CPfile << "#ifdef _VERBOSE" << endl;
#ifndef _CHRONO
      CPfile << "  //t_elapsed[\"" + title_ + "_workObj1_" + TheLabel + stm.str() + "\"] += difftime(time(NULL), t_start); // VERBOSE MODE" << endl; 
      CPfile << "  //t_elapsed[\"" + title_ + "_workObj1_" + TheLabel + stm.str() + "\"] += orz::wallclocktime() - t_start; // VERBOSE MODE" << endl; 
      CPfile << "  t_elapsed[\"" + title_ + "_workObj1_" + TheLabel + stm.str() + "\"] += ticktack() - t_start; // VERBOSE MODE" << endl; 
#else
      CPfile << "  boost::chrono::high_resolution_clock::time_point t_stop(boost::chrono::high_resolution_clock::now()); // VERBOSE MODE" << endl;
      CPfile << "  boost::chrono::duration<double> elapsed = t_stop - t_start;" << endl;
      CPfile << "  t_elapsed[\"" + title_ + "_workObj1_" + TheLabel + stm.str() + "\"] += elapsed.count(); // VERBOSE MODE" << endl; 
#endif
      CPfile << "  t_flops  [\"" + title_ + "_workObj1_" + TheLabel + stm.str() + "\"] += flops;                         // VERBOSE MODE" << endl; 
      CPfile << "#endif" << endl;
#endif
      
      //CPfile << "  } // End scope" << endl;

      CPfile << endl;
    } // End t

    if(ERIflag || D4Cflag) {

      if(LTensor_.get_name() == D4C_nameL())
	CPfile << "  orz::SaveBin(ctinp.dir()/(format(\"D4C_g[%d]\")%i_eri).str()) << retval; //   cout << \"made ok \" << format(\"D4C_g[%d]\")%i_eri << endl;" << endl << endl;

      CPfile << "//-@type(1).contraction(end)" << endl << endl;
      
      CPfile << "  } // End myrank" << endl;
      LoopEnd(Indents[0], CPfile);
      CPfile << "  orz::world().barrier();" << endl;
      
    } // End if

    //cout << "XXXXXXXXXXXXX " << endl; abort();   

    // Contract workObjs2 here
    if(workObjs2.size()){
      CPfile << endl;
      CPfile << "//-@type(2).contraction(begin)" << endl;
      CPfile << "  // -- Title : " << title_ << endl;      
      CPfile << "  //*-- Entering to take the type 2 contractions --*//" << endl;
    } // End if

    int numTerm2(0);
    for(auto t = workObjs2.begin();t != workObjs2.end();++t,++numTerm2){

      CPfile << endl;
      ostringstream stm;
      stm << (int)(t-workObjs2.begin());
      cout << "STM3 << " << stm.str() << endl;
      cout << SQcont<SQbinary>(*t) << endl;
      CPfile << "  { // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< " << endl;
      CPfile << boost::format("  // -- No. [%5d] -- ") % numTerm2 << endl;
      for(auto tt = t->begin();tt != t->end();++tt){
	CPfile << boost::format("  // |-- [%5d] --| ") % ((int)(tt-t->begin())) << *tt << endl;
      } // End tt
      CPfile << "  double flops(0); // Flop count  " << endl;
      makeContractions(LoopCount, CPfile, CHfile, F90file, stm.str(), "_type2_"+TheLabel, *t);
      CPfile << "  } // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> " << endl;
      cout << "STM4 >> " << stm.str() << endl;

    } // End t

    if(workObjs2.size()){ 
      //CPfile << "  } // End myrank" << endl;
      CPfile << "//-@type(2).contraction(end)" << endl; 
      CPfile << "  orz::world().barrier();" << endl;
    } // End if

#ifdef _VERBOSE_MODE
    CPfile << endl;
    CPfile << "#ifdef _VERBOSE" << endl;
    CPfile << "  // VERBOSE MODE" << endl;
    CPfile << "  double t_total(0);" << endl;
    CPfile << "  for(std::map<std::string, double>::const_iterator t = t_elapsed.begin(); t != t_elapsed.end();++t) {" << endl;
    CPfile << "    double flop(t->second ? t_flops[t->first.c_str()]/t->second : 0.0); " << endl;
    CPfile << "    t_total += t->second;" << endl; 
    CPfile << "    outfile << boost::format(\" * %5s : time=%14.9f, flops=%14.9e\") % t->first % t->second % flop << endl;                      " << endl;
    CPfile << "  } // End t " << endl;
    CPfile << "  outfile << \"Total time (s) \" << t_total << endl;" << endl;
    if(ERIflag){
      CPfile << "  outfile << \"Time (s) for loading  of ERI[" << TheLabel << "] \" << t_readERI << endl;" << endl;
    } // End if
    CPfile << "#endif" << endl;
#endif

    CPfile << endl;
    CPfile << "  } // End Femto" << endl;
    CPfile << "  //*-- FEMTO ends --//*" << endl << endl;;

    if(ERIflag || D4Cflag) {
      string sign;
      if     (nature.first == Femto::core) sign = "c";
      else if(nature.first == Femto::act ) sign = "a";
      else if(nature.first == Femto::virt) sign = "v";
      if (ERIflag) CPfile << "//-@loadERI(" << sign << ",end)" << endl << endl;
      else         CPfile << "//-@loadD4C(" << sign << ",end)" << endl << endl;
    } // End if
    
  }                                               

}} //End Femto
