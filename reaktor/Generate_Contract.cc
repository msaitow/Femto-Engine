//
//  Generate_Contract.cc
//  
//
//  Created by Masaaki Saitow on 12/10/29.
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
  void SQreaktor::generate_contract(ofstream &CPfile, ofstream &CHfile, ofstream &F90file, vector<vector<SQbinary> > &theBins, string TheLabel)
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
      cout << "Generate_Contract: Argument cannot be a spin-free unitary group generator ... " << endl;
      abort();
    }   

    bool ERIflag(false);
    bool D4Cflag(false);
    pair<char_state, string> nature;
    for(vector<vector<SQbinary> >::iterator t = theBins.begin();t != theBins.end();++t){
      for(vector<SQbinary>::iterator tt = t->begin();tt != t->end();++tt){
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

    vector<vector<SQbinary> > workObjs1; // :: Actual contraction groups (calculated within ERI loop)
    vector<vector<SQbinary> > workObjs2; // :: Actual contraction groups (calculated outside of ERI loop)
    // NOTE :: If there're no eri terms, all the terms are pushed back into workObjs1 and processed as usual.
    if(!ERIflag && !D4Cflag) workObjs1 = theBins; // In case of no-eri and no-d4c terms, it's easy
    // But if terms have eris or d4c, terms have to divided into two-groups, one of which is 
    // loops associated with eri ends after the first contraction, the other of which is 
    // sigma vector is summed up whithin the loop. 
    else{
      for(vector<vector<SQbinary> >::iterator t = theBins.begin();t != theBins.end();++t){
        // If the loading index of ERI, or D4C is contained in the X-tensor, first and second contractions are pushed back 
        // into workObjs1. In the another case, only the first contraction is treated as a member of workObjs1, 
        // and its counterpart is added to workObjs2.
        if(t->size() == 1) { workObjs1.push_back(*t); continue; }
        vector<SQindex*> x_inds(t->at(0).get_Ltensor().get_indices());
        SQindex* i_eri(NULL);
        for(size_t num = 0;num < t->at(0).get_Rtensors().size();++num)
          if     (t->at(0).get_Rtensors()[num].get_name() == name_h2_) i_eri = t->at(0).get_Rtensors()[num].get_indices()[exth2_]; 
          else if(is_D4C(t->at(0).get_Rtensors()[num].get_name())) i_eri = t->at(0).get_Rtensors()[num].get_indices()[extd4c_]; 
        if(find(x_inds.begin(), x_inds.end(), i_eri) != x_inds.end()) workObjs1.push_back(*t);
        else{
          t->at(0).set_Lindices(true);
          for(size_t num = 0;num < t->at(1).get_Rtensors().size();++num)
            //if(t->at(1).get_Rtensors()[num].get_name().at(0) == 'X') t->at(1).set_Rindices(num, true); 
            if(is_Interm(t->at(1).get_Rtensors()[num].get_name())) t->at(1).set_Rindices(num, true); 
          vector<SQbinary> cont1; cont1.push_back(t->at(0));
          vector<SQbinary> cont2; cont2.push_back(t->at(1));
          workObjs1.push_back(cont1);
          workObjs2.push_back(cont2);
	} // End else
      } // End t
    } // End else

    // Count number of the bachelors, which originates terms divided into workObjs1 and 2. 
    // And cut down connection in terms of internal and external indices
    pair<int, int> numBachelor;
    numBachelor.first = 0; numBachelor.second = 0;
    // Search in workObjs1
    for(vector<vector<SQbinary> >::iterator t = workObjs1.begin();t != workObjs1.end();++t){
      if(t->at(0).get_Lindices()){
        // Set all the indices internal, once.
        vector<SQindex*> summed(t->at(0).get_summedBody());
        for(vector<SQindex*>::iterator i = summed.begin();i != summed.end();++i) (*i)->switch_isExt(false);
        // And re-set indices as external, if it is necessary.
        for(size_t num = 0;num < t->at(0).get_Rtensors().size();++num){
          if     (t->at(0).get_Rtensors()[num].get_name() == name_h2_)    t->at(0).get_Rtensors()[num].get_indices()[exth2_]->switch_isExt(true);
          else if(t->at(0).get_Rtensors()[num].get_name() == name_amp_)   t->at(0).get_Rtensors()[num].get_indices()[extamp_]->switch_isExt(true);
          else if(is_D4C(t->at(0).get_Rtensors()[num].get_name()))        t->at(0).get_Rtensors()[num].get_indices()[extd4c_]->switch_isExt(true);
          else if(t->at(0).get_Rtensors()[num].get_name() == name_d4_)  { t->at(0).get_Rtensors()[num].get_indices()[0]->switch_isExt(true);
	                                                                  t->at(0).get_Rtensors()[num].get_indices()[1]->switch_isExt(true); }
	} // End num
        ++numBachelor.first;
      } // End if
    } // End t
    // Search in workObjs2
    for(vector<vector<SQbinary> >::iterator t = workObjs2.begin();t != workObjs2.end();++t){
      if(t->at(0).get_Rindices()[0] || t->at(0).get_Rindices()[1]){
        // Set all the indices internal, once.
        vector<SQindex*> summed(t->at(0).get_summedBody());
        for(vector<SQindex*>::iterator i = summed.begin();i != summed.end();++i) (*i)->switch_isExt(false);
        // And re-set indices as external, if it is necessary.
        for(size_t num = 0;num < t->at(0).get_Rtensors().size();++num){
          if      (t->at(0).get_Rtensors()[num].get_name() == name_h2_)    t->at(0).get_Rtensors()[num].get_indices()[exth2_]->switch_isExt(true);
          else if (t->at(0).get_Rtensors()[num].get_name() == name_amp_)   t->at(0).get_Rtensors()[num].get_indices()[extamp_]->switch_isExt(true);
          else if (is_D4C(t->at(0).get_Rtensors()[num].get_name()))        t->at(0).get_Rtensors()[num].get_indices()[extd4c_]->switch_isExt(true);
          else if (t->at(0).get_Rtensors()[num].get_name() == name_d4_)  { t->at(0).get_Rtensors()[num].get_indices()[0]->switch_isExt(true);
	                                                                   t->at(0).get_Rtensors()[num].get_indices()[1]->switch_isExt(true); }
	} // End num
        if(isBareLHS_) t->at(0).get_Ltensor().get_indices()[extamp_]->switch_isExt(true);
        ++numBachelor.second;
      } // End if
    } // End t
    
    //*OK* if(numBachelor.first != numBachelor.second){
    //*OK*   cout << "Generate_Conbtract: An error occured in setting external in workObja1 and 2" << endl;
    //*OK*   cout << " -- workObjs1 :: " << numBachelor.first  << endl;
    //*OK*   cout << " -- workObjs2 :: " << numBachelor.second << endl;
    //*OK*   abort();
    //*OK* } // End if

#ifdef _DEBUG_LVL1
    {
      int count(0);
      cout << "++ Contents of workObjs1 :: " << endl;
      for(vector<vector<SQbinary> >::iterator t = workObjs1.begin();t != workObjs1.end();++t,++count){
	cout << boost::format("[%5d] -- ") % count << t->at(0) << endl;
	if(t->size() == 2){
	  cout << "        -- " << t->at(1) << endl;
	} // End if
        cout << endl;
      } // End t
    }
    {
      int count(0);
      cout << "++ Contents of workObjs2 :: " << endl;
      for(vector<vector<SQbinary> >::iterator t = workObjs2.begin();t != workObjs2.end();++t,++count){
	cout << boost::format("[%5d] -- ") % count << t->at(0) << endl;
	if(t->size() == 2){
	  cout << "        -- " << t->at(1) << endl;
	} // End if
        cout << endl;
      } // End t
    }
#endif

    if(ERIflag || D4Cflag){
      string sign;
      if     (nature.first == Femto::core) sign = "c";
      else if(nature.first == Femto::act ) sign = "a";
      else if(nature.first == Femto::virt) sign = "v";
      if     (ERIflag) CPfile << "//-@loadERI(" << sign << ",begin)" << endl;
      else if(D4Cflag) CPfile << "//-@loadD4C(" << sign << ",begin)" << endl;
    }

    // Let's begin to generate!!!
    CPfile << "  //*-- FEMTO begins --//*" << endl;
    CPfile << "  // Label : " << TheLabel << endl;
    CPfile << "  {" << endl << endl;

#ifdef _VERBOSE_MODE
    {
      CPfile << "#ifdef _VERBOSE" << endl;
      { // * Timing each contraction
	CPfile << "  // Timing object for each term or contraction (VERBOSE MODE)" << endl;
	CPfile << "  std::map<std::string, double> t_elapsed;" << endl;
	int countObj1(0);
	int countObj2(0);
	for(vector<vector<SQbinary> >::const_iterator t = workObjs1.begin();t != workObjs1.end();++t){
	  ostringstream stm;
	  stm << countObj1++;
	  CPfile << "  t_elapsed[\"" + title_ + "_workObj1_" + TheLabel + stm.str() + "\"] = 0.0; // VERBOSE MODE" << endl;
	} // End t
	for(vector<vector<SQbinary> >::const_iterator t = workObjs2.begin();t != workObjs2.end();++t){
	  ostringstream stm;
	  stm << ++countObj2;
	  CPfile << "  t_elapsed[\"" + title_ + "_workObj2_" + TheLabel + stm.str() + "\"] = 0.0; // VERBOSE MODE" << endl;
	} // End t
	CPfile << endl;
      } // End scope
      { // * Flop counts
	CPfile << "  // Flop counts for each term or contraction (VERBOSE MODE)" << endl;
	CPfile << "  std::map<std::string, double> t_flops;" << endl;
	int countObj1(0);
	int countObj2(0);
	for(vector<vector<SQbinary> >::const_iterator t = workObjs1.begin();t != workObjs1.end();++t){
	  ostringstream stm;
	  stm << countObj1++;
	  CPfile << "  t_flops[\"" + title_ + "_workObj1_" + TheLabel + stm.str() + "\"] = 0.0; // VERBOSE MODE" << endl;
	} // End t
	for(vector<vector<SQbinary> >::const_iterator t = workObjs2.begin();t != workObjs2.end();++t){
	  ostringstream stm;
	  stm << ++countObj2;
	  CPfile << "  t_flops[\"" + title_ + "_workObj2_" + TheLabel + stm.str() + "\"] = 0.0; // VERBOSE MODE" << endl;
	} // End t
      } // End scope
      if(ERIflag) { // * Timing of loading of ERI if a bunch of binary contractions has it
	CPfile << endl;
	CPfile << "  // Time consumed for loading of ERI" << endl;
	CPfile << "  double t_readERI(0);" << endl;
      } // End if
      CPfile << "#endif" << endl;
      CPfile << endl;
    } // End scope
#endif

    size_t LoopCount(0); // :: Represents the depth of the loop of the current position
    string indent("");
    vector<string> Indents; Indents.reserve(100);
    for(int i = 0;i < 100;++i){
      Indents.push_back(indent);
      indent += "  ";
    } // End i 

    if(ERIflag || D4Cflag){
      vector<SQbinary*> tensors;
      for(vector<vector<SQbinary> >::iterator t = workObjs1.begin();t != workObjs1.end();++t){
        //if(t->size() == 1 && t->at(0).get_Ltensor().get_name().at(0) == 'X') tensors.push_back(&(t->at(0)));
        if(t->size() == 1 && is_Interm(t->at(0).get_Ltensor().get_name())) tensors.push_back(&(t->at(0)));
      } // End t

      if(tensors.size()) {
	CPfile << "//-@type(2).declaration(begin)" << endl;
	//CPfile << "// Label : " << TheLabel << endl;
	CPfile << "  // --  Title : " << title_ << endl;            
	CPfile << "  //  >> Intermediates for the type 2 contractions are defined here << " << endl;
      } // End if
      for(vector<SQbinary*>::iterator t = tensors.begin();t != tensors.end();++t){
        //int num_intind(0); // Number of internal indices
        int Ccount(0);
        int Ocount(0);
        int Vcount(0);
        if((*t)->get_Lindices()){
          for(size_t num = 0;num < (*t)->get_Ltensor().get_indices().size();++num)
            if     ((*t)->get_Ltensor().get_indices()[num]->get_char() == core) ++Ccount;
            else if((*t)->get_Ltensor().get_indices()[num]->get_char() == act ) ++Ocount;
            else if((*t)->get_Ltensor().get_indices()[num]->get_char() == virt) ++Vcount;
	} // End if
	else{
          cout << "Generate_Contract: Something is wrong in definition step of the intermediate tensor" << endl;
          abort();
	} // End else
        
	if(Ccount || Ocount || Vcount){
          vector<SQindex*> extI;
          string SymmLabel;
          if(!(*t)->get_Lindices())
	    for(size_t num_i = 0;num_i < (*t)->get_Ltensor().get_indices().size();++num_i){
	      if((*t)->get_Ltensor().get_indices()[num_i]->get_isExt())
		extI.push_back((*t)->get_Ltensor().get_indices()[num_i]);
	    } // End num_i
          if(!extI.size()) SymmLabel = "0";
          else{
            for(vector<SQindex*>::iterator i = extI.begin();i != extI.end();++i){
              SymmLabel += "s" + (*i)->get_index();
              if(*i != extI.back()) SymmLabel += "^";
	    } // End i
	  } // End else
          string DecSymm("  orz::DTensor " + (*t)->get_Ltensor().get_name());
          string Xname("");
          for(int c = 0;c < Ccount;++c) Xname += "c";
          for(int o = 0;o < Ocount;++o) Xname += "a";
          for(int v = 0;v < Vcount;++v) Xname += "v";
	  
	  //DecSymm += Xname + "(orz::mr::sizeof_sympack_X" + Xname + "(" + "symblockinfo, " + SymmLabel + "));"; //*MOD* 20130830
	  DecSymm += Xname + "_" + title_ +  "(orz::mr::sizeof_sympack_X" + Xname + "(" + "symblockinfo, " + SymmLabel + "));";
          CPfile << DecSymm << endl;
        } // End if
        else{
          //string DecX("  double " + (*t)->get_Ltensor().get_name() + " = 0;"); 
          string DecX("  double " + (*t)->get_Ltensor().get_name() + "_" + title_ + " = 0;"); 
          CPfile << DecX << endl;
        } // End else        
        
     } // End t
      if(tensors.size()) {
	CPfile << "//-@type(2).declaration(end)" << endl;
      } // End if

      SQindex i_eri("_eri", nature.first);
      CLoop(Indents[LoopCount], i_eri, CPfile);
      if(ERIflag){
	vector<SQindex*> inds(4, &i_eri);
	SQtensor ERI(name_h2_, inds, h2_symm());
	ReadERI(Indents[LoopCount], ERI, CPfile);

//*// 	string sign;
//*// 	if     (i_eri.get_char() == Femto::core) sign = "c";
//*// 	else if(i_eri.get_char() == Femto::act ) sign = "a";
//*// 	else if(i_eri.get_char() == Femto::virt) sign = "v";
//*// 	CPfile << "//-@loadERI(" << sign << ",begin)" << endl;

      }
      else if(D4Cflag){
	vector<SQindex*> inds(6, &i_eri);
        SQtensor D4C(D4C_name(), inds, d4c_symm());       
        MPIiproc(Indents[LoopCount], i_eri, CPfile);
	ReadD4C(Indents[LoopCount], D4C, CPfile);

//*// 	string sign;
//*// 	if     (i_eri.get_char() == Femto::core) sign = "c";
//*// 	else if(i_eri.get_char() == Femto::act ) sign = "a";
//*// 	else if(i_eri.get_char() == Femto::virt) sign = "v";
//*// 	CPfile << "//-@loadD4C(" << sign << ",begin)" << endl;

      }
      CPfile << endl;
    } // End if
//*OLD*//     else{
//*OLD*//       CPfile << "  if(myrank == 0){" << endl;
//*OLD*//     } // End else
    //*debug**// cout << "korekore1111" << endl; // *debug*
    //*debug**// cout << "Size :: " << theBins[0].size() << endl; // *debug
    //*debug**// cout << theBins[0][0] << endl; //*debug
    //*debug**// cout << theBins[0][1] << endl; //*debug
    if(theBins[0].size() > 2)
      if(theBins[0][1].get_Ltensor().get_name() == D4C_nameL())
	CPfile << "  orz::DTensor retval(orz::mr::sizeof_sympack_Xaaaaa(symblockinfo, s_eri));" << endl << endl;;

    // Contract workObjs1 here
    CPfile << "  //*-- Entering to take the type 1 contractions --*//" << endl;
    if(ERIflag || D4Cflag){
      CPfile << "//-@type(1).contraction(begin)" << endl;
      CPfile << "  // -- Title : " << title_ << endl;      
    }
    ////////cout << "korekore2222" << endl;
    int numTerm = 0;
    for(vector<vector<SQbinary> >::iterator t = workObjs1.begin();t != workObjs1.end();++t,++numTerm){

      // Base group of the indices
      vector<SQindex*> O1;
      // Even if don't re-set all the indices of terms divided into workObjs1 and 2, by modifing the ``if" sentence, it will work.
      // But it also will be rather less efficient.
      if(!(t->size() == 1 && is_Interm(t->at(0).get_Ltensor().get_name()))){
        for(size_t num = 0;num < t->at(0).get_Ltensor().get_indices().size();++num)
          if(t->at(0).get_Ltensor().get_indices()[num]->get_isExt() && find(O1.begin(), O1.end(), t->at(0).get_Ltensor().get_indices()[num]) == O1.end()) 
	    O1.push_back(t->at(0).get_Ltensor().get_indices()[num]);
      } // End if

      ///////////////////////////////////////////////////////////////////////
      // Flag for parallelization of the one-body integrals
      pair<size_t, bool> noIntflags(0, false);
      if(!ERIflag && !D4Cflag) noIntflags.second = true;
      if(!O1.size() && !ERIflag && !D4Cflag) { 
	CPfile << "  if(myrank == 0){" << endl;
	noIntflags.second           = false;
      } // End if
      ///////////////////////////////////////////////////////////////////////

      size_t LoopCount  =  0;
      size_t SigmaCount = -1;
      size_t D4Count    = -1;
      vector<SQindex*> Declared;
      vector<SQtensor> t_list;    // List of tensors already read from GA

      ostringstream stm;
      stm << numTerm;
      string title_1if("g_if_" + title_ + "_no0_x" + stm.str() + "_type1_" + TheLabel);
      string title_1  ("g_"    + title_ + "_no0_x" + stm.str() + "_type1_" + TheLabel);
      
      transform(title_1if.begin(), title_1if.end(), title_1if.begin(), (int(*)(int))tolower); 
      transform(title_1  .begin(), title_1  .end(), title_1  .begin(), (int(*)(int))tolower); 

      // Beginning of the Scope
      CPfile << "  { " << endl;
      CPfile << "  // No. " << stm.str() << ", [" << t->size() << "]" << endl;
      CPfile << "  // "     << t->at(0)  << endl;
      if(t->size() == 2)
      CPfile << "  // "     << t->at(1)  << endl;
      CPfile << "  double flops = 0; // Flop count" << endl;

#ifdef _VERBOSE_MODE
      CPfile << "#ifdef _VERBOSE" << endl;
#ifndef _CHRONO
      CPfile << "  //time_t t_start = time(NULL); // VERBOSE MODE" << endl;
      CPfile << "  //double t_start(orz::wallclocktime()); // VERBOSE MODE" << endl;
      CPfile << "  double t_start(ticktack()); // VERBOSE MODE" << endl;
#else
      CPfile << "  boost::chrono::high_resolution_clock::time_point t_start(boost::chrono::high_resolution_clock::now()); // VERBOSE MODE" << endl;
#endif
      CPfile << "#endif" << endl;
#endif

      // If there're ERIs (or D4Cs), catch the name of the loading index
      if(ERIflag || D4Cflag){
        SQindex* i_eri(NULL);
	for(size_t num = 0;num < t->at(0).get_Rtensors().size();++num)
	  if     (t->at(0).get_Rtensors()[num].get_name() == name_h2_)   i_eri = t->at(0).get_Rtensors()[num].get_indices()[exth2_];
	  else if(t->at(0).get_Rtensors()[num].get_name() == D4C_name()) i_eri = t->at(0).get_Rtensors()[num].get_indices()[extd4c_];
        if(i_eri == NULL){
          cout << "Generate_Contract: Can't find the loading index for contract " << stm.str() << endl;
          abort();
	} // End if
	CPfile << "  int s" + i_eri->get_index() + "(s_eri);" << endl;
	CPfile << "  int i" + i_eri->get_index() + "(i_eri);" << endl;
	vector<SQindex*>::iterator itr(find(O1.begin(), O1.end(), i_eri));
	if(itr != O1.end()) O1.erase(itr); // Eliminate the index of ERI if exists in O1
        Declared.push_back(i_eri);
	for(size_t num = 0;num < t->at(0).get_Rtensors().size();++num){
	  if(t->at(0).get_Rtensors()[num].get_name() == name_amp_ && t->at(0).get_Rtensors()[num].get_indices()[extamp_] == i_eri){
	    ReadAmp(Indents[LoopCount], t->at(0).get_Rtensors()[num], CPfile);
            t_list.push_back(t->at(0).get_Rtensors()[num]);
	  } // End if
	} // End num

        // In case of sigma (possibly) vector
	if(t->size() == 2){
	  if(isBareLHS_ && t->at(1).get_Ltensor().get_name() == LTensor_.get_name() && *(t->at(1).get_Ltensor().get_indices()[extamp_]) == *(i_eri)){
	    ReadRetval(Indents[LoopCount], t->at(1).get_Ltensor(), CPfile);
	    t_list.push_back(t->at(1).get_Ltensor());
            SigmaCount = LoopCount;
	  } // End if
	} // End if
	else{
	  if(isBareLHS_ && t->at(0).get_Ltensor().get_name() == LTensor_.get_name() && t->at(0).get_Ltensor().get_indices()[extamp_] == i_eri){
	    ReadRetval(Indents[LoopCount], t->at(0).get_Ltensor(), CPfile);
	    t_list.push_back(t->at(0).get_Ltensor());
            SigmaCount = 0;
	  } // End if
	} // End else

      } // End if

      // Prepare input map for each declaration method .... 
      vector<string> ExtInd;
      vector<string> NameTen;
      vector<string> Consts;
      //vector<string> Consts(patterns[optimal_num].get_Consts());
      for(size_t num_c = 0;num_c < t->at(0).get_Consts().size();++num_c)
        //if(patterns[optimal_num].get_Consts()[num_c] != "") Consts.push_back(patterns[optimal_num].get_Consts()[num_c]);
        if(t->at(0).get_Consts()[num_c] != "") Consts.push_back(t->at(0).get_Consts()[num_c]);

      contDecl DecFirst;
      for(size_t num_t = 0;num_t < t->at(0).get_Rtensors().size();++num_t){
        SQtensor tt = t->at(0).get_Rtensors()[num_t];
        if(!is_RDM(tt.get_name())        && tt.get_name() != kDelta_name() && tt.get_name() != "Fc1"     &&  
                   !is_C4(tt.get_name()) && !is_C6(tt.get_name())          && tt.get_name() != C2_name() &&
		     find(NameTen.begin(), NameTen.end(), tt.get_name()) == NameTen.end()){
          NameTen.push_back(tt.get_name());
	} // End if
      } // End num_t
      sort(NameTen.begin(), NameTen.end());
      NameTen.push_back(t->at(0).get_Ltensor().get_name());

      vector<SQindex*> temp_boddies(t->at(0).get_summedBody());
      for(vector<SQindex*>::iterator i = temp_boddies.begin(); i != temp_boddies.end();++i)
        if((*i)->get_isExt()) ExtInd.push_back((*i)->get_index());
      sort(ExtInd.begin(), ExtInd.end());

      DecFirst.insert(contDecl::value_type("DecInd",    ExtInd));  // Names of external indices
      DecFirst.insert(contDecl::value_type("DecConst",  Consts));  // All the constants
      DecFirst.insert(contDecl::value_type("DecTensor", NameTen)); // Names of the tensors 

      if(LTensor_.get_name() == t->at(0).get_Ltensor().get_name() && t->size() == 1){
	// PriConnt the calling section of the Fortran body
        makeCPP_header2(t->at(0), title_1if, CHfile, Indents[LoopCount], DecFirst, isBareLHS_);
        // Interface inbetween C++ and F90 codes
        makeF90_interface2(t->at(0), title_1if, F90file, Indents[LoopCount], DecFirst, isBareLHS_);
        // Body of the tensorial contraction
        if     (use_gemm_ && t->at(0).get_Rtensors().size() == 2 && !use_oldgemm_)
	  binary_contract3(t->at(0), title_1, F90file, Indents[LoopCount], DecFirst, isBareLHS_);
        else if(use_gemm_ && t->at(0).get_Rtensors().size() == 2 &&  use_oldgemm_)
	  binary_contract2(t->at(0), title_1, F90file, Indents[LoopCount], DecFirst, isBareLHS_);
        else
	  makeF90_contract2(t->at(0), title_1, F90file, Indents[LoopCount], DecFirst, isBareLHS_);
      } // End if
      else{
	// Print the calling section of the Fortran body
        makeCPP_header2(t->at(0), title_1if, CHfile, Indents[LoopCount], DecFirst, false);
        // Interface inbetween C++ and F90 codes
        makeF90_interface2(t->at(0), title_1if, F90file, Indents[LoopCount], DecFirst, false);
        // Body of the tensorial contraction
        if     (use_gemm_ && t->at(0).get_Rtensors().size() == 2 && !use_oldgemm_)
	  binary_contract3(t->at(0), title_1, F90file, Indents[LoopCount], DecFirst, false);
        else if(use_gemm_ && t->at(0).get_Rtensors().size() == 2 &&  use_oldgemm_)
	  binary_contract2(t->at(0), title_1, F90file, Indents[LoopCount], DecFirst, false);
        else
	  makeF90_contract2(t->at(0), title_1, F90file, Indents[LoopCount], DecFirst, false);
      } // End else

      //vector<SQtensor> t_list;  // List of tensors already read from GA
      vector<SQtensor> EarlyBird; // List of the second tensors already read from GA
      vector<SQtensor> theBody(t->at(0).get_Rtensors());

      // If there's D4 in the RHS of the second contraction, loading or construction of D4 should be appropriately upper place in 
      // the multiple loops. 
      SQtensor *D4_ptr(NULL); 
      vector<SQtensor> theBody2;
      vector<SQindex*> d4_inds;
      if(t->size() == 2) {
	theBody2 = t->at(1).get_Rtensors();
	for(vector<SQtensor>::iterator tt = theBody2.begin();tt != theBody2.end();++tt){
	  if(tt->get_name() == name_d4_){
	    d4_inds.push_back(tt->get_indices()[0]); // Take 0th index
	    d4_inds.push_back(tt->get_indices()[1]); // Take 1st index
	    D4_ptr = &(*tt);
	  } // End if
	} // End tt
      } // End if

      ///////////////////////////////////////////////////////////////////////////////////////
      // Loading, or construction of, 4-RDM is of the top priority
      pair<SQindex*, SQindex*> d4inds_ptr;
      if(D4_ptr != NULL) { 
	d4inds_ptr.first  = d4_inds[0];
	d4inds_ptr.second = d4_inds[1];
      } // End if
      else {
	bool found_d4(false);
	for(vector<SQtensor>::iterator tt = theBody.begin();tt != theBody.end();++tt){
          if(tt->get_name() == name_d4_) { 
            d4inds_ptr.first  = tt->get_indices()[0];
            d4inds_ptr.second = tt->get_indices()[1];
	    found_d4 = true;
	  } // End if 
	} // End tt
	if(!found_d4) { d4inds_ptr.first = NULL; d4inds_ptr.second = NULL; }
      } // End else
      if(d4inds_ptr.first != NULL && d4inds_ptr.second != NULL){
	vector<SQindex*> the_indices; // Loading index (or indices) of 4-RDM in O1
	for(vector<SQindex*>::iterator i = O1.begin();i != O1.end();++i){
	  if((*i)->get_index() == d4inds_ptr.first->get_index() || (*i)->get_index() == d4inds_ptr.second->get_index()) {
	    the_indices.push_back(*i);
	  } // End if
	} // End i
	for(vector<SQindex*>::iterator i = the_indices.begin();i != the_indices.end();++i){
	  vector<SQindex*>::iterator pos(find(O1.begin(), O1.end(), *i));
	  O1.erase(pos);
	} // End i
	if(the_indices.size()) O1.insert(O1.begin(), the_indices.begin(), the_indices.end());
      } // End if
      ///////////////////////////////////////////////////////////////////////////////////////

      for(vector<SQindex*>::iterator i = O1.begin();i != O1.end();++i){
        //*-- --*// string Label;
        //*-- --*// if     ((*i)->get_char() == (char_state)0) Label = "{core}";
        //*-- --*// else if((*i)->get_char() == (char_state)1) Label = "{occ}";
        //*-- --*// else if((*i)->get_char() == (char_state)2) Label = "{vir}";
        //*-- --*// cout << Indents[LoopCount] << "for " << (*i)->get_index() << " in " + Label << ":" << endl;
        CLoop(Indents[LoopCount], **i, CPfile);
        ++LoopCount;
        Declared.push_back(*i);
	if(noIntflags.first == 0 && noIntflags.second) { Cparallel(Indents[LoopCount], **i, CPfile); noIntflags.first = LoopCount; }
        for(vector<SQtensor>::iterator tt = theBody.begin();tt != theBody.end();++tt){
          // * In case of 4-RDM (in first contraction)
	  if(tt->get_name() == name_d4_ && find(Declared.begin(), Declared.end(), tt->get_indices()[0])!=Declared.end() 
	                                && find(Declared.begin(), Declared.end(), tt->get_indices()[1])!=Declared.end()){
	    //cout << Indents[LoopCount] << "Read " + t->get_name() + " from GA for " + t->get_indices()[0]->get_index() + ", " + t->get_indices()[1]->get_index() << endl;
	    if(!use_cumulant_)
	      ReadD4(Indents[LoopCount], *tt, CPfile);
	    else
	      ReadD4_Cumulant(Indents[LoopCount], *tt, CPfile);
	    D4Count = LoopCount;
	    t_list.push_back(*tt);
	  } // End if	
          // *In case of BareAmp
          else if(tt->get_name() == name_amp_ && (**i) == *(tt->get_indices()[extamp_])){
            //cout << Indents[LoopCount] << "Read " + t->get_name() + " from GA for " << (*i)->get_index() << endl;
            ReadAmp(Indents[LoopCount], *tt, CPfile);
            t_list.push_back(*tt);            
	  } // End if
	} // End tt
        // * If all the two contractions are of type 1
        if(isBareLHS_ && t->size() == 2){
          if((**i) == *(t->at(1).get_Ltensor().get_indices()[extamp_])){
            //cout << Indents[LoopCount] << "Read " + LTensor_.get_name() + " from GA for " << (*i)->get_index() << endl;
            ReadRetval(Indents[LoopCount], t->at(1).get_Ltensor(), CPfile);
            SigmaCount = LoopCount;
	  } // End if
	} // End if
        // * In case of this term is composed of only one contractions
	else if(isBareLHS_ && t->at(0).get_Ltensor().get_name() == LTensor_.get_name()){
	  if((**i) == *(t->at(0).get_Ltensor().get_indices()[extamp_])){
            //cout << Indents[LoopCount] << "Read " + LTensor_.get_name() + " from GA for " << (*i)->get_index() << endl;
            ReadRetval(Indents[LoopCount], t->at(0).get_Ltensor(), CPfile); //CPfile << "<-- Koko!!!!" << endl;//*TEST* 
            SigmaCount = LoopCount;
	  } // End if
	} // End if
	// * In case of 4-RDM (in second contraction)
	if(D4_ptr != NULL){
          pair<bool, bool> check;
	  check.first = false; check.second = false;
          for(vector<SQindex*>::iterator ii = Declared.begin();ii != Declared.end();++ii){
	    if(d4_inds[0]->get_index() == (*ii)->get_index()) check.first  = true;
	    if(d4_inds[1]->get_index() == (*ii)->get_index()) check.second = true;
	  } // End if
	  if(check.first && check.second && find(EarlyBird.begin(), EarlyBird.end(), *D4_ptr) == EarlyBird.end()) {
	    if(!use_cumulant_)
	      ReadD4(Indents[LoopCount], *D4_ptr, CPfile);
	    else
	      ReadD4_Cumulant(Indents[LoopCount], *D4_ptr, CPfile);
	    EarlyBird.push_back(*D4_ptr);
	  } // End if
	} // End if
      } // End i

      //*OBSOLETE?* // In case of 4-RDM
      //*OBSOLETE?* for(vector<SQtensor>::iterator tt = theBody.begin();tt != theBody.end();++tt){
      //*OBSOLETE?*   if(tt->get_name() == name_d4_ && find(O1.begin(), O1.end(), tt->get_indices()[0])!=O1.end() 
      //*OBSOLETE?* 	                              && find(O1.begin(), O1.end(), tt->get_indices()[1])!=O1.end()){
      //*OBSOLETE?*     //cout << Indents[LoopCount] << "Read " + t->get_name() + " from GA for " + t->get_indices()[0]->get_index() + ", " + t->get_indices()[1]->get_index() << endl;
      //*OBSOLETE?*     if(!use_cumulant_)
      //*OBSOLETE?*     ReadD4(Indents[LoopCount], *tt, CPfile);
      //*OBSOLETE?*     else
      //*OBSOLETE?*     ReadD4_Cumulant(Indents[LoopCount], *tt, CPfile);
      //*OBSOLETE?*     D4Count = LoopCount;
      //*OBSOLETE?*     t_list.push_back(*tt);
      //*OBSOLETE?* 	} // End if	
      //*OBSOLETE?* } // End t

      // Declaration of the intermediate (in case of the inTerm is composed of three tensors and only if Lindices == false)
      if(is_Interm(t->at(0).get_Ltensor().get_name()) && !(t->at(0).get_Lindices())){
        int Ccount(0);
        int Ocount(0);
        int Vcount(0);
        for(size_t num_i = 0;num_i < t->at(0).get_Ltensor().get_indices().size();++num_i){
	  SQindex* i = t->at(0).get_Ltensor().get_indices()[num_i];
          if     (!(i->get_isExt()) && i->get_char() == (char_state)0) ++Ccount;
          else if(!(i->get_isExt()) && i->get_char() == (char_state)1) ++Ocount;
          else if(!(i->get_isExt()) && i->get_char() == (char_state)2) ++Vcount;
        } // End num_i
        
        // Estimate the number of indices summed at the inner level
        int num_intind(0);
        //if(!X_indices_.size()){
        if(!t->at(0).get_Lindices()){
          for(size_t i = 0;i < t->at(0).get_Ltensor().get_indices().size();++i)
            if(!t->at(0).get_Ltensor().get_indices()[i]->get_isExt()) ++num_intind;
        } // End if
        else{
          num_intind = t->at(0).get_Ltensor().get_indices().size();
        } // End if

        if(num_intind){
          vector<SQindex*> extI;
          string SymmLabel;
          for(size_t num_i = 0;num_i < t->at(0).get_Ltensor().get_indices().size();++num_i){
            if(t->at(0).get_Ltensor().get_indices()[num_i]->get_isExt())
              extI.push_back(t->at(0).get_Ltensor().get_indices()[num_i]);
	  } // End num_i
          if(X_indices_.size()) extI.clear();
          if(!extI.size()) SymmLabel = "0";
          else{
            for(vector<SQindex*>::iterator i = extI.begin();i != extI.end();++i){
              SymmLabel += "s" + (*i)->get_index();
              if(*i != extI.back()) SymmLabel += "^";
	    } // End i
	  } // End else
          string DecSymm("  " + Indents[LoopCount] + "orz::DTensor " + t->at(0).get_Ltensor().get_name());
          string Xname;
          for(int c = 0;c < Ccount;++c) Xname += "c";
          for(int o = 0;o < Ocount;++o) Xname += "a";
          for(int v = 0;v < Vcount;++v) Xname += "v";
        
	  DecSymm += Xname + "_" + title_ + "(orz::mr::sizeof_sympack_X" + Xname + "(" + "symblockinfo, " + SymmLabel + "));";
          CPfile << DecSymm << endl;
        
        } // End if
        //else if(!(LTensor_.get_name() == intermediates[optimal_num].get_name())){
        else if(is_Interm(t->at(0).get_Ltensor().get_name())){
          string DecX("  " + Indents[LoopCount] + "double " + t->at(0).get_Ltensor().get_name() + "_" + title_ + " = 0;"); 
          CPfile << DecX << endl;
        } // End else

      } // End if

      // Declare the rest external indices
      vector<SQindex*> O2;
      for(vector<SQtensor>::iterator tt = theBody.begin();tt != theBody.end();++tt){
        //*-- --*// if(t->get_name() == name_h2_ && find(t_list.begin(), t_list.end(), *t) == t_list.end()){
        //*-- --*//   SQindex* i = t->get_indices()[exth2_];
	//*-- --*// 
        //*-- --*//   if(find(Declared.begin(), Declared.end(), i) == Declared.end()){
        //*-- --*//     O2.push_back(i);
        //*-- --*//     //*-- --*// string Label;
        //*-- --*//     //*-- --*// if     (i->get_char() == (char_state)0) Label = "{core}";  
        //*-- --*//     //*-- --*// else if(i->get_char() == (char_state)1) Label = "{occ}";  
        //*-- --*//     //*-- --*// else if(i->get_char() == (char_state)2) Label = "{vir}";
        //*-- --*//     //*-- --*// printVar += Indents[LoopCount] + "for " + i->get_index() + " in " + Label + ":\n";
        //*-- --*//     //*-- --*// CLoop(Indents[LoopCount], *i, CPfile);
        //*-- --*//     ++LoopCount;
        //*-- --*//     Declared.push_back(i);  
	//*-- --*//   } // End if
	//*-- --*// 
        //*-- --*//   //printVar += Indents[LoopCount] + "Read " + t->get_name() + " from GA for " + i->get_index() + "\n";
        //*-- --*//   ReadERI(Indents[LoopCount], *t, CPfile);
        //*-- --*//   ERICount = LoopCount;
        //*-- --*//   t_list.push_back(*t);
	//*-- --*// } // End if
	// * In case of T2 amplitude
        if(tt->get_name() == name_amp_ && find(t_list.begin(), t_list.end(), *tt) == t_list.end()){
          SQindex* i = tt->get_indices()[extamp_];

          if(find(Declared.begin(), Declared.end(), i) == Declared.end()){
            O2.push_back(i);
            //*-- --*// string Label;
            //*-- --*// if     (i->get_char() == (char_state)0) Label = "{core}";  
            //*-- --*// else if(i->get_char() == (char_state)1) Label = "{occ}";  
            //*-- --*// else if(i->get_char() == (char_state)2) Label = "{vir}";
            //*-- --*// printVar += Indents[LoopCount] + "for " + i->get_index() + " in " + Label + ":\n";
            CLoop(Indents[LoopCount], *i, CPfile);
            ++LoopCount;
	    if(noIntflags.first == 0 && noIntflags.second) { Cparallel(Indents[LoopCount], *i, CPfile); noIntflags.first = LoopCount; }
            Declared.push_back(i);  
	  } // End if

          //printVar += Indents[LoopCount] + "Read " + t->get_name() + " from GA for " + i->get_index() + "\n";
          ReadAmp(Indents[LoopCount], *tt, CPfile);
          t_list.push_back(*tt);
	} // End if
	// * In case of 4-RDM
        else if(tt->get_name() == name_d4_ && find(t_list.begin(), t_list.end(), *tt) == t_list.end()){
          SQindex* i1 = tt->get_indices()[0];
          SQindex* i2 = tt->get_indices()[1];
          vector<SQindex*> i;
          i.push_back(i1);
          i.push_back(i2); 

          for(vector<SQindex*>::iterator j = i.begin();j != i.end();++j){
            O2.push_back(*j);
            //*-- --*// string Label;
            //*-- --*// if     ((*j)->get_char() == (char_state)0) Label = "{core}";  
            //*-- --*// else if((*j)->get_char() == (char_state)1) Label = "{occ}";  
            //*-- --*// else if((*j)->get_char() == (char_state)2) Label = "{vir}";
            if(find(Declared.begin(), Declared.end(), *j) == Declared.end()){
              //printVar += Indents[LoopCount] + "for " + (*j)->get_index() + " in " + Label + ":\n";
              CLoop(Indents[LoopCount], **j, CPfile);
              ++LoopCount;
              Declared.push_back(*j);
            } // End if
            if(find(Declared.begin(), Declared.end(), i1)!=Declared.end() && find(Declared.begin(), Declared.end(), i2)!=Declared.end()
	       && find(t_list.begin(), t_list.end(), *tt) == t_list.end()){

              //printVar += Indents[LoopCount] + "Read " + t->get_name() + " from GA for " + i1->get_index() + "," + i2->get_index() + "\n";
              if(!use_cumulant_)
              ReadD4(Indents[LoopCount], *tt, CPfile);
              else
              ReadD4_Cumulant(Indents[LoopCount], *tt, CPfile);
              D4Count = LoopCount;
              t_list.push_back(*tt);

	    } // End if 
	  } // End j

	} // End if
      } // End tt

      // Print the calling section of the first contraction
      if(LTensor_.get_name() == t->at(0).get_Ltensor().get_name() && t->size() == 1){

        // Print the calling section of the Fortran body
        //makeCPP_body2(t->at(0), title_1if, CPfile, Indents[LoopCount], DecFirst, isBareLHS_); //MOD 20130830
        makeCPP_bodyType2(t->at(0), title_1if, CPfile, Indents[LoopCount], DecFirst, isBareLHS_);

        // Print D4Count post-processing
        int count(LoopCount);
        while(count > 0){
          if(isBareLHS_ && count == SigmaCount){
            SigmaCount = -1;
            AccAmp(Indents[count], t->at(0).get_Ltensor(), CPfile);
	  } // End if
          if(D4Count && count == D4Count){
            D4Count = 0;
            CPfile << "  " + Indents[count] + "FC_FUNC(g_if_unset_d4,G_IF_UNSET_D4)();" << endl;
	  } // End if
	  if(noIntflags.first != 0 && noIntflags.first == count && noIntflags.second) { CPfile << "  " + Indents[count-1] + "} // End para\n"; noIntflags.first = 0; }
          //*-- --*// if(ERICount && count == ERICount){
          //*-- --*//   ERICount = 0;
          //*-- --*//   CPfile << Indents[count+1] + "}// End myrank" << endl;
	  //*-- --*// } // End if        
          LoopEnd(Indents[count-1], CPfile);
          --count;
	} // End while
        if(SigmaCount == 0)
	  AccAmp(Indents[SigmaCount], t->at(0).get_Ltensor(), CPfile);
	//continue;
        //*-- --*// if(rank_guard){
        //*-- --*//   CPfile << "  } // End if" << endl;          
	//*-- --*// } // End if
	//*-- --*// if(!LTensor_.get_indices().size()){
	//*-- --*//   CPfile << "  double sum_" +  LTensor_.get_name() << " = 0.0;" << endl;
	//*-- --*//   CPfile << "  boost::mpi::reduce(orz::world()," + LTensor_.get_name() + ", sum_" + LTensor_.get_name() + ", std::plus<double>(), 0);" << endl;
	//*-- --*//   CPfile << "  " << LTensor_.get_name()  << " = sum_" + LTensor_.get_name() << ";" << endl;        
	//*-- --*// } // End if
	//*-- --*// 
        //*-- --*// if(do_timing_){
        //*-- --*//   CPfile << "  my_timer.push_back(boost::make_tuple(" << "\"" << title_ << "_no" << numTerm << "\", time.elapsed_cputime(), time.elapsed_wallclocktime(), flops/time.elapsed_cputime()));" << endl;
        //*-- --*// } // End if
        //*-- --*// CPfile << "  }" << endl;
	//*-- --*// CPfile << "  orz::world().barrier();" << endl;
        //*-- --*// X_indices_.clear(); // That's all right
        //*-- --*// continue;

      } // End if
      else{
       // Print the calling section of the Fortran body
	//makeCPP_body2(t->at(0), title_1if, CPfile, Indents[LoopCount], DecFirst, false); //MOD 20130830
	makeCPP_bodyType2(t->at(0), title_1if, CPfile, Indents[LoopCount], DecFirst, false);

	// Print D4Count post-processing
	int count(LoopCount);
	while(count > O1.size()){
	  if(isBareLHS_ && count == SigmaCount){
	    SigmaCount = 0;
	    AccAmp(Indents[count], t->at(0).get_Ltensor(), CPfile);
	  } // End if
	  if(noIntflags.first != 0 && noIntflags.first == count && noIntflags.second) { CPfile << "  " + Indents[count-1] + "} // End para\n"; noIntflags.first = 0; }
	  //*-- --*// if(ERICount && count == ERICount){
	  //*-- --*//   ERICount = 0;
	  //*-- --*//   CPfile << Indents[count+1] + "}// End my_imo" << endl;
	  //*-- --*// } // End if        
	  if(D4Count && count == D4Count){
	    D4Count = 0;
	    CPfile << "  " + Indents[count] + "FC_FUNC(g_if_unset_d4,G_IF_UNSET_D4)();" << endl;
	  } // End if        
	  
	  LoopEnd(Indents[count-1], CPfile);
	  --count;
	} // End while

      } // End else

      // Process the second element of workObjs1 if are there
      if(t->size() == 2){
	string title_2if("g_if_" + title_ + "_no1_x" + stm.str() + "_type1_" + TheLabel);
	string title_2  ("g_"    + title_ + "_no1_x" + stm.str() + "_type1_" + TheLabel);
	
	transform(title_2if.begin(), title_2if.end(), title_2if.begin(), (int(*)(int))tolower); 
	transform(title_2  .begin(), title_2  .end(), title_2  .begin(), (int(*)(int))tolower); 

	vector<SQindex*> theInd(t->at(1).get_summedBody());
	contDecl DecSecond;

	// Prepare the input for the code generation
	ExtInd.clear();
	NameTen.clear();
	Consts.clear();
	//Consts = thisTerm->get_Consts();
	for(size_t num_c = 0;num_c < t->at(1).get_Consts().size();++num_c)
	  if(t->at(1).get_Consts()[num_c] != "") Consts.push_back(t->at(1).get_Consts()[num_c]);
	
	for(size_t num_t = 0;num_t < t->at(1).get_Rtensors().size();++num_t){
	  if(!is_RDM(t->at(1).get_Rtensors()[num_t].get_name()) && t->at(1).get_Rtensors()[num_t].get_name() != kDelta_name()
             && t->at(1).get_Rtensors()[num_t].get_name() != "Fc1" && !is_C4(t->at(1).get_Rtensors()[num_t].get_name()) && 
             !is_C6(t->at(1).get_Rtensors()[num_t].get_name()) && t->at(1).get_Rtensors()[num_t].get_name() != C2_name()){
	    NameTen.push_back(t->at(1).get_Rtensors()[num_t].get_name());
	  } // End if
	} // End num_t
	sort(NameTen.begin(), NameTen.end());
	NameTen.push_back(LTensor_.get_name());
	
	for(vector<SQindex*>::iterator i = theInd.begin();i != theInd.end();++i){
	  if((*i)->get_isExt()) ExtInd.push_back((*i)->get_index());
	} // End i
	if(isBareLHS_ && find(ExtInd.begin(), ExtInd.end(), t->at(1).get_Ltensor().get_indices()[extamp_]->get_index()) == ExtInd.end()) 
	  ExtInd.push_back(t->at(1).get_Ltensor().get_indices()[extamp_]->get_index());
	for(size_t num_t = 0;num_t < t->at(1).get_Rtensors().size();++num_t){
	  if     (t->at(1).get_Rtensors()[num_t].get_name() == name_h2_ && 
		  find(ExtInd.begin(), ExtInd.end(),t->at(1).get_Rtensors()[num_t].get_indices()[exth2_]->get_index()) == ExtInd.end())
	    ExtInd.push_back(t->at(1).get_Rtensors()[num_t].get_indices()[exth2_]->get_index());
	  else if(t->at(1).get_Rtensors()[num_t].get_name() == name_amp_ && 
		  find(ExtInd.begin(), ExtInd.end(),t->at(1).get_Rtensors()[num_t].get_indices()[extamp_]->get_index()) == ExtInd.end())
	    ExtInd.push_back(t->at(1).get_Rtensors()[num_t].get_indices()[extamp_]->get_index());
	} // End num_t
	sort(ExtInd.begin(), ExtInd.end());
	
	DecSecond.insert(contDecl::value_type("DecInd",    ExtInd));  // Names of external indices
	DecSecond.insert(contDecl::value_type("DecConst",  Consts));  // All the constants
	DecSecond.insert(contDecl::value_type("DecTensor", NameTen)); // Names of the tensors 
	
	// Print the calling section of the Fortran body
	makeCPP_header2(t->at(1), title_2if, CHfile, Indents[LoopCount], DecSecond, isBareLHS_);
	// Print the calling section of the Fortran body
	makeF90_interface2(t->at(1), title_2if, F90file, Indents[LoopCount], DecSecond, isBareLHS_);
	// Print the body of the tensorial contraction
	if     (use_gemm_ && t->at(1).get_Rtensors().size() == 2 && !use_oldgemm_)
	  binary_contract3(t->at(1), title_2, F90file, Indents[LoopCount], DecSecond, isBareLHS_);
	else if(use_gemm_ && t->at(1).get_Rtensors().size() == 2 && use_oldgemm_)
	  binary_contract2(t->at(1), title_2, F90file, Indents[LoopCount], DecSecond, isBareLHS_);
	else
	  makeF90_contract2(t->at(1), title_2, F90file, Indents[LoopCount], DecSecond, isBareLHS_);
	
	// Convert O1 list to be better suited for thisTerm. Before that make sure, in the previous contraction,
        // index by which the ERI was loaded
        SQindex* i_eri(NULL);
	if(ERIflag || D4Cflag){
	  for(size_t num = 0;num < t->at(0).get_Rtensors().size();++num){
            if(t->at(0).get_Rtensors()[num].get_name() == name_h2_)
              i_eri = t->at(0).get_Rtensors()[num].get_indices()[exth2_];
            else if(is_D4C(t->at(0).get_Rtensors()[num].get_name()))
              i_eri = t->at(0).get_Rtensors()[num].get_indices()[extd4c_];
	  } // End num
          if(i_eri == NULL){
            cout << "Generate_Contract: Loading index of ERI, or D4C isn't found" << endl;
            abort();
	  } // End if
          O1.push_back(i_eri);
	}
	vector<SQindex*> newO1;
	for(vector<SQindex*>::iterator i = O1.begin();i != O1.end();++i){
	  for(vector<SQindex*>::iterator j = theInd.begin();j != theInd.end();++j){
	    if(**i == **j) newO1.push_back(*j);
	  } // End j
	} // End i
	O1 = newO1;
	
	// Starts the second contraction
	Declared.clear();
	Declared.insert(Declared.end(), O1.begin(), O1.end());
	
	LoopCount = (int)O1.size();
        // It's better not to count the loop over loading index of ERI 
        if(ERIflag || D4Cflag/*&& isBareLHS_ /*&& *(t->at(1).get_Ltensor().get_indices()[extamp_]) == *i_eri*/) --LoopCount;
	//printVar = "";
	string Label;

//*OLD* 2013/01/09 	if(isBareLHS_ && find(O1.begin(), O1.end(), t->at(1).get_Ltensor().get_indices()[extamp_]) == O1.end()){
//*OLD* 2013/01/09 	  //for (vector<SQindex*>::iterator i = O1.begin();i != O1.end();++i) cout << *i << endl; //*TEST* 
//*OLD* 2013/01/09 	  //*-- OLD --*// if     (LTensor_.get_indices()[extamp_]->get_char() == (char_state)0) Label = "{core}";
//*OLD* 2013/01/09 	  //*-- OLD --*// else if(LTensor_.get_indices()[extamp_]->get_char() == (char_state)1) Label = "{occ}";
//*OLD* 2013/01/09 	  //*-- OLD --*// else if(LTensor_.get_indices()[extamp_]->get_char() == (char_state)2) Label = "{vir}";
//*OLD* 2013/01/09 	  //*-- OLD --*// printVar += Indents[LoopCount] + "for " + LTensor_.get_indices()[extamp_]->get_index() + " in " + Label + ":\n";
//*OLD* 2013/01/09 	  CLoop(Indents[LoopCount], *(t->at(1).get_Ltensor().get_indices()[extamp_]), CPfile);
//*OLD* 2013/01/09 	  ++LoopCount;
//*OLD* 2013/01/09 	  Declared.push_back(t->at(1).get_Ltensor().get_indices()[extamp_]);
//*OLD* 2013/01/09 	  
//*OLD* 2013/01/09 	  //printVar += Indents[LoopCount] + "Read " + LTensor_.get_name() + " from GA for " + LTensor_.get_indices()[extamp_]->get_index();
//*OLD* 2013/01/09 	  ReadRetval(Indents[LoopCount], t->at(1).get_Ltensor(), CPfile);
//*OLD* 2013/01/09 	  //cout << printVar << endl;
//*OLD* 2013/01/09 	  SigmaCount = LoopCount;
//*OLD* 2013/01/09 	} // End if 
	
	vector<SQindex*> O4;
	t_list.clear();
	t_list.insert(t_list.end(), EarlyBird.begin(), EarlyBird.end()); // <-- NEW 2013/01/08
	//printVar = "";
	//for(vector<SQtensor>::iterator tt = t_list.begin();tt != t_list.end();++tt) CPfile << "t_list : " << *tt << endl;//*TEST*
	vector<SQtensor> tensors = t->at(1).get_Rtensors();

	///////////////////////////////////////////////////////////////////////////////////////
	// Loading of 4-RDM is of the top priority
	vector<SQtensor> vec_d4s;
	for(vector<SQtensor>::iterator tt = tensors.begin();tt != tensors.end();++tt){
          if(tt->get_name() == name_d4_) { vec_d4s.push_back(*tt); tensors.erase(tt); break; }
	} // End tt
	if(vec_d4s.size()) tensors.insert(tensors.begin(), vec_d4s.begin(), vec_d4s.end());
	///////////////////////////////////////////////////////////////////////////////////////

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// NEW 2013/01/09 (added instead of the *OLD* bunch of code seen above)
	bool ReadLT(false);
	if(isBareLHS_ && find(O1.begin(), O1.end(), t->at(1).get_Ltensor().get_indices()[extamp_]) == O1.end()) {
	  tensors.push_back(t->at(1).get_Ltensor());
	  ReadLT = true;
	} // End if
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

	for(vector<SQtensor>::iterator tt = tensors.begin();tt != tensors.end();++tt){
	  if     (tt->get_name() == name_h2_)  O4.push_back(tt->get_indices()[exth2_]);
	  else if(tt->get_name() == name_amp_) O4.push_back(tt->get_indices()[extamp_]);
	  else if(tt->get_name() == name_d4_){
	    O4.push_back(tt->get_indices()[0]);
	    O4.push_back(tt->get_indices()[1]);
	  } // End if
	  else if(tt->get_name() == LTensor_.get_name() && ReadLT && isBareLHS_) O4.push_back(tt->get_indices()[extamp_]);
	} // End tt
	
	for(vector<SQindex*>::iterator i = O4.begin();i != O4.end();++i){
	  //*-- OLD --*// if     ((*i)->get_char() == (char_state)0) Label = "{core}";
	  //*-- OLD --*// else if((*i)->get_char() == (char_state)1) Label = "{occ}";
	  //*-- OLD --*// else if((*i)->get_char() == (char_state)2) Label = "{vir}";
	  
	  if(find(Declared.begin(), Declared.end(), *i) == Declared.end()){
	    //cout << Indents[LoopCount] + "for " + (*i)->get_index() + " in " + Label << endl;
	    CLoop(Indents[LoopCount], **i, CPfile);
	    ++LoopCount;
	    Declared.push_back(*i);
	  } // End if

	  for(vector<SQtensor>::iterator tt = tensors.begin();tt != tensors.end();++tt){
	    //*-- OLD --*// if(t->get_name() == name_h2_ && find(t_list.begin(), t_list.end(), *t) == t_list.end()){
	    //*-- OLD --*//   printVar += Indents[LoopCount] + "Read " + t->get_name() + " from GA for " + t->get_indices()[exth2_]->get_index();
	    //*-- OLD --*//   ReadERI(Indents[LoopCount], *t, CPfile);
	    //*-- OLD --*//   ERICount = LoopCount;
	    //*-- OLD --*//   t_list.push_back(*t);
	    //*-- OLD --*//   cout << printVar << endl;
	    //*-- OLD --*// } // End if

	    //*TEST* ///////////////////////// *TEST* ////////////////////////////
	    //*TEST* CPfile << "khkhkhkhk " << EarlyBird.size(); //*TEST*
	    //*TEST* if(EarlyBird.size()) CPfile << " " << EarlyBird[0] << endl;
	    //*TEST* else                 CPfile << endl; 
	    //*TEST* ///////////////////////// *TEST* ////////////////////////////

	    // *In case of T2 amplitude
	    if(tt->get_name() == name_amp_ && find(t_list.begin(), t_list.end(), *tt) == t_list.end() && 
	       find(Declared.begin(), Declared.end(), tt->get_indices()[extamp_]) != Declared.end()){
	      //printVar += Indents[LoopCount] + "Read " + t->get_name() + " from GA for " + t->get_indices()[extamp_]->get_index();
	      ReadAmp(Indents[LoopCount], *tt, CPfile);
	      t_list.push_back(*tt);
	      //cout << printVar << endl;
	    } // End if
	    // *In case of 4-RDM
	    else if(tt->get_name() == name_d4_ && find(t_list.begin(), t_list.end(), *tt) == t_list.end() 
		    && find(Declared.begin(), Declared.end(), tt->get_indices()[0]) != Declared.end()
		    && find(Declared.begin(), Declared.end(), tt->get_indices()[1]) != Declared.end()){
	      //printVar += Indents[LoopCount] + "Read " + t->get_name() + " from GA for " + t->get_indices()[0]->get_index() + "," + t->get_indices()[1]->get_index();
	      if(!use_cumulant_)
		ReadD4(Indents[LoopCount], *tt, CPfile);
	      else
		ReadD4_Cumulant(Indents[LoopCount], *tt, CPfile);
	      D4Count = LoopCount;
	      t_list.push_back(*tt);
	      //cout << printVar << endl;
	    } // End if
	    // *In case of LHS amplitude
	    if(tt->get_name() == LTensor_.get_name() && ReadLT && isBareLHS_ && find(t_list.begin(), t_list.end(), *tt) == t_list.end() && 
	       find(Declared.begin(), Declared.end(), tt->get_indices()[extamp_]) != Declared.end()){
	      //printVar += Indents[LoopCount] + "Read " + t->get_name() + " from GA for " + t->get_indices()[extamp_]->get_index();
	      ReadRetval(Indents[LoopCount], *tt, CPfile);
	      SigmaCount = LoopCount;
	      t_list.push_back(*tt);
	      //cout << printVar << endl;
	    } // End if
	    //printVar = "";
	  } // End tt
	  
	} // End i
	  
	//*-- OLD --*// vector<SQindex*> Undeclared;
	//*-- OLD --*// vector<SQindex*> LTind(t->at(1).get_Ltensor().get_indices());
	//*-- OLD --*// for(vector<SQindex*>::iterator i = theInd.begin();i != theInd.end();++i){
	//*-- OLD --*//   if(find(Declared.begin(), Declared.end(), *i) == Declared.end() 
	//*-- OLD --*//      && find(LTind.begin(), LTind.end(), *i) == LTind.end()) Undeclared.push_back(*i);
	//*-- OLD --*// } // End i
		
	// Call the Fortran function
	//makeCPP_body2(t->at(1), title_2if, CPfile, Indents[LoopCount], DecSecond, isBareLHS_); //MOD 20130830
	makeCPP_bodyType2(t->at(1), title_2if, CPfile, Indents[LoopCount], DecSecond, isBareLHS_);
		
	int count(LoopCount);
	while(count > 0){
	  if(isBareLHS_ && count == SigmaCount) AccAmp(Indents[SigmaCount], t->at(1).get_Ltensor(), CPfile);
	  if(noIntflags.first != 0 && noIntflags.first == count && noIntflags.second) { CPfile << "  " + Indents[count] + "} // End para\n"; noIntflags.first = 0; }
	  //if(ERICount && count == ERICount)
	  //CPfile << Indents[count+1] << "}// End my_imo" << endl;
	  if(D4Count && count == D4Count)
	    CPfile << "  " + Indents[count] + "FC_FUNC(g_if_unset_d4,G_IF_UNSET_D4)();" << endl;
	  LoopEnd(Indents[count-1], CPfile);
	  --count;
	}// End while
        if(SigmaCount == 0)
	  AccAmp(Indents[SigmaCount], t->at(1).get_Ltensor(), CPfile);
	  
      } // End if

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
      
      CPfile << "  } // End scope" << endl;
      if(noIntflags.second) 
	CPfile << "  //orz::world().barrier();" << endl;
      else if(!ERIflag && !D4Cflag)
	CPfile << "  } // End myrank" << endl;

      CPfile << endl;
    } // End t

    if(ERIflag || D4Cflag) {

      if(LTensor_.get_name() == D4C_nameL())
	CPfile << "  orz::SaveBin(ctinp.dir()/(format(\"D4C_g[%d]\")%i_eri).str()) << retval; //   cout << \"made ok \" << format(\"D4C_g[%d]\")%i_eri << endl;" << endl << endl;

      CPfile << "//-@type(1).contraction(end)" << endl << endl;
      
      CPfile << "  } // End myrank" << endl;
      LoopEnd(Indents[0], CPfile);
      CPfile << "  orz::world().barrier();" << endl;
      
      //*// string sign;
      //*// if     (nature.first == Femto::core) sign = "c";
      //*// else if(nature.first == Femto::act ) sign = "a";
      //*// else if(nature.first == Femto::virt) sign = "v";
      //*// if (ERIflag) CPfile << "//-@loadERI(" << sign << ",end)" << endl;
      //*// else         CPfile << "//-@loadD4C(" << sign << ",end)" << endl;
    } // End if

//*BUG?*    // Before contract workObjs2, all_reduce the intermediate for that
//*BUG?*    for(vector<vector<SQbinary> >::iterator t = workObjs2.begin();t != workObjs2.end();++t,++numTerm){
//*BUG?*      for(size_t num = 0;num < t->at(0).get_Rtensors().size();++num){
//*BUG?*        if(is_Interm(t->at(0).get_Rtensors()[num].get_name())){
//*BUG?*          vector<SQindex*> ind(t->at(0).get_Rtensors()[num].get_indices());
//*BUG?*          string Xname(t->at(0).get_Rtensors()[num].get_name());
//*BUG?*          int numC(0);
//*BUG?*          int numO(0);
//*BUG?*          int numV(0);
//*BUG?*          if(t->at(0).get_Rindices()[num]){
//*BUG?*            for(vector<SQindex*>::iterator i = ind.begin();i != ind.end();++i){
//*BUG?*              if     ((*i)->get_char() == core) ++numC;
//*BUG?*              else if((*i)->get_char() == act ) ++numO;
//*BUG?*              else if((*i)->get_char() == virt) ++numV;
//*BUG?*	    } // End i
//*BUG?*	  } // End if
//*BUG?*	  else {
//*BUG?*            for(vector<SQindex*>::iterator i = ind.begin();i != ind.end();++i){
//*BUG?*              if     (!(*i)->get_isExt() && (*i)->get_char() == core) ++numC;
//*BUG?*              else if(!(*i)->get_isExt() && (*i)->get_char() == act ) ++numO;
//*BUG?*              else if(!(*i)->get_isExt() && (*i)->get_char() == virt) ++numV;
//*BUG?*	    } // End i
//*BUG?*	  } // End else
//*BUG?*          for(int n = 0;n < numC;++n) Xname += "c";
//*BUG?*          for(int n = 0;n < numO;++n) Xname += "a";
//*BUG?*          for(int n = 0;n < numV;++n) Xname += "v";
//*BUG?*          if(numC || numO || numV){
//*BUG?*            //CPfile << "  orz::DTensor sum_" + Xname + ";" << endl;
//*BUG?*            //CPfile << "  orz::tensor::mpi::broadcast(" + Xname + ", 0);";
//*BUG?*	    CPfile << "  boost::mpi::reduce(orz::world(), " + Xname + ".cptr(), " + Xname + ".size(), std::plus<double>(), 0);" << endl;
//*BUG?*	  } // End if
//*BUG?*          else{
//*BUG?*	    CPfile << "  double sum_" +  Xname + " = 0.0;" << endl;
//*BUG?*	    CPfile << "  boost::mpi::reduce(orz::world()," + Xname + ", sum_" + Xname + ", std::plus<double>(), 0);" << endl;
//*BUG?*	    CPfile << "  " + Xname  + " = sum_" + Xname + ";" << endl;        
//*BUG?*	  } // End else
//*BUG?*	} // End if
//*BUG?*      } // End num
//*BUG?*    } // End t

    // Contract workObjs2 here
    if(workObjs2.size()){
      CPfile << endl;
      CPfile << "//-@type(2).contraction(begin)" << endl;
      //CPfile << "// Label : " << TheLabel << endl;      
      CPfile << "  // -- Title : " << title_ << endl;      
      CPfile << "  //*-- Entering to take the type 2 contractions --*//" << endl;
      //CPfile << "  if(myrank == 0){" << endl;
    } // End if
    numTerm = 0;
    for(vector<vector<SQbinary> >::iterator t = workObjs2.begin();t != workObjs2.end();++t,++numTerm){
      size_t LoopCount  =  0;
      size_t SigmaCount = -1;
      size_t D4Count    = -1;
      vector<SQindex*> Declared;
      string indent("");
      bool is_Reduced(false); // The intermediate tensor is reduced to 0-th node, or not
      vector<string> Indents; Indents.reserve(100);
      for(int i = 0;i < 100;++i){
        Indents.push_back(indent);
        indent += "  ";
      } // End i 
      vector<SQindex*> thisBoddies(t->at(0).get_summedBody());

      ostringstream stm;
      stm << numTerm;
      string title_1if("g_if_" + title_ + "_no0_x" + stm.str() + "_type2_" + TheLabel);
      string title_1  ("g_"    + title_ + "_no0_x" + stm.str() + "_type2_" + TheLabel);
      
      transform(title_1if.begin(), title_1if.end(), title_1if.begin(), (int(*)(int))tolower); 
      transform(title_1  .begin(), title_1  .end(), title_1  .begin(), (int(*)(int))tolower); 

      // Beginning of the Scope
      CPfile << "  { " << endl;
      CPfile << "  // No. " << stm.str() << ", [" << t->size() << "]" << endl;
      CPfile << "  // "     << t->at(0)  << endl;
      CPfile << "  double flops = 0; // Flop count" << endl;

#ifdef _VERBOSE_MODE
      CPfile << "#ifdef _VERBOSE" << endl;
#ifndef _CHRONO 
      CPfile << "  //time_t t_start = time(NULL); // VERBOSE MODE" << endl;
      CPfile << "  //double t_start(orz::wallclocktime()); // VERBOSE MODE" << endl;
      CPfile << "  double t_start(ticktack()); // VERBOSE MODE" << endl;
#else
      CPfile << "  boost::chrono::high_resolution_clock::time_point t_start(boost::chrono::high_resolution_clock::now()); // VERBOSE MODE" << endl;
#endif
      CPfile << "#endif" << endl;
#endif

      // Search if the intermediate is of rank 0. If so, this should be reduced
      for(size_t num = 0;num < t->at(0).get_Rtensors().size();++num){
        if(is_Interm(t->at(0).get_Rtensors()[num].get_name()) && !(t->at(0).get_Rtensors()[num].get_indices().size())){
          string name(t->at(0).get_Rtensors()[num].get_name());
          CPfile << "  double sum_" + name << "(0.0);" << endl;
	  CPfile << "  boost::mpi::reduce(orz::world()," + name + "_" + title_ + ", sum_" + name + ", std::plus<double>(), 0);" << endl;
	  CPfile << "  " << name + "_" + title_ << " = sum_" + name << ";" << endl;
          CPfile << "  if(myrank == 0){" << endl;        
	  is_Reduced = true;
	} // End if
        if(is_Reduced) break;
      } // End num

      vector<SQindex*> O1(thisBoddies);
      for(vector<SQindex*>::iterator i = O1.begin();i != O1.end();)
        if(!((*i)->get_isExt())) i = O1.erase(i);
        else ++i;

      // Prepare input map for each declaration method .... 
      vector<string> ExtInd;
      vector<string> NameTen;
      vector<string> Consts;
      for(size_t num_c = 0;num_c < t->at(0).get_Consts().size();++num_c)
        if(t->at(0).get_Consts()[num_c] != "") Consts.push_back(t->at(0).get_Consts()[num_c]);

      contDecl DecFirst;
      for(size_t num_t = 0;num_t < t->at(0).get_Rtensors().size();++num_t){
        SQtensor tt(t->at(0).get_Rtensors()[num_t]);
        if(!is_RDM(tt.get_name()) && tt.get_name() != kDelta_name() && tt.get_name() != "Fc1" && tt.get_name() != C2_name() &&
		     find(NameTen.begin(), NameTen.end(), tt.get_name()) == NameTen.end()){
          NameTen.push_back(tt.get_name());
	} // End if
      } // End num_t
      sort(NameTen.begin(), NameTen.end());
      NameTen.push_back(t->at(0).get_Ltensor().get_name());

      for(vector<SQindex*>::iterator i = thisBoddies.begin(); i != thisBoddies.end();++i)
        if((*i)->get_isExt()) ExtInd.push_back((*i)->get_index());
      sort(ExtInd.begin(), ExtInd.end());

      DecFirst.insert(contDecl::value_type("DecInd",    ExtInd));  // Names of external indices
      DecFirst.insert(contDecl::value_type("DecConst",  Consts));  // All the constants
      DecFirst.insert(contDecl::value_type("DecTensor", NameTen)); // Names of the tensors 

      // Print the calling section of the Fortran body
      makeCPP_header2(t->at(0), title_1if, CHfile, Indents[LoopCount], DecFirst, isBareLHS_);
      // Interface inbetween C++ and F90 codes
      makeF90_interface2(t->at(0), title_1if, F90file, Indents[LoopCount], DecFirst, isBareLHS_);
      // Body of the tensorial contraction
      if     (use_gemm_ && t->at(0).get_Rtensors().size() == 2 && !use_oldgemm_)
	binary_contract3(t->at(0), title_1, F90file, Indents[LoopCount], DecFirst, isBareLHS_);
      else if(use_gemm_ && t->at(0).get_Rtensors().size() == 2 &&  use_oldgemm_)
	binary_contract2(t->at(0), title_1, F90file, Indents[LoopCount], DecFirst, isBareLHS_);
      else
	makeF90_contract2(t->at(0), title_1, F90file, Indents[LoopCount], DecFirst, isBareLHS_);

      vector<SQtensor> theBody(t->at(0).get_Rtensors());
      vector<SQtensor> t_list; // List of tensors already read from GA
      for(vector<SQindex*>::iterator i = O1.begin();i != O1.end();++i){
        //*-- --*// string Label;
        //*-- --*// if     ((*i)->get_char() == (char_state)0) Label = "{core}";
        //*-- --*// else if((*i)->get_char() == (char_state)1) Label = "{occ}";
        //*-- --*// else if((*i)->get_char() == (char_state)2) Label = "{vir}";
        //*-- --*// cout << Indents[LoopCount] << "for " << (*i)->get_index() << " in " + Label << ":" << endl;
        CLoop(Indents[LoopCount], **i, CPfile);
        ++LoopCount;
        Declared.push_back(*i);
        if(isBareLHS_)
          if((**i) == *(t->at(0).get_Ltensor().get_indices()[extamp_])){
            //*-- --*// cout << Indents[LoopCount] << "Read " + LTensor_.get_name() + " from GA for " << (*i)->get_index() << endl;
            ReadRetval(Indents[LoopCount], t->at(0).get_Ltensor(), CPfile);
            SigmaCount = LoopCount;
	  } // End if
        for(vector<SQtensor>::iterator tt = theBody.begin();tt != theBody.end();++tt){
          //*-- --*// if(tt->get_name() == name_h2_ && (**i) == *(tt->get_indices()[exth2_])){
          //*-- --*//   //*-- --*// cout << Indents[LoopCount] << "Read " + t->get_name() + " from GA for " << (*i)->get_index() << endl;
          //*-- --*//   ReadERI(Indents[LoopCount], *t, CPfile);
          //*-- --*//   ERICount = LoopCount;
          //*-- --*//   t_list.push_back(*t);
	  //*-- --*// } // End if
          if(tt->get_name() == name_amp_ && (**i) == *(tt->get_indices()[extamp_])){
            //*-- --*// cout << Indents[LoopCount] << "Read " + t->get_name() + " from GA for " << (*i)->get_index() << endl;
            ReadAmp(Indents[LoopCount], *tt, CPfile);
            t_list.push_back(*tt);            
	  } // End if
          else if(tt->get_name() == name_d4_ && find(Declared.begin(), Declared.end(), (tt->get_indices()[0]))!=Declared.end() 
		                             && find(Declared.begin(), Declared.end(), (tt->get_indices()[1]))!=Declared.end()
		                             && find(t_list.begin(), t_list.end(), *tt)==t_list.end()){
            //*-- --*// cout << Indents[LoopCount] << "Read " + t->get_name() + " from GA for " + t->get_indices()[0]->get_index() + ", " + t->get_indices()[1]->get_index() << endl;
            if(!use_cumulant_)
            ReadD4(Indents[LoopCount], *tt, CPfile);
            else
            ReadD4_Cumulant(Indents[LoopCount], *tt, CPfile);
            D4Count = LoopCount;
            t_list.push_back(*tt);
          } // End if

	} // End tt
      } // End i

      // Print the calling section of the Fortran body
      //makeCPP_body2(t->at(0), title_1if, CPfile, Indents[LoopCount], DecFirst, isBareLHS_); //MOD 20130830
      makeCPP_bodyType2(t->at(0), title_1if, CPfile, Indents[LoopCount], DecFirst, isBareLHS_);

      // Print D4Count post-processing
      while(LoopCount > 0){
        if(isBareLHS_ && LoopCount == SigmaCount){
          SigmaCount = 0;
          AccAmp(Indents[LoopCount], t->at(0).get_Ltensor(), CPfile);
	} // End if
        if(D4Count && LoopCount == D4Count){
          D4Count = 0;
          CPfile << "  " + Indents[LoopCount] + "FC_FUNC(g_if_unset_d4,G_IF_UNSET_D4)();" << endl;
	} // End if        
        LoopEnd(Indents[--LoopCount], CPfile);
      } // End while
      if(is_Reduced) CPfile << "  } // End if" << endl;

#ifdef _VERBOSE_MODE
      CPfile << "#ifdef _VERBOSE" << endl;
#ifndef _CHRONO
      CPfile << "  //t_elapsed[\"" + title_ + "_workObj2_" + TheLabel + stm.str() + "\"] += difftime(time(NULL), t_start); // VERBOSE MODE" << endl; 
      CPfile << "  //t_elapsed[\"" + title_ + "_workObj2_" + TheLabel + stm.str() + "\"] += orz::wallclocktime() - t_start; // VERBOSE MODE" << endl; 
      CPfile << "  t_elapsed[\"" + title_ + "_workObj2_" + TheLabel + stm.str() + "\"] += ticktack() - t_start; // VERBOSE MODE" << endl; 
#else
      CPfile << "  boost::chrono::high_resolution_clock::time_point t_stop(boost::chrono::high_resolution_clock::now()); // VERBOSE MODE" << endl;
      CPfile << "  boost::chrono::duration<double> elapsed = t_stop - t_start;" << endl;
      CPfile << "  t_elapsed[\"" + title_ + "_workObj1_" + TheLabel + stm.str() + "\"] += elapsed.count(); // VERBOSE MODE" << endl;
#endif 
      CPfile << "  t_flops  [\"" + title_ + "_workObj2_" + TheLabel + stm.str() + "\"] += flops;                         // VERBOSE MODE" << endl; 
      CPfile << "#endif" << endl;
#endif

      CPfile << "  } // End scope" << endl << endl;
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
