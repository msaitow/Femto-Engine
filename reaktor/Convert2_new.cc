//
//  Convert2_new.cc
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

//#define _CREATE_C4 

using namespace std;
using namespace Femto;
using namespace Femto::Core;

namespace Femto { namespace Reaktor {


  // *********************************************************
  //  This function is for the type2 inetrmediate
  // *********************************************************
  void SQreaktor::makeCPP_bodyType2_new(SQbinary &bin, string title, ofstream &f, string s, contDecl c, bool isBareLHS)
  {
    auto ExtInd (c["DecInd"]   );
    auto Consts (c["DecConst"] );
    auto NameTen(c["DecTensor"]);

    string FuncName("  " + s + "FC_FUNC(" + title + ",");
    transform(title.begin(), title.end(), title.begin(), (int(*)(int))toupper);
    FuncName += title + ")\n";
    string FuncBody("    " + s + "(");
    
    int count(0);
    for(auto i = ExtInd.begin();i != ExtInd.end();++i, ++count)
      FuncBody += "s" + *i + ", i" + *i + ", ";
    
    for(auto C = Consts.begin();C != Consts.end();++C, ++count)
      if(!is_Interm(*C)) FuncBody += "&" + *C + ", ";
      else               FuncBody += "&" + *C + "_" + title_ + ", ";
 
    for(auto t = NameTen.begin();t != NameTen.end();++t, ++count){
 
      /*- in case of ERI -*/
      if     (*t == name_h2_)  FuncBody += *t + "_sym.cptr(), ";
      /*- in case of D4C -*/
      else if(is_D4C(*t))      FuncBody += *t + ".cptr(), ";
      /*- in case of T2 amp -*/
      else if(*t == name_amp_) FuncBody += *t + "b.cptr(), ";
      /*- There is an intermedaite type tensor detected -*/
      else if(is_Interm(*t)){
        int Ccount(0);
        int Ocount(0);
        int Vcount(0);
        /*- in case of LHS isn't the intermediate -*/
        if(bin.get_Ltensor().get_name() != *t){
          for(size_t num_t = 0;num_t != bin.get_Rtensors().size();++num_t){
            /*- in case that RHS has intermediate of usual type -*/ 
	    //                      AND
            /*- in case that RHS has intermediate of the special type -*/
            if(bin.get_Rtensors()[num_t].get_name() == *t){
              auto inds(bin.get_RInnerIndices(num_t));
              for(auto i = inds.begin();i != inds.end();++i){
                if     ((*i)->get_char() == Femto::core) ++Ccount;
                else if((*i)->get_char() == Femto::act ) ++Ocount;
                else if((*i)->get_char() == Femto::virt) ++Vcount;
	      } // End i
	    } // End if
	  } // End num_t
	} // End if
        /*- LHS is an intermediate of the special type-*/
        else if(bin.get_Ltensor().get_name() == *t){
          auto inds(bin.get_LInnerIndices());
          for(auto i = inds.begin();i != inds.end();++i){
            if     ((*i)->get_char() == Femto::core) ++Ccount;
            else if((*i)->get_char() == Femto::act ) ++Ocount;
            else if((*i)->get_char() == Femto::virt) ++Vcount;
	  } // End i
	} // End else
        /*- In the other cases (LHS is an intermediate of the usual) -*/        
        else {
	  cout << "makeCPP_bodyType2_new: Can't handle this case >> " << *t << " << " << endl;
	  abort();
	} // End else
        if(Ccount || Ocount || Vcount){
          FuncBody += *t;
          for(int i = 0;i < Ccount;++i) FuncBody += "c";
          for(int i = 0;i < Ocount;++i) FuncBody += "a";
          for(int i = 0;i < Vcount;++i) FuncBody += "v";
	  FuncBody += "_" + title_ + ".cptr(), ";
	} // End if
        else{
          FuncBody += "&" + *t + "_" + title_ + ", ";
	} // End else
      } // End if
      /*- in case of RDM -*/
      else if(is_RDM(*t)) continue;
      /*- in case of one-body integral -*/
      else if(*t == name_h1_) FuncBody += "moint1_sym.cptr(), ";
      /*- Name of LHS is not of intermediate type -*/
      else if(*t == bin.get_Ltensor().get_name() && !is_Interm(bin.get_Ltensor().get_name())){
        /*- LHS is an bareamp type -*/
        if(isBareLHS) FuncBody += *t + "b.cptr(), ";
        /*- LHS is an scalar type -*/
        else if(!(bin.get_Ltensor().get_indices().size())) FuncBody += "&" + *t + ", ";
        /*- Otherwise case, usual DTensor type -*/
	else FuncBody += "retval.cptr(), ";
      } // End if
      /*- in case of core-Fock matrix -*/
      else if(t->size() >= 3 && t->at(0) == 'F' && t->at(1) == 'c') continue;
      /*- in case of T1 amplitude -*/
      else if(*t == T1_name()) FuncBody += *t + ".cptr(), ";
      /*- in case of CAS-Fock matrix -*/
      else if(*t == Fock_name()) FuncBody += "CFock.cptr(), ";
      /*- in case of 2-cumulant -*/
      else if(*t == C2_name()) continue;
      else {
        cout << "makeCPP_bodyType2_new: I cannot handle this .... " << *t << endl;
        abort();
      } // End else
      
    } // End t
    FuncBody += "nir, nsym, psym, &flops);\n";

    f << FuncName;
    f << FuncBody;

  }


  // *********************************************************
  // 
  // *********************************************************
  void SQreaktor::makeCPP_header2_new(SQbinary &bin, string title, ofstream &f, string s, contDecl c, bool isBareLHS)
  {
    auto ExtInd (c["DecInd"]   );
    auto Consts (c["DecConst"] );
    auto NameTen(c["DecTensor"]);

    string FuncName("void FC_FUNC(" + title + ",");
    transform(title.begin(), title.end(), title.begin(), (int(*)(int))toupper);
    FuncName += title + ")\n";
    string FuncBody("  (");

    int num(0);
    for(auto i = ExtInd.begin();i != ExtInd.end();++i, ++num){
      FuncBody += "const FC_INT &s" + *i + ", const FC_INT &i" + *i + ", ";
    } // End i
    if(ExtInd.size()) FuncBody += "\n   ";
    for(auto C = Consts.begin();C != Consts.end();++C, ++num){
      /*if(*C != "")*/ FuncBody += "const double * const " + *C + ", ";
    } // End c 
    if(Consts.size()) FuncBody += "\n   ";
    for(auto t = NameTen.begin();t != NameTen.end();++t, ++num){
      FuncBody += "const double * const " + *t + ", ";
    } // End t
    if(NameTen.size()) FuncBody += "\n   ";
    FuncBody += "const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym, const double * const);\n\n";

    f << FuncName;
    f << FuncBody;

  }


  // *********************************************************
  //  Now, *any* type of LTensor is not handlable
  // *********************************************************
  void SQreaktor::makeF90_interface2_new(SQbinary &bin, string title, ofstream &f, string s, contDecl c, bool isBareLHS)
  {
    auto ExtInd (c["DecInd"]   );
    auto Consts (c["DecConst"] );
    auto NameTen(c["DecTensor"]);

    string printGuard("\n\n");
    printGuard += "!                      >> makeF90_interface2_new <<                     \n";
    printGuard += "! **********************************************************************\n";
    printGuard += "!                                                                       \n";
    printGuard += "! **********************************************************************\n";

    string printName("subroutine " + title + " &\n  (");

    for(auto i = ExtInd.begin();i != ExtInd.end();++i)
      printName += "s" + *i + ", i" + *i + ", ";
    for(auto C = Consts.begin();C != Consts.end();++C)
      printName +=  *C + ", ";    
    for(auto t = NameTen.begin();t != NameTen.end();++t)
      printName += *t + ", ";
    printName += "nir, nsym, psym, flops)\n\n";

    string DecComm8("");
    DecComm8 += "use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6\n";
    DecComm8 += "use f_mr_if_module\n";
    DecComm8 += "implicit none\n";
    DecComm8 += "\n";

    f << printGuard;
    f << printName;
    f << DecComm8;

    // Represents intermediate construction is completely separated from that for the resultant tensor due to the priority of loading ERI from Disk.
    bool separate_Loops(false);
    int numInterm(0); // Number of the intermediate

    if(ExtInd.size()){
      string DecInd("integer, intent(inout) :: ");
      for(size_t num_i = 0;num_i < ExtInd.size();++num_i){
        string i(ExtInd[num_i]);
        DecInd += "s" + i + ", i" + i;
        if(num_i != ExtInd.size()-1) DecInd += ", ";
      } // End num_i
      DecInd += "\n";
      f << DecInd;
    } // End if
    if(Consts.size()){
      string DecConst("real(kind=8), intent(inout) :: ");
      for(size_t num_c = 0;num_c < Consts.size();++num_c){
        string c(Consts[num_c]);
        DecConst += c;
        if(num_c != Consts.size()-1) DecConst += ", ";
      } // End num_c
      DecConst += "\n";
      f << DecConst;
    } // End if
    if(NameTen.size()){
      vector<string> scalaNames;
      bool x_flag(false);
      if(is_Interm(bin.get_Ltensor().get_name())) { x_flag = true; separate_Loops = bin.get_Lindices(); }
      for(size_t num_t = 0;num_t < bin.get_Rtensors().size();++num_t){
        if(is_Interm(bin.get_Rtensors()[num_t].get_name())) { x_flag = true; separate_Loops = bin.get_Rindices()[num_t]; }
      } // End num_t
      if(x_flag && !separate_Loops){
        int num(bin.get_LInnerIndices().size());
        if(!num && is_Interm(bin.get_Ltensor().get_name())) scalaNames.push_back(bin.get_Ltensor().get_name());
        
        for(size_t num_t = 0;num_t < bin.get_Rtensors().size();++num_t){
          int num(bin.get_RInnerIndices(num_t).size());
          if(!num && is_Interm(bin.get_Rtensors()[num_t].get_name())) scalaNames.push_back(bin.get_Rtensors()[num_t].get_name());
        } // End num_t
      } // End if

      string DecTensor("real(kind=8), intent(inout) :: ");
      for(size_t num_t = 0;num_t < NameTen.size();++num_t){
        string t(NameTen[num_t]);
        if(t == bin.get_Ltensor().get_name() && !bin.get_Ltensor().get_indices().size() 
             || find(scalaNames.begin(), scalaNames.end(), t) != scalaNames.end())
          DecTensor += t;
        else
          DecTensor += t + "(*)";
        if(num_t != NameTen.size()-1) DecTensor += ", ";
      } // End num_t
      DecTensor += "\n";
      f << DecTensor;
    } // End if

    string DecSymm("");

    DecSymm += "! Information of the Irreps ....\n" ;
    DecSymm += "integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)\n";  
    DecSymm += "! Flop count\n";
    DecSymm += "real*8, intent(inout) :: flops \n";
    DecSymm += "! Some extra stuff\n";
    DecSymm += "integer :: sleft\n\n";
    f << DecSymm;

    map<string, string> extTen;
    auto tensors(bin.get_Rtensors());
    tensors.push_back(bin.get_Ltensor());
    for(auto t = tensors.begin();t != tensors.end();++t){
      /*- in case of intermediate of the usual type -*/
      if(is_Interm(t->get_name())){
        string SetX("");
        vector<string> exCore, inCore, exAct, inAct, exVir, inVir;
        auto inds(t->get_indices());
	int myNum(t-tensors.begin());
	vector<SQindex*> InnerInds;
	vector<SQindex*> OuterInds;
	if(myNum < 2){
	  InnerInds = bin.get_RInnerIndices((size_t)(t-tensors.begin()));
	  OuterInds = bin.get_ROuterIndices((size_t)(t-tensors.begin()));
	} // End if
	else{
	  InnerInds = bin.get_LInnerIndices();
	  OuterInds = bin.get_LOuterIndices();
	} // End else

	/////////////////////////////////////////////////////////////////////////////////
	// Inner indices 
        for(auto i = InnerInds.begin();i != InnerInds.end();++i){
          if     ((*i)->get_char() == Femto::core) inCore.push_back((*i)->get_index());
          else if((*i)->get_char() == Femto::act ) inAct. push_back((*i)->get_index());
          else if((*i)->get_char() == Femto::virt) inVir. push_back((*i)->get_index());
	} // End i
	/////////////////////////////////////////////////////////////////////////////////
	// Outer indices
        for(auto i = OuterInds.begin();i != OuterInds.end();++i){
          if     ((*i)->get_char() == Femto::core) exCore.push_back((*i)->get_index());
          else if((*i)->get_char() == Femto::act ) exAct. push_back((*i)->get_index());
          else if((*i)->get_char() == Femto::virt) exVir. push_back((*i)->get_index());
	} // End i
	/////////////////////////////////////////////////////////////////////////////////

        vector<string> temp(exCore);
        temp.insert(temp.end(), exAct.begin(), exAct.end());
        temp.insert(temp.end(), exVir.begin(), exVir.end());
	if     (temp.size() == 0) SetX += "sleft = 0";
        else if(temp.size() == 1) SetX += "sleft = s" + temp[0] + "\n";
        else if(temp.size() == 2) SetX += "sleft = IEOR(s" + temp[0] + ",s" + temp[1] + ")\n";
        else if(temp.size() == 3) 
          SetX += "sleft = IEOR(s" + temp[0] + ",IEOR(s" + temp[1] + ",s" + temp[2] + "))\n";
        else if(temp.size() == 4)
          SetX += "sleft = IEOR(IEOR(s" + temp[0] + ",s" + temp[1] + "),IEOR(s" + temp[2] + ",s" + temp[3] + "))";
        else if(temp.size() == 4)
          SetX += "sleft = IEOR(IEOR(s" + temp[0] + ",s" + temp[1] + "),IEOR(s" + temp[2] + ",s" + temp[3] + "))";
        else if(temp.size() == 5)
          SetX += "sleft = IEOR(IEOR(IEOR(s" + temp[0] + ",s" + temp[1] + "),IEOR(s" + temp[2] + ",s" + temp[3] + ")),s" + temp[4] + ")";
        else if(temp.size() == 6)
          SetX += "sleft = IEOR(IEOR(IEOR(s" + temp[0] + ",s" + temp[1] + "),IEOR(s" + temp[2] + ",s" + temp[3] + ")),IEOR(s" + temp[4] + ",s" + temp[5] + "))";
        else{
          cout << "makeF90_interfacce2: Number of external indices exceed the limit" << endl;
          abort();
	}
        SetX += "\n";
        f << SetX;

	ostringstream stm;
	if(numInterm) stm << numInterm;

        string nameX("X");
        if(inCore.size())
          for(size_t i = 0;i < inCore.size();++i) nameX += "c";
        if(inAct.size())
          for(size_t i = 0;i < inAct.size();++i)  nameX += "a";
        if(inVir.size())
          for(size_t i = 0;i < inVir.size();++i)  nameX += "v";
	if(numInterm) nameX += stm.str();
        SetX = "call set_symblock_" + nameX + "(sleft, " + t->get_name() + ", nir, nsym, psym) ! -> " + nameX + " (allocate) \n"; 
        if(inCore.size() || inAct.size() || inVir.size()){
	  f << SetX;
	  extTen.insert(map<string, string>::value_type(t->get_name(), nameX));
	} // End if
	++numInterm;          
      } // End if
      /*- in case of the T2 amp -*/
      else if(t->get_name() == name_amp_){
        string iname = t->get_indices()[extamp_]->get_index();
        string SetT("call set_symblock_av2(s" + iname + ", " + t->get_name() + ", nir, nsym, psym) ! -> av2_i (allocate)\n");
        f << SetT;
        extTen.insert(map<string,string>::value_type(name_amp_, "av2_i"));
      } // End if
      /*- in case of the one-body integral -*/
      else if(t->get_name() == name_h1_){
        string Seth1("call set_symblock_h1(" + t->get_name() + ", nir, nsym, psym) ! -> h1 (allocate)\n");
        f << Seth1;
        extTen.insert(map<string,string>::value_type(t->get_name(), "h1"));        
      } // End if
      /*- in case of the ERI -*/
      else if(t->get_name() == name_h2_){
        string iname(t->get_indices()[exth2_]->get_index());
        string Seth2("call set_symblock_h2(s" + iname + ", " + t->get_name() + ", nir, nsym, psym) ! -> h2_i (allocate)\n");
        f << Seth2;
        extTen.insert(map<string,string>::value_type(name_h2_, "h2_i"));
      } // End if
      /*- in case of T1 amplitude -*/
      else if(t->get_name() == T1_name()){ 
        string Seth1("call set_symblock_t1_1(" + t->get_name() + ", nir, nsym, psym) ! -> t1 (allocate)\n");
        f << Seth1;
        extTen.insert(map<string,string>::value_type(t->get_name(), "t1"));        
      } // End if
      /*- in case of D4C amplitude -*/
      else if(is_D4C(t->get_name())){ 
        string Seth1("call set_symblock_d4c(s"  + t->get_indices()[extd4c_]->get_index() + ", " + t->get_name() + ", nir, nsym, psym) ! -> d4cf (allocate)\n");
        f << Seth1;
        extTen.insert(map<string,string>::value_type(t->get_name(), "d4cf"));        
      } // End if
      /*- in case of CAS-Fock matrix -*/
      else if(t->get_name() == Fock_name()){ 
        string Seth1("call set_symblock_g1(" + t->get_name() + ", nir, nsym, psym) ! -> g1 (allocate)\n");
        f << Seth1;
        extTen.insert(map<string,string>::value_type(t->get_name(), "g1"));        
      } // End if
#ifdef _CREATE_C4
      /*- in case of C4 tensor -*/
      else if(t->get_name() == C4_name()){ 
        string Seth1("call set_symblock_c4f(" + t->get_name() + ", nir, nsym, psym) ! -> c4f (allocate)\n");
        f << Seth1;
        extTen.insert(map<string,string>::value_type(t->get_name(), "c4f"));        
      } // End if
#endif
    } // End t

    /*- if isBareLHS is true, LHS is treated as a bareamppack type symblock -*/ 
    if(isBareLHS){
      string SetLT("call set_symblock_av2_2(s" + bin.get_Ltensor().get_indices()[extamp_]->get_index() + ", " + bin.get_Ltensor().get_name() + ", nir, nsym, psym) ! -> av2_i2 (allocate)\n");
      f << SetLT;
      extTen.insert(map<string,string>::value_type(bin.get_Ltensor().get_name(), "av2_i2"));
    } // End if
    /*- if isBareLHS is false and LHS is S1 type tensor -*/
    else if(bin.get_Ltensor().get_name() == LTensor_.get_name() && bin.get_Ltensor().get_indices().size() && LTensor_.get_indices().size() == 2/*&& !isBareLHS_*/){
      bool isS1(true);
      for(map<string,string>::iterator s = extTen.begin();s != extTen.end();++s) 
        if(s->first == bin.get_Ltensor().get_name()) isS1 = false;
      if(isS1){
	string Seth1("call set_symblock_t1_2(" + bin.get_Ltensor().get_name() + ", nir, nsym, psym) ! -> s1 (allocate)\n");
	f << Seth1;
	extTen.insert(map<string,string>::value_type(bin.get_Ltensor().get_name(), "s1"));        
      } // End if
    } // End if
    /*- in case of LHS is a C5-tensor -*/
    else if(bin.get_Ltensor().get_name() == D4C_nameL()) {
      string Seth1("call set_symblock_d4c(s" + bin.get_Ltensor().get_indices()[5]->get_index() + ", " + bin.get_Ltensor().get_name() + ", nir, nsym, psym) ! -> d4cf (allocate)\n");
      f << Seth1;
      extTen.insert(map<string,string>::value_type(bin.get_Ltensor().get_name(), "d4cf"));        
    } // End if
    /*- in case of LHS is a C6-tensor -*/
    else if(bin.get_Ltensor().get_name() == C6_nameL()){
      string Seth1("call set_symblock_c6f(" + bin.get_Ltensor().get_name() + ", nir, nsym, psym) ! -> c6f (allocate)\n");
      f << Seth1;
      extTen.insert(map<string,string>::value_type(bin.get_Ltensor().get_name(), "c6f"));        
    } // End if
    /*- in case of LHS is a C4-tensor -*/
    else if(bin.get_Ltensor().get_name() == C4_nameL()){
      string Seth1("call set_symblock_c4f(" + bin.get_Ltensor().get_name() + ", nir, nsym, psym) ! -> c4f (allocate)\n");
      f << Seth1;
      extTen.insert(map<string,string>::value_type(bin.get_Ltensor().get_name(), "c4f"));        
    } // End if

    string FuncName(title);
    auto s_end(FuncName.begin());
    ++s_end; ++s_end; ++s_end; ++s_end;
    FuncName.erase(FuncName.begin(), s_end);
    string SetFunc("call g" + FuncName + " &\n  (");

    for(auto i = ExtInd.begin();i != ExtInd.end();++i)
      SetFunc += "s" + *i + ", i" + *i + ", ";
    for(auto C = Consts.begin();C != Consts.end();++C)
      SetFunc += *C + ", ";
    for(auto t = NameTen.begin();t != NameTen.end();++t){
      if(extTen.find(*t) != extTen.end())
        SetFunc += extTen[*t] + ", ";
      else
        SetFunc += *t + ", ";
    } // End t

    // List of names of RDM and other tensors pushed into the module
    vector<string> d_list;
    for(auto t = tensors.begin();t != tensors.end();++t) {
      if     (is_RDM(t->get_name()) || t->get_name() == "Fc1" || t->get_name() == C2_name()) 
                                                               { d_list.push_back(t->get_name()); }
#ifndef _CREATE_C4
      else if(is_C4(t->get_name()))                            { d_list.push_back("c4f");         }
#endif
      else if(is_C6(t->get_name()))                            { d_list.push_back("c6f");         }
    }
    for(auto s = d_list.begin();s != d_list.end();++s)
      if(*s == RDM_name() + "4") SetFunc += "d4_ij, ";
      else {
        transform(s->begin(), s->end(), s->begin(), (int(*)(int))tolower);
        SetFunc += *s + ", ";
      } // End else

    SetFunc += "nir, nsym, psym, flops)\n\n";
    f << SetFunc;

    for(auto k = extTen.begin();k != extTen.end();++k){
      if(is_Interm(k->second)) continue;
      string SetDealloc("deallocate(" + k->second + ")\n");
      f << SetDealloc;
    } // End s

    f << endl;

    f << "end subroutine " + title + "\n\n\n";

  }

  // *********************************************************
  // This method is now outdated. So if you're interested in use this,
  // some modification may be necessary (2012/12/20).
  // *********************************************************
  void SQreaktor::makeF90_contract2_new(SQbinary &bin, string title, ofstream &f, string s, contDecl c, bool isBareLHS)
  {

    if(bin.get_Lindices() || bin.get_Rindices()[0] || bin.get_Rindices()[1]){
      cout << "makeF90_contract: I can handle only the case of on-the-fly algorithm for now." << endl;
      abort();
    } // End if

    vector<string> ExtInd  = c["DecInd"];
    vector<string> Consts  = c["DecConst"];
    vector<string> NameTen = c["DecTensor"];

    f <<  "\n"                                                                              ;     
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
    vector<SQtensor> tensors(bin.get_Rtensors());
    vector<SQindex*> summedBody(bin.get_summedBody());
    vector<SQindex*> LTinds(bin.get_Ltensor().get_indices());
    for(vector<SQindex*>::iterator i = LTinds.begin();i != LTinds.end();++i){
      bool i_flag = true;
      for(vector<SQindex*>::iterator j = summedBody.begin();j != summedBody.end();++j)
        if(**i == **j) i_flag = false;
      if(i_flag) summedBody.push_back(*i);
    } // End i

    for(vector<SQindex*>::iterator i = summedBody.begin();i != summedBody.end();++i){
      if((*i)->get_isExt()) {
        DecExtInd += "integer, intent(in) :: i_" +  (*i)->get_index() + ", s_" + (*i)->get_index() + "\n";
        printName += "i_" + (*i)->get_index() + ", s_" + (*i)->get_index() + ", ";
      } // End if
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
    for(vector<string>::iterator s = Consts.begin();s != Consts.end();++s){
      DecConsts += "real(kind=8), intent(inout) :: " + *s;
      printName += *s + ", ";
    } // End s

    // Print all the names of tensors
    map<string, int> mapTen;
    vector<string> keys;
    if(isBareLHS){
      mapTen.insert( map<string, int>::value_type(bin.get_Ltensor().get_name(), 3) );
      keys.push_back(bin.get_Ltensor().get_name());
    } // End if
    //else if(!isBareLHS && LTensor.get_name() == "X"){
    else if(!isBareLHS && is_Interm(bin.get_Ltensor().get_name())){
      int NDcount = 0;
      vector<SQindex*> inds(bin.get_Ltensor().get_indices());
      for(vector<SQindex*>::iterator i = inds.begin();i != inds.end();++i)
        if(!(*i)->get_isExt()) ++NDcount;
      mapTen.insert( map<string, int>::value_type(bin.get_Ltensor().get_name(), NDcount) );
      keys.push_back(bin.get_Ltensor().get_name());
    } // End if
    else if(bin.get_Ltensor().get_name() == D4C_nameL()){
      mapTen.insert( map<string, int>::value_type(bin.get_Ltensor().get_name(), 5) );
      keys.push_back(bin.get_Ltensor().get_name());
    } // End if
    else{
      mapTen.insert( map<string, int>::value_type(bin.get_Ltensor().get_name(), bin.get_Ltensor().get_indices().size()) );
      keys.push_back(bin.get_Ltensor().get_name());
    } // End else
    for(vector<SQtensor>::iterator t = tensors.begin();t != tensors.end();++t){
      if(find(keys.begin(), keys.end(), t->get_name()) == keys.end()){
        if(t->get_name() == name_amp_ || t->get_name() == name_h2_){
          mapTen.insert( map<string, int>::value_type(t->get_name(), 3) );
          keys.push_back(t->get_name());
	} // End if
        else if(t->get_name() == RDM_name() + "4"){
          mapTen.insert( map<string, int>::value_type(t->get_name(), 6) );
          keys.push_back(t->get_name());      
	} // End if
        //else if(t->get_name() == "X"){
        else if(is_Interm(t->get_name())){
          int NDcount = 0;
          vector<SQindex*> inds(t->get_indices());
          for(vector<SQindex*>::iterator i = inds.begin();i != inds.end();++i)
            if(!(*i)->get_isExt()) ++NDcount;
          mapTen.insert( map<string, int>::value_type(t->get_name(), NDcount) );
          keys.push_back(t->get_name());
	} // End if
        else{
          mapTen.insert( map<string, int>::value_type(t->get_name(), t->get_indices().size()) );
          keys.push_back(t->get_name());      
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
    // This bunch of codes should be same to Convert_gemm2.cc.
    for(vector<string>::iterator t = keys.begin();t != keys.end();++t){
      if(is_RDM(*t) || *t == "Fc1"|| is_C4(*t) || is_C6(*t)) printName2 += *t + "_, ";
      size_t dimTen(mapTen[*t]);
      printName += *t + "_, ";
      if(!dimTen){
        DecTensor += "real(kind=8)                   :: " + *t + "_\n";
        continue;
      } // End if
      else{
        ostringstream stm;
        stm << dimTen;
        fraction = "type(symblock" + stm.str() + "), intent(inout) :: " + *t + "_(";
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

    //f << printData;
    //f << printName;
    f << printName2;
    f << DecComm8;
    f << DecExtInd;
    f << DecSymm;
    f << DecConsts;
    f << DecTensor;

    // Entering the body of the tensorial contraction part
    int LoopCount(0);
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
      if     (t->get_indices().size() == 1){ // In case of unity
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

    // OpenMP parallelization
    if(use_omp_){
      ostringstream stm;
      stm << NDindices.size();
      //f << "!$omp parallel do schedule(dynamic) collapse(" + stm.str() + ")";
      f << "!$omp parallel do";
      if(!bin.get_Ltensor().get_indices().size()) f << " reduction(+:" + bin.get_Ltensor().get_name() + "_)";
      f << endl;
    }
    // Write down the body of the tensorial contractions
    int count = 0;
    string flopcount("flops = flops + ");
    for(vector<SQindex*>::iterator i = NDindices.begin();i != NDindices.end();++i, ++count){
      string iname = "i_" + (*i)->get_index();
      string sname = "s_" + (*i)->get_index();
      string otype;
      if     ((*i)->get_char() == (char_state)0) otype = "I_C";
      else if((*i)->get_char() == (char_state)1) otype = "I_O";
      else if((*i)->get_char() == (char_state)2) otype = "I_V";
      f << "do " + iname + " = psym(I_BEGIN, " + otype + ", " + sname + "), psym(I_END, " + otype + ", " + sname+ ")" << endl;
      if(i != NDindices.begin()) flopcount += " * ";
      flopcount += "psym(I_LENGTH, " + otype + ", " + sname + ")";
      if(i != NDindices.end()-1) flopcount += " &";
      else                       flopcount += " * 2.0d+00";
      flopcount += "\n";
    } // End i
    if(count != LoopCount){
      cout << "Factorize: Something is wrong in the body of tensorial contraction part .... " << endl;
      abort();
    } // End if

    // Generate the left-hand side tensor
    vector<SQindex*> Lindlist;
    //if     (!isBareLHS && bin.get_Ltensor().get_name() != "X"){
    if     (!isBareLHS && !is_Interm(bin.get_Ltensor().get_name())){
      Lindlist = bin.get_Ltensor().get_indices();
    } // End if
    //else if(!isBareLHS && bin.get_Ltensor().get_name() == "X"){
    else if(!isBareLHS && is_Interm(bin.get_Ltensor().get_name())){
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

    if(Lindlist.size()){
      string nterm(bin.get_Ltensor().get_name() + "_(");
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
      f << nterm + " =  &\n    " + nterm + " &" << endl;
    } // End if
    else {
      f << bin.get_Ltensor().get_name() + "_ = " + bin.get_Ltensor().get_name() + "_ &" << endl;
    } // End else

    // Generate the numerical coefficients
    if(bin.get_numConst() > 0) {
      ostringstream stm;
      stm << (boost::format("%10.8f") % bin.get_numConst());
      f << "  + " + stm.str() + "d+00" + " & " << endl;
    } // end if
    else {
      ostringstream stm;
      stm << (boost::format("%10.8f") % fabs(bin.get_numConst()));
      f << "  - " + stm.str() + "d+00" + " & " << endl;
    } // end if

    // Print all the Coeffient
    for(size_t i = 0;i < bin.get_Consts().size();++i){
      if(bin.get_Consts()[i] != "") f << "  * " + bin.get_Consts()[i];
      if(bin.get_Consts()[i] != "" && tensors.size()) f << " & " << endl;
      else if(bin.get_Consts()[i] != "")              f << endl;
    } // End if

    // Generate all the tensors on the right-hand side
    for(vector<SQtensor>::iterator t = tensors.begin();t != tensors.end();++t){
      vector<SQindex*> Rindlist;
      if(t->get_name() == name_h2_){
        for(size_t num_i = 0;num_i < t->get_indices().size();++num_i){
          SQindex* i(t->get_indices()[num_i]);
          if(num_i != exth2_) Rindlist.push_back(i);
	} // End num_i
      } // End if
      else if(t->get_name() == name_amp_){
        for(size_t num_i = 0;num_i < t->get_indices().size();++num_i){
          SQindex* i(t->get_indices()[num_i]);
          if(num_i != extamp_) Rindlist.push_back(i);
	} // End num_i
      } // End if
      else if(t->get_name() == RDM_name() + "4"){
        int external_count = 0;
        for(size_t num_i = 0;num_i < t->get_indices().size();++num_i){
          SQindex* i(t->get_indices()[num_i]);
          if(num_i != 0 && num_i != 1) Rindlist.push_back(i);
          else ++external_count;
	} // End num_i
        if(external_count != 2) {
          cout << "Factorize: Something is wrong in treatment of 4-RDM .... " << endl;
          cout << "external_count: " << external_count << endl; 
          abort(); 
	} // End if
      } // End if
      //else if(t->get_name() == "X"){
      else if(is_Interm(t->get_name())){
        for(size_t num_i = 0;num_i < t->get_indices().size();++num_i){
          SQindex* i(t->get_indices()[num_i]);
          if(!i->get_isExt()) Rindlist.push_back(i);
	} // End num_i
      } // End if
      else{
        Rindlist = t->get_indices();
      } // End else

      string nterm("  * " + t->get_name() + "_");
      if(Rindlist.size()){
        nterm += "(";
        for(size_t num_i = 0;num_i < Rindlist.size();++num_i){
          nterm += "s_" + Rindlist[Rindlist.size()-(num_i+1)]->get_index();
          if(num_i != Rindlist.size()-1) nterm += ", ";
          else nterm += ")";
        } // End num_i

        nterm += "%array("; 
        for(size_t num_i = 0;num_i < Rindlist.size();++num_i){
          nterm += "i_" + Rindlist[Rindlist.size()-(num_i+1)]->get_index();
          if(num_i != Rindlist.size()-1) nterm += ", ";
          else nterm += ")";
        } // End num_i
      } // End if
      if(t != t_end) nterm += " & ";
      f << nterm << endl;

    } // End t

    f << endl;
    f << "! Flop count " << endl;
    f << flopcount << endl;

    for(int n = 0;n < LoopCount;++n) f << "end do ! Irrep Loop"   << endl;
    //if(use_omp_)                    f << "!$omp end parallel do" << endl;
                                     f << "end if ! Irrep Cond"   << endl;
    for(int n = 0;n < LoopCount;++n) f << "end do ! Orbital Loop" << endl;
    f << "! FEMTO END  ****************************************************************\n" << endl;
    f << printEnd;
  }


}} // Femto::
