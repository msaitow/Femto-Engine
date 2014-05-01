//
//  Convert.cc
//  
//
//  Created by Masaaki Saitow on 12/09/05.
//  Copyright (c) 2012 Masaaki Saitow. All rights reserved.
//

#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <SQreaktor.hpp>
#include <boost/format.hpp>

using namespace std;
using namespace Femto;
using namespace Femto::Core;

namespace Femto { namespace Reaktor {

  // *********************************************************
  // 
  // *********************************************************
  void SQreaktor::makeCPP_body(SQtensor &LTensor, SQterm &inTerm, string title, ofstream &f, string s, contDecl c, 
      bool isBareLHS)
  {
    vector<string> ExtInd  = c["DecInd"];
    vector<string> Consts  = c["DecConst"];
    vector<string> NameTen = c["DecTensor"];

    string FuncName("  " + s + "FC_FUNC(" + title + ", ");
    transform(title.begin(), title.end(), title.begin(), (int(*)(int))toupper);
    FuncName += title + ")\n";
    string FuncBody("    " + s + "(");
    
    int count = 0;
    for(vector<string>::iterator i = ExtInd.begin();i != ExtInd.end();++i, ++count)
      FuncBody += "s" + *i + ", i" + *i + ", ";
    
    for(vector<string>::iterator C = Consts.begin();C != Consts.end();++C, ++count)
      /*if(*C != "")*/ FuncBody += "&" + *C + ", ";
 
    for(vector<string>::iterator t = NameTen.begin();t != NameTen.end();++t, ++count){
      if     (*t == name_h2_)  FuncBody += *t + "_sym.cptr(), ";
      else if(*t == name_amp_) FuncBody += *t + "b.cptr(), ";
      else if(*t == "X"){
        int Ccount = 0;
        int Ocount = 0;
        int Vcount = 0;
        if(LTensor.get_name() != "X"){
          for(size_t num_t = 0;num_t != inTerm.get_tensors().size();++num_t){
            if(inTerm.get_tensors()[num_t].get_name() == "X" && !X_indices_.size()){
              vector<SQindex*> inds(inTerm.get_tensors()[num_t].get_indices());
              for(vector<SQindex*>::iterator i = inds.begin();i != inds.end();++i){
                if     ((*i)->get_char() == (char_state)0 && !(*i)->get_isExt()) ++Ccount;
                else if((*i)->get_char() == (char_state)1 && !(*i)->get_isExt()) ++Ocount;
                else if((*i)->get_char() == (char_state)2 && !(*i)->get_isExt()) ++Vcount;
	      } // End i
	    } // End if
            else if(inTerm.get_tensors()[num_t].get_name() == "X" && X_indices_.size()){
              vector<SQindex*> inds(X_indices_);
              for(vector<SQindex*>::iterator i = inds.begin();i != inds.end();++i){
                if     ((*i)->get_char() == (char_state)0) ++Ccount;
                else if((*i)->get_char() == (char_state)1) ++Ocount;
                else if((*i)->get_char() == (char_state)2) ++Vcount;
	      } // End i
	    } // End if
	  } // End num_t
	} // End if
        else if(LTensor.get_name() == "X" && X_indices_.size()){
          vector<SQindex*> inds(X_indices_);

          for(vector<SQindex*>::iterator i = inds.begin();i != inds.end();++i){
            if     ((*i)->get_char() == (char_state)0) ++Ccount;
            else if((*i)->get_char() == (char_state)1) ++Ocount;
            else if((*i)->get_char() == (char_state)2) ++Vcount;
	  } // End i
	} // End else        
        else {
          vector<SQindex*> inds(LTensor.get_indices());

          for(vector<SQindex*>::iterator i = inds.begin();i != inds.end();++i){
            if     ((*i)->get_char() == (char_state)0 && !(*i)->get_isExt()) ++Ccount;
            else if((*i)->get_char() == (char_state)1 && !(*i)->get_isExt()) ++Ocount;
            else if((*i)->get_char() == (char_state)2 && !(*i)->get_isExt()) ++Vcount;
	  } // End i
	} // End else
        if(Ccount || Ocount || Vcount){
          FuncBody += *t;
          for(int i = 0;i < Ccount;++i) FuncBody += "c";
          for(int i = 0;i < Ocount;++i) FuncBody += "a";
          for(int i = 0;i < Vcount;++i) FuncBody += "v";
	  FuncBody += ".cptr(), ";
	} // End if
        else{
          FuncBody += "&" + *t + ", ";
	} // End else
      } // End if
      else if((*t)[0] == 'Y' && LTensor.get_name().at(0) != 'Y'){
	FuncBody += *t + ".cptr(), ";
      } // End if
      else if(is_RDM(*t)) continue;
      else if(*t == name_h1_) FuncBody += "moint1_sym.cptr(), ";
      else if(*t == LTensor.get_name() && LTensor.get_name() != "X" || LTensor.get_name().at(0) == 'Y'){
        if(isBareLHS) FuncBody += *t + "b.cptr(), ";
        else if(!(LTensor.get_indices().size())) FuncBody += "&" + *t + ", ";
        else FuncBody += *t + ".cptr(), ";
      } // End if
      else {
        cout << "Factorize: I cannot handle this .... " << *t << endl;
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
  void SQreaktor::makeCPP_header(SQtensor &LTensor, SQterm &inTerm, string title, ofstream &f, string s, contDecl c, 
      bool isBareLHS)
  {
    vector<string> ExtInd  = c["DecInd"];
    vector<string> Consts  = c["DecConst"];
    vector<string> NameTen = c["DecTensor"];

    string FuncName("void FC_FUNC(" + title + ", ");
    transform(title.begin(), title.end(), title.begin(), (int(*)(int))toupper);
    FuncName += title + ")\n";
    string FuncBody("  (");

    int num = 0;
    for(vector<string>::iterator i = ExtInd.begin();i != ExtInd.end();++i, ++num){
      FuncBody += "const FC_INT &s" + *i + ", const FC_INT &i" + *i + ", ";
    } // End i
    if(ExtInd.size()) FuncBody += "\n   ";
    for(vector<string>::iterator C = Consts.begin();C != Consts.end();++C, ++num){
      /*if(*C != "")*/ FuncBody += "const double * const " + *C + ", ";
    } // End c 
    if(Consts.size()) FuncBody += "\n   ";
    for(vector<string>::iterator t = NameTen.begin();t != NameTen.end();++t, ++num){
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
  void SQreaktor::makeF90_interface(SQtensor &LTensor, SQterm &inTerm, string title, ofstream &f, string s, contDecl c, 
      bool isBareLHS)
  {
    vector<string> ExtInd  = c["DecInd"];
    vector<string> Consts  = c["DecConst"];
    vector<string> NameTen = c["DecTensor"];

    string printGuard("\n\n");
    printGuard += "! **********************************************************************\n";
    printGuard += "!                                                                       \n";
    printGuard += "! **********************************************************************\n";

    string printName("subroutine " + title + "(");

    for(vector<string>::iterator i = ExtInd.begin();i != ExtInd.end();++i)
      printName += "s" + *i + ", i" + *i + ", ";
    for(vector<string>::iterator C = Consts.begin();C != Consts.end();++C)
      printName +=  *C + ", ";    
    for(vector<string>::iterator t = NameTen.begin();t != NameTen.end();++t)
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
      bool x_flag = false;
      if(LTensor.get_name() == "X") x_flag = true;
      for(size_t num_t = 0;num_t < inTerm.get_tensors().size();++num_t){
        if(inTerm.get_tensors()[num_t].get_name() == "X") x_flag = true;
      } // End num_t
      if(x_flag && !X_indices_.size()){
        int num = 0;
        for(size_t num_i = 0;num_i < LTensor.get_indices().size();++num_i)
          if(!LTensor.get_indices()[num_i]->get_isExt()) ++num;
        if(!num && LTensor.get_name() == "X") scalaNames.push_back(LTensor.get_name());
        
        for(size_t num_t = 0;num_t < inTerm.get_tensors().size();++num_t){
          int num = 0;
          for(size_t num_i = 0;num_i < inTerm.get_tensors()[num_t].get_indices().size();++num_i)
            if(!inTerm.get_tensors()[num_t].get_indices()[num_i]->get_isExt()) ++num; 
          if(!num && inTerm.get_tensors()[num_t].get_name() == "X") 
            scalaNames.push_back(inTerm.get_tensors()[num_t].get_name());
        } // End num_t
      } // End if
      string DecTensor("real(kind=8), intent(inout) :: ");
      for(size_t num_t = 0;num_t < NameTen.size();++num_t){
        string t(NameTen[num_t]);
        if(t == LTensor.get_name() && !LTensor.get_indices().size() 
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
    DecSymm += "integer :: sleft\n";
    DecSymm += "integer :: sleft2\n\n";
    f << DecSymm;

    map<string, string> extTen;
    vector<string> keys;
    vector<SQtensor> tensors(inTerm.get_tensors());
    tensors.push_back(LTensor);
    for(vector<SQtensor>::iterator t = tensors.begin();t != tensors.end();++t){
      if(t->get_name() == "X" && !X_indices_.size()){
        string SetX("");
        vector<string> exCore, inCore, exAct, inAct, exVir, inVir;
        vector<SQindex*> inds(t->get_indices());
        for(vector<SQindex*>::iterator i = inds.begin();i != inds.end();++i){
          if     ( (*i)->get_isExt() && (*i)->get_char() == (char_state)0) exCore.push_back((*i)->get_index());
          else if( (*i)->get_isExt() && (*i)->get_char() == (char_state)1) exAct. push_back((*i)->get_index());
          else if( (*i)->get_isExt() && (*i)->get_char() == (char_state)2) exVir. push_back((*i)->get_index());
          else if(!(*i)->get_isExt() && (*i)->get_char() == (char_state)0) inCore.push_back((*i)->get_index());
          else if(!(*i)->get_isExt() && (*i)->get_char() == (char_state)1) inAct. push_back((*i)->get_index());
          else if(!(*i)->get_isExt() && (*i)->get_char() == (char_state)2) inVir. push_back((*i)->get_index());
	} // End i
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
          cout << "Factorize: Number of external indices exceed the limit" << endl;
          abort();
	}
        SetX += "\n";
        f << SetX;

        string nameX("x");
        if(inCore.size())
          for(size_t i = 0;i < inCore.size();++i) nameX += "c";
        if(inAct.size())
          for(size_t i = 0;i < inAct.size();++i)  nameX += "a";
        if(inVir.size())
          for(size_t i = 0;i < inVir.size();++i)  nameX += "v";
        SetX = "call set_symblock_" + nameX + "(sleft, x, nir, nsym, psym) ! -> " + nameX + " (allocate) \n"; 
        if(inCore.size() || inAct.size() || inVir.size())
	  f << SetX;
	extTen.insert(map<string, string>::value_type("X", nameX));
        keys.push_back("X");
          
      } // End if
      else if(t->get_name() == "X" && X_indices_.size()){
        string SetX("");
        vector<string> inCore, inAct, inVir;
        for(vector<SQindex*>::iterator i = X_indices_.begin();i != X_indices_.end();++i){
          if     ((*i)->get_char() == (char_state)0) inCore.push_back((*i)->get_index());
          else if((*i)->get_char() == (char_state)1) inAct. push_back((*i)->get_index());
          else if((*i)->get_char() == (char_state)2) inVir. push_back((*i)->get_index());          
	} // End if
	SetX += "sleft = 0\n";
        f << SetX;

        string nameX("x");
        if(inCore.size())
          for(size_t i = 0;i < inCore.size();++i) nameX += "c";
        if(inAct.size())
          for(size_t i = 0;i < inAct.size();++i)  nameX += "a";
        if(inVir.size())
          for(size_t i = 0;i < inVir.size();++i)  nameX += "v";
        SetX = "call set_symblock_" + nameX + "(sleft, x, nir, nsym, psym) ! -> " + nameX + " (allocate) \n"; 
        if(inCore.size() || inAct.size() || inVir.size())
	  f << SetX;
	extTen.insert(map<string, string>::value_type("X", nameX));
        keys.push_back("X");
          
      } // End if
      else if(t->get_name()[0] == 'Y' && t->get_indices().size()){
        vector<string> Core, Act, Vir;
        vector<SQindex*> inds(t->get_indices());
        for(vector<SQindex*>::iterator i = inds.begin();i != inds.end();++i){
          if     ((*i)->get_char() == (char_state)0) Core.push_back((*i)->get_index());
          else if((*i)->get_char() == (char_state)1) Act. push_back((*i)->get_index());
          else if((*i)->get_char() == (char_state)2) Vir. push_back((*i)->get_index());
	} // End i

        string nameY("Y");
        if(Core.size())
          for(size_t i = 0;i < Core.size();++i) nameY += "c";
        if(Act.size())
          for(size_t i = 0;i < Act.size();++i)  nameY += "a";
        if(Vir.size())
          for(size_t i = 0;i < Vir.size();++i)  nameY += "v";
        string SetY("call set_symblock_" + nameY + "(sleft2, " + t->get_name() + ", nir, nsym, psym) ! -> " + nameY + " (allocate) \n"); 
        f << "sleft2 = 0" << endl;
	f << SetY;
	extTen.insert(map<string, string>::value_type(t->get_name(), nameY));
        keys.push_back(t->get_name());
          
      } // End if
      else if(t->get_name() == name_amp_){
        string iname = t->get_indices()[extamp_]->get_index();
        string SetT("call set_symblock_av2(s" + iname + ", " + t->get_name() + ", nir, nsym, psym) ! -> av2_i (allocate)\n");
        f << SetT;
        extTen.insert(map<string,string>::value_type(name_amp_, "av2_i"));
        keys.push_back(name_amp_);
      } // End if
      else if(t->get_name() == name_h1_){
        string Seth1("call set_symblock_h1(" + t->get_name() + ", nir, nsym, psym) ! -> h1 (allocate)\n");
        f << Seth1;
        extTen.insert(map<string,string>::value_type(t->get_name(), "h1"));        
        keys.push_back(t->get_name());
      } // End if
      else if(t->get_name() == name_h2_){
        string iname(t->get_indices()[exth2_]->get_index());
        string Seth2("call set_symblock_h2(s" + iname + ", " + t->get_name() + ", nir, nsym, psym) ! -> h2_i (allocate)\n");
        f << Seth2;
        extTen.insert(map<string,string>::value_type(name_h2_, "h2_i"));
        keys.push_back(name_h2_);
      } // End if
    } // End t

    if(isBareLHS){
      string SetLT("call set_symblock_av2_2(s" + LTensor.get_indices()[extamp_]->get_index() + ", " + LTensor.get_name() + ", nir, nsym, psym) ! -> av2_i2 (allocate)\n");
      f << SetLT;
      extTen.insert(map<string,string>::value_type(LTensor.get_name(), "av2_i2"));
      keys.push_back(LTensor.get_name()); 
    } // End if

    string FuncName(title);
    string::iterator s_end = FuncName.begin();
    ++s_end; ++s_end; ++s_end; ++s_end;
    FuncName.erase(FuncName.begin(), s_end);
    string SetFunc("call g" + FuncName + "(");

    for(vector<string>::iterator i = ExtInd.begin();i != ExtInd.end();++i)
      SetFunc += "s" + *i + ", i" + *i + ", ";
    for(vector<string>::iterator C = Consts.begin();C != Consts.end();++C)
      SetFunc += *C + ", ";
    for(vector<string>::iterator t = NameTen.begin();t != NameTen.end();++t){
      if(find(keys.begin(), keys.end(), *t) != keys.end())
        SetFunc += extTen[*t] + ", ";
      else
        SetFunc += *t + ", ";
    } // End t

    vector<string> d_list;
    for(vector<SQtensor>::iterator t = tensors.begin();t != tensors.end();++t) {
      //f << t->get_name() << endl; //*TEST*
      if(is_RDM(t->get_name())) {d_list.push_back(t->get_name()); /*f<<"OK " << t->get_name();*/}
    }
    for(vector<string>::iterator s = d_list.begin();s != d_list.end();++s)
      if(*s == RDM_name() + "4") SetFunc += "d4_ij, ";
      else {
        transform(s->begin(), s->end(), s->begin(), (int(*)(int))tolower);
        SetFunc += *s + ", ";
      } // End else

    SetFunc += "nir, nsym, psym, flops)\n\n";
    f << SetFunc;

    for(vector<string>::iterator s = keys.begin();s != keys.end();++s){
      if(extTen[*s] == "X" || extTen[*s] == "x") continue;
      string SetDealloc("deallocate(" + extTen[*s] + ")\n");
      f << SetDealloc;
    } // End s
    f << endl;

    f << "end subroutine " + title + "\n\n\n";

  }


  // *********************************************************
  // 
  // *********************************************************
  void SQreaktor::makeF90_contract(SQtensor &LTensor, SQterm &inTerm, string title, ofstream &f, string s, contDecl c, 
      bool isBareLHS)
  {

    if(X_indices_.size()){
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
    string printName2("subroutine " + title + "(");
    
    string printEnd("end subroutine " + title + "\n\n");

    string DecComm8;
    DecComm8  = "! FEMTO BEGIN  **************************************************************\n";
    DecComm8 += "use f_mr_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6\n";
    DecComm8 += "\n";
    if(use_omp_)
    DecComm8 += "use omp_lib\n";
    DecComm8 += "implicit none\n\n";
    
    string DecExtInd("");
    vector<SQtensor> tensors(inTerm.get_tensors());
    vector<SQindex*> summedBody(inTerm.get_summedBody());
    vector<SQindex*> LTinds(LTensor.get_indices());
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
      mapTen.insert( map<string, int>::value_type(LTensor.get_name(), 3) );
      keys.push_back(LTensor.get_name());
    } // End if
    else if(!isBareLHS && LTensor.get_name() == "X"){
      int NDcount = 0;
      vector<SQindex*> inds(LTensor.get_indices());
      for(vector<SQindex*>::iterator i = inds.begin();i != inds.end();++i)
        if(!(*i)->get_isExt()) ++NDcount;
      mapTen.insert( map<string, int>::value_type(LTensor.get_name(), NDcount) );
      keys.push_back(LTensor.get_name());
    } // End if
    else{
      mapTen.insert( map<string, int>::value_type(LTensor.get_name(), LTensor.get_indices().size()) );
      keys.push_back(LTensor.get_name());
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
        else if(t->get_name() == "X"){
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
    for(vector<string>::iterator t = keys.begin();t != keys.end();++t){
      if(is_RDM(*t)) printName2 += *t + "_, ";
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
    f << "! " << LTensor << " <-- \n";
    f << "! " << inTerm << endl;

    for(vector<SQindex*>::iterator i = NDindices.begin();i != NDindices.end();++i){
      f << "do s_" + (*i)->get_index() + " = 0, nir-1" + "\n";
      ++LoopCount;
    } // End i

    // Symmetry constraints
    f << "if( &" << endl;
    if     (LTensor.get_indices().size() == 1){ // In case of unity
      vector<SQindex*> inds(LTensor.get_indices());
      f << "s_" + inds[0]->get_index() + " == 0";
    } // End if
    else if(LTensor.get_indices().size() == 2){ // In case of two
      vector<SQindex*> inds(LTensor.get_indices());
      f << "IEOR(s_" + inds[0]->get_index() + ",s_" + inds[1]->get_index() + ") == 0";
    } // End if
    else if(LTensor.get_indices().size() == 3){ // In case of three
      vector<SQindex*> inds(LTensor.get_indices());
      string cstr;
      cstr  = "IEOR(IEOR(s_" + inds[0]->get_index() + ",s_" + inds[1]->get_index() + "),";
      cstr += "s_" + inds[2]->get_index() + ") == 0";
      f << cstr;
    } // End if
    else if(LTensor.get_indices().size() == 4){ // In case of four
      vector<SQindex*> inds(LTensor.get_indices());
      string cstr;
      cstr  = "IEOR(s_" + inds[0]->get_index() + ",s_" + inds[1]->get_index() + ") == ";
      cstr += "IEOR(s_" + inds[2]->get_index() + ",s_" + inds[3]->get_index() + ")";
      f << cstr;
    } // End if
    else if(LTensor.get_indices().size() == 5){ // In case of five
      vector<SQindex*> inds(LTensor.get_indices());
      string cstr;
      cstr  = "IEOR(IEOR(s_" + inds[0]->get_index() + ",s_" + inds[1]->get_index() + "),s_" + inds[2]->get_index() + ") == ";
      cstr += "IEOR(s_" + inds[3]->get_index() + ",s_" + inds[4]->get_index() + ")";
      f << cstr;
    } // End if
    else if(LTensor.get_indices().size() == 6){ // In case of six
      vector<SQindex*> inds(LTensor.get_indices());
      string cstr;
      cstr  = "IEOR(IEOR(s_" + inds[0]->get_index() + ",s_" + inds[1]->get_index() + "),s_" + inds[2]->get_index() + ") == ";
      cstr += "IEOR(IEOR(s_" + inds[3]->get_index() + ",s_" + inds[4]->get_index() + "),s_" + inds[5]->get_index() + ")";
      f << cstr;
    } // End if
    else if(LTensor.get_indices().size() == 7){ // In case of seven
      vector<SQindex*> inds(LTensor.get_indices());
      string cstr;
      cstr  = "IEOR(IEOR(s_" + inds[0]->get_index() + ",s_" + inds[1]->get_index() + "),IEOR(s_" + inds[2]->get_index() + ",s_" + inds[3]->get_index() + ")) == ";
      cstr += "IEOR(IEOR(s_" + inds[4]->get_index() + ",s_" + inds[5]->get_index() + "),s_" + inds[6]->get_index()+ ")";
      f << cstr;
    } // End if
    else if(LTensor.get_indices().size() == 8){ // In case of eight
      vector<SQindex*> inds(LTensor.get_indices());
      string cstr;
      cstr  = "IEOR(IEOR(s_" + inds[0]->get_index() + ",s_" + inds[1]->get_index() + "),IEOR(s_" + inds[2]->get_index() + ",s_" + inds[3]->get_index() + ")) == ";
      cstr += "IEOR(IEOR(s_" + inds[4]->get_index() + ",s_" + inds[5]->get_index() + "),IEOR(s_" + inds[6]->get_index() + ",s_" + inds[7]->get_index() + "))";
      f << cstr;
    } // End if
    else if(LTensor.get_indices().size() != 0){
      cout << "Factorize: Number of indices in tensors must be 0 - 8" << endl;
      abort();
    } // End if
    if (LTensor.get_indices().size() && tensors.size())
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
      if(!LTensor.get_indices().size()) f << " reduction(+:" + LTensor.get_name() + "_)";
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
    if     (!isBareLHS && LTensor.get_name() != "X"){
      Lindlist = LTensor.get_indices();
    } // End if
    else if(!isBareLHS && LTensor.get_name() == "X"){
      for(size_t num_i = 0;num_i < LTensor.get_indices().size();++num_i){
        if(!LTensor.get_indices()[num_i]->get_isExt()) 
          Lindlist.push_back(LTensor.get_indices()[num_i]);
      } // End num_i
    } // End if
    else{
      for(size_t num_i = 0;num_i < LTensor.get_indices().size();++num_i){
        if(num_i != extamp_) 
          Lindlist.push_back(LTensor.get_indices()[num_i]);
      } // End num_i      
    } // End 

    if(Lindlist.size()){
      string nterm(LTensor.get_name() + "_(");
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
      f << LTensor.get_name() + "_ = " + LTensor.get_name() + "_ &" << endl;
    } // End else

    // Generate the numerical coefficients
    if(inTerm.get_numConst() > 0) {
      ostringstream stm;
      stm << (boost::format("%10.8f") % inTerm.get_numConst());
      f << "  + " + stm.str() + "d+00" + " & " << endl;
    } // end if
    else {
      ostringstream stm;
      stm << (boost::format("%10.8f") % fabs(inTerm.get_numConst()));
      f << "  - " + stm.str() + "d+00" + " & " << endl;
    } // end if

    // Print all the Coeffient
    for(size_t i = 0;i < inTerm.get_Consts().size();++i){
      if(inTerm.get_Consts()[i] != "") f << "  * " + inTerm.get_Consts()[i];
      if(inTerm.get_Consts()[i] != "" && tensors.size()) f << " & " << endl;
      else if(inTerm.get_Consts()[i] != "")              f << endl;
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
      else if(t->get_name() == "X"){
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
