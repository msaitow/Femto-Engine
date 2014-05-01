//
//  Convert_gemm2.cc
//  
//
//  Created by Masaaki Saitow on 12/10/30.
//  Copyright (c) 2012 Masaaki Saitow. All rights reserved.
//

#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <SQreaktor.hpp>
#include <boost/format.hpp>

#define _DEBUG_BLAS // Debug option for use of blas
//#define _USE_DAXPY // Debug option for complementary use of daxpy 

using namespace std;
using namespace Femto;
using namespace Femto::Core;

namespace Femto { namespace Reaktor {

  // *********************************************************
  // 
  // *********************************************************
  void SQreaktor::binary_contract2(SQbinary &bin, string title, ofstream &f, string s, contDecl c, bool isBareLHS)
  {

    if(bin.get_Rtensors().size() != 2){
      cout << "Number of terms pushed into binary_contract should be 2." << endl;
      abort();
    } 

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

    // Set-up the X_indices, which is true associated indices to intermediate array.
    vector<SQindex*> X_indices;
    if(bin.get_Lindices()) X_indices = bin.get_Ltensor().get_indices();
    else {
      for(size_t i = 0;i < bin.get_Rindices().size();++i)
        if(bin.get_Rindices()[i]) X_indices = bin.get_Rtensors()[i].get_indices();
    } // End else

    // Print all the names of tensors
    map<string, int> mapTen; // Correspondence between tensor name and number of the associated internal indices
    vector<string> keys;
    if(isBareLHS){
      mapTen.insert( map<string, int>::value_type(bin.get_Ltensor().get_name(), 3) );
      keys.push_back(bin.get_Ltensor().get_name());
    } // End if
    //else if(!isBareLHS && bin.get_Ltensor().get_name() == "X"){
    else if(!isBareLHS && is_Interm(bin.get_Ltensor().get_name())){
      int NDcount = 0;
      if(!X_indices.size()){
        vector<SQindex*> inds(bin.get_Ltensor().get_indices());
        for(vector<SQindex*>::iterator i = inds.begin();i != inds.end();++i)
          if(!(*i)->get_isExt()) ++NDcount;
      } // End if
      else NDcount = X_indices.size();
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
          if(!X_indices.size()){
          vector<SQindex*> inds(t->get_indices());
          for(vector<SQindex*>::iterator i = inds.begin();i != inds.end();++i)
            if(!(*i)->get_isExt()) ++NDcount;
	  } // End if
          else NDcount = X_indices.size();
          mapTen.insert( map<string, int>::value_type(t->get_name(), NDcount) );
          keys.push_back(t->get_name());
	} // End if
        else if(is_D4C(t->get_name())){
          mapTen.insert( map<string, int>::value_type(t->get_name(), 5) );
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
    // This bunch of codes should be same to Convert2.cc.
    for(vector<string>::iterator t = keys.begin();t != keys.end();++t){
      if(is_RDM(*t) || *t == "Fc1" || is_C4(*t) || is_C6(*t)) printName2 += *t + "_, ";
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

    // Determine the order of the tensors in dgemm
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
    } // End if
    else{
      vector<SQindex*> temp_i(bin.get_Rtensors()[0].get_indices());
      if(find(temp_i.begin(), temp_i.end(), L_nd_ind[0]) != temp_i.end()){
        the_tensors.push_back(bin.get_Rtensors()[0]);
        the_tensors.push_back(bin.get_Rtensors()[1]);
      } // End if
      else{
        the_tensors.push_back(bin.get_Rtensors()[1]);
        the_tensors.push_back(bin.get_Rtensors()[0]);
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
    //vector<SQindex*> I3(bin.get_Ltensor().get_indices());
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

    // Variables that define which intermediate is allocated 
    vector<bool> is_allocated(3);

    bool use_daxpy = false; // If dim of Z2 becomes zero, daxpy is used instead of dgemm
    int count = 1;
    string DecInterm("! Intermediate arrays                \n");
    for(vector<SQtensor>::iterator t = the_tensors.begin();t != the_tensors.end();++t, ++count){
      ostringstream stm;
      stm << count;
      size_t id_count = (count == 1 ? rowInd.size() : colInd.size()) + summedInd.size();
      if(!id_count) {
#ifdef _USE_DAXPY
        use_daxpy = true; //if(use_daxpy) f << "****** count : " << count << endl; //*TEST*
#endif
        if(count == 1 && use_daxpy){
          cout << "binary_contract: !!!WARNING!!! Weird contraction form detected." << endl; //*TEST* 
          //*TEST*cout << "binary_contract: An error occured" << endl;
          //*TEST*abort();
	}
        //*TEST*break;
      } // End if
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

    // Allocate the first tensor
    f << "! Z1 <-- " << the_tensors[0] << endl;
    if(is_allocated[0]) f << "allocate(Z1_(";
    vector<SQindex*> Z1inds(rowInd);
    Z1inds.insert(Z1inds.end(), summedInd.begin(), summedInd.end());
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
    else if(is_Interm(the_tensors[0].get_name()) && !X_indices.size()){
      for(size_t num_i = 0;num_i < the_tensors[0].get_indices().size();++num_i){
        SQindex* i(the_tensors[0].get_indices()[num_i]);
        if(!i->get_isExt()) Rindlist.push_back(i);
      } // End num_i
    } // End if
    //else if(the_tensors[0].get_name() == "X" &&  X_indices.size()){
    else if(is_Interm(the_tensors[0].get_name()) &&  X_indices.size()){
      Rindlist.insert(Rindlist.end(), X_indices.begin(), X_indices.end());
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

    if(!use_daxpy){        
    // Allocate the second tensor
    f << "! Z2 <-- " << the_tensors[1] << endl;
    if(is_allocated[1]) f << "allocate(Z2_(";
    vector<SQindex*> Z2inds(summedInd);
    Z2inds.insert(Z2inds.end(), colInd.begin(), colInd.end());
    for(vector<SQindex*>::iterator i = Z2inds.begin();i != Z2inds.end();++i){
      string sname("s_" + (*i)->get_index());
      string otype;
      if     ((*i)->get_char() == (char_state)0) otype = "I_C";
      else if((*i)->get_char() == (char_state)1) otype = "I_O";
      else if((*i)->get_char() == (char_state)2) otype = "I_V";

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
      if     ((*i)->get_char() == (char_state)0) otype = "I_C";
      else if((*i)->get_char() == (char_state)1) otype = "I_O";
      else if((*i)->get_char() == (char_state)2) otype = "I_V";
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
    else if(is_Interm(the_tensors[1].get_name()) && !X_indices.size()){
      for(size_t num_i = 0;num_i < the_tensors[1].get_indices().size();++num_i){
        SQindex* i(the_tensors[1].get_indices()[num_i]);
        if(!i->get_isExt()) Rindlist.push_back(i);
      } // End num_i
    } // End if
    //else if(the_tensors[1].get_name() == "X" &&  X_indices.size()){
    else if(is_Interm(the_tensors[1].get_name()) &&  X_indices.size()){
      Rindlist.insert(Rindlist.end(), X_indices.begin(), X_indices.end());
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
    } // End if (!use_daxpy)

    vector<SQindex*> Z3inds(rowInd);
    Z3inds.insert(Z3inds.end(), colInd.begin(), colInd.end());
    if(Z3inds.size()){
    // Allocate the resultant tensor
    f << "! Z3 <-- " << bin.get_Ltensor() << endl;
    f << "allocate(Z3_(";
    for(vector<SQindex*>::iterator i = Z3inds.begin();i != Z3inds.end();++i){
      string sname("s_" + (*i)->get_index());
      string otype;
      if     ((*i)->get_char() == (char_state)0) otype = "I_C";
      else if((*i)->get_char() == (char_state)1) otype = "I_O";
      else if((*i)->get_char() == (char_state)2) otype = "I_V";
      
      if(i == Z3inds.begin())
      f << "psym(I_BEGIN," + otype + ", " + sname + "):psym(I_END," + otype + ", " + sname + ")";
      else
      f << "             psym(I_BEGIN," + otype + ", " + sname + "):psym(I_END," + otype + ", " + sname + ")";

      if(*i == Z3inds.back()) f << "))\n";
      else                    f << ", &\n";
    } // End i
    // Zero out Z3_, if neccessary
    if(use_daxpy){
      size_t size = colInd.size() + rowInd.size();
      f << "Z3_(";
      for(size_t i = 0;i < size;++i){
        if(i != size-1) f << ":,";
        else            f << ":) = 0.0d+00" << endl;;
      } // End i
    } // End if
    f << endl;
    } // End if

    // Flop count
    string flopcount("flops = flops + ");

    if(!use_daxpy){
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
    f << "                     Z1_,&"               << endl;
    f << "                     " << m_label << ",&" << endl;
    f << "                     Z2_,&"               << endl;
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
    } // End if (use_daxpy)
    else{
    // In case of daxpy
    f << "call daxpy(";

    string length("");         
    for(vector<SQindex*>::iterator i = Z1inds.begin();i != Z1inds.end();++i){
      string sname("s_" + (*i)->get_index());
      string otype;
      if     ((*i)->get_char() == (char_state)0) otype = "I_C";
      else if((*i)->get_char() == (char_state)1) otype = "I_O";
      else if((*i)->get_char() == (char_state)2) otype = "I_V";
      
      length += "psym(I_LENGTH," + otype + ", " + sname + ")";

      if(*i != Z1inds.back()) length += "*"; 
    } // End i
    if(!Z1inds.size()) length += "1";
    f << length << ", &" << endl;
    flopcount += length + "\n";

    string alpha("           ");    
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
    //f << alpha;

    if(bin.get_Rtensors().size() == 2) {
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
    else if(is_Interm(the_tensors[1].get_name()) && !X_indices.size()){
      for(size_t num_i = 0;num_i < the_tensors[1].get_indices().size();++num_i){
        SQindex* i(the_tensors[1].get_indices()[num_i]);
        if(!i->get_isExt()) Rindlist.push_back(i);
      } // End num_i
    } // End if
    //else if(the_tensors[1].get_name() == "X" &&  X_indices.size()){
    else if(is_Interm(the_tensors[1].get_name()) &&  X_indices.size()){
      Rindlist.insert(Rindlist.end(), X_indices.begin(), X_indices.end());
    } // End if
    else{
      Rindlist = the_tensors[1].get_indices();
    } // End else
    
    for(vector<SQindex*>::iterator i = Rindlist.begin();i != Rindlist.end();++i){
      if(!(*i)->get_isExt() && !X_indices.size()){
        cout << "Something wrong in binary_contract: in case of daxpy .... " << **i << endl;
        //*TEST*abort();
      } // End if
    } // End i

    string nterm2(the_tensors[1].get_name() + "_");
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
    alpha += "*" + nterm2;
    } // End if

    f << alpha << ", &"<< endl;
    f << "           Z1_, &" << endl;
    f << "           1, &" << endl;
    f << "           Z3_, &" << endl;
    f << "           1)" << endl << endl;
    }

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
    if     (!isBareLHS && !is_Interm(bin.get_Ltensor().get_name())){
      Lindlist = bin.get_Ltensor().get_indices();
    } // End if
    else if(!isBareLHS && is_Interm(bin.get_Ltensor().get_name())){
      if(!X_indices.size()){
        for(size_t num_i = 0;num_i < bin.get_Ltensor().get_indices().size();++num_i){
          if(!bin.get_Ltensor().get_indices()[num_i]->get_isExt()) 
            Lindlist.push_back(bin.get_Ltensor().get_indices()[num_i]);
        } // End num_i
      } // End if
      else{
        Lindlist.insert(Lindlist.end(), X_indices.begin(), X_indices.end());
      } // End else
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
    if(is_allocated[0])               Z_names.push_back("Z1_");
    if(is_allocated[1] && !use_daxpy) Z_names.push_back("Z2_");
    if(is_allocated[2])               Z_names.push_back("Z3_");
    if(Z_names.size()) f << "deallocate(";
    for(vector<string>::iterator s = Z_names.begin();s != Z_names.end();++s)
      if(s != Z_names.end()-1) f << *s << ", ";
      else                     f << *s << ")\n";
//*OLD*     if(use_daxpy)
//*OLD*       f << "deallocate(Z1_, Z3_)\n";
//*OLD*     else if(Z3inds.size())
//*OLD*       f << "deallocate(Z1_, Z2_, Z3_)\n";
//*OLD*     else
//*OLD*       f << "deallocate(Z1_, Z2_)\n";

    if(extInd.size())                f << "end do ! Orbital Loop" << endl << endl;
    if(dimCount)                     f << "end if ! Dim Const"    << endl;
    f << endl;
                                     f << "end if ! Irrep Cond"   << endl;
    for(int n = 0;n < LoopCount;++n) f << "end do ! Irrep Loop"   << endl;
    f << "! FEMTO END  ****************************************************************\n" << endl;
    f << printEnd;
    f.flush();
    //X_indices.clear();
  }


}} // Femto::

