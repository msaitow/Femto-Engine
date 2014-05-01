//
//  Orz.hpp
//
//  Some facilities to convert the tensorial terms gerenerated by femto to
//  working code pf Orz suite.  
//
//
//  Created by Masaaki Saitow on 12/03/26.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#pragma once

//#include <femto.hpp>
#include <SQterm.hpp>
#include <SQtensor.hpp>
#include <SQindex.hpp>
#include <map>

using namespace std;

namespace femto {

  typedef map<string, vector<string> > contDecl; 

  // detects 2J-K typed terms (/orz)
  void find_JK_pairs(vector<SQterm> &inTerm, vector<SQterm> &outTerms, string name_h1="h", string name_h2="V2");

  // detects 2J-K typed terms (/orz)
  vector<SQterm>& find_JK_pairs(vector<SQterm> &inTerm, string name_h1="h", string name_h2="V2");

  // Convert tensorial terms into factorized Orz code (/Orz) .... 
  void factorize(SQtensor &LTensor, vector<SQterm> &inTerms, string title, bool isBareLHS=false, string name_h2="V2", string name_amp="T2");

  // Convert tensorial terms into Orz code, in which tensorial contraction is carried out in simple loops (/Orz) .... 
  void simpleloops(SQtensor &LTensor, vector<SQterm> &inTerms, string title, bool isBareLHS=false, string name_h2="V2", string name_amp="T2");

  // Loop in C++ style (/orz)
  void CLoop(const string s, const SQindex &i, ofstream &f);

  // Read BareAmp from GA (/orz)
  void ReadAmp(const string s, const SQtensor &t, ofstream &f);

  // Read retval (/orz)
  void ReadRetval(const string s, const SQtensor &t, ofstream &f);

  // Read ERIs from GA (/orz)
  void ReadERI(const string s, const SQtensor &t, ofstream &f);

  // Read D4 from GA (/orz)
  void ReadD4(const string s, const SQtensor &t, ofstream &f);

  // Generate D4 by cumulant expansion (/orz)
  void ReadD4_Cumulant(const string s, const SQtensor &t, ofstream &f);

  // Accumulate BareAmp (/orz)
  void AccAmp(const string s, const SQtensor &t, ofstream &f);

  // Print end of the loop (/orz)
  void LoopEnd(string s, ofstream &f);

  // Print the calling section of the Fortran routine in C++ side (/orz)
  void makeCPP_body(SQtensor &LTensor, SQterm &inTerm, string title, ofstream &f, string s, contDecl c, 
      bool isBareLHS, string name_h2, string name_amp);

  // Print header of the Fortran routine in C++ side (/orz)
  void makeCPP_header(SQtensor &LTensor, SQterm &inTerm, string title, ofstream &f, string s, contDecl c, 
      bool isBareLHS, string name_h2, string name_amp);

  // Print interface of the F90 subroutine (/orz)
  void makeF90_interface(SQtensor &LTensor, SQterm &inTerm, string title, ofstream &f, string s, contDecl c, 
      bool isBareLHS, string name_h2, string name_amp);

  // Print the body of tensorial contraction (/orz)
  void makeF90_contract(SQtensor &LTensor, SQterm &inTerm, string title, ofstream &f, string s, contDecl c, 
      bool isBareLHS, string name_h2, string name_amp);

  // Print the body of tensorial contraction processed with dgemm(/orz)
  void binary_contract(SQtensor &LTensor, SQterm &inTerm, string title, ofstream &f, string s, contDecl c, 
      bool isBareLHS, string name_h2, string name_amp);


} // End femto

