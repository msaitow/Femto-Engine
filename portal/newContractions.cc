//
//  SQportal.cc
//  
//
//  Created by Masaaki Saitow on 14/02/21.
//  Copyright (c) 2012 Masaaki Saitow. All rights reserved.
//

#include <SQportal.hpp>
#include <boost/tokenizer.hpp>

#define _DEBUG1

using namespace std;

using namespace Femto;
using namespace Femto::Core;
using namespace Femto::Reaktor;

namespace Femto { namespace Portal {

       
    // *********************************************************
    // 
    // *********************************************************
    SQbinary 
    SQportal::readContractionNew(string myContraction) const
    {

      cout << " >>----------------------------------------------------<< " << endl;
      cout << " >> readContractionNew is called                       << " << endl;
      cout << " >>----------------------------------------------------<< " << endl;

      SQcont<SQcont<SQindex> > myIndices;
      SQcont<SQtensor> myBodies; 
      SQcont<string> myTensors; 
      
      typedef boost::char_separator<char> char_separator;
      typedef boost::tokenizer<char_separator> tokenizer;
      
      //SQcont<string> myTensors;
      string myFactor;
      string myLiteral;
      char_separator sep(" ", "", boost::keep_empty_tokens);
      tokenizer tokens(myContraction, sep);
      
      // Extract components of the contraction
      int myNum(0);
      for (tokenizer::iterator tok_iter = tokens.begin(); tok_iter != tokens.end(); ++tok_iter, ++myNum) 
	if(myNum == 0 || myNum == 4 || myNum == 5)
	  myTensors <= *tok_iter;
	else if(myNum == 2)
	  myFactor = *tok_iter;
	else if(myNum == 3)
	  myLiteral = *tok_iter;
      
      // Collect all the indices
      { // Scope begins
	
	for(auto t = myTensors.begin();t != myTensors.end();++t){
	  SQcont<SQindex> batch;
	  size_t myBegin(t->find("(", 0));
	  string myCont(t->substr(myBegin+1));
	  myCont.erase(--myCont.end());
	  cout << " -- My contraction >> " << myCont << endl;
	  
	  char_separator my_sep(",", "", boost::keep_empty_tokens);
	  tokenizer myTokens(myCont, my_sep);
	  
	  // Read the indices
	  for (tokenizer::iterator tok_iter = myTokens.begin(); tok_iter != myTokens.end(); ++tok_iter) {
	    batch <= readIndex(*tok_iter);
	    cout << " --    >> ind " << " - "<< *tok_iter << "- " << endl;
	    cout << " --    >> body " << batch.back() << endl; 
	  } // End tok_iter
	  myIndices <= batch;
	  
	} // End t
      } // Scope ends
      
      // Collect the tensors
      { // Scope begins
	
	vector<vector<SQindex*> > myIndices_ptr;
	for(auto i = myIndices.begin();i != myIndices.end();++i){
	  vector<SQindex*> temp;
	  for(auto I = i->begin();I != i->end();++I) temp.push_back(&(*I));
	  myIndices_ptr.push_back(temp);
	} // End i
	
	for(auto t = myTensors.begin();t != myTensors.end();++t){
	  string copied(*t);
	  size_t myBegin(t->find("@", 0));
	  size_t myEnd(t->find("(", 0));
	  string name(copied.substr(myBegin, calc_distance(myBegin, myEnd)+1));
	  cout << " -- My name is : " << name << endl;
	  
	  size_t myCount((size_t)(t-myTensors.begin()));
	  // In case of the contraction is given in the Dirac notation
	  if(myNotation_ == Femto::Dirac){
	    if     (  t_t2Amps().count(name)) myBodies <= Femto::Core::SQtensor("T2",  myIndices_ptr[myCount], Femto::t2_symm());
	    else if(  t_s2Amps().count(name)) myBodies <= Femto::Core::SQtensor("S2",  myIndices_ptr[myCount], Femto::t2_symm());
	    else if(      t_h1().count(name)) myBodies <= Femto::Core::SQtensor("h" ,  myIndices_ptr[myCount], Femto::h1_symm());
	    else if(      t_h2().count(name)) myBodies <= Femto::Core::SQtensor("V2",  myIndices_ptr[myCount], Femto::h2_symm());
	    else if(    t_Rdms().count(name)) myBodies <= Femto::Core::RDM(            myIndices_ptr[myCount]                  );	
	    else if(     t_Cum().count(name)) myBodies <=         Cumulant(            myIndices_ptr[myCount]                  );
	    else if( t_CasFock().count(name)) myBodies <= Femto::Core::SQtensor("P1",  myIndices_ptr[myCount], Femto::h1_symm());
	    else if(t_CoreFock().count(name)) myBodies <= Femto::Core::SQtensor("Fc1", myIndices_ptr[myCount], Femto::h1_symm());
	    else if(ltensorName_ == name){
	      string temp(ltensorName_);
	      temp.erase(temp.begin());
	      myBodies <= Femto::Core::SQtensor(temp, myIndices_ptr[myCount], Femto::uni_symm(myIndices_ptr[myCount].size()));
	    } // End if
	    else if(     is_Interm(name)){
	      name.erase(name.begin());
	      myBodies <= Femto::Core::SQtensor(name, myIndices_ptr[myCount], Femto::uni_symm(myIndices_ptr[myCount].size()));
	    } // End if
	    else if(  t_Null().count(name)) { cout << " =>> A dummy tensor detected " << endl; }
	    else{
	      cout << "readContraction: Can't read this >> " << *t << endl;
	      abort();
	    } // Else
	  } // End if(notation)
	  // In case of the contraction is given in the Mulliken notation
	  if(myNotation_ == Femto::Mulliken){
	    if     (  t_t2Amps().count(name)) myBodies <= Femto::Core::SQtensor ("T2",  myIndices_ptr[myCount], Femto::t2_symm());
	    else if(  t_s2Amps().count(name)) myBodies <= Femto::Core::SQtensor ("S2",  myIndices_ptr[myCount], Femto::t2_symm());
	    else if(      t_h1().count(name)) myBodies <= Femto::Portal::Mtensor("h" ,  myIndices_ptr[myCount], Femto::h1_symm());
	    else if(      t_h2().count(name)) myBodies <= Femto::Portal::Mtensor("V2",  myIndices_ptr[myCount], Femto::h2_symm());
	    else if(    t_Rdms().count(name)) myBodies <= Femto::Portal::MRDM(          myIndices_ptr[myCount]                  ); 
	    else if(     t_Cum().count(name)) myBodies <=        MCumulant(             myIndices_ptr[myCount]                  );
	    else if( t_CasFock().count(name)) myBodies <= Femto::Portal::Mtensor("P1",  myIndices_ptr[myCount], Femto::h1_symm());
	    else if(t_CoreFock().count(name)) myBodies <= Femto::Portal::Mtensor("Fc1", myIndices_ptr[myCount], Femto::h1_symm());
	    else if(ltensorName_ == name){
	      string temp(ltensorName_);
	      temp.erase(temp.begin());
	      myBodies <= Femto::Core::SQtensor(temp, myIndices_ptr[myCount], Femto::uni_symm(myIndices_ptr[myCount].size()));
	    } // End if
	    else if(     is_Interm(name)){
	      name.erase(name.begin());
	      myBodies <= Femto::Core::SQtensor(name, myIndices_ptr[myCount], Femto::uni_symm(myIndices_ptr[myCount].size()));
	    } // End if
	    else if(  t_Null().count(name)) { cout << " =>> A dummy tensor detected " << endl; }
	    else{
	      cout << "readContraction: Can't read this >> " << *t << endl;
	      abort();
	    } // Else
	  } // End if(notation)
	  
	} // End t
	
      } // Scope ends
      
      // Construct the SQbinary object and returns the result
      { // Scope begins
	SQtensor ltensor(myBodies[0]);
	vector<SQtensor> rtensors;
	myBodies.erase(myBodies.begin());
	for(auto t = myBodies.begin();t != myBodies.end();++t) rtensors.push_back(*t);
	double factor(atof(myFactor.c_str()));
	vector<string> coeffs;
	if (!c_Unit().count(myLiteral)) 
	  if(c_Ecore().count(myLiteral)) coeffs.push_back("Fc0"    );
	  else                           coeffs.push_back(myLiteral);
	cout << "myLiteral : >>" << myLiteral << "<< " << endl; //*TEST* 
	SQbinary myTerm(factor, coeffs, ltensor, rtensors);
	
#ifdef _DEBUG1
	cout << "Interpreted: " << myTerm << endl;
	myTerm.print_summedBody();
#endif
	return myTerm;
      } // Scope ends
      
    }
      

    // *********************************************************
    // 
    // *********************************************************
    SQterm
    SQportal::readTensor(string myTensor) const
    {

      cout << " >>----------------------------------------------------<< " << endl;
      cout << " >> readTensor is called                               << " << endl;
      cout << " >>----------------------------------------------------<< " << endl << endl;
    
      string tensorBody("");
      { // Begin elimination
	cout << " >> Eliminating the redundant white spaces ~~~~~~~ ";
	for(auto i = myTensor.begin();i != myTensor.end();++i) if(*i != ' ') tensorBody += *i;
	cout << "    -->> White space eliminated" << endl << endl;
      } // End scope

      string cutBody("");
      { // Begin to cut the redundant parts in the tensor body
	bool isAt   (false); // Is there "@"
	bool isBegin(false); // Is there "("
	bool isEnd  (false); // Is there ")"
	for(auto i = tensorBody.begin();i != tensorBody.end();++i){
	  // >> In case of "@"
	  if(*i == '@'){
	    if(isAt    == true) { cout << "readTensor: Something is wrong <jfdiw92ue92u>" << endl; abort(); }
	    if(isBegin == true) { cout << "readTensor: Something is wrong <2e280>"        << endl; abort(); }
	    if(isEnd   == true) { cout << "readTensor: Something is wrong <i329r398r>"    << endl; abort(); }
	    isAt = true;
	  } // End if
	  // >> In case of "("
	  if(*i == '('){
	    if(isAt    == false) { cout << "readTensor: Something is wrong <jr20h32>"     << endl; abort(); }
	    if(isBegin == true)  { cout << "readTensor: Something is wrong <r037817>"     << endl; abort(); }
	    if(isEnd   == true)  { cout << "readTensor: Something is wrong <3rr3>"        << endl; abort(); }
	    isBegin = true;
	  } // End if
	  // >> In case of ")"
	  if(*i == ')'){
	    if(isAt    == false) { cout << "readTensor: Something is wrong <17e>"         << endl; abort(); }
	    if(isBegin == false) { cout << "readTensor: Something is wrong <hg3e8g38>"    << endl; abort(); }
	    if(isEnd   == true)  { cout << "readTensor: Something is wrong <3r230>"       << endl; abort(); }
	    isEnd = true;
	  } // End if

	  /////////////////////////////////////////////////////////////
	  cutBody += *i;
	  if(isAt && isBegin && isEnd) break;
	  /////////////////////////////////////////////////////////////

	} // End i
      } // End scope

      { // Check the correctness of the input
	//if(cutBody.find("@", 0) != 0){ cout << "readTensor: Something is wrong <eh1e1hy383>" << endl; abort(); } 
	if(cutBody.find("@", 0) == string::npos){ cout << "readTensor: Something is wrong <eh1e1hy383>" << endl; abort(); } 
      } // End if

      typedef boost::char_separator<char> char_separator;
      typedef boost::tokenizer<char_separator> tokenizer;

      SQcont<SQindex> myIndices;
      // Read the associated indices
      { // Begin scope
	string copied(cutBody);
	size_t myBegin(copied.find("(", 0));
	string myCont(copied.substr(myBegin+1));
	myCont.erase(--myCont.end());
	cout << " -- My contraction >> " << myCont << endl;
	
	char_separator my_sep(",", "", boost::keep_empty_tokens);
	tokenizer myTokens(myCont, my_sep);
	  
	// Read the indices
	for (tokenizer::iterator tok_iter = myTokens.begin(); tok_iter != myTokens.end(); ++tok_iter) {
	  myIndices <= readIndex(*tok_iter);
	  cout << " --    >> ind " << " --- "<< *tok_iter << " --- " << endl;
	  cout << " --    >> body " << myIndices.back() << endl; 
	} // End tok_iter
	
      } // End scope

      // Prepare the indices_ptr to create the SQtensor object
      vector<SQindex*> myIndices_ptr;
      for(auto i = myIndices.begin();i != myIndices.end();++i) myIndices_ptr.push_back(&(*i));

      vector<SQtensor> myBodies;
      { 
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	//--*                                            Create SQtensor                                                     *--//
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	string copied(cutBody);
	size_t myBegin(cutBody.find("@", 0));
	size_t myEnd(cutBody.find("(", 0));
	string name(copied.substr(myBegin, calc_distance(myBegin, myEnd)+1));
	cout << " -- My name is : " << name << endl;
	
	// In case of the contraction is given in the Dirac notation
	if(myNotation_ == Femto::Dirac){
	  if     (  t_t2Amps().count(name)) myBodies.push_back(Femto::Core::SQtensor("T2",  myIndices_ptr, Femto::t2_symm()                     ));
	  else if(  t_s2Amps().count(name)) myBodies.push_back(Femto::Core::SQtensor("S2",  myIndices_ptr, Femto::t2_symm()                     ));
	  else if(      t_h1().count(name)) myBodies.push_back(Femto::Core::SQtensor("h" ,  myIndices_ptr, Femto::h1_symm()                     ));
	  else if(      t_h2().count(name)) myBodies.push_back(Femto::Core::SQtensor("V2",  myIndices_ptr, Femto::h2_symm()                     ));

	  else if(      t_C5().count(name)) myBodies.push_back(Femto::Core::SQtensor("C5",  myIndices_ptr, Femto::uni_symm(myIndices_ptr.size()))); // Symmetry is disabled
	  else if(      t_C4().count(name)) myBodies.push_back(Femto::Core::SQtensor("C4",  myIndices_ptr, Femto::uni_symm(myIndices_ptr.size()))); // Symmetry is disabled
	  else if(      t_C6().count(name)) myBodies.push_back(Femto::Core::SQtensor("C6",  myIndices_ptr, Femto::uni_symm(myIndices_ptr.size()))); // Symmetry is disabled

	  else if(    t_JInt().count(name)) myBodies.push_back(Femto::Core::SQtensor("h4_int",  myIndices_ptr, Femto::uni_symm(myIndices_ptr.size()))); // Symmetry is disabled
	  else if(    t_KInt().count(name)) myBodies.push_back(Femto::Core::SQtensor("h5_int",  myIndices_ptr, Femto::uni_symm(myIndices_ptr.size()))); // Symmetry is disabled

	  else if(    t_Rdms().count(name)) myBodies.push_back(     Femto::Core::RDM(       myIndices_ptr                                       ));	
	  else if(     t_Cum().count(name)) myBodies.push_back(             Cumulant(       myIndices_ptr                                       ));
	  else if( t_CasFock().count(name)) myBodies.push_back(Femto::Core::SQtensor("P1",  myIndices_ptr, Femto::h1_symm()                     ));
	  else if(t_CoreFock().count(name)) myBodies.push_back(Femto::Core::SQtensor("Fc1", myIndices_ptr, Femto::h1_symm()                     ));
	  else if(ltensorName_ == name){
	    string temp(ltensorName_);
	    temp.erase(temp.begin());
	    myBodies.push_back(Femto::Core::SQtensor(temp, myIndices_ptr, Femto::uni_symm(myIndices_ptr.size())));
	  } // End if
	  else if(     is_Interm(name)){
	    name.erase(name.begin());
	    myBodies.push_back(Femto::Core::SQtensor(name, myIndices_ptr, Femto::uni_symm(myIndices_ptr.size())));
	  } // End if
	  else{
	    cout << "readContraction: Can't read this >> " << cutBody << endl;
	    abort();
	  } // Else
	} // End if(notation)
	// In case of the contraction is given in the Mulliken notation
	if(myNotation_ == Femto::Mulliken){
	  if     (  t_t2Amps().count(name)) myBodies.push_back(Femto::Core::SQtensor ("T2",  myIndices_ptr, Femto::t2_symm()));
	  else if(  t_s2Amps().count(name)) myBodies.push_back(Femto::Core::SQtensor ("S2",  myIndices_ptr, Femto::t2_symm()));
	  else if(      t_h1().count(name)) myBodies.push_back(Femto::Portal::Mtensor("h" ,  myIndices_ptr, Femto::h1_symm()));
	  else if(      t_h2().count(name)) myBodies.push_back(Femto::Portal::Mtensor("V2",  myIndices_ptr, Femto::h2_symm()));

	  else if(      t_C5().count(name)) myBodies.push_back(Femto::Core::SQtensor ("C5",  myIndices_ptr, Femto::uni_symm(myIndices_ptr.size()))); // Symmetry is disabled
	  else if(      t_C4().count(name)) myBodies.push_back(Femto::Core::SQtensor ("C4",  myIndices_ptr, Femto::uni_symm(myIndices_ptr.size()))); // Symmetry is disabled
	  else if(      t_C6().count(name)) myBodies.push_back(Femto::Core::SQtensor ("C6",  myIndices_ptr, Femto::uni_symm(myIndices_ptr.size()))); // Symmetry is disabled

	  else if(    t_JInt().count(name)) myBodies.push_back(Femto::Core::SQtensor("h4_int",  myIndices_ptr, Femto::uni_symm(myIndices_ptr.size()))); // Symmetry is disabled
	  else if(    t_KInt().count(name)) myBodies.push_back(Femto::Core::SQtensor("h5_int",  myIndices_ptr, Femto::uni_symm(myIndices_ptr.size()))); // Symmetry is disabled

	  else if(    t_Rdms().count(name)) myBodies.push_back(Femto::Portal::MRDM   (       myIndices_ptr                  )); 
	  else if(     t_Cum().count(name)) myBodies.push_back(             MCumulant(       myIndices_ptr                  ));
	  else if( t_CasFock().count(name)) myBodies.push_back(Femto::Portal::Mtensor("P1",  myIndices_ptr, Femto::h1_symm()));
	  else if(t_CoreFock().count(name)) myBodies.push_back(Femto::Portal::Mtensor("Fc1", myIndices_ptr, Femto::h1_symm()));
	  else if(ltensorName_ == name){
	    string temp(ltensorName_);
	    temp.erase(temp.begin());
	    myBodies.push_back(Femto::Core::SQtensor(temp, myIndices_ptr, Femto::uni_symm(myIndices_ptr.size())));
	  } // End if
	  else if(     is_Interm(name)){
	    name.erase(name.begin());
	    myBodies.push_back(Femto::Core::SQtensor(name, myIndices_ptr, Femto::uni_symm(myIndices_ptr.size())));
	  } // End if
	  else{
	    cout << "readTensor: Can't read this >> " << cutBody << endl;
	    abort();
	  } // Else
	} // End if(notation)
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
      } // End scope

      vector<string> coeff; coeff.push_back("");
      return SQterm(1.0, coeff, myBodies);

    }


     // *********************************************************
     // 
     // *********************************************************
     SQbinary 
     SQportal::readContractionNew2(string myContraction) const
     {

       cout << " >>----------------------------------------------------<< " << endl;
       cout << " >> readContractionNew2 is called                      << " << endl;
       cout << " >>----------------------------------------------------<< " << endl;

       pair<int, vector<size_t> > numTerms; // <number of terms, position of '@' in myContraction> //
       numTerms.first = 0;
       for(auto i = myContraction.begin();i != myContraction.end();++i) 
	 if(*i == '@') 
	   {
	     numTerms.second.push_back((size_t)(i-myContraction.begin()));
	     ++numTerms.first;
	   } // End if
       cout << " >> Number of tensors detected ... " << numTerms.first << endl;

       SQcont<string> myTensors;
       for(int n = 0;n < numTerms.first;++n){
	 string temp("");
	 bool foundEnd(false);
	 for(size_t pos = numTerms.second[n];pos < myContraction.length();++pos){
	   if(myContraction[pos] == ')') foundEnd = true;
	   temp += myContraction[pos];
	   if(foundEnd) break;
	 } // End pos
	 if(!foundEnd){ cout << "readContractionsNew2: Somwthingis wrong <2d62e12>" << endl; abort(); } // End if
	 myTensors <= temp;
       } // End n

       SQcont<SQterm> tensorData;
       for(auto t = myTensors.begin();t != myTensors.end();++t)
	 tensorData <= readTensor(*t);

       double myCoeff(0.0);
       SQcont<string> literals;
       { // Read the numerical and literal coefficient
	 
	 typedef boost::char_separator<char> char_separator;
	 typedef boost::tokenizer<char_separator> tokenizer;
	 
	 //SQcont<string> myTensors;
	 string myFactor;
	 string myLiteral;
	 char_separator sep(" ", "", boost::keep_empty_tokens);
	 tokenizer tokens(myContraction, sep);
	 
	 // Extract components of the contraction
	 bool foundArrow(false);
	 SQcont<string> myLiterals;
	 for (tokenizer::iterator tok_iter = tokens.begin(); tok_iter != tokens.end(); ++tok_iter) {
	   if("<<=" == *tok_iter){ foundArrow = true; }
	   if(foundArrow){
	     if(tok_iter->find('@', 0) != string::npos) break;
	     if("<<=" != *tok_iter) myLiterals <= *tok_iter;
	   } // End if
	 } // End tok_iter
	 if(!foundArrow) { cout << "readContractionsNew2: \"<<==\" can't be found" << endl; abort(); }

	 // Separate the numerical and lietral coefficients
	 for(auto s = myLiterals.begin();s != myLiterals.end();++s){
	   double myVal(atof(s->c_str()));
	   if(myVal != 0.0) myCoeff   = myVal;
	   else if(c_Ecore().count(*s)) literals <= "Fc0";
	   else if(c_Ecas().count(*s))  literals <= "h6_int";
	   else                         literals <= *s;
	 } // End s
	   
       } // End scope

       ostringstream printTensors;
       for(auto t = tensorData.begin();t != tensorData.end();++t) printTensors << t->get_tensors()[0];

       string literalCoeff("");
       for(auto c = literals.begin();c != literals.end();++c){
	 string temp(*c); temp += "  ";
	 literalCoeff += temp;
       } // End c

       cout << endl;
       cout << " >>~~~~~~ summary of the cotraction ~~~~~~<< " << endl;
       cout << " >> Numerical factor : " << (boost::format("%10.8f") % myCoeff) << endl;
       cout << " >> Literal factor   : " << (literals.size() ? literalCoeff : " None ") << endl;
       cout << " >> Tensor data      : " << printTensors.str()  << endl;
       cout << " >>~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~<<" << endl << endl;

       SQtensor         ltensor;
       vector<SQtensor> tensorBody;
       for(auto t = tensorData.begin();t != tensorData.end();++t) 
	 if(t == tensorData.begin()) ltensor = t->get_tensors()[0]; 
	 else tensorBody.push_back(t->get_tensors()[0]);

       return SQbinary(myCoeff, literals.p(), ltensor, tensorBody);

     }
  
}} //Femto::portal
