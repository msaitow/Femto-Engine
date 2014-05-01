//
//  SQportal.cc
//  
//
//  Created by Masaaki Saitow on 13/10/17.
//  Copyright (c) 2012 Masaaki Saitow. All rights reserved.
//

#include <SQportal.hpp>
#include <boost/tokenizer.hpp>

#define _DEBUG1

using namespace std;

using namespace Femto;
using namespace Femto::Core;
using namespace Femto::Reaktor;
//using namespace Femto::Portal;

namespace Femto { namespace Portal {

    
    // *********************************************************
    // 
    // *********************************************************
    SQportal::SQportal(const string &file_name, const string &ltensorName, Femto::notation myNotation)
      : Faf_name_(file_name),
	ltensorName_(ltensorName),
	myNotation_(myNotation)
    {
      if(myNotation != Femto::Dirac && myNotation != Femto::Mulliken)
	{ cout << "SQportal: Notation unkown" << endl; abort(); }
    }
    
    // *********************************************************
    // 
    // *********************************************************
    void 
    SQportal::transportContractions(vector<vector<SQbinary> > &theBins, Femto::Portal::version myVersion)
    {
      ifstream ifs(Faf_name_.c_str(), ios::in );
      if(ifs.fail()){
	cout << "transportContractions: Can't open file that contains the information of header files >> " << Faf_name_ << " << " << endl;
	//return;    // Toriaezu
	abort(); // Toriaezu
      } // End if
      
      ///////////////////////////////////
      // Symbol for the instruction
      const string instKey ("--|");
      const string beginKey("Begin");
      const string endKey  ("End");
      const string bigKey  ("AllContras");
      const string smallKey("Contras");
      const string titleKey("Title");
      ///////////////////////////////////
      
      string line;
      int numContractions(0);      // Number of the contraction terms
      int numSmallContractions(0); // Number of the binary contractions for each term
      //int totalCount(0);
      //int smallCount(0);
      
      bool nowInAllContra(false);
      bool nowInContra(false);
      
      vector<SQbinary> batch;
      //while(!ifs.eof()){
      //  ifs >> line;
      while(getline(ifs, line)){
	cout << "MyLine: " << line << endl;
	size_t myBegin = line.find(instKey, 0);
	if(myBegin != string::npos){
	  
	  // ==> In case of some instructions
	  string instruction(line.substr(myBegin));
	  size_t begin = instruction.find(beginKey, 0);
	  size_t end   = instruction.find(endKey, 0);
	  size_t title = instruction.find(titleKey, 0);
	  
	  // ==>> In case of begin
	  if     (begin != string::npos && end == string::npos && title == string::npos){
	    string contents(instruction.substr(begin+beginKey.size()+1));
	    size_t bigBegin   = contents.find(bigKey, 0);
	    size_t smallBegin = contents.find(smallKey, 0);
	    
	    // ==>> In case of Begin.AllContras
	    if     (bigBegin != string::npos && smallBegin != string::npos){
	      // Read the number of the contraction terms
	      size_t pos_forward (contents.find("(", 0));
	      size_t pos_backward(contents.find(")", 0));
	      size_t distance(calc_distance(pos_forward, pos_backward));
	      std::string myNum(contents.substr(pos_forward+1, distance));
	      if(numContractions != 0) { 
		cout << "transportContractions: Algorithmic Error[1]" << endl; 
		cout << "    At: " << line << endl;
		cout << "        Where: " << contents << endl;
		abort(); 
	      } // End if
	      else{
		// Set counter for the Begin.AllContras
		numContractions = atoi(myNum.c_str());
		nowInAllContra = true;
		theBins.reserve(numContractions);
	      } // End else
	    } // End if
	    
	    // ==>> In case of Begin.Contras
	    else if(bigBegin == string::npos && smallBegin != string::npos){
	      // Read the number of the binary contractions
	      size_t pos_forward (contents.find("(", 0));
	      size_t pos_backward(contents.find(")", 0));
	      size_t distance(calc_distance(pos_forward, pos_backward));
	      std::string myNum(contents.substr(pos_forward+1, distance));
	      if(numSmallContractions != 0 || batch.size() != 0) { 
		cout << "transportContractions: Algorithmic Error[2]" << endl; 
		cout << "    At: " << line << endl;
		cout << "        Where: " << contents << endl;
		abort(); 
	      } // End if
	      else{
		// Set counter for the Begin.Contras
		numSmallContractions = atoi(myNum.c_str());
		nowInContra          = true;
		batch.reserve(numSmallContractions);
	      } // End else
	    } // End if
	    else {
	      cout << "transportContractions: Algorithmic Error[3]" << endl;
	      cout << "    At: " << line << endl;
	      cout << "        Where: " << contents << endl;
	      abort();
	    } // End else
	    
	  } // End if (begin)
	  
	  // ==>> In case of end
	  else if(begin == string::npos && end != string::npos && title == string::npos){
	    string contents(instruction.substr(begin+beginKey.size()+1));
	    size_t bigBegin   = contents.find(bigKey, 0);
	    size_t smallBegin = contents.find(smallKey, 0);
	    // ==>> In case of End.AllContras
	    if     (bigBegin != string::npos && smallBegin != string::npos){
	      if (theBins.size() != numContractions || !nowInAllContra){ 
		cout << "transportContractions: Number of contraction terms mismatched " << endl; 
		cout << "    At: " << line << endl;
		abort(); 
	      } // End if
	      else {
		nowInAllContra  = false;
		numContractions = 0;
		//totalCount      = 0;
	      } // End else
	    } // End if	  
	    // ==>> In case of End.Contras
	    else if(bigBegin == string::npos && smallBegin != string::npos){
	      if (batch.size() != numSmallContractions || !nowInContra){ 
		cout << "transportContractions: Number of binary contractions mismatched " << endl; 
		cout << "    At: " << line << endl;	      
		abort(); 
	      } // End if
	      else {
		nowInContra          = false;
		numSmallContractions = 0;
		//smallCount           = 0;
		theBins.push_back(batch);
		batch.clear();
	      } // End else
	      
	    } // End if
	    else {
	      cout << "transportContractions: Algorithmic Error[4]"  << endl;
	      cout << "    At: " << line << endl;
	      abort();
	    } // End else
	    
	  } // End if (end)
	  else if(title == string::npos){
	    cout << "transportContractions: Algorithmic Error[5]" << endl;
	    cout << "    At: " << line << endl;
	    abort();
	  } // End else
	  
	} // End if
	
	// ---------------------------------- //
	// -- Read the binary contractions -- //
	// ---------------------------------- //
	else if(nowInAllContra && nowInContra){
	  size_t myBegin(line.find("@", 0));
	  string myContraction(line.substr(myBegin));
	  cout << " >> Reading ... " << myContraction << endl;
	  if      (myVersion == Femto::Portal::OLD)
	    batch.push_back(readContraction(myContraction));     // Reading [old]
	  else if (myVersion == Femto::Portal::NEW)
	    batch.push_back(readContractionNew(myContraction));  // Reading [new]
	  else if (myVersion == Femto::Portal::NEW2)
	    batch.push_back(readContractionNew2(myContraction)); // Reading [new2]
	  //{ cout << "Not yet implemented <3r93rh>" << endl; abort(); }
	  else {
	    cout <<  "transportContractions: Can't handle this type" << endl;
	    abort();
	  } // End else
	} // End if
	
      } // End while
      
    }
    
    
    // *********************************************************
    // 
    // *********************************************************
    SQbinary 
    SQportal::readContraction(string myContraction) const
    {
      SQcont<SQcont<SQindex> > myIndices;
      SQcont<SQtensor> myBodies; 
      SQcont<string> myTensors; 
      
      typedef boost::char_separator<char> char_separator;
      typedef boost::tokenizer<char_separator> tokenizer;
      
      string myFactor;
      char_separator sep(" ", "", boost::keep_empty_tokens);
      tokenizer tokens(myContraction, sep);
      
      // Extract components of the contraction
      int myNum(0);
      for (tokenizer::iterator tok_iter = tokens.begin(); tok_iter != tokens.end(); ++tok_iter, ++myNum) 
	if(myNum == 0 || myNum == 3 || myNum == 4)
	  myTensors <= *tok_iter;
	else if(myNum == 2)
	  myFactor = *tok_iter;
      
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
	rtensors.push_back(myBodies[1]);
	rtensors.push_back(myBodies[2]);
	double factor(atof(myFactor.c_str()));
	vector<string> empty;
	SQbinary myTerm(factor, empty, ltensor, rtensors);
	
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
    SQindex 
    SQportal::readIndex(string myIndex) const
    {
      char_state myState;
      bool isDummy;
      
      // >> Core and non-dummy
      if     (coreNonDummy().count(myIndex))
	{ myState = Femto::core; isDummy = false; }
      // >> Active and non-dummy
      else if( actNonDummy().count(myIndex))
	{ myState = Femto::act; isDummy = false; }
      // >> Virtual and non-dummy
      else if(virtNonDummy().count(myIndex))
	{ myState = Femto::virt; isDummy = false; }
      
      // >> Core and dummy
      else if(coreDummy().count(myIndex))
	{ myState = Femto::core; isDummy = true; }
      // >> Active and dummy
      else if( actDummy().count(myIndex))
	{ myState = Femto::act; isDummy = true; }
      // >> Virtual and dummy
      else if(virtDummy().count(myIndex))
	{ myState = Femto::virt; isDummy = true; }
      else{
	cout << "readIndex: Can't handle this >> " << myIndex << " << "<< endl; abort();
      } // End else 
      
      return (SQindex(myIndex, myState, isDummy, false));
    }
    
    // *********************************************************
    // 
    // *********************************************************
    void 
    SQportal::transformD2M_INT_RDM(vector<vector<SQbinary> > &theBins, bool doMask)
    {
      for(auto ts = theBins.begin();ts != theBins.end();++ts)
	for(auto t = ts->begin();t != ts->end();++t){
	  if (doMask) t->masquerade();
	  vector<SQtensor*> temp(t->get_Rtensors_ptr());
	  for(auto t2 = temp.begin();t2 != temp.end();++t2)
	    if((*t2)->get_name() == name_h2() || (*t2)->get_name() == name_h1()  || is_RDM((*t2)->get_name()) || 
	       (*t2)->get_name() == C2_name() || (*t2)->get_name() == name_Fock()) (*t2)->convertD2M(); 
	} // End t
    }
    
}} //Femto::portal
