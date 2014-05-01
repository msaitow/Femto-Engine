//
//  GenUtils.cc
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

using namespace std;
using namespace Femto;
using namespace Femto::Core;

namespace Femto { namespace Reaktor {


  // *********************************************************
  // 
  // *********************************************************
  void
  SQreaktor::declareInterm(int myIndent, ofstream &CPfile, SQtensor ten, SQcont<SQindex*> decInds, decType myType)
  {
    string Xname   ("");
    string SymLabel("");
    string indents ("");
    SQcont<SQindex> extInds;
    if(myType == External){
      for(size_t num = 0;num < ten.get_indices().size();++num)
	if(decInds.count(ten.get_indices()[num]) == 0) extInds <= *ten.get_indices()[num];
 
    } // End if
    else {
      for(size_t num = 0;num < ten.get_indices().size();++num)
	if (ten.get_indices()[num]->get_isExt()) extInds <= *ten.get_indices()[num];
    } // End else

    if(extInds.size()){
      for(auto i = extInds.begin();i != extInds.end();++i){
	SymLabel += "s" + i->get_index();
	if(*i != extInds.back()) SymLabel += "^";
      } // End if
    } // End if
    else SymLabel = "0";

    for(int c = 0;c < myIndent;++c) indents += "  ";

    int Ccount(0);
    int Ocount(0);
    int Vcount(0);
    for(auto i = decInds.begin();i != decInds.end();++i)
      if     ((*i)->get_char() == Femto::core) ++Ccount;
      else if((*i)->get_char() == Femto::act ) ++Ocount;
      else if((*i)->get_char() == Femto::virt) ++Vcount;

    if(Ccount || Ocount || Vcount){
      for(int i = 0;i < Ccount;++i) Xname += "c";
      for(int i = 0;i < Ocount;++i) Xname += "a";
      for(int i = 0;i < Vcount;++i) Xname += "v";

      string DecInterm(indents + "  orz::DTensor " + ten.get_name() + Xname + "_" + title_);
      DecInterm += "(orz::mr::sizeof_sympack_X" + Xname + "(" + "symblockinfo, " + SymLabel + "));";
      CPfile << DecInterm << endl;
    } // End if
    else{
      string DecInterm(indents + "  double " + ten.get_name() + "_" + title_ + "(0);");
      CPfile << DecInterm << endl; 
    } // End else

  }


  // *********************************************************
  // 
  // *********************************************************
  SQindex
  SQreaktor::findERIIndex(vector<SQbinary> &contras) const
  {
    SQindex retVal;
    for(auto t = contras.cbegin();t != contras.cend();++t){
      auto tensors(t->get_Rtensors());
      for(auto ten = tensors.cbegin();ten != tensors.cend();++ten)
	if(ten->get_name() == name_h2_) retVal = *ten->get_indices()[exth2_];
    } // End t
    return retVal;
  }


  // *********************************************************
  // 
  // *********************************************************
  SQindex
  SQreaktor::findD4CIndex(vector<SQbinary> &contras) const
  {
    SQindex retVal;
    for(auto t = contras.cbegin();t != contras.cend();++t){
      auto tensors(t->get_Rtensors());
      for(auto ten = tensors.cbegin();ten != tensors.cend();++ten)
	if(is_D4C(ten->get_name())) retVal = *ten->get_indices()[extd4c_];
    } // End t
    return retVal;
  }


  // *********************************************************
  // 
  // *********************************************************
  SQcont<SQindex>
  SQreaktor::returnsExtIndices(SQbinary &contra) const
  {
    SQcont<SQindex> retVal;
    auto inds(contra.get_summedBody());
    for(auto i = inds.begin();i != inds.end();++i) if ((*i)->get_isExt()) retVal <= **i;
    return eliminate_duplicate(retVal);
  }


  // *********************************************************
  // 
  // *********************************************************
  void
  SQreaktor::makeContractions(int myIndent, ofstream &CPfile, ofstream &CHfile, ofstream &F90file, string myLabel, string theLabel, vector<SQbinary> &myBins, Femto::Reaktor::ParaFlag myFlag)
  {

    int depthScope(myIndent);
    string indent("");
    SQcont<string> Indents;
    for(int i = 0;i < 100;++i){
      Indents <= indent;
      indent += "  ";
    } // End i 

    //////////////////////////////////
    SQindex* ERIInd(NULL);
    SQindex* AMPInd(NULL);
    SQindex* SIGInd(NULL); 
    SQtensor Amp;
    SQtensor Sig;
    int posAmp(-1); 
    SQcont<SQindex> alreadyLoaded;
    SQcont<SQindex>  extIndices;
    SQcont<SQtensor> loadedInterms;
    //////////////////////////////////
    for(auto c = myBins.begin(); c != myBins.end();++c){
      auto tempR(c->get_Rtensors_ptr());
      for(auto t = tempR.begin();t != tempR.end();++t){
	if((*t)->get_name() == name_amp_){
	  Amp         = **t;
	  extIndices <= *(*t)->get_indices()[extamp_];
	  AMPInd      =  (*t)->get_indices()[extamp_];
	  posAmp      = (int)(c-myBins.begin()); // Position of the contraction that has the T2-amp

	  if(!AMPInd->get_isExt()){
	    cout << "makeContractions: Loading index of Amp isn't set as the external one >> " << *AMPInd << " << " << endl;
	    cout << SQcont<SQbinary>(myBins) << endl;
	    CPfile.flush();
	    abort();
	  } // End if

	} // End if
	if ((*t)->get_name() == name_h2_) ERIInd = &(*(*t)->get_indices()[exth2_ ]);
	if (is_D4C((*t)->get_name())    ) ERIInd = &(*(*t)->get_indices()[extd4c_]);
      } // End if
      if(c->get_Ltensor().get_name() == LTensor_.get_name() && isBareLHS_) {
	Sig         = c->get_Ltensor();
	extIndices <= *(c->get_Ltensor().get_indices()[extamp_]);
	SIGInd      =   c->get_Ltensor().get_indices()[extamp_];

	if(!SIGInd->get_isExt()){
	  cout << "makeContractions: Loading index of Sigma isn't set as the external one >> " << *SIGInd << " << " << endl;
	  cout << SQcont<SQbinary>(myBins) << endl;
	  CPfile.flush();
	  abort();
	} // End if

      } // End if
    } // End c

///////    ////////////////////////////////////////////////////
///////    // Index specifications
///////    ////////////////////////////////////////////////////
///////    if(ERIInd != NULL) alreadyLoaded <= *ERIInd;
///////    if(AMPInd != NULL) if(alreadyLoaded.count(*AMPInd)) ReadAmp(Indents[depthScope], Amp, CPfile);
///////    if(SIGInd != NULL) if(alreadyLoaded.count(*SIGInd)) ReadAmp(Indents[depthScope], Sig, CPfile);

    // Determine which intermediate needs special cares
    SQcont<SQcont<SQindex> > outerIndices; // vector<OuterIndices>;
    SQcont<pair<int, int> >  LandR;        // Numbers of intermediate that appears on the LHS and RHS
    for(auto contra = myBins.begin();contra != myBins.end();++contra){
      int numL(contra-myBins.begin());
      int numR;
      auto myLtensor(contra->get_Ltensor());
      SQbinary *thisOne(NULL);
      // Find where the contraction that has myLtensor on the RHS is
      for(auto contra2 = contra; contra2 < myBins.end();++contra2){
	auto Rtensors(contra2->get_Rtensors());
	for(auto t = Rtensors.begin();t != Rtensors.end();++t) 
	  if(t->get_name() == myLtensor.get_name()) {
	    thisOne = &(*contra2);
	    numR = (int)(contra2-myBins.begin());
	    break; 
	  } // End 
      } // End contra2
      if(thisOne != NULL){
	auto LInds(returnsExtIndices(*contra ));
	auto RInds(returnsExtIndices(*thisOne));
	outerIndices <= make_intersection(LInds, RInds); // Only the intersections have to be the outer indices
	LandR        <= make_pair(numL, numR);
      } // End if
      else{
	SQcont<SQindex> temp;
	outerIndices <= temp;
	LandR        <= make_pair(-2, -1);
      } // End else
    } // End if
    if(outerIndices.size() != myBins.size()) { cout << "makeContractions: Number of outerIndices mismatched." << endl; abort(); }
    if(LandR.size()        != myBins.size()) { cout << "makeContractions: Number of LandR mismatched." << endl; abort(); }

    // Now we've got relationship between the indices and contractions
    size_t numOuter(0);
    for(auto i = outerIndices.begin();i != outerIndices.end();++i,++numOuter){

      // ---------------------------------------------------------------------------------------
      // [1] In case that the external indices associated with the contractions are different
      // ---------------------------------------------------------------------------------------
      if(i->size()){
	myBins[numOuter].set_LOuterIndices(*i);
	{ // Verify whether all the indices are external
	  auto lIndices(myBins[numOuter].get_summedBody());
	  SQcont<SQindex> extInds;
	  for(auto j = lIndices.begin();j != lIndices.end();++j) if((*j)->get_isExt()) extInds <= **j;
	  if(extInds.size() != i->size()) {
	    CPfile << "// ==> makeContractions: [1] Special care is taken for the contraction " << numOuter << " <" << extInds.size() << ", " << i->size() << "> " << endl; //*TEST*//
	  } // End if
	} // End verification
	auto Rtensor(myBins[LandR[numOuter].second].get_Rtensors());
	int numTen(-1);
	for(auto n = 0;n < Rtensor.size();++n) if(myBins[numOuter].get_Ltensor().get_name() == Rtensor[n].get_name()) numTen = n;
	if(numTen == -1) { 
	  cout << "makeContractions: Algorithmic Error <12353> \n" << SQcont<SQbinary>(myBins) << endl; 
	  abort();
	} // End if
	myBins[LandR[numOuter].second].set_ROuterIndices(numTen, *i);
      } // End if

      // ---------------------------------------------------------------------------------------
      // [2] In case that the L- and R-contractions are not contiguous
      // ---------------------------------------------------------------------------------------
      if(LandR[numOuter].first+1 != LandR[numOuter].second){
	if(LandR[numOuter].second - LandR[numOuter].first != 2) { cout << "makeContractions: Not implemented yet <47364836>" << endl; abort(); }
	auto lInds(returnsExtIndices(myBins[numOuter  ]));
	auto rInds(returnsExtIndices(myBins[numOuter+1]));
	auto outers(make_intersection(lInds, rInds));
	myBins[numOuter].set_LOuterIndices(outers);
	CPfile << "// ==> makeContractions: [2] L- and R-contractions are not contiguous >> " << numOuter << endl; //*TEST*//
	auto Rtensor(myBins[LandR[numOuter].second].get_Rtensors());
	int numTen(-1);
	for(auto n = 0;n < Rtensor.size();++n) if(myBins[numOuter].get_Ltensor().get_name() == Rtensor[n].get_name()) numTen = n;
	if(numTen == -1) { 
	  cout << "makeContractions: Algorithmic Error <757b398> \n" << SQcont<SQbinary>(myBins) << endl; 
	  abort();
	} //End if
	myBins[LandR[numOuter].second].set_ROuterIndices(numTen, outers);
      } // End if

    } // End i
    ////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////
    // Care about allocation of a batch of the sigma-vector
    /////////////////////////////////////////////////////////
//new//bug?     pair<size_t, bool> allocSigma; // <Position of contraction that has the sigma-vector, where sigma is allocated>
//new//bug?     allocSigma.first  = -1;
//new//bug?     allocSigma.second = false;
//new//bug?     if(myBins.back().get_Ltensor().get_name() == LTensor_.get_name() && isBareLHS_){
//new//bug?       for(auto contra = myBins.rbegin();contra != myBins.rend();++contra){
//new//bug? 	auto currentExtInds(returnsExtIndices(*contra));
//new//bug? 	if(currentExtInds.count(*SIGInd)) allocSigma.first = myBins.size() - (size_t)(contra-myBins.rbegin()) - 1;
//new//bug? 	else break;
//new//bug? 	CPfile << "// |------> " << allocSigma.first << " (alloc)" << endl; //TEST
//new//bug?       } // End contra
//new//bug?     } // End if
//new//bug?     if(allocSigma.first != -1){
//new//bug?       CPfile << "// allocSigma.first = " << allocSigma.first << endl; //TEST
//new//bug?     } // End if 

    pair<size_t, bool> allocSigma; // <Position of contraction that has the sigma-vector, where sigma is allocated>
    allocSigma.first  = -1;
    allocSigma.second = false;
    if(myBins.back().get_Ltensor().get_name() == LTensor_.get_name() && isBareLHS_){
      auto sigExtInds(returnsExtIndices(myBins.back()));
      for(auto contra = myBins.rbegin();contra != myBins.rend();++contra){
	auto currentExtInds(returnsExtIndices(*contra));
	if(currentExtInds == sigExtInds) allocSigma.first = myBins.size() - (size_t)(contra-myBins.rbegin()) - 1;
	else break;
	CPfile << "// |------> " << allocSigma.first << " (alloc)" << endl; //TEST
      } // End contra
    } // End if
    if(allocSigma.first != -1){
      CPfile << "// allocSigma.first = " << allocSigma.first << endl; //TEST
    } // End if 
    /////////////////////////////////////////////////////////    

    /////////////////////////////////////////////////////////
    // Care about the loop order of the loading indices of 
    // the T2-amplitude and sigma-vector
    /////////////////////////////////////////////////////////
    enum Preference {t2first, s2first, none};
    Preference ampPref(none);
    if(AMPInd != NULL && SIGInd != NULL){
      SQcont<size_t> contrasHaveSigInd;
      SQcont<size_t> contrasHaveAmpInd;
      for(auto c = myBins.begin();c != myBins.end();++c){
	auto thisExtIndices(returnsExtIndices(*c));
	if(thisExtIndices.count(*AMPInd)) contrasHaveAmpInd <= (size_t)(c-myBins.begin());
	if(thisExtIndices.count(*SIGInd)) contrasHaveSigInd <= (size_t)(c-myBins.begin());
      } // End c    
      if     (contrasHaveAmpInd.size() > contrasHaveSigInd.size()){
	auto intersec(make_intersection(contrasHaveAmpInd, contrasHaveSigInd));
	sort(intersec.begin(), intersec.end());
	if(intersec == contrasHaveSigInd) ampPref = t2first;
      } // End if
      else if(contrasHaveAmpInd.size() < contrasHaveSigInd.size()){
	auto intersec(make_intersection(contrasHaveAmpInd, contrasHaveSigInd));
	sort(intersec.begin(), intersec.end());
	if(intersec == contrasHaveAmpInd) ampPref = s2first;
      } // End if
      else ampPref = none;
    } // End if
    CPfile << "  // Pref: " << ampPref << endl; //TEST
    /////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////
    // Care about incontiguous external indices
    ////////////////////////////////////////////////////////
    if(ERIInd != NULL && AMPInd != NULL && SIGInd != NULL){
      enum anomality {anormal, normal};
      SQcont<anomality>        incontiguity1;
      SQcont<anomality>        incontiguity2;
      SQcont<SQcont<SQindex> > overlaps;
      SQcont<SQcont<SQindex> > extIndices;
      for(auto c = myBins.begin();c != myBins.end();++c) {
	incontiguity1 <= normal;
	incontiguity2 <= normal;
	extIndices   <= returnsExtIndices(*c);
      } // End c
      for(size_t num = 0;num < extIndices.size();++num) sort(extIndices[num].begin(), extIndices[num].end());
      { // -- Model 1 -- //
	SQcont<SQcont<SQindex> > model1;
	model1 <= SQcont<SQindex>(*ERIInd, *AMPInd);
	model1 <= SQcont<SQindex>(*ERIInd, *AMPInd, *SIGInd);
	model1 <= SQcont<SQindex>(*ERIInd, *SIGInd);
	sort(model1[0].begin(), model1[0].end());
	sort(model1[1].begin(), model1[1].end());
	sort(model1[2].begin(), model1[2].end());
	for(auto i = extIndices.begin();i != extIndices.end();++i){
	  size_t myPos(i-extIndices.begin());
	  if(myPos < extIndices.size()-1){
	    if(*i == model1[0]) {
	      auto next1(i); ++next1;
	      auto next2(i); ++next2; ++next2;	      
	      if(*next1 == model1[1] && *next2 == model1[2]) incontiguity1[myPos] = anormal;
	    } // End if
	  } // End if
	} // End i
      } // End scope
      { // -- Model 2 -- //
	SQcont<SQcont<SQindex> > model2;
	model2 <= SQcont<SQindex>(*ERIInd, *SIGInd);
	model2 <= SQcont<SQindex>(*ERIInd, *AMPInd, *SIGInd);
	model2 <= SQcont<SQindex>(*ERIInd, *AMPInd);
	sort(model2[0].begin(), model2[0].end());
	sort(model2[1].begin(), model2[1].end());
	sort(model2[2].begin(), model2[2].end());
	for(auto i = extIndices.begin();i != extIndices.end();++i){
	  size_t myPos(i-extIndices.begin());
	  if(myPos < extIndices.size()-1){
	    if(*i == model2[0]) {
	      auto next1(i); ++next1;
	      auto next2(i); ++next2; ++next2;	      
	      if(*next1 == model2[1] && *next2 == model2[2]) incontiguity2[myPos] = anormal;
	    } // End if
	  } // End if
	} // End i
      } // End scope
      for(auto c = incontiguity1.begin(); c != incontiguity1.end();++c)
        if(*c == anormal) {
	  size_t myPos(c-incontiguity1.begin());
	  CPfile << "  // Type1 anormality found for external indices of the contraction >> " << (c-incontiguity1.begin()) << " <<" << endl;
	  if(LandR[myPos].first+1 != LandR[myPos].second){
	    cout << "makeContractions: Algorithmic Error <8723hfqu3> \n" << SQcont<SQbinary>(myBins) << endl;
	    abort();
	  } // End if
	  // Set the loading index of the t2-amplitude as outer
	  myBins[LandR[myPos].first].set_LOuterIndices(SQcont<SQindex>(*ERIInd));
	  // Make correspondence to the Rtensor
	  auto Rtensor(myBins[LandR[myPos].second].get_Rtensors());
	  int numTen(-1);
	  for(auto n = 0;n < Rtensor.size();++n) if(myBins[myPos].get_Ltensor().get_name() == Rtensor[n].get_name()) numTen = n;
	  if(numTen == -1) { 
	    cout << "makeContractions: Algorithmic Error <uh39h4t3h4> \n" << SQcont<SQbinary>(myBins) << endl; 
	    abort();
	  } //End if
	  myBins[LandR[myPos].second].set_ROuterIndices(numTen, SQcont<SQindex>(*ERIInd));
	  // Set the laoding index of the t2-amplitude as internal in the intermediate contraction
	  auto dummies(myBins[LandR[myPos].second].get_summedBody());
	  for(auto i = dummies.begin();i != dummies.end();++i) if(**i == *AMPInd) (*i)->switch_isExt(false); 
	} // End if
      for(auto c = incontiguity2.begin(); c != incontiguity2.end();++c)
        if(*c == anormal) {
	  CPfile << "  // Type2 anormality found for external indices of the contraction >> " << (c-incontiguity2.begin()) << " <<" << endl;
	  abort();
	} // End if
    } // End scope
    ////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////
    // Index specifications
    ////////////////////////////////////////////////////
    if(ERIInd != NULL) alreadyLoaded <= *ERIInd;
    if(AMPInd != NULL) if(alreadyLoaded.count(*AMPInd)) ReadAmp(Indents[depthScope], Amp, CPfile);
    if(SIGInd != NULL) if(alreadyLoaded.count(*SIGInd)) { ReadAmp(Indents[depthScope], Sig, CPfile); allocSigma.second = true; } 

    // Let's begin to generate the code!
    for(auto contra = myBins.begin();contra != myBins.end();++contra){

      ///////////////////////////////////////////////////////
      // Allocation of the S2-tensor [1]
      ///////////////////////////////////////////////////////
      if( allocSigma.first == (size_t)(contra-myBins.begin()) && alreadyLoaded.count(*SIGInd) && !allocSigma.second) {
	ReadAmp(Indents[depthScope], Sig, CPfile);
	allocSigma.second = true;
      } // End if

      ////////////////////////////////////////////////////////
      // Allocate the intermediate of the external type [1]
      ////////////////////////////////////////////////////////
      for(auto contra2 = myBins.begin(); contra2 < myBins.end();++contra2){
	bool allLoaded(true);
	auto outers(contra2->get_LOuterIndices());
	for(auto i = outers.begin();i != outers.end();++i) if(alreadyLoaded.count(**i) == 0) allLoaded = false;
	if(!contra2->get_Lindices() && contra2->get_Ltensor().get_name() != LTensor_.get_name() && allLoaded && loadedInterms.count(contra2->get_Ltensor()) == 0){ 
	  declareInterm(depthScope, CPfile, contra2->get_Ltensor(), SQcont<SQindex*>(contra2->get_LInnerIndices()), Femto::Reaktor::External);
	  loadedInterms <= contra2->get_Ltensor();
	} // End if
      } // End scope
      ////////////////////////////////////////////////////////

      //////////////////////////////////////////////
      // Deepen the depth of the loop if necessary
      //////////////////////////////////////////////
      if(ampPref == t2first || ampPref == none){
	// ------------------------------
	// [1] In case of the T2-amp
	// ------------------------------
	if(AMPInd != NULL){
	  SQcont<SQindex*> inds(contra->get_LOuterIndices());
	  for(size_t n = 0;n < contra->get_Rtensors().size();++n) inds += SQcont<SQindex*>(contra->get_ROuterIndices(n));
	  bool hasInd(false);
	  for(auto i = inds.begin();i != inds.end();++i) if(**i == *AMPInd) hasInd = true;
	  if(hasInd && alreadyLoaded.count(*AMPInd) == 0){
	    CLoop(Indents[depthScope], *AMPInd, CPfile);
	    ++depthScope;
	    ReadAmp(Indents[depthScope], Amp, CPfile);
	    if(alreadyLoaded.count(*AMPInd) == 0) alreadyLoaded <= *AMPInd;
	  } // End if
	} // End if

	////////////////////////////////////////////////////////
	// Allocate the intermediate of the external type [1]
	////////////////////////////////////////////////////////
	for(auto contra2 = myBins.begin(); contra2 < myBins.end();++contra2){
	  bool allLoaded(true);
	  auto outers(contra2->get_LOuterIndices());
	  for(auto i = outers.begin();i != outers.end();++i) if(alreadyLoaded.count(**i) == 0) allLoaded = false;
	  if(!contra2->get_Lindices() && contra2->get_Ltensor().get_name() != LTensor_.get_name() && allLoaded && loadedInterms.count(contra2->get_Ltensor()) == 0){ 
	    declareInterm(depthScope, CPfile, contra2->get_Ltensor(), SQcont<SQindex*>(contra2->get_LInnerIndices()), Femto::Reaktor::External);
	    loadedInterms <= contra2->get_Ltensor();
	  } // End if
	} // End scope
	////////////////////////////////////////////////////////

	// ------------------------------
	// [2] In case of the S2-amp
	// ------------------------------
	if(SIGInd != NULL){
	  SQcont<SQindex*> inds(contra->get_LOuterIndices());
	  for(size_t n = 0;n < contra->get_Rtensors().size();++n) inds += SQcont<SQindex*>(contra->get_ROuterIndices(n));
	  bool hasInd(false);
	  for(auto i = inds.begin();i != inds.end();++i) if(**i == *SIGInd) hasInd = true;
	  if(hasInd && alreadyLoaded.count(*SIGInd) == 0){
	    CLoop(Indents[depthScope], *SIGInd, CPfile);
	    ++depthScope;
	    if(alreadyLoaded.count(*SIGInd) == 0) alreadyLoaded <= *SIGInd;
	  } // End if
	} // End if
      } // End if
      else if(ampPref == s2first){
	// ------------------------------
	// [1] In case of the S2-amp
	// ------------------------------
	if(SIGInd != NULL){
	  SQcont<SQindex*> inds(contra->get_LOuterIndices());
	  for(size_t n = 0;n < contra->get_Rtensors().size();++n) inds += SQcont<SQindex*>(contra->get_ROuterIndices(n));
	  bool hasInd(false);
	  for(auto i = inds.begin();i != inds.end();++i) if(**i == *SIGInd) hasInd = true;
	  if(hasInd && alreadyLoaded.count(*SIGInd) == 0){
	    CLoop(Indents[depthScope], *SIGInd, CPfile);
	    ++depthScope;
	    if(alreadyLoaded.count(*SIGInd) == 0) alreadyLoaded <= *SIGInd;
	  } // End if
	} // End if

	////////////////////////////////////////////////////////
	// Allocate the intermediate of the external type [1]
	////////////////////////////////////////////////////////
	for(auto contra2 = myBins.begin(); contra2 < myBins.end();++contra2){
	  bool allLoaded(true);
	  auto outers(contra2->get_LOuterIndices());
	  for(auto i = outers.begin();i != outers.end();++i) if(alreadyLoaded.count(**i) == 0) allLoaded = false;
	  if(!contra2->get_Lindices() && contra2->get_Ltensor().get_name() != LTensor_.get_name() && allLoaded && loadedInterms.count(contra2->get_Ltensor()) == 0){ 
	    declareInterm(depthScope, CPfile, contra2->get_Ltensor(), SQcont<SQindex*>(contra2->get_LInnerIndices()), Femto::Reaktor::External);
	    loadedInterms <= contra2->get_Ltensor();
	  } // End if
	} // End scope
	////////////////////////////////////////////////////////

	// ------------------------------
	// [2] In case of the T2-amp
	// ------------------------------
	if(AMPInd != NULL){
	  SQcont<SQindex*> inds(contra->get_LOuterIndices());
	  for(size_t n = 0;n < contra->get_Rtensors().size();++n) inds += SQcont<SQindex*>(contra->get_ROuterIndices(n));
	  bool hasInd(false);
	  for(auto i = inds.begin();i != inds.end();++i) if(**i == *AMPInd) hasInd = true;
	  if(hasInd && alreadyLoaded.count(*AMPInd) == 0){
	    CLoop(Indents[depthScope], *AMPInd, CPfile);
	    ++depthScope;
	    ReadAmp(Indents[depthScope], Amp, CPfile);
	    if(alreadyLoaded.count(*AMPInd) == 0) alreadyLoaded <= *AMPInd;
	  } // End if
	} // End if
      } // End if

      ///////////////////////////////////////////////////////
      // Allocation of the S2-tensor [2]
      ///////////////////////////////////////////////////////
      if( allocSigma.first == (size_t)(contra-myBins.begin()) && alreadyLoaded.count(*SIGInd) && !allocSigma.second) {
	ReadAmp(Indents[depthScope], Sig, CPfile);
	allocSigma.second = true;
      } // End if
      
      ////////////////////////////////////////////////////////
      // Allocate the intermediate of the external type [2]
      ////////////////////////////////////////////////////////
      for(auto contra2 = myBins.begin(); contra2 < myBins.end();++contra2){
	bool allLoaded(true);
	auto outers(contra2->get_LOuterIndices());
	for(auto i = outers.begin();i != outers.end();++i) if(alreadyLoaded.count(**i) == 0) allLoaded = false;
	if(!contra2->get_Lindices() && contra2->get_Ltensor().get_name() != LTensor_.get_name() && allLoaded && loadedInterms.count(contra2->get_Ltensor()) == 0){ 
	  declareInterm(depthScope, CPfile, contra2->get_Ltensor(), SQcont<SQindex*>(contra2->get_LInnerIndices()), Femto::Reaktor::External);
	  loadedInterms <= contra2->get_Ltensor();
	} // End if
      } // End scope
      ////////////////////////////////////////////////////////

      ////////////////////////////////////////////
      // Generate the binary contraction
      ////////////////////////////////////////////
      vector<string> ExtInd;
      vector<string> NameTen;
      vector<string> Consts;
      contDecl DecFirst;

      // ----------------------------
      // [1] Character constants
      // ----------------------------
      for(size_t numc = 0;numc < contra->get_Consts().size();++numc)
	if(contra->get_Consts()[numc] != "") Consts.push_back(contra->get_Consts()[numc]);

      // ----------------------------
      // [2] Tensors
      // ----------------------------
      for(size_t numt = 0;numt < contra->get_Rtensors().size();++numt){
	SQtensor tt(contra->get_Rtensors()[numt]);
        if(!is_RDM(tt.get_name())        && tt.get_name() != kDelta_name() && tt.get_name() != "Fc1"     &&  
	   !is_C4(tt.get_name()) && !is_C6(tt.get_name())          && tt.get_name() != C2_name() &&
	   find(NameTen.begin(), NameTen.end(), tt.get_name()) == NameTen.end()){
          NameTen.push_back(tt.get_name());
	} // End if
      } // End numt
      sort(NameTen.begin(), NameTen.end());
      NameTen.push_back(contra->get_Ltensor().get_name());

      // -----------------------------
      // [3] External indices
      // -----------------------------
      for(auto i = alreadyLoaded.begin();i != alreadyLoaded.end();++i) ExtInd.push_back(i->get_index());
      sort(ExtInd.begin(), ExtInd.end());

      DecFirst.insert(contDecl::value_type("DecInd",    ExtInd));  // Names of external indices
      DecFirst.insert(contDecl::value_type("DecConst",  Consts));  // All the constants
      DecFirst.insert(contDecl::value_type("DecTensor", NameTen)); // Names of the tensors 

      //////////////////////////////////////////
      // Generate each part ////////////////////
      //////////////////////////////////////////
      ostringstream stm1;
      stm1 << (int)(contra-myBins.begin());
      string myTitle_if("g_if_" + title_ + "_no" + myLabel + "_x" + stm1.str() + theLabel);
      string myTitle   ("g_"    + title_ + "_no" + myLabel + "_x" + stm1.str() + theLabel);
      transform(myTitle_if.begin(), myTitle_if.end(), myTitle_if.begin(), (int(*)(int))tolower);
      transform(myTitle   .begin(), myTitle   .end(), myTitle   .begin(), (int(*)(int))tolower);

      bool myCondition(contra->get_Ltensor().get_name() == LTensor_.get_name() && isBareLHS_);
      makeCPP_header2(*contra, myTitle_if, CHfile, Indents[depthScope], DecFirst, myCondition);
      makeF90_interface2_new(*contra, myTitle_if, F90file, Indents[depthScope], DecFirst, myCondition);
      makeCPP_bodyType2_new(*contra, myTitle_if, CPfile, Indents[depthScope], DecFirst, myCondition);
      if(use_gemm_ && contra->get_Rtensors().size() == 2)
	binary_contract3_new(*contra, myTitle, F90file, Indents[depthScope], DecFirst, myCondition);
      else
	makeF90_contract2(*contra, myTitle, F90file, Indents[depthScope], DecFirst, myCondition);
      //{ cout << "Not implemented yet <119087>" << endl; abort(); }

      ///////////////////////////////////////////
      // Close the loop structure
      ///////////////////////////////////////////
      SQcont<SQindex> stillNeeded; 
      if(ERIInd != NULL) if(alreadyLoaded.count(*ERIInd) != 0 && stillNeeded.count(*ERIInd) == 0) stillNeeded <= *ERIInd;     
      size_t currentContra(contra-myBins.begin());
      if(currentContra != myBins.size()-1){
	size_t next(currentContra+1);
	SQcont<SQindex> temp(returnsExtIndices(myBins[next]));
	for(auto i = alreadyLoaded.begin();i != alreadyLoaded.end();++i)
	  if(temp.count(*i) && stillNeeded.count(*i) == 0) stillNeeded <= *i;

      } // End if

      SQcont<SQindex> notNeeded;
      for(auto i = alreadyLoaded.rbegin();i != alreadyLoaded.rend();++i)
	if(stillNeeded.count(*i) == 0) {
	  notNeeded <= *i;
	}

      if(alreadyLoaded.size() != stillNeeded.size() + notNeeded.size()){ 
	cout << "makeContractions: Algorithmic error ocurred >> " << *contra <<  " << " << endl;
	cout << " << Loaded indices >> " << endl;
	for(auto i = alreadyLoaded.begin();i != alreadyLoaded.end();++i) cout << " >> " << *i << endl;
	cout << " << Needed indices >> " << endl;
	for(auto i = stillNeeded.begin();i != stillNeeded.end();++i) cout << " >> " << *i << endl;
	cout << " << NotNeeded indices >> " << endl;
	for(auto i = notNeeded.begin();i != notNeeded.end();++i) cout << " >> " << *i << endl;
	CPfile.flush();
	abort(); 
      } // End if
      else alreadyLoaded = stillNeeded;

    CPfile.flush();
    // Close the loop
    for(auto i = notNeeded.begin();i != notNeeded.end();++i){
      CPfile << "  // --> " << *i << " [notNeeded]" << endl;

    cout << "depth  >> " << depthScope << endl;
    cout << "length >> " << notNeeded.size() << endl; 
    cout << "       >> " << *i << endl;
        if(SIGInd != NULL)
	  if(*SIGInd == *i && allocSigma.second) {
	    int myScope(depthScope == 0 ? depthScope+1 : depthScope);
	    AccAmp(Indents[myScope], Sig, CPfile);
	  } // End i
	LoopEnd_i(Indents[depthScope-1], *i, CPfile);
	--depthScope;

      } // End i

    } // End contra
    
    
  }

}} // Femto::Reaktor
