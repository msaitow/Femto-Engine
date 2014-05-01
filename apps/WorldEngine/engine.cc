
#include <map>
#include <iostream>
#include <vector>
#include <fstream>
#include <boost/format.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/tokenizer.hpp>
#include <Femto.hpp>
#include <SQterm.hpp>
#include <SQindex.hpp>
#include <SQtensor.hpp>
#include <SQreaktor.hpp>
#include <SQportal.hpp>

//using namespace std;
using std::cout;
using std::endl;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
class contractionReader{

public:

  contractionReader(const Femto::SQcont<std::pair<std::string,std::string> > tensor2elem, const std::string mySynonym, const std::string myFormat);
  std::pair<std::string,std::string> readElements(const std::string fileName);
  bool isTensorFile(const std::string fileName) const;
  
private:
  
  std::string theFormat_;  // Format for the Fafnir tensor file
  std::string theSynonym_; // Synonym for the sigma contarctions
  std::map<std::string, std::string> tensor2elem_;
}; // End class

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
contractionReader::contractionReader(const Femto::SQcont<std::pair<std::string,std::string> > tensorLabels, const std::string theSynonym, const std::string theFormat)
  : theFormat_ (theFormat),
    theSynonym_(theSynonym)
{ for(auto t = tensorLabels.cbegin(); t != tensorLabels.cend();++t) tensor2elem_[t->first] = t->second; }

//////////////////////////////////////////////////////////////////////////////////////
std::pair<std::string,std::string> 
contractionReader::readElements(const std::string fileName)
{
  
  if(theSynonym_.size() + fileName.size() < theFormat_.size()) { cout << "readElements: Algorithmic error <r203r92>" << endl; abort(); }
  else{
    std::string mySynonym("");
    for(int i = 0;i < theSynonym_.size();++i) mySynonym += fileName[i];
    if(mySynonym != theSynonym_) { cout << "readElements: Algorithmic error <ry23y30reu> >> " << mySynonym << " << " << endl; abort(); }
    else {
      size_t formatStart(fileName.size()-theFormat_.size());
      std::string myFormat(fileName.substr(formatStart));
      if(myFormat != theFormat_) { cout << "readElements: Algorithmic error <3hr92h4r9uwe>" << endl; abort(); }
      else {
	
	typedef boost::char_separator<char> char_separator;
	typedef boost::tokenizer<char_separator> tokenizer;
	char_separator sep("_", "", boost::keep_empty_tokens);
	tokenizer tokens(fileName, sep);
	int pos(0);
	std::string lBasis(""); // Left space
	std::string rBasis(""); // Right space
	for(auto tok_iter = tokens.begin();tok_iter != tokens.end();++tok_iter,++pos){
	  if(pos == 1) lBasis = *tok_iter;
	  if(pos == 2) rBasis = *tok_iter;
	} // End tok_iter
	if(lBasis == "" || rBasis == "") { cout << "readElements: Algorithmic error <rh283xgr28yr>" << endl; abort(); }
	std::string rBasisTrimmed("");
	size_t myPos(rBasis.size()-theFormat_.size());
	for(int i = 0;i < myPos;++i) rBasisTrimmed += rBasis[i];
	
	return std::pair<std::string, std::string>(tensor2elem_[lBasis], tensor2elem_[rBasisTrimmed]);
	
      } // End else
    } // End else
  } // End else
  
} // End func

//////////////////////////////////////////////////////////////////////////////////////
bool contractionReader::isTensorFile(const std::string fileName) const
{
  std::string mySynonym("");
  for(int i = 0;i < theSynonym_.size();++i) mySynonym += fileName[i];
  if(theSynonym_ != mySynonym) return false;
  
  if(fileName.size() < theFormat_.size()) return false;
  else{
    size_t formatStart(fileName.size()-theFormat_.size());
    std::string myFormat(fileName.substr(formatStart));
    if(myFormat == theFormat_) return true;
    else                       return false;
  } // End else
  
} // End func

///////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[]){
  
  if(argc == 1){
    cout << "worldEngine: No target specified" << endl;
    abort(); 
  } // End if
  else if(argc > 2){
    cout << "worldEngine: Only one target can be specified at once" << endl;
    abort();
  } // End if

  const std::string myTag(argv[1]);

  ///// Directory for the tensor files
  boost::filesystem::path tensorDir("./TensorFiles");
  Femto::SQcont<std::pair<std::string, std::string> > LabelAndTensor; // <Label, Tensor file name>
  
  LabelAndTensor <= std::make_pair("ref",  "g"   );
  LabelAndTensor <= std::make_pair("aavv", "oovv");
  LabelAndTensor <= std::make_pair("aaav", "ooov");
  LabelAndTensor <= std::make_pair("ccvv", "ccvv");
  LabelAndTensor <= std::make_pair("cavv", "covv");
  LabelAndTensor <= std::make_pair("ccaa", "ccoo");
  LabelAndTensor <= std::make_pair("ccav", "ccov");
  LabelAndTensor <= std::make_pair("caaa", "cooo");
  LabelAndTensor <= std::make_pair("caav", "coov");
  LabelAndTensor <= std::make_pair("acav", "ocov");

  contractionReader myReader(LabelAndTensor, myTag, ".tensor");

  //////////// Tensor files //////////////////////////////////////////////////////
  Femto::SQcont<std::pair<std::string,std::pair<std::string,std::string> > > sigmaAndElems; // <Name of left-hand side tensor, <Name of tensor file, Name of the element> >
  { // Extract the tensor file name and the corresponding matrix block name
    boost::filesystem::directory_iterator end;
    for(boost::filesystem::directory_iterator f(tensorDir); f != end;++f){
      auto myPath(f->path().string());
      auto myName(myPath.substr(tensorDir.string().size()+1));
      if(myReader.isTensorFile(myName) && !boost::filesystem::is_directory(*f)){
        std::pair<std::string, std::string> myLabel(myReader.readElements(myName));
	std::ostringstream myBlock;
	myBlock << myTag << "_" << myLabel.first << "_" << myLabel.second;
	std::string sigmaName("");
	if(myLabel.first == "g") sigmaName = "@S0";
	else                     sigmaName = "@S2";
	sigmaAndElems <= std::make_pair(sigmaName, std::pair<std::string, std::string>(myName,myBlock.str()));
	cout << "Tensor file detected => " << myName << " ..." << endl;
      } // End if
    } // End f
  } // End scope
  ////////////////////////////////////////////////////////////////////////////////

  std::string myHeaders("");
  for(auto n = sigmaAndElems.begin();n != sigmaAndElems.end();++n){
    std::string sigName(n->first);
    std::string tenName(n->second.first);
    std::string elemName(n->second.second);
    bool   isBareLHS(sigName != "@S0");
    cout << elemName << " ---> " << sigName << endl;
    std::vector<std::vector<Femto::Core::SQbinary> > myBins;

    std::ostringstream tensorFile;
    tensorFile << tensorDir.string() << "/" << tenName;
 
    Femto::Portal::SQportal myPortal(tensorFile.str(), sigName, Femto::Mulliken);
    
    myPortal.transportContractions(myBins, Femto::Portal::NEW2);

    cout << endl << endl;
    cout << "----------------------------------------------------------" << endl;
    cout << "----------------------------------------------------------" << endl;
    for(auto ts = myBins.begin();ts != myBins.end();++ts){
      cout << boost::format(">> No. %d\n") % ((int)(ts-myBins.begin())+1);
      Femto::SQcont<Femto::Core::SQbinary> temp(*ts);
      cout << temp << endl << endl; 
    } // End ts
    
    Femto::Reaktor::SQreaktor gen(myBins, elemName, isBareLHS, "V2", "T2"); // for aavv/aavv
    gen.init_header("mrci_header");
    myHeaders += gen.generate(Femto::Reaktor::Orz, Femto::Reaktor::Fafnir);
  } // End n
  
  std::ofstream headerFractions(myTag + ".headerFractions.hpp");
  headerFractions << myHeaders;

}
