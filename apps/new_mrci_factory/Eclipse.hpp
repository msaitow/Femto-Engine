
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <boost/format.hpp>
#include <Femto.hpp>
#include <SQterm.hpp>
#include <SQindex.hpp>
#include <SQtensor.hpp>
#include <SQreaktor.hpp>

// MM""""""""`M          dP oo                            
// MM  mmmmmmmM          88                               
// M`      MMMM .d8888b. 88 dP 88d888b. .d8888b. .d8888b. 
// MM  MMMMMMMM 88'  `"" 88 88 88'  `88 Y8ooooo. 88ooood8 
// MM  MMMMMMMM 88.  ... 88 88 88.  .88       88 88.  ... 
// MM        .M `88888P' dP dP 88Y888P' `88888P' `88888P' 
// MMMMMMMMMMMM                88                         

// Code to generate all the tensorial contraction codes for fic-MRCI

#define _NEW_CON
#define _CODE_GENE
#define _MAKE_GENE

// This enables the quasi-commutator approach, <E_{L} F E_{R}> ~= <[E_{L}, [H, E_{R}]] > + Eres <E_{L} E_{R}>
//#define _NEW_ENO

class Eclipse{

public:
  Eclipse();
  void generate_overlap();
  void generate_diag();
  void generate_g_e();
  void generate_e_g();
  void generate_e_e();
  void generate_fock_vector();
  void generate_fock_diag();
  ~Eclipse();

private:
  std::pair<int, int> nums_rdm4;
  std::map<std::string, std::vector<Femto::Core::SQterm> > d4_terms;
  std::map<std::string, std::vector<Femto::Core::SQterm> > d4_terms_truely;
  int total_num_terms;
  std::string status;
  std::ofstream formulas;
  std::ofstream print_eqns;
  std::ofstream make_fock_vector;
  std::ofstream header_fock_vector;
  std::ofstream gaxpy_fock_vector;
  std::vector<std::string> names;
  std::vector<std::vector<Femto::Core::SQindex*> > EL_indices;
  std::vector<std::vector<Femto::Core::SQindex*> > ER_indices;
  std::vector<std::vector<Femto::Core::SQindex*> > ET_indices;
  std::vector<Femto::Core::SQindex> Ps;   
  std::vector<Femto::Core::SQindex> Qs;   
  std::vector<Femto::Core::SQindex> Rs;   
  std::vector<Femto::Core::SQindex> Ss;   

//*//   ///////////////////////////////////////////////////////
//*//   // * Utilities that stands for where orz is
//*//   ///////////////////////////////////////////////////////
//*//   inline std::string f_source()
//*//   { return "tensors/"; }

  ///////////////////////////////////////////////////////
  // * Basic indices used in the code generation process
  ///////////////////////////////////////////////////////
  // Virtual
  Femto::Core::SQindex A;
  Femto::Core::SQindex B; // dummy
  Femto::Core::SQindex C;
  Femto::Core::SQindex D; // dummy
		
  // Active  
  Femto::Core::SQindex I;
  Femto::Core::SQindex J; // dummy
  Femto::Core::SQindex K;
  Femto::Core::SQindex L; // dummy
  Femto::Core::SQindex M;
  Femto::Core::SQindex N; // dummy
		
  // Core  
  Femto::Core::SQindex W;
  Femto::Core::SQindex X; // dummy
  Femto::Core::SQindex Y;
  Femto::Core::SQindex Z; // dummy

};

