#include"Eclipse.hpp"

// MM""""""""`M          dP oo                            
// MM  mmmmmmmM          88                               
// M`      MMMM .d8888b. 88 dP 88d888b. .d8888b. .d8888b. 
// MM  MMMMMMMM 88'  `"" 88 88 88'  `88 Y8ooooo. 88ooood8 
// MM  MMMMMMMM 88.  ... 88 88 88.  .88       88 88.  ... 
// MM        .M `88888P' dP dP 88Y888P' `88888P' `88888P' 
// MMMMMMMMMMMM                88                         

using namespace std;

Eclipse::Eclipse()
{
  
  vector<Femto::Core::SQindex*> EL;
  vector<Femto::Core::SQindex*> ER;
  vector<Femto::Core::SQindex*> ET;

  A = Femto::Core::SQindex("a", Femto::virt);
  B = Femto::Core::SQindex("b", Femto::virt, true);
  C = Femto::Core::SQindex("c", Femto::virt);
  D = Femto::Core::SQindex("d", Femto::virt, true);

  I = Femto::Core::SQindex("i", Femto::act);
  J = Femto::Core::SQindex("j", Femto::act, true);
  K = Femto::Core::SQindex("k", Femto::act);
  L = Femto::Core::SQindex("l", Femto::act, true);
  M = Femto::Core::SQindex("m", Femto::act);
  N = Femto::Core::SQindex("n", Femto::act, true);

  W = Femto::Core::SQindex("w", Femto::core);
  X = Femto::Core::SQindex("x", Femto::core, true);
  Y = Femto::Core::SQindex("y", Femto::core);
  Z = Femto::Core::SQindex("z", Femto::core, true);

  for(size_t space = 0;space < 3;++space){
    Ps.push_back(Femto::Core::SQindex("p", (Femto::char_state)space, true));
    Qs.push_back(Femto::Core::SQindex("q", (Femto::char_state)space, true));
    Rs.push_back(Femto::Core::SQindex("r", (Femto::char_state)space, true));
    Ss.push_back(Femto::Core::SQindex("s", (Femto::char_state)space, true));
  }

  //////////////////////////////////////////////////////////////////
  // ***** Definition of the left-side excitation operators ***** 
  //////////////////////////////////////////////////////////////////  
  // 1) (o,o)->(v,v) operator (left-side)
  EL.push_back(&I);
  EL.push_back(&K);  
  EL.push_back(&A);
  EL.push_back(&C);  
  EL_indices.push_back(EL);
  EL.clear();
  
  // 2) (o,o)->(o,v) operator (left-side)
  EL.push_back(&I);
  EL.push_back(&K);  
  EL.push_back(&M);
  EL.push_back(&A);  
  EL_indices.push_back(EL);
  EL.clear();

  // 3) (c,c)->(v,v) operator (left-side)
  EL.push_back(&W);
  EL.push_back(&Y);  
  EL.push_back(&A);
  EL.push_back(&C);  
  EL_indices.push_back(EL);
  EL.clear();

  // 4) (c,a)->(v,v) operator (left-side)
  EL.push_back(&W);
  EL.push_back(&I);  
  EL.push_back(&A);
  EL.push_back(&C);  
  EL_indices.push_back(EL);
  EL.clear();

  // 5) (c,c)->(a,a) operator (left-side)
  EL.push_back(&W);
  EL.push_back(&Y);  
  EL.push_back(&I);
  EL.push_back(&K);  
  EL_indices.push_back(EL);
  EL.clear();

  // 6) (c,c)->(a,v) operator (left-side)
  EL.push_back(&W);
  EL.push_back(&Y);  
  EL.push_back(&I);
  EL.push_back(&A);  
  EL_indices.push_back(EL);
  EL.clear();

  // 7) (c,a)->(a,a) operator (left-side)
  EL.push_back(&W);
  EL.push_back(&I);  
  EL.push_back(&K);
  EL.push_back(&M);  
  EL_indices.push_back(EL);
  EL.clear();

  // 8) (c,a)->(a,v) operator (left-side)
  EL.push_back(&W);
  EL.push_back(&I);  
  EL.push_back(&K);
  EL.push_back(&A);  
  EL_indices.push_back(EL);
  EL.clear();

  // 9) (a,c)->(a,v) operator (left-side)
  EL.push_back(&I);
  EL.push_back(&W);  
  EL.push_back(&K);
  EL.push_back(&A);  
  EL_indices.push_back(EL);
  EL.clear();


  //////////////////////////////////////////////////////////////////
  // ***** Definition of the right-side excitation operators ***** 
  //////////////////////////////////////////////////////////////////  
  // 1) (o,o)->(v,v) operator (right-side)
  ER.push_back(&B);
  ER.push_back(&D);
  ER.push_back(&J);
  ER.push_back(&L);
  ER_indices.push_back(ER);
  ER.clear();

  ET.push_back(&J);
  ET.push_back(&L);
  ET.push_back(&B);
  ET.push_back(&D);
  ET_indices.push_back(ET);
  ET.clear();

  // 2) (o,o)->(o,v) operator (right-side)
  ER.push_back(&N);
  ER.push_back(&B);
  ER.push_back(&L);
  ER.push_back(&J);
  ER_indices.push_back(ER);
  ER.clear();

  ET.push_back(&L);
  ET.push_back(&J);
  ET.push_back(&N);
  ET.push_back(&B);
  ET_indices.push_back(ET);
  ET.clear();

  // 3) (c,c)->(v,v) operator (right-side)
  ER.push_back(&B);
  ER.push_back(&D);
  ER.push_back(&X);
  ER.push_back(&Z);
  ER_indices.push_back(ER);
  ER.clear();

  ET.push_back(&X);
  ET.push_back(&Z);
  ET.push_back(&B);
  ET.push_back(&D);
  ET_indices.push_back(ET);
  ET.clear();

  // 4) (c,a)->(v,v) operator (right-side)
  ER.push_back(&B);
  ER.push_back(&D);
  ER.push_back(&X);
  ER.push_back(&J);
  ER_indices.push_back(ER);
  ER.clear();

  ET.push_back(&X);
  ET.push_back(&J);
  ET.push_back(&B);
  ET.push_back(&D);
  ET_indices.push_back(ET);
  ET.clear();

  // 5) (c,c)->(a,a) operator (right-side)
  ER.push_back(&J);
  ER.push_back(&L);
  ER.push_back(&X);
  ER.push_back(&Z);
  ER_indices.push_back(ER);
  ER.clear();

  ET.push_back(&X);
  ET.push_back(&Z);
  ET.push_back(&J);
  ET.push_back(&L);
  ET_indices.push_back(ET);
  ET.clear();

  // 6) (c,c)->(a,v) operator (right-side)
  ER.push_back(&J);
  ER.push_back(&B);
  ER.push_back(&X);
  ER.push_back(&Z);
  ER_indices.push_back(ER);
  ER.clear();

  ET.push_back(&X);
  ET.push_back(&Z);
  ET.push_back(&J);
  ET.push_back(&B);
  ET_indices.push_back(ET);
  ET.clear();

  // 7) (c,a)->(a,a) operator (right-side)
  ER.push_back(&J);
  ER.push_back(&L);
  ER.push_back(&X);
  ER.push_back(&N);
  ER_indices.push_back(ER);
  ER.clear();

  ET.push_back(&X);
  ET.push_back(&N);
  ET.push_back(&J);
  ET.push_back(&L);
  ET_indices.push_back(ET);
  ET.clear();

  // 8) (c,a)->(a,v) operator (right-side)
  ER.push_back(&J);
  ER.push_back(&B);
  ER.push_back(&X);
  ER.push_back(&N);
  ER_indices.push_back(ER);
  ER.clear();

  ET.push_back(&X);
  ET.push_back(&N);
  ET.push_back(&J);
  ET.push_back(&B);
  ET_indices.push_back(ET);
  ET.clear();

  // 9) (a,c)->(a,v) operator (right-side)
  ER.push_back(&J);
  ER.push_back(&B);
  ER.push_back(&N);
  ER.push_back(&X);
  ER_indices.push_back(ER);
  ER.clear();

  ET.push_back(&N);
  ET.push_back(&X);
  ET.push_back(&J);
  ET.push_back(&B);
  ET_indices.push_back(ET);
  ET.clear();


  //////////////////////////////////////////////////////////////////
  // ***** Definition of the name of the excitation operators ***** 
  //////////////////////////////////////////////////////////////////  
  names.push_back("oovv"); // 1) oovv
  names.push_back("ooov"); // 2) ooov
  names.push_back("ccvv"); // 3) ccvv
  names.push_back("covv"); // 4) covv
  names.push_back("ccoo"); // 5) ccoo
  names.push_back("ccov"); // 6) ccov
  names.push_back("cooo"); // 7) cooo
  names.push_back("coov"); // 8) coov
  names.push_back("ocov"); // 9) ocov

  // File to write the sigma equations (LaTeX)
  formulas.open("formulas.tex");

  // File to write the names of source file for non-zero block of Fock hamiltonian matrix
  make_fock_vector.open("make_fock_v"); 

  // File to write the declaration added to the header
  header_fock_vector.open("header_fock_v");

  // File to write the calling section in C++ source
  gaxpy_fock_vector.open("gaxpy_fock_v");

  // File to flush all the terms
  print_eqns.open("fic-mrci_terms.out");
  print_eqns << Femto::Femto_logo("");
  print_eqns << endl;
  print_eqns << "           " << Femto::Femto_date();
  print_eqns << endl;

  status = "\n### SUMMARY ###\n";
  total_num_terms = 0;
  nums_rdm4.first  = 0;
  nums_rdm4.second = 0;
}

Eclipse::~Eclipse() { 
  { // For standard output 
  cout << status << endl; 
  cout << " ++ Number of 4-RDMs in Hdiag : " << nums_rdm4.first  << endl; 
  cout << " ++ Number of 4-RDMs in sigma : " << nums_rdm4.second << endl; 
  cout << " ++ Total number of terms     : " << total_num_terms  << endl; 
  cout << " ++ Terms with 4-RDMs :: " << endl;
  for(std::map<std::string, std::vector<Femto::Core::SQterm> >::iterator symbol = d4_terms.begin();symbol != d4_terms.end();++symbol){
    cout << "<" << symbol->first << ">" << endl;
    int count(0);
    for(std::vector<Femto::Core::SQterm>::iterator t = symbol->second.begin();t != symbol->second.end();++t)
      cout << "    -- "  << count++ << " : " << *t << endl;
  } // End symbol
  cout << endl;
  cout << " ++ Terms with 4-RDMs (actually survived) :: " << endl;
  for(std::map<std::string, std::vector<Femto::Core::SQterm> >::iterator symbol = d4_terms_truely.begin();symbol != d4_terms_truely.end();++symbol){
    cout << "<" << symbol->first << ">" << endl;
    int count(0);
    for(std::vector<Femto::Core::SQterm>::iterator t = symbol->second.begin();t != symbol->second.end();++t)
      cout << "    -- "  << count++ << " : " << *t << endl;
  } // End symbol
  cout << endl;
  }
  { // For the file
  print_eqns << status << endl; 
  print_eqns << " ++ Number of 4-RDMs in Hdiag : " << nums_rdm4.first  << endl; 
  print_eqns << " ++ Number of 4-RDMs in sigma : " << nums_rdm4.second << endl; 
  print_eqns << " ++ Total number of terms     : " << total_num_terms  << endl; 
  print_eqns << " ++ Terms with 4-RDMs :: " << endl;
  for(std::map<std::string, std::vector<Femto::Core::SQterm> >::iterator symbol = d4_terms.begin();symbol != d4_terms.end();++symbol){
    print_eqns << "<" << symbol->first << ">" << endl;
    int count(0);
    for(std::vector<Femto::Core::SQterm>::iterator t = symbol->second.begin();t != symbol->second.end();++t)
      print_eqns << "    -- "  << count++ << " : " << *t << endl;
  } // End symbol
  print_eqns << endl;
  print_eqns << " ++ Terms with 4-RDMs (actually survived) :: " << endl;
  for(std::map<std::string, std::vector<Femto::Core::SQterm> >::iterator symbol = d4_terms_truely.begin();symbol != d4_terms_truely.end();++symbol){
    print_eqns << "<" << symbol->first << ">" << endl;
    int count(0);
    for(std::vector<Femto::Core::SQterm>::iterator t = symbol->second.begin();t != symbol->second.end();++t)
      print_eqns << "    -- "  << count++ << " : " << *t << endl;
  } // End symbol
  print_eqns << endl;
  }
}
