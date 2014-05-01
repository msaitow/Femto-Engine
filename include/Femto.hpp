//
//  Femto.hpp
//  
//  FEMTO :: An Integrated Toolset for Automated Tensor Generation 
//
//  Created by Masaaki Saitow on 12/06/11.
//  Copyright (c) 2012 Masaaki Saitow. All rights reserved.
//

//#ifndef FEMTO_HPP
//#define FEMTO_HPP
#pragma once

#include <vector>
#include <utility>
#include <string>
#include <stdlib.h>
#include <SQcont.hpp>
//using namespace std;

namespace Femto{

  // :: typedef for simple int vectors
  typedef std::vector<int> Ivector;
  typedef std::vector<std::vector<int> > IIvector;

  // :: Symmetry < Pattern, Factors > 
  typedef std::pair<IIvector, Ivector> Symmetry;

  // ::::::: Charactar of one-electronic state :::::::
  // :: Orbital groups      :: spin        :: Index
  // ---------------------------------------------------
  //  0 <- core             :: spin-free   :: SQindex
  //  1 <- active           :: spin-free   :: SQindex
  //  2 <- virtual          :: spin-free   :: SQindex
  //  3 <- auxiliary basis  :: spin-free   :: SQdf_aux
  //       (used in DF/RI) 
  //  4 <- core             :: alpha       :: SQindex
  //  5 <- active           :: alpha       :: SQindex
  //  6 <- virtual          :: alpha       :: SQindex
  //  7 <- core             :: beta        :: SQindex
  //  8 <- active           :: beta        :: SQindex
  //  9 <- virtual          :: beta        :: SQindex
  // 10 <- AO basis         :: --          :: SQindex 
  enum char_state {core  =0, act  =1, virt  =2, aux=3, 
		   core_a=4, act_a=5, virt_a=6,
		   core_b=7, act_b=8, virt_b=9, ao=10};  

  // :: Types of tensor notations
  // 0 <- Dirac
  // 1 <- Mulliken
  // 2 <- None for second quantized operator
  enum notation{Dirac=0, Mulliken=1, None=2}; 

  // Utilities for handling permutations and n-tuples
  IIvector makePermutations(int n);
  IIvector makeTuples1  (int n,            Ivector &inList);
  IIvector makeTuples2  (int n,            Ivector &inList);
  IIvector makeTuples3  (int n, int order, Ivector &inList);

  int fact(const int n);
  Symmetry h1_symm();
  Symmetry h2_symm();
  Symmetry h2_symmM(); // In Mulliken notation
  Symmetry u2_symm(); 
  Symmetry t2_symm(); 
  Symmetry u4_symm();
  Symmetry uni_symm(int length);

  // Auxiliary structs for sorting object
  struct SecGreat { bool operator()(const std::pair<int,int> &i, const std::pair<int,int> &j) const 
    { return (i.second < j.second ? true : false); } };

  // Get number of permutations inbetween two Ivectors. But such the object must be 
  // composed of comtiguous integers. If not, it may behaves crazy (taken from Eric's SQA).
  int get_num_perms(std::vector<int> &ti, std::vector<int> &bi);

  // Print date and time
  std::string Femto_date();
  // Print logo of Femto
  std::string Femto_logo(const std::string s);

  // *********************************************************
  // Maximum nNumber of terms appear in the normal ordering 
  // *********************************************************
  inline const int Nterms()
  { return 2000; }

  // *********************************************************
  // Definition of name of kDelta
  // *********************************************************
  inline std::string kDelta_name()
  { return "kDelta"; }

  // *********************************************************
  // Definition of name of sfGen
  // *********************************************************
  inline std::string sfGen_name()
  { return "E"; }

  // *********************************************************
  // Definition of name of RDM
  // *********************************************************
  inline std::string RDM_name()
  { return "D"; }

  // *********************************************************
  // Definition of name of creation operator
  // *********************************************************
  inline std::string aCre_name()
  { return "aC"; }

  // *********************************************************
  // Definition of name of destruction operator
  // *********************************************************
  inline std::string aDes_name()
  { return "aD"; }

  // *********************************************************
  // Whether name is equal to that of sfGen
  // *********************************************************
  inline bool is_sfGen(const std::string name)
  { 
    std::string num(name.substr(1));
    if(name.at(0) == sfGen_name()[0] && atoi(num.c_str())) return true;
    else                                                   return false;
  }

  // *********************************************************
  // Whether name is equal to that of RDM
  // *********************************************************
  inline bool is_RDM(const std::string name)
  { 
    std::string num(name.substr(1));
    if(name.at(0) == RDM_name()[0] && atoi(num.c_str())) return true;
    else                                                 return false;
  }

} //Femto::

//#endif
