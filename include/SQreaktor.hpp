//
//  SQreaktor.hpp
//  
//  Class for the automated code generator (the reaktor)
//
//  Created by Masaaki Saitow on 12/09/04.
//  Copyright (c) 2012 Masaaki Saitow. All rights reserved.
//

#pragma once

#include <cstring>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <fstream>
#include <SQterm.hpp>
#include <SQtensor.hpp>
#include <SQindex.hpp>
#include <SQcontract.hpp>
#include <SQbinary.hpp>
#include <map>

//#define _FOR_PUBLICATION

//#define _VERBOSE_MODE
#define _CHRONO // Turn on coost::chrono

namespace Femto { namespace Reaktor {
  
  // Definition of the intermediate tensor
  inline std::string Interm_name()
  { return "W"; }

  // Whether name given is of intermediate tensor
  inline bool is_Interm(const std::string name)
  {
    if(name.at(0) == Interm_name()[0]) return true;
    else                               return false;
  }

  // Definition of the T1 amplitude
  inline std::string T1_name()
  { return "T1"; }

  // Whether name given is of T1 amplitude
  inline bool is_T1(const std::string name)
  {
    if(name == T1_name()) return true;
    else                  return false;
  }

  // Definition of the intermediate tensor constructed from D4 and V2 of six indices
  // D4C(*,o,o,o,o,o) <-- D4(o,o,o,o,o,o1,o2,o3) V(*,o1,o2,o3)
  inline std::string D4C_name()
  { return "C5"; } 

  // Definition of the intermediate tensor constructed from D4 and V2 of six indices (on LHS)
  // D4C(*,o,o,o,o,o) <-- D4(o,o,o,o,o,o1,o2,o3) V(*,o1,o2,o3)
  inline std::string D4C_nameL()
  { return "C4D_5"; } 

  // Definition of the intermediate tensor constructed from D4 and V2 of six indices (on LHS)
  // C4(o,o,o,o) <-- D4(o,o,o,o,o1,o2,o3,o4) V(o1,o2,o3,04)
  inline std::string C4_nameL()
  { return "C4D_4"; } 

  // Definition of the intermediate tensor constructed from D4 and V2 of six indices (on LHS)
  // C6(o,o,o,o,o,o) <-- D4(o,o,o,o,o,o,o1,o2) Fock(o1,o2)
  inline std::string C6_nameL()
  { return "C4D_6"; } 

  // Definition of the 2-body cumulant
  inline std::string C2_name()
  { return "C2"; } 

  // Whether name given is same to that of D4C
  inline bool is_D4C(std::string name)
  { 
    if(name == D4C_name()) return true;
    else                   return false; 
  }

  // Definition of the inetermediate tenor constructed from D4 and V2 of 4 indices
  inline std::string C4_name()
  { return "C4"; }

  // Whether name given is same to that of C4
  inline bool is_C4(std::string name)
  { 
    if(name == C4_name()) return true;
    else                  return false; 
  }

  // Definition of the inetermediate tenor constructed from D4 and Fock matrix of 6 indices
  inline std::string C6_name()
  { return "C6"; }

  // Whether name given is same to that of C6
  inline bool is_C6(std::string name)
  { 
    if(name == C6_name()) return true;
    else                  return false; 
  }

  // Definition of the CAS-Fock tensor
  inline std::string Fock_name()
  { return "P1"; }

  // Whether the internediate is the external type
  enum decType {Internal, External};

  // Priority flag
  enum priority {polynomial=0, memory=1};
  // Code generation paradigms of generator
  enum Modes {SimpleLoops=0, Factorize=1, Cas=2, Fafnir=3};
  // Target program
  enum Codes {Orz=0, Psi=1};
  // Whether the generateContract cares about the MPI parallelisms
  enum ParaFlag {NonPara=0, Para=1};

  // Container of information of declaration
  typedef std::map<std::string, std::vector<std::string> > contDecl; 

  ////////////////////////////////////////////////////////////////////////////////
  // SQreaktor: A tensor reactor class that produces a stream of optimized code //
  //            (ver.0.2, Mark2)                                                //
  ////////////////////////////////////////////////////////////////////////////////
  
  ///////////////////////////////////////////////////////////////////////////////////////
  // NOTE: In the factorize method, it is assumed that all the indices on the LHS      //
  //       appear also on the RHS and such indices are of non-dummy type. In case that // 
  //       these conditions are not valid (such as preconditionor), use of simpleloops //
  //       is strongly recommended.                                                    //
  ///////////////////////////////////////////////////////////////////////////////////////
  class SQreaktor{
      
  public:
    SQreaktor(const Femto::Core::SQtensor &LTensor, const std::vector<Femto::Core::SQterm> &inTerms, const std::string title="default_scheme", 
              const bool isBareLHS=false, const std::string name_h2 = "V2", const std::string name_amp="T2");

    SQreaktor(const std::vector<std::vector<Femto::Core::SQbinary> > &inContras, const std::string title="default_scheme", 
	        const bool isBareLHS=false, const std::string name_h2 = "V2", const std::string name_amp="T2");
    
    //~SQreaktor();

    // Process kDeltas to produce kDelta-free stream of Femto::Core::SQterms, and the information of kDeltas
    // are stored in kDeltas_
    void process_kDeltas();

    // Relpace the un-linked terms with the core-Fock matrix eleemnts
    void replace_Fock();
    void replace_Fock2();

    // Contract D4 and V2 to form D4C tensor
    void construct_D4C();

    // Regulate the index matching 
    void regulate_indices();

    // Set and regulate the index matching 
    void set_and_regulate_indices(std::vector<std::vector<Femto::Core::SQbinary> > &theBins);

    // Decompose a stream of Femto::Core::SQterms into a group of Femto::Core::SQbinary
    void binary_decomposition(std::vector<Femto::Core::SQterm> &theTerms, std::vector<Femto::Core::SQterm> &theLTensors);

    // Generate Orz code
    std::string generate(const Codes code_name=Orz, const Modes mode_name=SimpleLoops);

    // Analysis tool that returns the terms with 4-RDMs
    std::vector<Femto::Core::SQterm> get_terms_d4() const;

    std::vector<Femto::Core::SQterm> get_inTerms() const;
    size_t num_inTerms() const;

    // Contents of the header file
    std::string cont_header_;

    // Whether print guard not to calculate redundant contraction if there are no orbitals in specific orbital groups
    bool guard_core_; 
    bool guard_act_;

    ///////// Interfaces to set up parameters /////////
    void set_oldgemm     (const bool flag);
    void set_priority    (const priority flag);
    void set_h2_prior    (const bool flag);
    void set_isOF        (const bool flag);
    void set_use_cumulant(const bool flag);
    void set_use_omp     (const bool flag);
    void set_use_gemm    (const bool flag);
    void set_do_timing   (const bool flag);
    void set_name_h1     (const std::string &name);
    void set_V_D         (const bool flag);
    void set_T_D         (const bool flag); // <-- NEW 2013/01/08
    void set_D4C         (const bool flag); // <-- NEW 2012/12/11

    void set_guard_core  (const bool flag); // <-- NEW 2013/01/10
    void set_guard_act   (const bool flag); // <-- NEW 2013/01/10

    // Initialize the content of header files included from the C++ source
    void init_header(std::string file_name);
 
    void set_numCore(const int num);
    void set_numOcc (const int num);
    void set_numVirt(const int num);
    ///////////////////////////////////////////////////

#ifdef CASYANAI
    // CAS (yanai) //////////////////////////////////////////////
    void CAS_generate(std::map<int, std::vector<std::vector<Femto::Core::SQbinary> > > &bins,
                      std::map<int, std::vector<std::vector<int> > > &scl1,
                      std::map<int, std::vector<std::vector<int> > > &scl2,
                      std::map<int, std::vector<std::vector<int> > > &scl3);
    void CAS_factorize_new(std::map<int, std::vector<std::vector<Femto::Core::SQbinary> > > &bins,
                           std::map<int, std::vector<std::vector<int> > > &scl1,
                           std::map<int, std::vector<std::vector<int> > > &scl2,
                           std::map<int, std::vector<std::vector<int> > > &scl3);
    void CAS_regulate_indices();
    void CAS_decomp(std::vector<Femto::Core::SQterm> &theTerms,
                    std::map<int, std::vector<std::vector<Femto::Core::SQbinary> > > &bins,
                    std::map<int, std::vector<std::vector<int> > > &scl1,
                    std::map<int, std::vector<std::vector<int> > > &scl2,
                    std::map<int, std::vector<std::vector<int> > > &scl3);
    // END OF CAS (yanai) ///////////////////////////////////////
#endif

  friend std::ostream &operator<<(std::ostream &os, const SQreaktor &r);
      
  private:

    /////////////////////// Body of generating function ////////////////////  
    // Convert tensorial terms into factorized Orz code ....
    // OLD verison 
    void factorize();
    // Femto::Core::SQbinary is utilized throughout this version to carry binary contractions
    std::string factorize_new();
    // Femto::Core::SQbinary is utilized throughout this version to carry binary contractions
    std::string factorize_fafnir();
    // In this version, input tensor (such as T1 amplitude) are assumed to be DTensor type object
    std::string factorize_cas();

    // Main function to generate comutational code for Orz
    void generate_contract(std::ofstream &CPfile, std::ofstream &CHfile, std::ofstream &F90file, std::vector<std::vector<Femto::Core::SQbinary> > &theBins, std::string TheLabel);

    // Main function to generate comutational code for Orz
    void generate_contract_new(std::ofstream &CPfile, std::ofstream &CHfile, std::ofstream &F90file, std::vector<std::vector<Femto::Core::SQbinary> > &theBins, std::string TheLabel);

    // Small function just to count and print the order of the contractions
    void count_order(std::vector<Femto::Core::SQbinary> &bin);

    // Convert tensorial terms into Orz code, in which tensorial contraction is carried out in simple loops .... 
    std::string simpleloops();
    ////////////////////////////////////////////////////////////////////////
    //////////////////////// Small utilities ///////////////////////////////
    // Returns symmetry for contracted D4 tensor (D4C)
    Femto::Symmetry d4c_symm();

    // Returns symmetry for contracted D4 tensor (C6)
    Femto::Symmetry c6_symmM();

    // Loop in C++ style (/orz)
    void CLoop(const std::string s, const Femto::Core::SQindex &i, std::ofstream &f);

    // Parallelization statement in C++ style (/orz)
    void Cparallel(const std::string s, const Femto::Core::SQindex &i, std::ofstream &f);

    // MPIipric in C++ style (/orz)
    void MPIiproc(const std::string s, const Femto::Core::SQindex &i, std::ofstream &f);
  
    // Read BareAmp from GA (/orz)
    void ReadAmp(const std::string s, const Femto::Core::SQtensor &t, std::ofstream &f);
  
    // Read retval (/orz)
    void ReadRetval(const std::string s, const Femto::Core::SQtensor &t, std::ofstream &f);
  
    // Read ERIs from GA (/orz)
    void ReadERI(const std::string s, const Femto::Core::SQtensor &t, std::ofstream &f);

    // Read D4Cs from GA (/orz)
    void ReadD4C(const std::string s, const Femto::Core::SQtensor &t, std::ofstream &f);
  
    // Read D4 from GA (/orz)
    void ReadD4(const std::string s, const Femto::Core::SQtensor &t, std::ofstream &f);
  
    // Generate D4 by cumulant expansion (/orz)
    void ReadD4_Cumulant(const std::string s, const Femto::Core::SQtensor &t, std::ofstream &f);
  
    // Accumulate BareAmp (/orz)
    void AccAmp(const std::string s, const Femto::Core::SQtensor &t, std::ofstream &f);
  
    // Print end of the loop (/orz)
    void LoopEnd(std::string s, std::ofstream &f);

    // Print end of the loop (/orz)
    void LoopEnd_i(std::string s, const Femto::Core::SQindex &i, std::ofstream &f);
    ////////////////////////////////////////////////////////////////////////
    //////////////////////// Worker functions //////////////////////////////
    // Print the calling section of the Fortran routine in C++ side (/orz)
    void makeCPP_body(Femto::Core::SQtensor &LTensor, Femto::Core::SQterm &inTerm, std::string title, std::ofstream &f, std::string s, contDecl c, bool isBareLHS);

    // Modified version based on Femto::Core::SQbinary object
    void makeCPP_body2(Femto::Core::SQbinary &bin, std::string title, std::ofstream &f, std::string s, contDecl c, bool isBareLHS);

    // // Modified version based on Femto::Core::SQbinary object
    // void makeCPP_body2_new(Femto::Core::SQbinary &bin, std::string title, std::ofstream &f, std::string s, contDecl c, bool isBareLHS);

    // Modified version for the ClipSource.hs for minimizing the loading ERIs from disk
    void makeCPP_bodyType2(Femto::Core::SQbinary &bin, std::string title, std::ofstream &f, std::string s, contDecl c, bool isBareLHS);

    // Modified version for the ClipSource.hs for minimizing the loading ERIs from disk
    void makeCPP_bodyType2_new(Femto::Core::SQbinary &bin, std::string title, std::ofstream &f, std::string s, contDecl c, bool isBareLHS);
  
    // Print header of the Fortran routine in C++ side (/orz)
    void makeCPP_header(Femto::Core::SQtensor &LTensor, Femto::Core::SQterm &inTerm, std::string title, std::ofstream &f, std::string s, contDecl c, bool isBareLHS);

    // Modified version based on Femto::Core::SQbinary object
    void makeCPP_header2(Femto::Core::SQbinary &bin, std::string title, std::ofstream &f, std::string s, contDecl c, bool isBareLHS);

    // Modified version based on Femto::Core::SQbinary object
    void makeCPP_header2_new(Femto::Core::SQbinary &bin, std::string title, std::ofstream &f, std::string s, contDecl c, bool isBareLHS);
  
    // Print interface of the F90 subroutine (/orz)
    void makeF90_interface(Femto::Core::SQtensor &LTensor, Femto::Core::SQterm &inTerm, std::string title, std::ofstream &f, std::string s, contDecl c, bool isBareLHS);

    // Modified version based on Femto::Core::SQbinary object
    void makeF90_interface2(Femto::Core::SQbinary &bin, std::string title, std::ofstream &f, std::string s, contDecl c, bool isBareLHS);

    // Modified version based on Femto::Core::SQbinary object
    void makeF90_interface2_new(Femto::Core::SQbinary &bin, std::string title, std::ofstream &f, std::string s, contDecl c, bool isBareLHS);

    // Print the body of tensorial contraction (/orz)
    void makeF90_contract(Femto::Core::SQtensor &LTensor, Femto::Core::SQterm &inTerm, std::string title, std::ofstream &f, std::string s, contDecl c, bool isBareLHS);

    // Modified version based on Femto::Core::SQbinary object
    void makeF90_contract2(Femto::Core::SQbinary &bin, std::string title, std::ofstream &f, std::string s, contDecl c, bool isBareLHS);

    // Modified version based on Femto::Core::SQbinary object
    void makeF90_contract2_new(Femto::Core::SQbinary &bin, std::string title, std::ofstream &f, std::string s, contDecl c, bool isBareLHS);
   
    // Print the body of tensorial contraction processed with dgemm(/orz)
    void binary_contract(Femto::Core::SQtensor &LTensor, Femto::Core::SQterm &inTerm, std::string title, std::ofstream &f, std::string s, contDecl c, bool isBareLHS);

    // Modified version based on Femto::Core::SQbinary object
    void binary_contract2(Femto::Core::SQbinary &bin, std::string title, std::ofstream &f, std::string s, contDecl c, bool isBareLHS);

    // Memory swapping is modified
    void binary_contract3(Femto::Core::SQbinary &bin, std::string title, std::ofstream &f, std::string s, contDecl c, bool isBareLHS);

    // Memory swapping is modified
    void binary_contract3_new(Femto::Core::SQbinary &bin, std::string title, std::ofstream &f, std::string s, contDecl c, bool isBareLHS);

    ////////////////////////////////////////////////////////////////////////
    ///////////////////////// Small utilities //////////////////////////////
    void declareInterm(int myIndent, std::ofstream &CPfile, Femto::Core::SQtensor ten, Femto::SQcont<Femto::Core::SQindex*> decInds, decType myType=Internal);
    void makeContractions(int myIndent, std::ofstream &CPfile, std::ofstream &CHfile, std::ofstream &F90file, std::string myLabel, std::string theLabel, std::vector<Femto::Core::SQbinary> &myBins, Femto::Reaktor::ParaFlag myFlag=Femto::Reaktor::NonPara);
    Femto::Core::SQindex findERIIndex(std::vector<Femto::Core::SQbinary> &contras) const;
    Femto::Core::SQindex findD4CIndex(std::vector<Femto::Core::SQbinary> &contras) const;
    Femto::SQcont<Femto::Core::SQindex> returnsExtIndices(Femto::Core::SQbinary &contra) const;
    //////////////////////////////////////////////////////////////////////// 

    // Whether to use old dgemm generator or not
    bool use_oldgemm_;

    priority is_priority_;

    // Whether ERI is processed with the top priority
    bool h2_prior_;

    // Whether intermediate is processed in ad hoc fashoion in every situation
    bool isOF_;

    // Number of core orbitals in the model systems 
    int numCore_;

    // Number of active  orbitals in the model systems 
    int numOcc_;
    
    // Number of virtual orbitals in the model systems 
    int numVirt_;

    // Whether cumulant decomposition is turned on
    bool use_cumulant_;

    // Whether openmp is turned on
    bool use_omp_;

    // Whether dgemm is utilized
    bool use_gemm_;

    // Whether time the tensor contraction
    bool do_timing_;

    // Whether the legs of D4 and ERI are swap to match to each other
    bool Set_V_D_;

    // Whether the legs of T2 and ERI are swap to match to each other
    bool Set_T_D_;

    // Whether to contract D4 with V2 before the tensor contraction steps
    bool Set_D4C_;

    // Maximum length of the tensor that can be processed in ad hoc fashion
    unsigned long int MaxSize_;

    // External index of ERI
    const int exth2_;

    // External index of BareAmp
    const int extamp_;

    // External index of D4c
    const int extd4c_;

    // Whether LHS is a bareamp ot npt
    bool isBareLHS_;

    // Title of this scheme
    std::string title_;

    std::string name_h1_;  // Name of one-body integral 
    std::string name_h2_;  // Name of two-body integral
    std::string name_amp_; // Name of the anplitude
    std::string name_d4_;  // Name of 4-RDM 
 
    Femto::Core::SQtensor LTensor_;             // Tensor on the left-hand side
    std::vector<Femto::Core::SQterm> inTerms_ ; // Terms on the right-hand side
    std::vector<Femto::Core::SQterm> LTensors_; // List of LTensor, each member may be same to LTensor_ in case of sigma, or overlap vector, but in case of
                                   // diagonal preconditionor, this may be different in indices due to the existence of Kronecker's delta
    std::vector<Femto::Core::SQindex*> X_indices_; // *OBSOLETE* Indices of X intermediate for the time it is not processed on the fly (now used only in Factorize.cc, Convert.cc)

    // Keys :: 
    //   ++ "noeri" : terms without ERIs
    //   ++ "eri_c" : terms with ERI loaded by core index
    //   ++ "eri_o" : terms with ERI loaded by active index
    //   ++ "eri_v" : terms with ERI loaded by virtual index
    //   ++ "d4c_c" : terms with D4C loaded by core index
    //   ++ "d4c_o" : terms with D4C loaded by active index
    //   ++ "d4c_v" : terms with D4C loaded by virtual index
    std::map<std::string, std::vector<Femto::Core::SQterm> > type_terms_;    // inTerms_ classified by type of loading index of ERI
    std::map<std::string, std::vector<Femto::Core::SQterm> > type_LTensors_; // LTensors_ classified by type

    std::vector<std::vector<Femto::Core::SQbinary> > binaries_; // Divided form of inTerms into group of the binary contractions

    // Keys :: 
    //   ++ "noeri" : terms without ERIs
    //   ++ "eri_c" : terms with ERI loaded by core index
    //   ++ "eri_o" : terms with ERI loaded by active index
    //   ++ "eri_v" : terms with ERI loaded by virtual index
    //   ++ "d4c_c" : terms with D4C loaded by core index
    //   ++ "d4c_o" : terms with D4C loaded by active index
    //   ++ "d4c_v" : terms with D4C loaded by virtual index
    std::map<std::string, std::vector<std::vector<Femto::Core::SQbinary> > > FafBinaries_; // Binary contractions generated by Fafnir module

    // Some stuffs for serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int ver){
      ar & boost::serialization::make_nvp("guard_core",   guard_core_);
      ar & boost::serialization::make_nvp("guard_act",    guard_act_);
      ar & boost::serialization::make_nvp("is_priority",  is_priority_);
      ar & boost::serialization::make_nvp("h2_prior",     h2_prior_);
      ar & boost::serialization::make_nvp("isOF",         isOF_);
      ar & boost::serialization::make_nvp("numCore",      numCore_);
      ar & boost::serialization::make_nvp("numOcc",       numOcc_);
      ar & boost::serialization::make_nvp("numVirt",      numVirt_);
      ar & boost::serialization::make_nvp("use_cumulant", use_cumulant_);
      ar & boost::serialization::make_nvp("use_omp",      use_omp_);
      ar & boost::serialization::make_nvp("use_gemm",     use_gemm_);
      ar & boost::serialization::make_nvp("do_timing",    do_timing_);
      ar & boost::serialization::make_nvp("Set_V_D",      Set_V_D_);
      ar & boost::serialization::make_nvp("Set_T_D",      Set_T_D_);
      ar & boost::serialization::make_nvp("MaxSize",      MaxSize_);
      ar & boost::serialization::make_nvp("exth2",        exth2_);
      ar & boost::serialization::make_nvp("extamp",       extamp_);
      ar & boost::serialization::make_nvp("extd4c",       extd4c_);
      ar & boost::serialization::make_nvp("isBareLHS",    isBareLHS_);
      ar & boost::serialization::make_nvp("title",        title_);
      ar & boost::serialization::make_nvp("name_h1",      name_h1_);
      ar & boost::serialization::make_nvp("name_h2",      name_h2_);
      ar & boost::serialization::make_nvp("name_amp",     name_amp_);
      ar & boost::serialization::make_nvp("name_d4",      name_d4_);
      ar & boost::serialization::make_nvp("LTensor",      LTensor_);
      ar & boost::serialization::make_nvp("inTerms",      inTerms_);
      ar & boost::serialization::make_nvp("X_indices",    X_indices_);
    }

  };
  
}} // Femto::reaktor

