//
//  SQinput.hpp
//  
//  Class for the information container as input for SQreaktor object
//
//  Created by Masaaki Saitow on 12/11/13.
//  Copyright (c) 2012 Masaaki Saitow. All rights reserved.
//

#pragma once

#include <cstring>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/nvp.hpp>
#include <map>
#include <Femto.hpp>

namespace Femto { namespace Reaktor2 {
  
  // Priority flag
  enum priority {polynomial=0, memory=1};
  // Code generation paradigms of generator
  enum Modes {SimpleLoops=0, Factorize=1};
  // Target program
  enum Codes {Orz=0, Psi=1};
  // Container of information of declaration
  typedef std::map<std::string, std::vector<std::string> > contDecl; 

  class SQinput{
      
  public:
    SQinput(const std::string title="default_scheme", const bool isBareLHS=false);
    ~SQinput();

    ///////// Interfaces to set up parameters /////////
    // Paramaters for code generation
    void set_priority       (const priority flag);
    void set_use_cumulant   (const bool flag);
    void set_use_omp        (const bool flag);
    void set_use_gemm       (const bool flag);
    void set_do_timing      (const bool flag);
    void set_V_D            (const bool flag);

    priority get_priority       () const;
    bool get_use_cumulant   () const;
    bool get_use_omp        () const;
    bool get_use_gemm       () const;
    bool get_do_timing      () const;
    bool get_V_D            () const;

    // Names of various tensor
    void set_name_h1        (const std::string &name);
    void set_name_h2        (const std::string &name);
    void set_name_t1        (const std::string &name);
    void set_name_t2        (const std::string &name);

    std::string get_name_h1 () const;
    std::string get_name_h2 () const;
    std::string get_name_t1 () const;
    std::string get_name_t2 () const;

    // Loading index for ERI and BareAamp
    void set_exth2          (const int i);
    void set_extamp         (const int i);

    int get_exth2           () const;
    int get_extamp          () const;

    // Numbers that define the size of target system
    void set_numCore        (const int num);
    void set_numOcc         (const int num);
    void set_numVirt        (const int num);

    int get_numCore         () const;
    int get_numOcc          () const;
    int get_numVirt         () const;
    ///////////////////////////////////////////////////

  friend std::ostream &operator<<(std::ostream &os, const SQinput &r);
      
  private:

    priority is_priority_;

    // Whether cumulant decomposition is turned on
    bool use_cumulant_;

    // Whether openmp is turned on
    bool use_omp_;

    // Whether dgemm is utilized
    bool use_gemm_;

    // Whether time the tensor contraction
    bool do_timing_;

    // Whether the legs of D4 and ERI are swap to match to each other
    bool set_V_D_;

    // Whether LHS is a bareamp ot npt
    bool isBareLHS_;

    // Title of this scheme
    std::string title_;

    std::string name_h1_;  // Name of one-body integral 
    std::string name_h2_;  // Name of two-body integral 
    std::string name_t1_;  // Name of t1 amplitude 
    std::string name_t2_;  // Name of t2 amplitude (treated as BareAmp)
    std::string name_amp_; // Name of the anplitude
    std::string name_d4_;  // Name of 4-RDM 

    // External index of ERI
    int exth2_;

    // External index of BareAmp
    int extamp_;

    // Number of core orbitals in the model systems 
    int numCore_;

    // Number of active  orbitals in the model systems 
    int numOcc_;
    
    // Number of virtual orbitals in the model systems 
    int numVirt_;

    // Maximum length of the tensor that can be processed in ad hoc fashion
    unsigned long int MaxSize_;
 
    // Some stuffs for serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int ver){
      ar & boost::serialization::make_nvp("numCore",      numCore_);
      ar & boost::serialization::make_nvp("numOcc",       numOcc_);
      ar & boost::serialization::make_nvp("numVirt",      numVirt_);
      ar & boost::serialization::make_nvp("use_cumulant", use_cumulant_);
      ar & boost::serialization::make_nvp("use_omp",      use_omp_);
      ar & boost::serialization::make_nvp("use_gemm",     use_gemm_);
      ar & boost::serialization::make_nvp("do_timing",    do_timing_);
      ar & boost::serialization::make_nvp("Set_V_D",      set_V_D_);
      ar & boost::serialization::make_nvp("MaxSize",      MaxSize_);
      ar & boost::serialization::make_nvp("exth2",        exth2_);
      ar & boost::serialization::make_nvp("extamp",       extamp_);
      ar & boost::serialization::make_nvp("isBareLHS",    isBareLHS_);
      ar & boost::serialization::make_nvp("title",        title_);
      ar & boost::serialization::make_nvp("name_h1",      name_h1_);
      ar & boost::serialization::make_nvp("name_h2",      name_h2_);
      ar & boost::serialization::make_nvp("name_t1",      name_t1_);
      ar & boost::serialization::make_nvp("name_t2",      name_t2_);
      ar & boost::serialization::make_nvp("name_amp",     name_amp_);
      ar & boost::serialization::make_nvp("name_d4",      name_d4_);
    }

  };
  
}} // Femto::input

