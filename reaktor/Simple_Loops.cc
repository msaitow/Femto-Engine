//
//  Simple_Loops.cc
//  
//
//  Created by Masaaki Saitow on 12/09/04.
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
  string SQreaktor::simpleloops()
  {
    if(LTensor_.get_indices().size() != 4 && isBareLHS_){
      cout << "Argument 0 cannnot be treated as a bareampack .... " << endl;
      abort();
    } // End if

    vector<SQindex*> Linds = LTensor_.get_indices();
    for(vector<SQindex*>::iterator i = Linds.begin();i != Linds.end();++i)
      if((*i)->get_isSummed()){
        cout << *i << " in 0th argument is a kind of dummy index." << endl;
        abort(); 
      } // End if

    for(vector<SQterm>::iterator t = inTerms_.begin();t != inTerms_.end();++t){
      if(!t->get_isInCanonical()){
        cout << "Term, " << *t << " is not in canonical form .... " << endl;
        abort();
      } // End if
    } // End t

    if(is_sfGen(LTensor_.get_name())){
      cout << "Factorize: 1st argument cannot be a spin-free unitary group generator ... " << endl;
      abort();
    }
    for(vector<SQterm>::iterator t = inTerms_.begin();t != inTerms_.end();++t){
      for(size_t num_t = 0;num_t < t->get_tensors().size();++num_t){
        if(is_sfGen(t->get_tensors()[num_t].get_name())){
          cout << "Factorize: A spin-free unitary group generator is detected in the 2nd argument .... " << endl;
          abort();
	} // End if
      } // End num_t
    } // End t

    if(!inTerms_.size()) return "";

    // Header file for C++ 
    string CHname("c_" + title_ + ".h");
    string CHguard("c_" + title_ + "_h");
    // It's a kind of trick, though .... 
    transform(CHguard.begin(), CHguard.end(), CHguard.begin(), (int(*)(int))toupper); 
    ofstream CHfile( CHname.c_str() );
    CHfile << "\n" << endl;
    CHfile << "#ifndef " << CHguard << endl;
    CHfile << "#define " << CHguard << endl;
    CHfile << "" << endl;
    CHfile << "// #include <tensor/tensor.h>                  " << endl;  
    CHfile << "// #include <sci/hint/hintmo/hintmo.h>         " << endl;  
    CHfile << "// #include <sci/ctnew2/ctclass_input.h>       " << endl;   
    CHfile << "// #include <sci/ctnew2/ctclass_symblock.h>    " << endl;   
    CHfile << "// #include <sci/ctnew2/ctclass_rdmpack.h>     " << endl;   
    CHfile << "// #include <sci/ctnew2/ctclass_bareamppack.h> " << endl;       
    CHfile << "                                               " << endl;            
    CHfile << "extern \"C\"{                                  " << endl;      
    CHfile << "                                               " << endl;           
    CHfile << "                                               " << endl;            
    CHfile  << Femto_logo("// ") << endl; 

    string CHend("");
    CHend += "      \n ";
    CHend += "}     \n ";       
    CHend += "      \n ";  
    CHend += "      \n ";      
    CHend += "#endif\n ";       
    CHend += "      \n ";        
    CHend += "      \n ";       

    // Output for tensorial contractions .... 
    string Fname = "f_" + title_ + ".F90";
    ofstream F90file( Fname.c_str() );
    F90file << "#include <sci/icmr/fsrc/f_mr.fh>\n\n" << endl;    
    F90file << Femto_logo("! ")              << endl;

    // Output for the body ob tensorial contractions .... 
    bool Bareflag = false;
    bool D4flag   = false;
    for(vector<SQterm>::iterator t = inTerms_.begin();t != inTerms_.end();++t){
      for(size_t num_t = 0;num_t < t->get_tensors().size();++num_t){
        if(t->get_tensors()[num_t].get_name() == name_amp_) Bareflag = true; 
        if(t->get_tensors()[num_t].get_name() == name_d4_)  D4flag   = true; 
      } // End t    
    } // End term

    // inTerms_ have constants or not
    bool Constsflag = false;
    vector<string> Consts;
    for(vector<SQterm>::iterator t = inTerms_.begin();t != inTerms_.end();++t){
      for(size_t num_c = 0;num_c < t->get_Consts().size();++num_c)
        if(find(Consts.begin(), Consts.end(), t->get_Consts()[num_c]) == Consts.end() && t->get_Consts()[num_c] != "")
          Consts.push_back(t->get_Consts()[num_c]);
    } // End t
    if(Consts.size()) Constsflag = true;

    string CPname = "c_" + title_ + ".cpp";
    ofstream CPfile( CPname.c_str() );
    CPfile << "                                                                                " << endl;      
    CPfile << "#include <orz/orz.h>                                                            " << endl;    
    CPfile << "#include <orz/openmp.h>                                                         " << endl;     
    CPfile << "#include <orz/cblas.h>                                                          " << endl;      
    CPfile << "#include <orz/clapack.h>                                                        " << endl;   
    CPfile << "#include <tensor/tensor.h>                                                      " << endl;     
    CPfile << "#include <sci/hint/para_disttools.h>                                            " << endl; 
    CPfile << "#include <sci/ctnew2/ct.h>                                                      " << endl; 
    CPfile << "#include <sci/ctnew2/ct_f.h>                                                    " << endl;  
    CPfile << "#include <sci/ctnew2/ctclass_input.h>                                           " << endl;     
    CPfile << "#include <sci/ctnew2/ctclass_symblock.h>                                        " << endl;     
    CPfile << "#include <sci/ctnew2/ctclass_hintmo.h>                                          " << endl;  
    CPfile << "#include <sci/ctnew2/ctclass_rdmpack.h>                                         " << endl;     
    CPfile << "#include <sci/ctnew2/ctclass_bareamppack.h>                                     " << endl;  
    CPfile << "#include <sci/ctnew2/ctclass_orthamppack.h>                                     " << endl;      
    CPfile << "#include <sci/ctnew2/diaghessian.h>                                             " << endl;      
    CPfile << "#include <sci/ctnew2/symamp2.h>                                                 " << endl;    
    CPfile << "#include <sci/ctnew2/mrci.h>                                                    " << endl;     
    CPfile << "#include <sci/ctnew2/" + CHname + ">                                            " << endl;
    CPfile << "                                                                                " << endl;        
    CPfile << "using std::cout;                                                                " << endl;      
    CPfile << "using std::endl;                                                                " << endl;  
    CPfile << "                                                                                " << endl;      
    CPfile << "#define FLOPCOUNT                                                               " << endl;  
    CPfile << "                                                                                " << endl;   
    CPfile << Femto_logo("// ")                                                                  << endl;
    CPfile << "                                                                                " << endl;
    CPfile << "// ***************************************************************************  " << endl;           
    CPfile << "// orz::mr::mrci                                                                " << endl;      
    CPfile << "// ***************************************************************************  " << endl;          
    CPfile << "									    /*!        " << endl; 
    CPfile << "   @brief CT input                                                              " << endl;    
    CPfile << "									     */        " << endl;    
    CPfile << "                                                                                " << endl;     

    if(isBareLHS_)
    CPfile << "orz::mr::BareAmpPack orz::mr::" + title_ + "(const orz::mr::Input &ctinp,                                    " << endl;
    else if(!LTensor_.get_indices().size() && Bareflag)
    CPfile << "double orz::mr::" + title_ + "(const orz::mr::Input &ctinp,                                                  " << endl;
    CPfile << "					 const orz::mr::SymBlockInfo &symblockinfo,                                " << endl;
    CPfile << "					 const orz::mr::HintMO &hintmo,                                            " << endl;
    if(isBareLHS_)
    CPfile << "                                  const int alloc_type,                                            " << endl;
    if(D4flag){
    CPfile << "					 const orz::mr::RdmPack &rdmPack,                                          " << endl;
    CPfile << "					 const orz::DTensor &rdm4,                                                 " << endl;
    }
    if(!LTensor_.get_indices().size() && !isBareLHS_)
    CPfile << "      			         const double init_value,                                                  " << endl;
    if(Bareflag)
    CPfile << "					 const orz::mr::BareAmpPack &" + name_amp_ + ",                             " << endl;
    CPfile << "					 const int num_sigma";
    if(Consts.size()) CPfile << "," << endl;

    for(vector<string>::iterator c = Consts.begin();c != Consts.end();++c){
    CPfile << "					 const double " + *c;
    if(*c != Consts.back()) CPfile << "," << endl;     
    } // End c
    CPfile << ") {" << endl;     

    //////////////////////////////////////////////////////////////////////////////////////////////
    ostringstream retVal("");
    { // Construction of the header file
      if(isBareLHS_)
      retVal << "orz::mr::BareAmpPack " + title_ + "(const orz::mr::Input &ctinp,                                    " << endl;
      else if(!LTensor_.get_indices().size() && Bareflag)
      retVal << "double " + title_ + "(const orz::mr::Input &ctinp,                                                  " << endl;
      else if(LTensor_.get_name() == D4C_nameL())
      retVal << "void " + title_ +  "(const orz::mr::Input &ctinp,                                                     " << endl;
      else if(!isBareLHS_ && LTensor_.get_indices().size())
      retVal << "orz::DTensor " + title_ +  "(const orz::mr::Input &ctinp,                                                     " << endl;
      else{
	cout << "Cannot detect the output type .... " << endl;
	abort();
      }
      retVal << "                                  const orz::mr::SymBlockInfo &symblockinfo,                                 " << endl;
      retVal << "                                  const orz::mr::HintMO &hintmo,                                             " << endl;
      if(isBareLHS_)
      retVal << "                                  const int alloc_type,                                                      " << endl;
      if(D4flag){
      retVal << "                                  const orz::mr::RdmPack &rdmPack,                                           " << endl;
      retVal << "                                  const orz::DTensor &rdm4,                                                  " << endl;
      }
      if(!LTensor_.get_indices().size() && !isBareLHS_)
      retVal << "                                  const double init_value,                                                   " << endl;
      if(Bareflag)
      retVal << "                                  const orz::mr::BareAmpPack &" + name_amp_ + ",                             " << endl;
      retVal << "                                  const int num_sigma";
      if(Consts.size()) retVal << "," << endl;
      
      for(vector<string>::iterator c = Consts.begin();c != Consts.end();++c){
      retVal << "                                  const double " + *c;
      if(*c != Consts.back()) retVal << "," << endl;     
      } // End c
      retVal << ");" << endl << endl;           
    } // End scope
    //////////////////////////////////////////////////////////////////////////////////////////////

    CPfile << "                                                                                                                " << endl;         
    CPfile << "                                                                                                                " << endl;
    CPfile << "  // set up nmo nclosed, nocc                                                                                   " << endl;    
    CPfile << "  const FC_INT nclosed = ctinp.nclosed();                                                                       " << endl;         
    CPfile << "  const FC_INT nocc    = ctinp.nocc();                                                                          " << endl;                
    CPfile << "  const FC_INT nvir    = ctinp.nvir();                                                                          " << endl;                  
    CPfile << "  const FC_INT nmo     = nclosed + nocc + nvir;                                                                 " << endl;               
    CPfile << "  const FC_INT nir     = symblockinfo.nir();                                                                    " << endl;               
    CPfile << "  const FC_INT * const nsym    = symblockinfo.nsym().cptr();                                                    " << endl;                       
    CPfile << "  const FC_INT * const psym    = symblockinfo.psym().cptr();                                                    " << endl;                   
    CPfile << "  const FC_INT * const amo2imo = symblockinfo.amo2imo().cptr();                                                 " << endl;               
    CPfile << "                                                                                                                " << endl;
    if(isBareLHS_){
    CPfile << "  std::ostringstream stm;                                                                                       " << endl;                       
    CPfile << "  stm << num_sigma;                                                                                             " << endl;                        
    CPfile << "  std::string name_of_sigma = \"" + LTensor_.get_name() +  "\" + stm.str() + \"]\"; // Name of the Sigma vector " << endl;                              
    CPfile << "  orz::mr::BareAmpPack retval                                                                                   " << endl;                         
    CPfile << "    = orz::mr::BareAmpPack(ctinp, symblockinfo, name_of_sigma, alloc_type); // Sigma(a, a', e, e') tensor                   " << endl;                 
    CPfile << "                                                                                                                " << endl;                           
    CPfile << "  orz::DTensor " + LTensor_.get_name() + "b; // Container of S2_aae,[b] tensor                                  " << endl;                       
    CPfile << "                                                                                                                " << endl;
    } // End if
    else if(!LTensor_.get_indices().size() && Bareflag)
    CPfile << "  double " + LTensor_.get_name() + " = init_value;                                                              " << endl;
    if(Bareflag)
    CPfile << "  orz::DTensor " + name_amp_ + "b; // Container of T2_aae,[b] tensor                                            " << endl;                       
    if(D4flag){
    CPfile << "  orz::DTensor rdm4_sym;                                                                                        " << endl;                       
    CPfile << "  orz::DTensor rdm4_ij_sliced(ctinp.use_d4cum_of() ? nocc*nocc*nocc*nocc*nocc*nocc : 0);                        " << endl;
    }

    // Print OpenMPI parameters
    CPfile << "  // set nproc, myrank                      " << endl;                           
    CPfile << "  const int nproc = orz::world().size();    " << endl;                    
    CPfile << "  const int myrank = orz::world().rank();   " << endl;                     
    CPfile << endl;
    CPfile << "  // ==== set imopair_distlist, nmopair_distlist, and iproc_havingpairs ====                      " << endl;     
    CPfile << "  ivector imo_distlist;                                                                           " << endl;        
    CPfile << "  ivector nmo_distlist;                                                                           " << endl;         
    CPfile << "  ivector iproc_havingimo; // [i,j] and [j,i] are stored on the same node                         " << endl;               
    CPfile << "  boost::tie(imo_distlist, nmo_distlist, iproc_havingimo) = orz::hint::make_distlist(nmo, nproc); " << endl;                
    CPfile << "                                                                                                  " << endl;       
    CPfile << "  // my_imopair and my_nmopair                                                                    " << endl;       
    CPfile << "  const int my_imo = imo_distlist[myrank];                                                        " << endl;           
    CPfile << "  const int my_nmo = nmo_distlist[myrank];                                                        " << endl;            
    CPfile << "  orz::DTensor moint1 = hintmo.int1(); // Setting up one-body integrals                                         " << endl;                        
    CPfile << "  const orz::DTensor moint1_sym = (myrank == 0) ? orz::mr::sympack_int1(symblockinfo, moint1) : orz::DTensor(); // moint1=(IR-COV index)" << endl;         
    CPfile << "  orz::DTensor " + name_h2_ + "(nmo,nmo,nmo);                                                                    " << endl;                             
    CPfile << "  double * const " + name_h2_ + "_ptr = " + name_h2_ + ".cptr();                                                  " << endl;                          

    // Print logos
    cout << endl;
    cout << Femto_logo("");
    cout << endl;

    vector<string> Flags;
    vector<SQterm> Yeffs;
    int num_replaced = 0;
    for(vector<SQterm>::iterator t = inTerms_.begin();t != inTerms_.end();++t){
      vector<SQindex*>  indices(t->get_summedBody());
      vector<SQtensor>  tensors(t->get_tensors());
      vector<SQindex*>  summed;
      vector<int>       counters;
      vector<SQtensor*> holders;
      for(vector<SQindex*>::iterator i = indices.begin();i != indices.end();++i){
        if((*i)->get_isSummed()) summed.push_back(*i);
      } // End i
      // Count how many tensors share a dummy index i
      for(vector<SQindex*>::iterator i = summed.begin();i != summed.end();++i){
        int counter = 0;
        SQtensor* holder;
        for(vector<SQtensor>::iterator ten = tensors.begin();ten != tensors.end();++ten){
          vector<SQindex*> inds(ten->get_indices());
          if(find(inds.begin(), inds.end(), *i) != inds.end()) { ++counter; holder = &(*ten); }
	} // End ten
        counters.push_back(counter);
        holders.push_back(holder);
      } // End i
      for(size_t num_i = 0;num_i < summed.size();++num_i){
        if(counters[num_i] == 1){

          SQindex* i = summed[num_i];
          vector<SQtensor> ts;
          ts.push_back(*holders[num_i]);
          vector<string> coeff;
	  coeff.push_back("");
          Yeffs.push_back(SQterm(1.0, coeff, ts));
          vector<SQtensor> new_ten;
          for(vector<SQtensor>::iterator ten = tensors.begin();ten != tensors.end();++ten){
            if(!(*ten == *holders[num_i])) new_ten.push_back(*ten);  
	  } // End ten
          vector<SQindex*> new_inds(holders[num_i]->get_indices());
          for(vector<SQindex*>::iterator j = new_inds.begin();j != new_inds.end();){
            vector<SQindex*>::iterator j_ptr = find(summed.begin(), summed.end(), *j);
            
            if(j_ptr != summed.end()){
              if(counters[(size_t)(j_ptr-summed.begin())] == 1) 
                j = new_inds.erase(j);
              else ++j;
	    } // End if
            else ++j;
	  } // End j

          ostringstream stm;
          stm << num_replaced++;
          // Currently the permutational symmetry of the reduced tensor is not considered 
          if(new_inds.size()){
            Symmetry my_symm;
            if     (new_inds.size() == 2) my_symm = u2_symm();
            else if(new_inds.size() == 4) my_symm = u4_symm();
            new_ten.push_back(SQtensor("Y"+stm.str(), new_inds, my_symm));
	  } // End is
          else{
            vector<string> my_coeff(t->get_Consts());
            my_coeff.push_back("Y" + stm.str());
            t->set_Consts(my_coeff);
	  } // End else
          t->set_tensors(new_ten);
          break;

	} // End if
      } // End num_i      
    } // End t

    cout << endl;
    cout << "! * " << num_replaced << " terms are replaced in the linking process .... " << endl << endl;;
    if(num_replaced){
      cout << endl;
      int count = 0;
      cout << "The linked formulas .... " << endl;
      for(vector<SQterm>::iterator t = inTerms_.begin();t != inTerms_.end();++t, ++count){
        cout << count << " : " << *t << endl;
      } // End t
      cout << endl;
      cout << "The content of each effective tensor .... " << endl;
      int num = 0;
      for(vector<SQterm>::iterator t = Yeffs.begin();t != Yeffs.end();++t, ++num)
        cout << "Y" << num << " <-- " << *t << endl;
      cout << endl << endl;
    } // End if

    // Before entering to the generation step, pick up the indices of LTensor_ as a reference
    vector<SQindex*> ref_indices = LTensor_.get_indices();
    vector<SQindex> ref_bodies; ref_bodies.reserve(ref_indices.size());
    int cnt = 0;
    for(vector<SQindex*>::iterator i = ref_indices.begin();i != ref_indices.end();++i,++cnt){
      ref_bodies.push_back(**i);
      ref_indices[cnt] = &ref_bodies[cnt];
    } // End i

    const SQtensor LTdummy(LTensor_);

    // Convert term by term ....
    int numTerm = 0; 
    for(vector<SQterm>::iterator thisTerm = inTerms_.begin();thisTerm != inTerms_.end();++thisTerm, ++numTerm){
    // Starts from L.1646 in sqaConvert.py ... 
      cout << "! No."  << numTerm << endl;
      cout << "! " << LTdummy << " <-- " << endl;
      cout << "! " << *thisTerm << endl;

      // Retain the consistency between LTensor_ and thisTerm
      vector<SQindex*> Linds(ref_indices);
      vector<SQindex*> thisBoddies(thisTerm->get_summedBody());
      for(vector<SQindex*>::iterator i = Linds.begin();i != Linds.end();++i){
        //if(find(thisBoddies.begin(), thisBoddies.end(), *i) == thisBoddies.end())
        bool i_flag = true;
        for(vector<SQindex*>::iterator Bi = thisBoddies.begin();Bi != thisBoddies.end();++Bi)
          if(**i == **Bi) i_flag = false;
        if(i_flag) thisBoddies.push_back(*i);
      } // End i
      for(size_t num_i = 0;num_i < Linds.size();++num_i){
        SQindex* Li = Linds[num_i];
        for(vector<SQindex*>::iterator Bi = thisBoddies.begin();Bi != thisBoddies.end();++Bi){
          //cout << **Bi << endl;
          if(*Li == *(*Bi)) LTensor_.put_indices(num_i, *Bi);
	} // End Bi
        //cout << endl;
      } // End Li

      // Break Kronecker's delta
      vector<SQtensor> new_t;
      vector<SQindex> newLTind;
      for(size_t i = 0;i < thisTerm->get_tensors().size();++i){
        SQtensor t = thisTerm->get_tensors()[i];
        if(t.get_name() != kDelta_name()) new_t.push_back(t);
        else{
          SQindex* i1 = t.get_indices()[0];
          SQindex* i2 = t.get_indices()[1];
          *i2 = *i1;
          newLTind.clear();
          for(size_t j = 0;j < LTensor_.get_indices().size();++j) newLTind.push_back(*LTensor_.get_indices()[j]);
	}
      } // End i

      //cout << LTensor_ << " <=== " << *thisTerm << endl;

      vector<SQindex> LTbodies; LTbodies.reserve(LTensor_.get_indices().size());
      vector<SQindex*> LT_ptrs(LTensor_.get_indices());
      int num_i = 0;
      for(vector<SQindex*>::iterator i = LT_ptrs.begin();i != LT_ptrs.end();++i, ++num_i){
        LTbodies.push_back(**i);
        LT_ptrs[num_i] = &LTbodies[num_i];
      } // End num_i

      if(new_t.size() != thisTerm->get_tensors().size()){
        thisTerm->set_tensors(new_t);
        thisBoddies = thisTerm->get_summedBody();

        // If so, consistency between LTensor_ and thisTerm has to be retain again
        for(vector<SQindex*>::iterator i = LT_ptrs.begin();i != LT_ptrs.end();++i){
          //if(find(thisBoddies.begin(), thisBoddies.end(), *i) == thisBoddies.end())
          bool i_flag = true;
          for(vector<SQindex*>::iterator Bi = thisBoddies.begin();Bi != thisBoddies.end();++Bi)
            if(**i == **Bi) i_flag = false;
          if(i_flag) thisBoddies.push_back(*i);
        } // End i
        for(size_t num_i = 0;num_i < LT_ptrs.size();++num_i){
          SQindex* Li = LT_ptrs[num_i];
          for(vector<SQindex*>::iterator Bi = thisBoddies.begin();Bi != thisBoddies.end();++Bi){
            //cout << **Bi << endl;
            if(*Li == *(*Bi)) LTensor_.put_indices(num_i, *Bi);
	  } // End Bi
        } // End Li
        cout << "! Kronecker's delta removed." << endl;
      } // End if

      // Achieve the best matching
      typedef pair<bool, SQtensor> conTen;
      conTen h2f, ampf, d4f;
      h2f.first  = false;
      ampf.first = false;
      d4f.first  = false;
      for(size_t num_t = 0;num_t < thisTerm->get_tensors().size();++num_t){
        if     (thisTerm->get_tensors()[num_t].get_name() == name_h2_){
          h2f.first  = true;
          h2f.second = (thisTerm->get_tensors()[num_t]);
	} // End if
        else if(thisTerm->get_tensors()[num_t].get_name() == name_amp_){
          ampf.first  = true;
          ampf.second = (thisTerm->get_tensors()[num_t]); 
          //cout << thisTerm->get_tensors()[num_t].get_indices()[extamp_];
	} // End if
        else if(thisTerm->get_tensors()[num_t].get_name() == name_d4_){
          d4f.first  = true;
          d4f.second = (thisTerm->get_tensors()[num_t]);
	} // End if
      } // End num_t

      // Match the external indices of the bareamp and eri, if possible.
      // If impossible, consider such the possibility of d4 and eri.
      bool SetT_V = false;
      bool SetL_V = false;
      bool SetL_T = false;
      bool SetV_D = false;
      bool SetV_v = false;
      // LHS v.s. T2
      if(ampf.first && isBareLHS_){
        vector<SQindex*> ampind(ampf.second.get_indices());
        vector<SQindex*> LTind(LTensor_.get_indices());
        for(vector<SQindex*>::iterator i = ampind.begin();i != ampind.end();++i){
          for(vector<SQindex*>::iterator j = LTind.begin();j != LTind.end();++j){
            if(**i == **j){
              IIvector ampperm(ampf.second.get_perms());
              IIvector LTperm(LTensor_.get_perms());
              typedef pair<bool, size_t> found;
              found amp, LT;
              amp.first = false;
              LT.first = false;
              amp.second = 0;
              LT.second = 0;
              for(size_t k = 0;k < ampperm.size();++k){
                if(*ampind[ampperm[k][extamp_]] == **i) { amp.first = true; amp.second = k; break; } 
                if(amp.first) break;
	      } //End k
              for(size_t k = 0;k < LTperm.size();++k){
                if(*LTind[LTperm[k][extamp_]] == **j) { LT.first = true; LT.second = k; break; } 
                if(LT.first) break;
	      } //End k

              if(amp.first && LT.first) {
                vector<SQtensor*> ten_ptr(thisTerm->get_tensors_ptr()); 
                for(vector<SQtensor*>::iterator t = ten_ptr.begin();t != ten_ptr.end();++t){
                  if     ((*t)->get_name() == name_amp_) 
                    (*t)->rotateIndices(amp.second);
	        } // End t
                LTensor_.rotateIndices(LT.second);
                cout << "! Indices of BareAmp are rotated to match with LHS. " << endl;
                SetL_T = true;
	      } // End if
              if(SetL_T) break;
	    } // End if
	  } // End j
          if(SetL_T) break;
	} // End i
      } // End if
      // LHS v.s. V2
      if(h2f.first && isBareLHS_){
        vector<SQindex*> h2ind(h2f.second.get_indices());
        vector<SQindex*> LTind(LTensor_.get_indices());
        for(vector<SQindex*>::iterator i = h2ind.begin();i != h2ind.end();++i){
          for(vector<SQindex*>::iterator j = LTind.begin();j != LTind.end();++j){
            if(**i == **j){
              IIvector h2perm(h2f.second.get_perms());
              IIvector LTperm(LTensor_.get_perms());
              typedef pair<bool, size_t> found;
              found h2, LT;
              h2.first = false;
              LT.first = false;
              h2.second = 0;
              LT.second = 0;
              for(size_t k = 0;k < h2perm.size();++k){
                if(*h2ind[h2perm[k][exth2_]] == **i) { h2.first = true; h2.second = k; break; } 
                if(h2.first) break;
	      } //End k
              for(size_t k = 0;k < LTperm.size();++k){
                if(*LTind[LTperm[k][extamp_]] == **j) { LT.first = true; LT.second = k; break; } 
                if(LT.first) break;
	      } //End k

              if(h2.first && LT.first) {
                vector<SQtensor*> ten_ptr(thisTerm->get_tensors_ptr()); 
                for(vector<SQtensor*>::iterator t = ten_ptr.begin();t != ten_ptr.end();++t){
                  if     ((*t)->get_name() == name_h2_) 
                    (*t)->rotateIndices(h2.second);
	        } // End t
                LTensor_.rotateIndices(LT.second);
                cout << "! Indices of ERI are rotated to match with LHS. " << endl;
                SetL_V = true;
	      } // End if
              if(SetL_V) break;
	    } // End if
	  } // End j
          if(SetL_V) break;
	} // End i
      } // End if
      // T2 v.s. V2
      if(h2f.first && ampf.first && !SetL_V || SetL_T){
        vector<SQindex*> h2ind(h2f.second.get_indices());
        vector<SQindex*> ampind(ampf.second.get_indices());
        for(vector<SQindex*>::iterator i = h2ind.begin();i != h2ind.end();++i){
          for(vector<SQindex*>::iterator j = ampind.begin();j != ampind.end();++j){
            if(**i == **j){
              IIvector h2perm(h2f.second.get_perms());
              IIvector ampperm(ampf.second.get_perms());
              typedef pair<bool, size_t> found;
              found h2, amp;
              h2.first = false;
              amp.first = false;
              h2.second = 0;
              amp.second = 0;
              for(size_t k = 0;k < h2perm.size();++k){
                if(*h2ind[h2perm[k][exth2_]] == **i) { h2.first = true; h2.second = k; break; } 
                if(h2.first) break;
	      } //End k
              for(size_t k = 0;k < ampperm.size();++k){
                if(*ampind[ampperm[k][extamp_]] == **j) { amp.first = true; amp.second = k; break; } 
                if(amp.first) break;
	      } //End k

              if(h2.first && amp.first) {
                vector<SQtensor*> ten_ptr(thisTerm->get_tensors_ptr()); 
                for(vector<SQtensor*>::iterator t = ten_ptr.begin();t != ten_ptr.end();++t){
                  if     ((*t)->get_name() == name_h2_) 
                    (*t)->rotateIndices(h2.second);
                  else if((*t)->get_name() == name_amp_) 
                    (*t)->rotateIndices(amp.second);
	        } // End t
                cout << "! Indices of ERI are rotated to match with Bareamp. " << endl;
                SetT_V = true;
	      } // End if
              if(SetT_V) break;
	    } // End if
	  } // End j
          if(SetT_V) break;
	} // End i
      } // End if
      // V2 v.s. D4
      if(h2f.first && d4f.first && !SetT_V && !SetL_V){
        vector<SQindex*> h2ind(h2f.second.get_indices());
        vector<SQindex*> d4ind(d4f.second.get_indices());
        //SQindex* the_first = h2ind[exth2_];
        for(vector<SQindex*>::iterator i = h2ind.begin();i != h2ind.end();++i){
          for(vector<SQindex*>::iterator j = d4ind.begin();j != d4ind.end();++j){
            if(**i == **j){
              IIvector h2perm(h2f.second.get_perms());
              IIvector d4perm(d4f.second.get_perms());
              typedef pair<bool, size_t> found;
              found h2, d4;
              h2.first = false;
              d4.first = false;
              h2.second = 0;
              d4.second = 0;
              for(size_t k = 0;k < h2perm.size();++k){
                if(*h2ind[h2perm[k][exth2_]] == **i) { h2.first = true; h2.second = k; break; } 
                if(h2.first) break;
	      } //End k
              for(size_t k = 0;k < d4perm.size();++k){
                if(*d4ind[d4perm[k][0]] == **j) { d4.first = true; d4.second = k; break; } 
                if(d4.first) break;
	      } //End k

              if(h2.first && d4.first) {
                vector<SQtensor*> ten_ptr(thisTerm->get_tensors_ptr()); 
                for(vector<SQtensor*>::iterator t = ten_ptr.begin();t != ten_ptr.end();++t){
                  if     ((*t)->get_name() == name_h2_) 
                    (*t)->rotateIndices(h2.second);
                  else if((*t)->get_name() == name_d4_) 
                    (*t)->rotateIndices(d4.second);
	        } // End t
                cout << "! Indices of ERI and D4 are rotated to match with each other. " << endl;
                //cout << "H2: " << h2.second << " D4: " << d4.second << endl;
                //cout << "*TEST* " << *thisTerm << endl; 
                SetV_D = true;
	      } // End if
	    } // End if
            if(SetV_D) break;
	  } // End j
          if(SetV_D) break;
	} // End i
      } // End if
      if(!SetT_V && !SetL_V && !SetV_D && h2f.first){
        vector<SQindex*> h2ind(h2f.second.get_indices());
        for(size_t num_i = 0;num_i < h2f.second.get_indices().size();++num_i){
          SQindex* i = h2f.second.get_indices()[num_i];
          if(i->get_char() == (char_state)2){
            IIvector h2perm(h2f.second.get_perms());
            typedef pair<bool, size_t> found;
            found h2;
            h2.first = false;
            h2.second = 0;
            for(size_t k = 0;k < h2perm.size();++k){
              if(*h2ind[h2perm[k][exth2_]] == *i) { h2.first = true; h2.second = k; break; } 
	    } // End k

            if(h2.first){
              vector<SQtensor*> ten_ptr(thisTerm->get_tensors_ptr()); 
              for(vector<SQtensor*>::iterator t = ten_ptr.begin();t != ten_ptr.end();++t){
                if((*t)->get_name() == name_h2_) (*t)->rotateIndices(h2.second);
	      } // End t
              cout << "! The indices of ERI are rotated to became virtual. " << endl;
              SetV_v = true;              
	    } // End if
	  } // End if
          if(SetV_v) break;
	} // End num_i
      } // End if

      // Set all the external indices 
      vector<SQtensor> tensors = thisTerm->get_tensors();
      for(vector<SQtensor>::iterator t = tensors.begin();t != tensors.end();++t){
        // In case of ERI
        if(t->get_name() == name_h2_) {
          vector<SQindex*> inds = t->get_indices();
          inds[exth2_]->switch_isExt(true);
	} // End if
        // In case of Bareamp
        else if(t->get_name() == name_amp_) {
          vector<SQindex*> inds = t->get_indices();
          inds[extamp_]->switch_isExt(true);
	} // End if
        else if(t->get_name() == name_d4_) {
          vector<SQindex*> inds = t->get_indices();
          inds[0]->switch_isExt(true);        
          inds[1]->switch_isExt(true);        
	}
      } // End t
      if(isBareLHS_){
        vector<SQindex*> inds = LTensor_.get_indices();
        inds[extamp_]->switch_isExt(true);
      } // End if

      CPfile << "\n\n  {\n  // No." << numTerm << endl;
      CPfile << "  double flops = 0; // Flop count" << endl;;
      CPfile << "  //* " << LTensor_ << " <--  " << *thisTerm << endl;

      bool rank_guard = false; 
      bool dec_guard  = false;
      cout << "! * Begin scaling analysis .... *" << endl << endl;

      // Evaluate the reduced tensor if this scheme contains one
      {
        pair<bool, string> Yflag;
        Yflag.first = false;
        const string Yname("Y");
        for(size_t num_t = 0;num_t < thisTerm->get_tensors().size();++num_t){
          string ten_name(thisTerm->get_tensors()[num_t].get_name());
          if(ten_name.at(0) == Yname.c_str()[0]){ 
            Yflag.first  = true; 
            Yflag.second = ten_name;
            break;
	  } // End if 
	} // End num_t
        if(!Yflag.first && thisTerm->get_Consts().size()){
          vector<string> Consts(thisTerm->get_Consts());
          for(vector<string>::iterator s = Consts.begin();s != Consts.end();++s){
            if(!s->size()) continue;
            if(s->at(0) == Yname.c_str()[0]){
              Yflag.first  = true; 
              Yflag.second = *s;
              break;
	    } // End if
	  } // End s
	} // End if
        if(!Yflag.first){
          for(size_t num_t = 0;num_t < thisTerm->get_tensors().size();++num_t){
            string ten_name(thisTerm->get_tensors()[num_t].get_name());
            if(ten_name.at(0) == Yname.c_str()[0]){ 
              Yflag.first  = true; 
              Yflag.second = ten_name;
              break;
  	    } // End if 
  	  } // End num_t  
  	} // End if
        if(!Yflag.first && thisTerm->get_Consts().size()){
          vector<string> Consts(thisTerm->get_Consts());
          for(vector<string>::iterator s = Consts.begin();s != Consts.end();++s){
            if(!s->size()) continue;
            if(s->at(0) == Yname.c_str()[0]){
              Yflag.first  = true; 
              Yflag.second = *s;
              break;
	    } // End if
	  } // End s
	} // End if

        // Now we got the name of Y if there is.
        if(Yflag.first){
          Yflag.second.erase(Yflag.second.begin());
          int Ynum(atoi(Yflag.second.c_str()));
          SQterm Ypattern(Yeffs[(size_t)Ynum]);

          // If Y tensor is not composed of ERI ....
          if(Ypattern.get_tensors()[0].get_name() != name_h2_){
            CPfile << "  if(myrank == 0){" << endl;
            rank_guard = true;
            dec_guard  = true;
	  } // End if

          vector<SQindex*> Yinds(Ypattern.get_summedBody());
          vector<SQindex*> new_inds;
          vector<SQindex*> left_inds;          
          //cout << Yname + Yflag.second << endl;
          //Ypattern.print_summedBody();
          vector<int> counters;
          for(size_t num_i = 0;num_i < Ypattern.get_tensors()[0].get_indices().size();++num_i){
            SQindex* i = Ypattern.get_tensors()[0].get_indices()[num_i];
            int count = 0;
            for(size_t num_j = 0;num_j < Ypattern.get_tensors()[0].get_indices().size();++num_j){
              SQindex* j = Ypattern.get_tensors()[0].get_indices()[num_j];
              if(i == j) ++count; 
	    } // End num_j
            counters.push_back(count);
            if(count == 1 || !i->get_isSummed()) new_inds.push_back(i);
            else           left_inds.push_back(i);
	  } // End num_i
          Symmetry symm;
          // Now all indices in Yten is inevitably linked to summedIndices in Ypattern
          SQtensor Yten(Yname + Yflag.second, new_inds, symm);
          //cout << Yten << endl;

          // Declare Y first
          CPfile << "  // The effective tensor is detected .... " << endl;
          cout << "Declare " + Yten.get_name() + " as a " + (Yten.get_indices().size() ? "tensor " : "scalar ") << endl;
          if(Yten.get_indices().size()){
            int Ccount = 0;
            int Ocount = 0;
            int Vcount = 0;
            string Yname("");
            vector<string> Var;
            for(size_t num_i = 0;num_i < Yten.get_indices().size();++num_i){
              if     (Yten.get_indices()[num_i]->get_char() == (char_state)0) ++Ccount;
              else if(Yten.get_indices()[num_i]->get_char() == (char_state)1) ++Ocount;
              else if(Yten.get_indices()[num_i]->get_char() == (char_state)2) ++Vcount;
	    } // End num_i
            for(int c = 0;c < Ccount;++c) { Yname += "c"; Var.push_back("nclosed"); }
            for(int o = 0;o < Ocount;++o) { Yname += "a"; Var.push_back("nocc");    }
            for(int v = 0;v < Vcount;++v) { Yname += "v"; Var.push_back("nvir");    }
            string DecY("  orz::DTensor Y(");
            for(size_t num_s = 0;num_s < Var.size();++num_s){
              string s = Var[num_s];
              DecY += s;
	      if(num_s != Var.size()-1) DecY += ", ";
              else DecY += ");";
	    } // End num_s

            string DecSymm("  orz::DTensor " + Yten.get_name()  + " = orz::mr::sympack_X" + Yname + "(symblockinfo, 0, Y);");
	    //cout << "Y" + Yname << endl;
            CPfile << DecY    << endl;
            CPfile << DecSymm << endl;
	  } // End if
          else{
            CPfile << "  double " + Yten.get_name() + " = 0;" << endl;;
	  } // End else

          size_t LoopCount = 0;
          string indent("");
          vector<string> Indents; Indents.reserve(100);
          for(int i = 0;i < 10;++i){
            Indents.push_back(indent);
            indent += "  ";
          } // End i 

          // Read ERI is if it has.
          size_t ERICount   = 0;          
          if(Ypattern.get_tensors()[0].get_name() == name_h2_)
            Ypattern.get_summedBody()[exth2_]->switch_isExt(true);
          for(size_t num_i = 0;num_i < Ypattern.get_summedBody().size();++num_i){
            SQindex* i = Ypattern.get_summedBody()[num_i];
            if(i->get_isExt()){
              string Label;
              if     (i->get_char() == (char_state)0) Label = "{core}";  
              else if(i->get_char() == (char_state)1) Label = "{occ}";  
              else if(i->get_char() == (char_state)2) Label = "{vir}";              
              cout << "for " + i->get_index() + " in " + Label + ":\n";
              CLoop("", *i, CPfile);
              ++LoopCount;
              cout << Indents[LoopCount] + "Read " + Ypattern.get_tensors()[0].get_name() + " from GA for " + i->get_index() + "\n";
              ReadERI(Indents[LoopCount], Ypattern.get_tensors()[0], CPfile);
              ERICount = LoopCount;
              break; 
	    } // End if
	  } // End num_i
          
          // Prepare input for the small functions
          vector<string> ExtInd;
          vector<string> NameTen;
	  vector<string> Consts;
          contDecl DecZero;
          for(size_t num_c = 0;num_c < Ypattern.get_Consts().size();++num_c)
            if(Ypattern.get_Consts()[num_c] != "") Consts.push_back(Ypattern.get_Consts()[num_c]);
          NameTen.push_back(Ypattern.get_tensors()[0].get_name());
          NameTen.push_back(Yten.get_name());
          for(vector<SQindex*>::iterator i = Yinds.begin();i != Yinds.end();++i)
            if((*i)->get_isExt()) ExtInd.push_back((*i)->get_index());

          DecZero.insert(contDecl::value_type("DecInd",    ExtInd));  // Names of external indices
          DecZero.insert(contDecl::value_type("DecConst",  Consts));  // All the constants
          DecZero.insert(contDecl::value_type("DecTensor", NameTen)); // Names of the tensors 

          ostringstream stm;
          stm << numTerm;
          string title_0if("g_if_" + title_ + "_y" + stm.str());
          string title_0  ("g_"    + title_ + "_y" + stm.str());

	  // Print the calling section of the Fortran body
          makeCPP_header(Yten, Ypattern, title_0if, CHfile, Indents[LoopCount], DecZero, false);
          // Interface inbetween C++ and F90 codes
          makeF90_interface(Yten, Ypattern, title_0if, F90file, Indents[LoopCount], DecZero, false);
          // Body of the tensorial contraction
          makeF90_contract(Yten, Ypattern, title_0, F90file, Indents[LoopCount], DecZero, false);
          // Print the calling section of the Fortran body
          makeCPP_body(Yten, Ypattern, title_0if, CPfile, Indents[LoopCount], DecZero, false);

          while(LoopCount > 0){
            if(ERICount && LoopCount == ERICount)
	      CPfile << Indents[LoopCount+1] + "} // End my_imo" << endl;
            LoopEnd(Indents[--LoopCount], CPfile);
            //--LoopCount;
  	  } // End while

	  cout << endl;
	} // End if
      } // End scoping

      typedef pair<bool, string> find_y;
      find_y y_scala;  y_scala.first  = false;
      find_y y_tensor; y_tensor.first = false;
      vector<string> coeffs(thisTerm->get_Consts());
      for(vector<string>::iterator c = coeffs.begin();c != coeffs.end();++c)
        if(c->at(0) == 'Y') { y_scala.first = true; y_scala.second = *c; }
      for(size_t num_t = 0;num_t < thisTerm->get_tensors().size();++num_t){
        if(thisTerm->get_tensors()[num_t].get_name().at(0) == 'Y') { 
	  y_tensor.first = true;
          y_tensor.second = thisTerm->get_tensors()[num_t].get_name();
	} // End num_t
      } // End num_t
      if(y_scala.first && !rank_guard){
        CPfile << "  double sum_" + y_scala.second << " = 0.0;" << endl;
        CPfile << "  boost::mpi::reduce(orz::world()," + y_scala.second + ", sum_" + y_scala.second + ", std::plus<double>(), 0);" << endl;
        CPfile << "  " << y_scala.second  << " = sum_" + y_scala.second << ";" << endl;        
      } // End if
      if(y_tensor.first && !rank_guard){
        CPfile << "  " << y_tensor.second << " = orz::tensor::mpi::reduce_plus(" + y_tensor.second + ", 0);" << endl;
      } // End if

      // If there is no ERIs .... 
      if(!h2f.first && !dec_guard){
        CPfile << "  if(myrank == 0){" << endl;
        rank_guard = true;
        dec_guard  = true;
      } // End if

      size_t LoopCount  = 0;
      size_t SigmaCount = 0;
      size_t D4Count    = 0;
      size_t ERICount   = 0;
      vector<SQindex*> Declared;
      string indent("");
      vector<string> Indents; Indents.reserve(100);
      for(int i = 0;i < 100;++i){
        Indents.push_back(indent);
        indent += "  ";
      } // End i 

      // ********** First Contraction ***********
      ostringstream stm;
      stm << numTerm;
      string title_1if("g_if_" + title_ + "_no0_x" + stm.str());
      string title_1  ("g_"    + title_ + "_no0_x" + stm.str());

      transform(title_1if.begin(), title_1if.end(), title_1if.begin(), (int(*)(int))tolower); 
      transform(title_1  .begin(), title_1  .end(), title_1  .begin(), (int(*)(int))tolower); 

      vector<SQindex*> O1(thisBoddies);
      for(vector<SQindex*>::iterator i = O1.begin();i != O1.end();)
        if(!((*i)->get_isExt())) i = O1.erase(i);
        else ++i;

      // Prepare input map for each declaration method .... 
      vector<string> ExtInd;
      vector<string> NameTen;
      vector<string> Consts;
      for(size_t num_c = 0;num_c < thisTerm->get_Consts().size();++num_c)
        if(thisTerm->get_Consts()[num_c] != "") Consts.push_back(thisTerm->get_Consts()[num_c]);

      contDecl DecFirst;
      for(size_t num_t = 0;num_t < thisTerm->get_tensors().size();++num_t){
        SQtensor t = thisTerm->get_tensors()[num_t];
        if(!is_RDM(t.get_name()) && t.get_name() != kDelta_name() && 
		     find(NameTen.begin(), NameTen.end(), t.get_name()) == NameTen.end()){
          NameTen.push_back(t.get_name());
	} // End if
      } // End num_t
      sort(NameTen.begin(), NameTen.end());
      NameTen.push_back(LTensor_.get_name());

      for(vector<SQindex*>::iterator i = thisBoddies.begin(); i != thisBoddies.end();++i)
        if((*i)->get_isExt()) ExtInd.push_back((*i)->get_index());
      sort(ExtInd.begin(), ExtInd.end());

      DecFirst.insert(contDecl::value_type("DecInd",    ExtInd));  // Names of external indices
      DecFirst.insert(contDecl::value_type("DecConst",  Consts));  // All the constants
      DecFirst.insert(contDecl::value_type("DecTensor", NameTen)); // Names of the tensors 

      // Print the calling section of the Fortran body
      makeCPP_header(LTensor_, *thisTerm, title_1if, CHfile, Indents[LoopCount], DecFirst, isBareLHS_);
      // Interface inbetween C++ and F90 codes
      makeF90_interface(LTensor_, *thisTerm, title_1if, F90file, Indents[LoopCount], DecFirst, isBareLHS_);
      // Body of the tensorial contraction
      if(use_gemm_ && thisTerm->get_tensors().size() == 2)
      binary_contract(LTensor_, *thisTerm, title_1, F90file, Indents[LoopCount], DecFirst, isBareLHS_);
      else
      makeF90_contract(LTensor_, *thisTerm, title_1, F90file, Indents[LoopCount], DecFirst, isBareLHS_);

//*OLD ?*       // If O1 has both irdm1 and irdm2, irdm1 must be prior to irdm2
//*OLD ?*       tensors = thisTerm->get_tensors();
//*OLD ?*       vector<SQtensor> theBody(thisTerm->get_tensors());
//*OLD ?*       for(vector<SQtensor>::iterator t = theBody.begin();t != theBody.end();++t){
//*OLD ?*         if(t->get_name()==name_d4_) {
//*OLD ?*           bool irdm1flag = false;
//*OLD ?*           vector<SQindex*>::iterator irdm1;
//*OLD ?*           bool irdm2flag = false;
//*OLD ?*           vector<SQindex*>::iterator irdm2;
//*OLD ?*           for(vector<SQindex*>::iterator i = O1.begin();i != O1.end();++i){
//*OLD ?*             if(**i == *(t->get_indices()[0])) { irdm1flag = true; irdm1 = i; }
//*OLD ?*             if(**i == *(t->get_indices()[1])) { irdm2flag = true; irdm2 = i; }
//*OLD ?* 	  } // End i
//*OLD ?*           if(irdm1flag && irdm2flag){
//*OLD ?*             if((irdm1-O1.begin()) > (irdm2-O1.begin())){
//*OLD ?*               SQindex* temp_i = *irdm2;
//*OLD ?*               O1.erase(irdm2);
//*OLD ?*               O1.push_back(temp_i);
//*OLD ?*             } // End if
//*OLD ?* 	  } // End if
//*OLD ?* 	} // End if
//*OLD ?*       } // End t

      // If O1 has both irdm1 and irdm2, irdm1 must be prior to irdm2
      tensors = thisTerm->get_tensors();
      vector<SQtensor> theBody(thisTerm->get_tensors());
      for(vector<SQtensor>::iterator t = theBody.begin();t != theBody.end();++t){
        if(t->get_name()==name_d4_) {
          bool irdm1flag = false;
          vector<SQindex*>::iterator irdm1;
          bool irdm2flag = false;
          vector<SQindex*>::iterator irdm2;
          for(vector<SQindex*>::iterator i = O1.begin();i != O1.end();++i){
            if(**i == *(t->get_indices()[0])) { irdm1flag = true; irdm1 = i; }
            if(**i == *(t->get_indices()[1])) { irdm2flag = true; irdm2 = i; }
	  } // End i
          if(irdm1flag && irdm2flag){
            if((irdm1-O1.begin()) > (irdm2-O1.begin())){
              SQindex* temp_i1 = *irdm1;
              SQindex* temp_i2 = *irdm2;
              O1.erase(irdm1);
              O1.erase(irdm2);
              //O1.push_back(temp_i);
              O1.insert(O1.begin(), temp_i2);
              O1.insert(O1.begin(), temp_i1);
	    } // End if
	  } // End if
	} // End if
      } // End t

      // ERI should be the top!
      for(vector<SQtensor>::iterator t = theBody.begin();t != theBody.end();++t){
        if(t->get_name() == name_h2_){
          SQtensor temp = *t;
          t = theBody.erase(t);
          theBody.insert(theBody.begin(), temp);
	}
      } // End t

      // Order of the ERI must be top!
      for(vector<SQtensor>::iterator t = theBody.begin();t != theBody.end();++t){
        if(t->get_name() == name_h2_){
          vector<SQindex*>::iterator i_eri = find(O1.begin(), O1.end(), t->get_indices()[exth2_]);
          if(i_eri != O1.end()) {
            SQindex* temp_i = *i_eri;
            O1.erase(i_eri);
            O1.insert(O1.begin(), temp_i);
	  } // End if
	} // End if
      } // End t

      vector<SQtensor> t_list; // List of tensors already read from GA
      for(vector<SQindex*>::iterator i = O1.begin();i != O1.end();++i){
        string Label;
        if((*i)->get_char() == (char_state)0) Label = "{core}";
        else if((*i)->get_char() == (char_state)1) Label = "{occ}";
        else if((*i)->get_char() == (char_state)2) Label = "{vir}";
        cout << Indents[LoopCount] << "for " << (*i)->get_index() << " in " + Label << ":" << endl;
        CLoop(Indents[LoopCount], **i, CPfile);
        ++LoopCount;
        Declared.push_back(*i);
        if(isBareLHS_)
          if((**i) == *(LTensor_.get_indices()[extamp_])){
            cout << Indents[LoopCount] << "Read " + LTensor_.get_name() + " from GA for " << (*i)->get_index() << endl;
            ReadRetval(Indents[LoopCount], LTensor_, CPfile);
            SigmaCount = LoopCount;
	  } // End if
        for(vector<SQtensor>::iterator t = theBody.begin();t != theBody.end();++t){
          if(t->get_name() == name_h2_ && (**i) == *(t->get_indices()[exth2_])){
            cout << Indents[LoopCount] << "Read " + t->get_name() + " from GA for " << (*i)->get_index() << endl;
            ReadERI(Indents[LoopCount], *t, CPfile);
            ERICount = LoopCount;
            t_list.push_back(*t);
	  } // End if
          else if(t->get_name() == name_amp_ && (**i) == *(t->get_indices()[extamp_])){
            cout << Indents[LoopCount] << "Read " + t->get_name() + " from GA for " << (*i)->get_index() << endl;
            ReadAmp(Indents[LoopCount], *t, CPfile);
            t_list.push_back(*t);            
	  } // End if
          else if(t->get_name() == name_d4_ && find(Declared.begin(), Declared.end(), (t->get_indices()[0]))!=Declared.end() 
		                            && find(Declared.begin(), Declared.end(), (t->get_indices()[1]))!=Declared.end()
		                            && find(t_list.begin(), t_list.end(), *t)==t_list.end()){
            cout << Indents[LoopCount] << "Read " + t->get_name() + " from GA for " + t->get_indices()[0]->get_index() + ", " + t->get_indices()[1]->get_index() << endl;
            if(!use_cumulant_)
            ReadD4(Indents[LoopCount], *t, CPfile);
            else
            ReadD4_Cumulant(Indents[LoopCount], *t, CPfile);
            D4Count = LoopCount;
            t_list.push_back(*t);
          } // End if

	} // End t
      } // End i

//*OLD ?*       // In case of 4-RDM
//*OLD ?*       for(vector<SQtensor>::iterator t = theBody.begin();t != theBody.end();++t){
//*OLD ?*         if(t->get_name() == name_d4_ && find(O1.begin(), O1.end(), t->get_indices()[0])!=O1.end() 
//*OLD ?* 	                             && find(O1.begin(), O1.end(), t->get_indices()[1])!=O1.end()){
//*OLD ?*           cout << Indents[LoopCount] << "Read " + t->get_name() + " from GA for " + t->get_indices()[0]->get_index() + ", " + t->get_indices()[1]->get_index() << endl;
//*OLD ?*           if(!use_cumulant_)
//*OLD ?*           ReadD4(Indents[LoopCount], *t, CPfile);
//*OLD ?*           else
//*OLD ?*           ReadD4_Cumulant(Indents[LoopCount], *t, CPfile);
//*OLD ?*           D4Count = LoopCount;
//*OLD ?*           t_list.push_back(*t);
//*OLD ?* 	} // End if
//*OLD ?*       } // End t

      int iCount = 0;
      string T1name = LTensor_.get_name() + "_(";
      if(!isBareLHS_){
        for(size_t num_i = 0;num_i < O1.size();++num_i){
          SQindex* i = O1[num_i];
          T1name += i->get_index();
          if(num_i != O1.size()-1) T1name += ",";
	} // End num_i
        T1name += ")(";
        for(size_t num_i = 0;num_i < LTensor_.get_indices().size();++num_i){
          SQindex* i = LTensor_.get_indices()[num_i];
          if(find(O1.begin(), O1.end(), i) == O1.end()){
            T1name += i->get_index();
            ++iCount;
            if(iCount != LTensor_.get_indices().size()-O1.size()) T1name += ",";
            else T1name += ")";
	  } // End if
	} // End num_i
      } // End if
      else{
        T1name += LTensor_.get_indices()[extamp_]->get_index() + ")(" + LTensor_.get_indices()[0]->get_index() + ", " + LTensor_.get_indices()[1]->get_index() + ", "
	  + LTensor_.get_indices()[2]->get_index() + ")";
      } // End else

      string sumname = "sum(";
      vector<SQindex*> O3;
      vector<SQindex*> XInd(LTensor_.get_indices());
      for(vector<SQtensor>::iterator t = theBody.begin();t != theBody.end();++t){
        for(size_t num_i = 0;num_i < t->get_indices().size();++num_i){
          SQindex* i = t->get_indices()[num_i];
          if(find(O1.begin(), O1.end(), i)==O1.end() && find(XInd.begin(), XInd.end(), i)==XInd.end()
	     && find(O3.begin(), O3.end(), i)==O3.end()) O3.push_back(i);
	} // End num_i
      } // End t

      for(size_t num_i = 0;num_i < O3.size();++num_i){
        SQindex* i = O3[num_i];
        sumname += i->get_index();
        if(num_i != O3.size()-1) sumname += ",";
      } // End num_i
      sumname += ")";

      string T2name(" ");
      for(vector<SQtensor>::iterator t = theBody.begin();t != theBody.end();++t){
        T2name += t->get_name() + "(";
        for(size_t num_i = 0;num_i < t->get_indices().size();++num_i){
          SQindex* i = t->get_indices()[num_i];
          if     (t->get_name() == name_h2_  && num_i != exth2_)         T2name += i->get_index();
          else if(t->get_name() == name_d4_  && num_i != 0 && num_i != 1) T2name += i->get_index();
          else if(t->get_name() == name_amp_ && num_i != extamp_)        T2name += i->get_index();
          else if(t->get_name() != name_h2_ && t->get_name() != name_amp_ && t->get_name() != name_d4_)
            T2name += i->get_index();

          if(num_i != t->get_indices().size()-1) T2name += ",";
	} // End num_I
        T2name += ")";
        if(!(*t == theBody.back())) T2name += " * ";
      } // End t
      
      string Constnames("");
      ostringstream stm2; 
      if(fabs(thisTerm->get_numConst()) != 1) {
        stm2 << thisTerm->get_numConst();
        Constnames += stm2.str() + " ";
      } // End if
      for(size_t num = 0;num < thisTerm->get_Consts().size();++num){
        Constnames += thisTerm->get_Consts()[num];
      } // End num
      Constnames += " ";
      cout << Indents[LoopCount] + T1name + " += " + Constnames + sumname + T2name << endl;

      // Print the calling section of the Fortran body
      makeCPP_body(LTensor_, *thisTerm, title_1if, CPfile, Indents[LoopCount], DecFirst, isBareLHS_);

      if(isBareLHS_){
        string LTname(LTensor_.get_name() + "_(" + LTensor_.get_indices()[extamp_]->get_index());
        LTname += ")(";
        for(size_t num = 0;num < LTensor_.get_indices().size()-1;++num){
          SQindex* i = LTensor_.get_indices()[num];
          LTname += i->get_index();
          if(num != LTensor_.get_indices().size()-2) LTname += ",";
	} // nd num 
        LTname += ")";

        cout << Indents[SigmaCount] << "Accumulate " + LTname + " for " + LTensor_.get_indices()[extamp_]->get_index() << endl;
        //AccAmp(Indents[LoopCount], LTensor_, CPfile);
      } // End if

      cout << ""                                                                         << endl;
      cout << "! ----------------------------------------------------------------------" << endl;
      cout << "! ----------------------------------------------------------------------" << endl;
      cout << ""                                                                         << endl;

      // Print D4Count post-processing
      while(LoopCount > 0){
        if(ERICount && LoopCount == ERICount){
          ERICount = 0;
          CPfile << Indents[LoopCount+1] + "}// End my_imo" << endl;
	} // End if        
        if(isBareLHS_ && LoopCount == SigmaCount){
          SigmaCount = 0;
          AccAmp(Indents[LoopCount], LTensor_, CPfile);
	} // End if
        if(D4Count && LoopCount == D4Count){
          D4Count = 0;
          CPfile << "  " + Indents[LoopCount] + "FC_FUNC(g_if_unset_d4,G_IF_UNSET_D4)();" << endl;
	} // End if        
        LoopEnd(Indents[--LoopCount], CPfile);
      } // End while

      if(rank_guard){
        CPfile << "  } // End if" << endl;
      } // End if
      if(!LTensor_.get_indices().size()){
	CPfile << "  double sum_" +  LTensor_.get_name() << " = 0.0;" << endl;
	CPfile << "  boost::mpi::reduce(orz::world()," + LTensor_.get_name() + ", sum_" + LTensor_.get_name() + ", std::plus<double>(), 0);" << endl;
	CPfile << "  " << LTensor_.get_name()  << " = sum_" + LTensor_.get_name() << ";"<< endl;        
      } // End if

      CPfile << "  }" << endl;
      CPfile << "  orz::world().barrier();" << endl;
      for(size_t num_i = 0;num_i < LTensor_.get_indices().size();++num_i) 
        LTensor_.put_indices(num_i, ref_indices[num_i]);

    } // End thisTerm

    if(!LTensor_.get_indices().size() && Bareflag){
      CPfile << endl;
      CPfile << "  return  " << LTensor_.get_name()  << ";" << endl;
      CPfile << "} " << endl;
    } // End if
    else{
      if(isBareLHS_ && LTensor_.get_indices()[0]->get_char() == (char_state)1 && LTensor_.get_indices()[1]->get_char() == (char_state)1
	           && LTensor_.get_indices()[2]->get_char() == (char_state)2 && LTensor_.get_indices()[3]->get_char() == (char_state)2){

        CPfile << endl;
        CPfile << "  if(myrank == 0){                                                                                                          " << endl;
        CPfile << "  for(int ssig = 0;ssig < nir;++ssig){                                                                                      " << endl;        
        CPfile << "  for(int isig = symblockinfo.psym()(ssig,I_V,I_BEGIN);isig <= symblockinfo.psym()(ssig,I_V,I_END);++isig){                 " << endl;    
        CPfile << "    " + LTensor_.get_name() + "b = retval.get_amp2(isig);                                                                   " << endl;      
        CPfile << "    FC_FUNC(g_if_sigma_oovv_scale,G_IF_SIGMA_OOVV_SCALE) // to make hamiltonian symmetric (scale by some factor of overlap) " << endl;       
        CPfile << "      (ssig, isig, " + LTensor_.get_name() + "b.cptr(), nir, nsym, psym);                                                   " << endl;              
        CPfile << "    retval.put_amp2(isig, " + LTensor_.get_name() + "b); // S2ija, [b] <<-- Sb                                              " << endl;                   
        CPfile << "  } // End isig                                                                                                             " << endl;            
        CPfile << "  } // End ssig                                                                                                             " << endl;     
        CPfile << "  } // End if                                                                                                               " << endl;     
	CPfile << "  orz::world().barrier();" << endl;
      } // End if
      if(isBareLHS_ && LTensor_.get_indices()[0]->get_char() == (char_state)0 && LTensor_.get_indices()[1]->get_char() == (char_state)0
	           && LTensor_.get_indices()[2]->get_char() == (char_state)2 && LTensor_.get_indices()[3]->get_char() == (char_state)2){

        CPfile << endl;
        CPfile << "  if(myrank == 0){                                                                                                          " << endl;
        CPfile << "  for(int ssig = 0;ssig < nir;++ssig){                                                                                      " << endl;        
        CPfile << "  for(int isig = symblockinfo.psym()(ssig,I_V,I_BEGIN);isig <= symblockinfo.psym()(ssig,I_V,I_END);++isig){                 " << endl;    
        CPfile << "    " + LTensor_.get_name() + "b = retval.get_amp2(isig);                                                                   " << endl;      
        CPfile << "    FC_FUNC(g_if_sigma_ccvv_scale,G_IF_SIGMA_CCVV_SCALE) // to make hamiltonian symmetric (scale by some factor of overlap) " << endl;       
        CPfile << "      (ssig, isig, " + LTensor_.get_name() + "b.cptr(), nir, nsym, psym);                                                   " << endl;              
        CPfile << "    retval.put_amp2(isig, " + LTensor_.get_name() + "b); // S2ija, [b] <<-- Sb                                              " << endl;                   
        CPfile << "  } // End isig                                                                                                             " << endl;            
        CPfile << "  } // End ssig                                                                                                             " << endl;     
        CPfile << "  } // End if                                                                                                               " << endl;     
	CPfile << "  orz::world().barrier();" << endl;
      } // End if
      if(isBareLHS_ && LTensor_.get_indices()[0]->get_char() == (char_state)0 && LTensor_.get_indices()[1]->get_char() == (char_state)0
	           && LTensor_.get_indices()[2]->get_char() == (char_state)1 && LTensor_.get_indices()[3]->get_char() == (char_state)1){

        CPfile << endl;
        CPfile << "  if(myrank == 0){                                                                                                          " << endl;
        CPfile << "  for(int ssig = 0;ssig < nir;++ssig){                                                                                      " << endl;        
        CPfile << "  for(int isig = symblockinfo.psym()(ssig,I_O,I_BEGIN);isig <= symblockinfo.psym()(ssig,I_O,I_END);++isig){                 " << endl;    
        CPfile << "    " + LTensor_.get_name() + "b = retval.get_amp2(isig);                                                                   " << endl;      
        CPfile << "    FC_FUNC(g_if_sigma_ccoo_scale,G_IF_SIGMA_CCOO_SCALE) // to make hamiltonian symmetric (scale by some factor of overlap) " << endl;       
        CPfile << "      (ssig, isig, " + LTensor_.get_name() + "b.cptr(), nir, nsym, psym);                                                   " << endl;              
        CPfile << "    retval.put_amp2(isig, " + LTensor_.get_name() + "b); // S2ija, [b] <<-- Sb                                              " << endl;                   
        CPfile << "  } // End isig                                                                                                             " << endl;            
        CPfile << "  } // End ssig                                                                                                             " << endl;     
        CPfile << "  } // End if                                                                                                               " << endl;     
	CPfile << "  orz::world().barrier();" << endl;
      } // End if

      CPfile << endl;
      CPfile << "  return retval; " << endl;
      CPfile << "} " << endl;
    } // End else

    CHfile << CHend;

    return retVal.str();    
  }                                               

}} //End Femto
