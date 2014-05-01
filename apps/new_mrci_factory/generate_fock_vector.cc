
#include "Eclipse.hpp"

// MM""""""""`M          dP oo                            
// MM  mmmmmmmM          88                               
// M`      MMMM .d8888b. 88 dP 88d888b. .d8888b. .d8888b. 
// MM  MMMMMMMM 88'  `"" 88 88 88'  `88 Y8ooooo. 88ooood8 
// MM  MMMMMMMM 88.  ... 88 88 88.  .88       88 88.  ... 
// MM        .M `88888P' dP dP 88Y888P' `88888P' `88888P' 
// MMMMMMMMMMMM                88                         

//#define _DIAG_ONLY
//#define _NEW_ENO

using namespace std;

void
Eclipse::generate_fock_vector()
{

  // Construction of Sigma_{aa'}^{ee'} += <Psi0|EL H TR ER |Psi0>
  int num_block = 0;
  int num_non_zero = 0;
  string name("< Construction of F_{aa'}^{ee'} += <Psi0|EL F TR ER|Psi0> > \n\n");
  cout << name;
  status += name;

  print_eqns << " ++ Sigma like vector of Fock Hamiltonian ::" << endl << endl;

  for(auto Li = EL_indices.begin();Li != EL_indices.end();++Li){
    for(auto Ri = ER_indices.begin();Ri != ER_indices.end();++Ri){
      Femto::Core::sfGen EL(*Li);
      Femto::Core::sfGen ER(*Ri);

      // T2 amplitude
      size_t num_l = (size_t)(Li-EL_indices.begin());
      size_t num_r = (size_t)(Ri-ER_indices.begin());
      Femto::Symmetry symm(Femto::u4_symm());
      {
	// If the the amplitude is of T_{XX'}^{YY'} type excitation, allow permutational symmetry
        if(ET_indices[num_r].at(0)->get_char() == ET_indices[num_r].at(1)->get_char() && 
           ET_indices[num_r].at(2)->get_char() == ET_indices[num_r].at(3)->get_char()) symm = Femto::t2_symm();
      }       
      Femto::Core::SQtensor Tamp("T2", ET_indices[num_r], symm);

#ifdef _DIAG_ONLY    
      if(num_l != num_r) { cout << "DIAG ONLY" << endl; continue; } 
#endif

      ostringstream stm_l;
      ostringstream stm_r;
      ostringstream stm1;
      stm_l << num_l;
      stm_r << num_r;
      stm1  << num_block++;
      string title("* " + stm1.str() + " <" + names[num_l] + "/" + names[num_r] + ">  ");
      cout << "* " << num_block << " <" + names[num_l] + "/" + names[num_r] + ">" << endl << endl;
      status += title;

      // Unit coefficient
      vector<string> coeff1;
    
      vector<Femto::Core::SQterm> result; result.reserve(Femto::Nterms());
      //vector<Femto::Core::SQterm> CoreTerms; CoreTerms.reserve(Femto::Nterms());
    
      // One-body part ....
      for(vector<Femto::Core::SQindex>::iterator p = Ps.begin();p != Ps.end();++p){
        for(vector<Femto::Core::SQindex>::iterator q = Qs.begin();q != Qs.end();++q){
    
          std::vector<Femto::Core::SQindex*> EH_indices;
          EH_indices.push_back(&(*p));
          EH_indices.push_back(&(*q));
    
          Femto::Core::SQtensor h("P1", EH_indices, Femto::h1_symm());
          Femto::Core::sfGen EH(EH_indices);
    
          // < Psi | EL f EH T2 ER | Psi >
          vector<Femto::Core::SQtensor> ten1;
          ten1.push_back(EL);
          ten1.push_back(h);
          ten1.push_back(EH);
          ten1.push_back(Tamp);
          ten1.push_back(ER);

          double numcoeff = 1.0;
#ifdef _NEW_ENO
          if(names[num_l] == names[num_r] && (names[num_r] == "ooov" || names[num_l] == "cooo"))
            numcoeff = 0.5;           
#endif
          Femto::Core::SQterm term1(numcoeff, coeff1, ten1);
          vector<Femto::Core::SQterm> batch; batch.push_back(term1);

#ifdef _NEW_ENO
          if(names[num_l] == names[num_r] && (names[num_r] == "ooov" || names[num_l] == "cooo")){
	    Femto::Core::normalOrderComm(&batch); 
	  } else{
	    Femto::Core::normalOrder(&batch); 
	  }
#else
	  Femto::Core::normalOrder(&batch); 
#endif

          for(vector<Femto::Core::SQterm>::iterator b_it = batch.begin();b_it != batch.end();++b_it){ 
            b_it->contractkDeltas(); // Burst Kronecker's deltas 
            b_it->decomposeRDMVirt();    // Decompose RDM, or sfGen 
          } 
          vector<Femto::Core::SQterm> temp;
          Femto::Core::screenTerms(batch, &temp); // Screen terms with negligible factor
    
          for(vector<Femto::Core::SQterm>::iterator b_it = temp.begin();b_it != temp.end();++b_it) 
            b_it->transform2RDM(); // Transform sfGen to RDM (only if isInCanonical is true)
    
          result.insert(result.end(), temp.begin(), temp.end()); //*SLOW* 
    
#ifdef _NEW_ENO
	  if(names[num_l] == names[num_r] && (names[num_r] == "ooov" || names[num_l] == "cooo")){
            // EL_tilde = ER^{\dagger}
            vector<Femto::Core::SQindex*> Lid;
            Lid.push_back(ER.get_indices().at(2));
            Lid.push_back(ER.get_indices().at(3));
            Lid.push_back(ER.get_indices().at(0));
            Lid.push_back(ER.get_indices().at(1));
            // ER_tilde = EL^{\dagger}
            vector<Femto::Core::SQindex*> Rid;
            Rid.push_back(EL.get_indices().at(2));
            Rid.push_back(EL.get_indices().at(3));
            Rid.push_back(EL.get_indices().at(0));
            Rid.push_back(EL.get_indices().at(1));

	    Femto::Core::sfGen EL_tilde(Lid);
	    Femto::Core::sfGen ER_tilde(Rid);

            // +1/2 <Psi| EL_tilde f Eh T ER_tilde |Psi>
            vector<Femto::Core::SQtensor> ten4;
            ten4.push_back(EL_tilde);
            ten4.push_back(h);
            ten4.push_back(EH);
            ten4.push_back(Tamp);
            ten4.push_back(ER_tilde);
    
            Femto::Core::SQterm term4(+0.5, coeff1, ten4);
            vector<Femto::Core::SQterm> batch4; batch4.push_back(term4);

            Femto::Core::normalOrderComm(&batch4); 

            for(vector<Femto::Core::SQterm>::iterator b_it = batch4.begin();b_it != batch4.end();++b_it){ 
              b_it->contractkDeltas(); // Burst Kronecker's deltas 
              b_it->decomposeRDMVirt();    // Decompose RDM, or sfGen 
            } 
            vector<Femto::Core::SQterm> temp4;
            Femto::Core::screenTerms(batch4, &temp4); // Screen terms with negligible factor
    
            for(vector<Femto::Core::SQterm>::iterator b_it = temp4.begin();b_it != temp4.end();++b_it) 
              b_it->transform2RDM(); // Transform sfGen to RDM (only if isInCanonical is true)
    
            result.insert(result.end(), temp4.begin(), temp4.end()); //*SLOW* 
	  } // End if
#endif

        } // End q
      } // End p

#ifdef _NEW_ENO
      if(names[num_l] == names[num_r] && (names[num_l] == "ooov" || names[num_l] == "cooo")){
      // Zeroth-body part
      // Eres < Psi | EL T2 ER | Psi >
        vector<Femto::Core::SQtensor> ten;
        ten.push_back(EL);
        ten.push_back(Tamp);
        ten.push_back(ER);
  
        vector<string> Eres;
        Eres.push_back("Eres");
  
        Femto::Core::SQterm term(1.0, Eres, ten);
        vector<Femto::Core::SQterm> batch; batch.push_back(term);
        Femto::Core::normalOrder(&batch);
  
        for(vector<Femto::Core::SQterm>::iterator b_it = batch.begin();b_it != batch.end();++b_it){ 
          b_it->contractkDeltas(); // Burst Kronecker's deltas 
          b_it->decomposeRDMVirt();    // Decompose RDM, or sfGen 
          b_it->masquerade();
        }
  
        vector<Femto::Core::SQterm> temp;
        Femto::Core::screenTerms(batch, &temp); // Screen terms with negligible factor
  
        for(vector<Femto::Core::SQterm>::iterator b_it = temp.begin();b_it != temp.end();++b_it) 
          { cout << *b_it << endl; b_it->transform2RDM(); }// Transform sfGen to RDM (only if isInCanonical is true)
        result.insert(result.end(), temp.begin(), temp.end());
      } // End if
#endif
        
      cout << endl << "Decompose RDMs ..... " << endl;
      int cnt = 0;
      vector<Femto::Core::SQterm> result2;
      for(vector<Femto::Core::SQterm>::iterator t = result.begin();t != result.end();++t){
        vector<Femto::Core::SQterm> batch;
        //cout << boost::format("%5d : ") % cnt++ << *t << endl;
        Femto::Core::decomposeRDMCore(*t, &batch);
        vector<Femto::Core::SQterm>::iterator new_t = batch.begin();
        for(;new_t != batch.end();++new_t){
          new_t->contractkDeltas();
          new_t->transform2RDM();
        } // End new_t
        result2.insert(result2.end(), batch.begin(), batch.end());
      } // End t
    
      vector<Femto::Core::SQterm> combined_result2;
      Femto::Core::combineTerms2(result2, &combined_result2);
    
      cout << endl;
      cnt = 0;
          
      // TeX conversion
      Femto::Core::SQtensor Sig("S2", *Li, Femto::u4_symm());
      formulas << "<" + names[num_l] + "\\" + names[num_r] + ">" << endl; 
      formulas << Sig.convert2LaTeX() << "&=& ";  
      for(vector<Femto::Core::SQterm>::iterator t = combined_result2.begin();t != combined_result2.end();++t, ++cnt){
        formulas << t->convert2LaTeX();
        if((cnt+1)%6==0) formulas << "\\\\" << endl << "& &";
      } // End t
      formulas << endl << endl;
        
      int count2 = 0;
      cout << endl;
      cout << "< RESULT >" << endl;
      for(vector<Femto::Core::SQterm>::iterator t = combined_result2.begin();t != combined_result2.end();++t){
        cout << boost::format("%5d : ") % count2 << *t << endl;
        ++count2;
      }
      cout << endl;
    
      // ****** Conversion of tensor into the code starts ********
      count2 = 0;
      cout << endl;
      cout << "< RESULT2 >" << endl;
      // Convert Dirac 2 Mulliken, do not do twice!
      for(vector<Femto::Core::SQterm>::iterator t = combined_result2.begin();t != combined_result2.end();++t){
        for(size_t i = 0;i < t->get_tensors().size();++i){
          if(t->get_tensors()[i].get_name() == "V2") t->get_tensors_ptr()[i]->convertD2M();
          if(Femto::is_RDM(t->get_tensors()[i].get_name())) t->get_tensors_ptr()[i]->convertD2M();
        } // End i
        cout << boost::format("%5d : ") % count2 << *t << endl;
        ++count2;
      }
        
//       count2 = 0;
//       cout << endl;
//       vector<Femto::Core::SQterm> res;
//       Femto::find_JK_pairs(combined_result2, res);
//       cout << "< RESULT3 >" << endl;
//       for(vector<Femto::Core::SQterm>::iterator t = res.begin();t != res.end();++t){
//         cout << boost::format("%5d : ") % count2 << *t << endl;
//         ++count2;
//       }
    
      int count = 0;
      cout << endl;
      cout << "< RESULT >" << endl;
      for(vector<Femto::Core::SQterm>::iterator t = combined_result2.begin();t != combined_result2.end();++t){
        //cout << count << " : " << *t << endl;
        cout << boost::format("%5d : ") % count << *t << endl;
        ++count;
      }

      ostringstream stm2;
      //*OLD* total_num_terms += count;
      //*OLD* stm2 << count;
      //*OLD* status += stm2.str() + " terms are generated .... \n\n";

#ifdef _CODE_GENE    
      // ********* Construct Sigma *************
      string thisname("fock_v_" + names[num_l] + "_" + names[num_r]);
      Femto::Reaktor::SQreaktor gen(Sig, combined_result2, thisname, true, "V2", "T2");
      gen.set_use_cumulant(true);
      //gen.set_priority(Femto::memory);
      if(names[num_l] != names[num_r]) //*This is OK, not the test*
	//if(names[num_l] == "ooov" || names[num_r] == "ooov") //*TEST* 
        gen.set_V_D(false);
      //gen.set_isOF(true); //*DEBUG* 
      //gen.set_use_gemm(false); //*DEBUG* 
      //gen.replace_Fock();
      gen.set_D4C(true);
      if(names[num_l] == "ooov" && names[num_r] == "cooo") gen.set_T_D(true);
      if(names[num_l].find("c") != string::npos || names[num_r].find("c") != string::npos) gen.set_guard_core(true);
      if(names[num_l].find("o") != string::npos || names[num_r].find("o") != string::npos) gen.set_guard_act(true);
      gen.init_header("mrci_header");
      gen.generate(Femto::Reaktor::Orz, Femto::Reaktor::Factorize);

      // Print all the terms
      if(combined_result2.size()) {
	std::vector<Femto::Core::SQterm> terms(gen.get_inTerms());
	print_eqns << " :::: " << names[num_l] << "/" << names[num_r] << " element :::: " << endl;
	print_eqns << "   " << Sig << " += " << endl;
	for(std::vector<Femto::Core::SQterm>::const_iterator t = terms.begin();t != terms.end();++t) 
	  print_eqns << "     " << *t << endl;
	print_eqns << endl;
      } // End scope

      total_num_terms += gen.num_inTerms();
      stm2 << gen.num_inTerms();
      status += stm2.str() + " terms are generated .... \n\n";
      // ********* Construct Sigma *************

#ifdef _MAKE_GENE
      // Write down name of the file is EL and ER interacts with each other
      if(gen.num_inTerms()){
        // Should be aded to Makefile.am
        //*// make_fock_vector << "c_" << thisname << ".cpp \\" << endl;
        //*// make_fock_vector << f_source() << "f_" << thisname << ".F90 \\" << endl;
        make_fock_vector << "Femto/elems/" << "c_" << thisname << ".cpp \\" << endl;
        make_fock_vector << "Femto/tensors/" << "f_" << thisname << ".F90 \\" << endl;

        // Should be added to C++ source
        if(num_l == num_r && (names[num_l] == "ooov" || names[num_l] == "cooo"))
	  gaxpy_fock_vector << "retval_.gaxpy(+1.0, " + thisname + "(ctinp, symblockinfo, hintmo, CFock, alloc_type, rdmPack, rdm4, T2, num_sigma), +1.0); // No. " << num_non_zero++ << endl;
	else if(names[num_l] == "covv")
	  gaxpy_fock_vector << "ret_sym.gaxpy(+1.0, " + thisname + "(ctinp, symblockinfo, hintmo, CFock, alloc_type,                T2, num_sigma), +1.0); // No. " << num_non_zero++ << endl;
	else
	  gaxpy_fock_vector << "retval_.gaxpy(+1.0, " + thisname + "(ctinp, symblockinfo, hintmo, CFock, alloc_type,                T2, num_sigma), +1.0); // No. " << num_non_zero++ << endl;

        // Should be added to mrci.h
        if(num_l == num_r && (names[num_l] == "ooov" || names[num_l] == "cooo")){
	  header_fock_vector << "orz::ct::BareAmpPack " + thisname + "(const orz::ct::Input &ctinp,        " << endl;                             
	  header_fock_vector << "				       const orz::ct::SymBlockInfo &symblockinfo, " << endl;                      
	  header_fock_vector << "				       const orz::ct::HintMO &hintmo,             " << endl;                      
	  header_fock_vector << "				       const orz::DTensor &CFock,                 " << endl;                      
	  header_fock_vector << "				       const int alloc_type,                      " << endl;                      
	  header_fock_vector << "				       const orz::ct::RdmPack &rdmPack,           " << endl;                      
	  header_fock_vector << "				       const orz::DTensor &rdm4,                  " << endl;                      
	  header_fock_vector << "				       const orz::ct::BareAmpPack &T2,            " << endl;                      
	  header_fock_vector << "				       const int num_sigma);                      " << endl << endl;;                      
	} // End if
	else{
	  header_fock_vector << "orz::ct::BareAmpPack " + thisname + "(const orz::ct::Input &ctinp,        " << endl; 
	  header_fock_vector << "				       const orz::ct::SymBlockInfo &symblockinfo, " << endl;
	  header_fock_vector << "				       const orz::ct::HintMO &hintmo,             " << endl;
	  header_fock_vector << "				       const orz::DTensor &CFock,                 " << endl;     
	  header_fock_vector << "				       const int alloc_type,                      " << endl;
	  header_fock_vector << "				       const orz::ct::BareAmpPack &T2,            " << endl;   
	  header_fock_vector << "				       const int num_sigma);                      " << endl << endl;         
	} // End else

      } // End if
#endif

#endif
    } // End Li
  } // End Ri

}


