
#include "Eclipse.hpp"

// MM""""""""`M          dP oo                            
// MM  mmmmmmmM          88                               
// M`      MMMM .d8888b. 88 dP 88d888b. .d8888b. .d8888b. 
// MM  MMMMMMMM 88'  `"" 88 88 88'  `88 Y8ooooo. 88ooood8 
// MM  MMMMMMMM 88.  ... 88 88 88.  .88       88 88.  ... 
// MM        .M `88888P' dP dP 88Y888P' `88888P' `88888P' 
// MMMMMMMMMMMM                88                         

using namespace std;

void
Eclipse::generate_diag()
{

  // Construction of the diagonal elements of the Hamiltonian
  int num_ham = 0;
  string name("< Calculation of the diagonal elements starts > \n\n");
  cout << name;
  status += name;

  print_eqns << " ++ Diagonal elements of fic-MRCI ::" << endl << endl;

  for(auto Li = EL_indices.begin();Li != EL_indices.end();++Li){

    vector<Femto::Core::SQindex*> Lid;
    Lid.push_back(Li->at(2));
    Lid.push_back(Li->at(3));
    Lid.push_back(Li->at(0));
    Lid.push_back(Li->at(1));

    Femto::Core::sfGen EL(*Li);
    Femto::Core::sfGen ER(Lid);

    // T2 amplitude
    size_t num_l = (size_t)(Li-EL_indices.begin());
    size_t num_r = (size_t)(Li-EL_indices.begin());
    Femto::Core::SQtensor Tamp("T2", ET_indices[num_r], Femto::u4_symm());
  
    //if(names[num_l] != "oovv") { cout << "TEST" << endl; continue; } //*TEST*

    ostringstream stm1;
    stm1 << num_ham++;    
    string title("* " + stm1.str() + " L:" + names[num_l] + "/R:" + names[num_r] + "  ");
    cout << "* " << num_ham << " L:" + names[num_l] + "/R:" + names[num_r] << endl << endl;
    status += title;
  
    // Unit coefficient
    vector<string> coeff1;
  
    vector<Femto::Core::SQterm> result; result.reserve(Femto::Nterms());
    //vector<Femto::Core::SQterm> CoreTerms; CoreTerms.reserve(Femto::Nterms());
  
    // One-body part ....
    for(auto p = Ps.begin();p != Ps.end();++p){
      for(auto q = Qs.begin();q != Qs.end();++q){
  
        std::vector<Femto::Core::SQindex*> EH_indices;
        EH_indices.push_back(&(*p));
        EH_indices.push_back(&(*q));
  
        Femto::Core::SQtensor h("h", EH_indices, Femto::h1_symm());
        Femto::Core::sfGen EH(EH_indices);
  
        // < Psi | EL h EH T2 ER | Psi >
        vector<Femto::Core::SQtensor> ten1;
        ten1.push_back(EL);
        ten1.push_back(h);
        ten1.push_back(EH);
        ten1.push_back(ER);
  
        Femto::Core::SQterm term1(1.0, coeff1, ten1);
        vector<Femto::Core::SQterm> batch; batch.push_back(term1);
        Femto::Core::normalOrder(&batch); 
  
        for(vector<Femto::Core::SQterm>::iterator b_it = batch.begin();b_it != batch.end();++b_it){ 
          b_it->contractkDeltas(); // Burst Kronecker's deltas 
          b_it->decomposeRDMVirt();    // Decompose RDM, or sfGen 
        } 
        vector<Femto::Core::SQterm> temp;
        Femto::Core::screenTerms(batch, &temp); // Screen terms with negligible factor
  
        for(vector<Femto::Core::SQterm>::iterator b_it = temp.begin();b_it != temp.end();++b_it) 
          b_it->transform2RDM(); // Transform sfGen to RDM (only if isInCanonical is true)
  
        result.insert(result.end(), temp.begin(), temp.end()); //*SLOW* 

        // Use of the commutator
        if(names[num_l] == "ooov" || names[num_l] == "cooo"){
          // < Psi | EL ER EH | Psi >
          vector<Femto::Core::SQtensor> ten2;
          ten2.push_back(EL);
          ten2.push_back(ER);
          ten2.push_back(h);
          ten2.push_back(EH);
  
          Femto::Core::SQterm term2(-1.0, coeff1, ten2);
          vector<Femto::Core::SQterm> batch2; batch2.push_back(term2);
          Femto::Core::normalOrder(&batch2); 
  
          for(vector<Femto::Core::SQterm>::iterator b_it = batch2.begin();b_it != batch2.end();++b_it){ 
            b_it->contractkDeltas(); // Burst Kronecker's deltas 
            b_it->decomposeRDMVirt();    // Decompose RDM, or sfGen 
          } 
          vector<Femto::Core::SQterm> temp2;
          Femto::Core::screenTerms(batch2, &temp2); // Screen terms with negligible factor
  
          for(vector<Femto::Core::SQterm>::iterator b_it = temp2.begin();b_it != temp2.end();++b_it) 
            b_it->transform2RDM(); // Transform sfGen to RDM (only if isInCanonical is true)
  
          result.insert(result.end(), temp2.begin(), temp2.end()); //*SLOW* 
        } // End if

      } // End q
    } // End p
  
    // Two-body part ....
    int c = 0;
    for(auto p = Ps.begin();p != Ps.end();++p){
      for(auto q = Qs.begin();q != Qs.end();++q){
	for(auto r = Rs.begin();r != Rs.end();++r){
	  for(auto s = Ss.begin();s != Ss.end();++s){
  
            std::vector<Femto::Core::SQindex*> EH_indices;
            EH_indices.push_back(&(*p));
            EH_indices.push_back(&(*q));
            EH_indices.push_back(&(*r));
            EH_indices.push_back(&(*s));
  
            Femto::Core::SQtensor V("V2", EH_indices, Femto::h2_symm());
            Femto::Core::sfGen EH(EH_indices);
  
            // < Psi | EL V EH ER | Psi >
            vector<Femto::Core::SQtensor> ten1;
            ten1.push_back(EL);
            ten1.push_back(V);
            ten1.push_back(EH);
            ten1.push_back(ER);
                  
            Femto::Core::SQterm term1(0.5, coeff1, ten1);
  
            vector<Femto::Core::SQterm> batch; batch.push_back(term1);
            Femto::Core::normalOrder(&batch); 
  
            for(vector<Femto::Core::SQterm>::iterator b_it = batch.begin();b_it != batch.end();++b_it){ 
              b_it->contractkDeltas(); // Burst Kronecker's deltas 
              b_it->decomposeRDMVirt();    // Decompose RDM, or sfGen 
            }
  
            vector<Femto::Core::SQterm> temp;
            Femto::Core::screenTerms(batch, &temp); // Screen terms with negligible factor
  
            for(vector<Femto::Core::SQterm>::iterator b_it = temp.begin();b_it != temp.end();++b_it) 
  	        { b_it->masquerade(); cout << *b_it << endl; b_it->transform2RDM(); }// Transform sfGen to RDM (only if isInCanonical is true)
            result.insert(result.end(), temp.begin(), temp.end());

            // Use of the commutator
            if(names[num_l] == "ooov" || names[num_l] == "cooo"){
              // < Psi | EL ER V EH | Psi >
              vector<Femto::Core::SQtensor> ten2;
              ten2.push_back(EL);
              ten2.push_back(ER);
              ten2.push_back(V);
              ten2.push_back(EH);
                  
              Femto::Core::SQterm term2(-0.5, coeff1, ten2);
  	        //cout << "this " << term2 << endl;
              vector<Femto::Core::SQterm> batch2; batch2.push_back(term2);
              Femto::Core::normalOrder(&batch2); 
  
              for(vector<Femto::Core::SQterm>::iterator b_it = batch2.begin();b_it != batch2.end();++b_it){ 
                b_it->contractkDeltas(); // Burst Kronecker's deltas 
                b_it->decomposeRDMVirt();    // Decompose RDM, or sfGen 
              }
  
              vector<Femto::Core::SQterm> temp2;
              Femto::Core::screenTerms(batch2, &temp2); // Screen terms with negligible factor
  
       	      for(vector<Femto::Core::SQterm>::iterator b_it = temp2.begin();b_it != temp2.end();++b_it)
              // Transform sfGen to RDM (only if isInCanonical is true) 
  	        { b_it->masquerade(); cout << *b_it << endl; b_it->transform2RDM(); }
  
              result.insert(result.end(), temp2.begin(), temp2.end());
  	    } // End if

          } // End s
        } // End r
      } // End q
    } // End p

    if(names[num_l] == "ooov"  || names[num_l] == "cooo"){
    // Zeroth-body part
    // Ecas < Psi | EL T2 ER | Psi >
      vector<Femto::Core::SQtensor> ten;
      ten.push_back(EL);
      ten.push_back(ER);

      vector<string> Ecas;
      Ecas.push_back("Ecas");

      Femto::Core::SQterm term(1.0, Ecas, ten);
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

    //vector<Femto::Core::SQterm> combined_result;
    //Femto::combineTerms(result, &combined_result);
    
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
    Femto::Core::SQtensor Hdiag("Hdiag", *Li, Femto::u4_symm());
    formulas << "<" + names[num_l] + "\\" + names[num_r] + ">" << endl; 
    formulas << Hdiag.convert2LaTeX() << "&=& ";  
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
      
//     count2 = 0;
//     cout << endl;
//     vector<Femto::Core::SQterm> res;
//     Femto::find_JK_pairs(combined_result2, res);
//     cout << "< RESULT3 >" << endl;
//     for(vector<Femto::Core::SQterm>::iterator t = res.begin();t != res.end();++t){
//       cout << boost::format("%5d : ") % count2 << *t << endl;
//       ++count2;
//     }
  
    int count = 0;
    cout << endl;
    cout << "< RESULT >" << endl;
    for(vector<Femto::Core::SQterm>::iterator t = combined_result2.begin();t != combined_result2.end();++t){
      //cout << count << " : " << *t << endl;
      cout << boost::format("%5d : ") % count << *t << endl;
      ++count;
    }

    // Count the number of the 4-RDM
    const std::string name_d4(Femto::RDM_name()+"4");
    for(std::vector<Femto::Core::SQterm>::iterator t = combined_result2.begin();t != combined_result2.end();++t){
      for(size_t num_t = 0;num_t < t->get_tensors().size();++num_t){
        if(t->get_tensors()[num_t].get_name() == name_d4) ++nums_rdm4.first;
      } // End num_t
    } // End t

    //*OLD* ostringstream stm2;
    //*OLD* total_num_terms += count;
    //*OLD* stm2 << count;
    //*OLD* status += stm2.str() + " terms are generated .... \n\n";

#ifdef _CODE_GENE  
    // ********* Construct Sigma *************
    string thisname("diag_" + names[num_l]);
    Femto::Reaktor::SQreaktor gen(Hdiag, combined_result2, thisname, true, "V2", "T2");
    gen.set_use_cumulant(true);
    //gen.set_use_gemm(false); //*DEBUG* 
    if(names[num_l] == "ooov" && names[num_r] == "cooo") gen.set_T_D(true);
    if(names[num_l].find("c") != string::npos || names[num_r].find("c") != string::npos) gen.set_guard_core(true);
    if(names[num_l].find("o") != string::npos || names[num_r].find("o") != string::npos) gen.set_guard_act(true);
    gen.init_header("mrci_header");
    //gen.generate(Femto::Orz, Femto::SimpleLoops);    
    gen.generate(Femto::Reaktor::Orz, Femto::Reaktor::Factorize);    

    // Print all the terms
    if(combined_result2.size()) {
      std::vector<Femto::Core::SQterm> terms(gen.get_inTerms());
      print_eqns << " :::: " << names[num_l] << " element :::: " << endl;
      print_eqns << "   " << Hdiag << " += " << endl;
      for(std::vector<Femto::Core::SQterm>::const_iterator t = terms.begin();t != terms.end();++t) 
	print_eqns << "     " << *t << endl;
      print_eqns << endl;
    } // End scope

    ostringstream stm2;
    total_num_terms += gen.num_inTerms();
    stm2 << gen.num_inTerms();
    status += stm2.str() + " terms are generated .... \n\n";    

    // ********* Construct Sigma *************
#endif
  } // End Li  

}
