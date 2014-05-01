
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
Eclipse::generate_g_e()
{

  int num_ham = 0;
  string name("< Calculation of the Sigma_{0} += <Psi0| H T2 ER |Psi0> > \n\n");
  cout << name;
  status += name;

  print_eqns << " ++ Sigma vector of fic-MRCI (g-e) ::" << endl << endl;

  for(auto Ri = ER_indices.begin();Ri != ER_indices.end();++Ri){
    Femto::Core::sfGen ER(*Ri);

    // T2 amplitude
    size_t num_r = (size_t)(Ri-ER_indices.begin());
    Femto::Core::SQtensor Tamp("T2", ET_indices[num_r], Femto::t2_symm());

    ostringstream stm1;
    stm1 << num_ham++;  
    string title("* " + stm1.str() + " <g/" + names[num_r] + ">  ");
    cout << "* " << num_ham << " <g/" + names[num_r] + ">" << endl << endl;
    status += title;  

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
  
        // < Psi | h EH T2 ER | Psi >
        vector<Femto::Core::SQtensor> ten1;
        ten1.push_back(h);
        ten1.push_back(EH);
        ten1.push_back(Tamp);
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
  
            // < Psi | V EH T2 ER | Psi >
            vector<Femto::Core::SQtensor> ten1;
            ten1.push_back(V);
            ten1.push_back(EH);
            ten1.push_back(Tamp);
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

          } // End s
        } // End r
      } // End q
    } // End p
  
    vector<Femto::Core::SQterm> combined_result;
    Femto::Core::combineTerms(result, &combined_result);
    
    cout << endl << "Decompose RDMs ..... " << endl;
    int cnt = 0;
    vector<Femto::Core::SQterm> result2;
    for(vector<Femto::Core::SQterm>::iterator t = combined_result.begin();t != combined_result.end();++t){
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
    Femto::Core::combineTerms(result2, &combined_result2);
  
    cout << endl;
    cnt = 0;
        
    // TeX conversion
    vector<Femto::Core::SQindex*> null;
    Femto::Symmetry null_symm;
    Femto::Core::SQtensor Sig("S0", null, null_symm);
    formulas << "<g\\" + names[num_r] + ">" << endl; 
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

    //*OLD* ostringstream stm2;
    //*OLD* total_num_terms += count;
    //*OLD* stm2 << count;
    //*OLD* status += stm2.str() + " terms are generated .... \n\n";

#ifdef _CODE_GENE  
    // ********* Construct Sigma *************
    string thisname("sigma_g_" + names[num_r]);
    Femto::Reaktor::SQreaktor gen(Sig, combined_result2, thisname, false, "V2", "T2");
    if(names[num_r].find("c") != string::npos) gen.set_guard_core(true);
    if(names[num_r].find("o") != string::npos) gen.set_guard_act(true);
    gen.init_header("mrci_header");
    gen.generate(Femto::Reaktor::Orz, Femto::Reaktor::Factorize);

    // Print all the terms
    if(combined_result2.size()){
      std::vector<Femto::Core::SQterm> terms(gen.get_inTerms());
      print_eqns << " :::: " << "g/" << names[num_r] << " element :::: " << endl;
      print_eqns << "   " << Sig << " += " << endl;
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
  } // End Ri  

}
