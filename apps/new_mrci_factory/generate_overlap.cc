
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
Eclipse::generate_overlap()
{

  // Generation of the diagonal elements for test
  string name("< Calculation of the overlap vector > \n\n");
  cout << name;
  status += name;
  int count = 0;
  for(size_t num_l = 0;num_l < EL_indices.size();++num_l){
  for(size_t num_r = 0;num_r < ER_indices.size();++num_r){

    if(!((num_l == 8 && num_r == 7) || (num_l == 7 && num_r == 8)) && num_l != num_r) continue;

    Femto::Core::sfGen EL(EL_indices[num_l]);
    Femto::Core::sfGen ER(ER_indices[num_r]); 
  
    ostringstream stm1;
    stm1 << count++;
    string title("* " + stm1.str() + " <" + names[num_l] + "/" + names[num_r] + ">  ");
    cout << title << endl;
    status += title;
  
    // T2 amplitude
    Femto::Core::SQtensor Tamp("T2", ET_indices[num_r], Femto::u4_symm());
    vector<Femto::Core::SQterm> result; result.reserve(Femto::Nterms());
    //vector<Femto::Core::SQterm> CoreTerms; CoreTerms.reserve(Femto::Nterms());
    
    vector<Femto::Core::SQtensor> ten;
    ten.push_back(EL);
    ten.push_back(Tamp);
    ten.push_back(ER);

    vector<string> coeff;
    Femto::Core::SQterm term(1.0, coeff, ten);
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
    Femto::Core::combineTerms(result2, &combined_result2);
  
    // ****** Conversion of tensor into the code starts ********
    int count2 = 0;
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

    ostringstream stm2;
    total_num_terms += count2;
    stm2 << count2;
    status += stm2.str() + " terms are generated .... \n\n";

#ifdef _CODE_GENE      
    // ********* Construct Sigma *************
    string thisname("overlap_" + names[num_l]);
    if(num_l != num_r) thisname += "_" + names[num_r];
    Femto::Core::SQtensor O2("O2", EL_indices[num_l], Femto::u4_symm());
    Femto::Reaktor::SQreaktor gen(O2, combined_result2, thisname, true, "V2", "T2");
    gen.generate(Femto::Reaktor::Orz, Femto::Reaktor::Factorize);
    // ********* Construct Sigma *************
#endif
  } // End num_l
  } // End num_r

}
