
#include <iostream>
#include <vector>
#include <fstream>
#include <boost/format.hpp>
#include <femto.hpp>
#include <SQterm.hpp>
#include <SQindex.hpp>
#include <SQtensor.hpp>
#include <SQreaktor.hpp>

using namespace std;

// Test code to derive the single-reference CID equations based on the contravariant formulations
// <tilde{Phi}^{wy}_{ac}| = [2<Phi^{wy}_{ac}| + <Phi^{yw}_{ac}|] / 6
// sigma_{wy}^{ac} = <tilde{Phi}^{wy}_{ac}| H |Phi_{xz}^{bd}> T_{xz}^{bd}

int main(){

  std::vector<femto::Core::SQindex*> EL_indices1;
  std::vector<femto::Core::SQindex*> EL_indices2;
  std::vector<femto::Core::SQindex*> ET_indices;
  std::vector<femto::Core::SQindex*> ER_indices;

  femto::Core::SQindex A = femto::Core::SQindex("a", femto::virt);
  femto::Core::SQindex B = femto::Core::SQindex("b", femto::virt, true);
  femto::Core::SQindex C = femto::Core::SQindex("c", femto::virt);
  femto::Core::SQindex D = femto::Core::SQindex("d", femto::virt, true);

  femto::Core::SQindex I = femto::Core::SQindex("w", femto::core);
  femto::Core::SQindex J = femto::Core::SQindex("x", femto::core, true);
  femto::Core::SQindex K = femto::Core::SQindex("y", femto::core);
  femto::Core::SQindex L = femto::Core::SQindex("z", femto::core, true);

  vector<femto::Core::SQindex> Ps;
  vector<femto::Core::SQindex> Qs;
  vector<femto::Core::SQindex> Rs;
  vector<femto::Core::SQindex> Ss;
  for(size_t space = 0;space < 3;++space){
    if((femto::char_state)space == femto::act) continue;
    Ps.push_back(femto::Core::SQindex("p", (femto::char_state)space, true));
    Qs.push_back(femto::Core::SQindex("q", (femto::char_state)space, true));
    Rs.push_back(femto::Core::SQindex("r", (femto::char_state)space, true));
    Ss.push_back(femto::Core::SQindex("s", (femto::char_state)space, true));
  }

  EL_indices1.push_back(&I);
  EL_indices1.push_back(&K);  
  EL_indices1.push_back(&A);
  EL_indices1.push_back(&C);  

  EL_indices2.push_back(&K);
  EL_indices2.push_back(&I);  
  EL_indices2.push_back(&A);
  EL_indices2.push_back(&C);  

  ER_indices.push_back(&B);
  ER_indices.push_back(&D);
  ER_indices.push_back(&J);
  ER_indices.push_back(&L);

  ET_indices.push_back(&J);
  ET_indices.push_back(&L);
  ET_indices.push_back(&B);
  ET_indices.push_back(&D);

  femto::Core::sfGen EL1(EL_indices1);
  femto::Core::sfGen EL2(EL_indices2);
  femto::Core::sfGen ER(ER_indices);

  // Unit coefficient
  vector<string> coeff1;
  //  coeff1.push_back("");

  // T2 amplitude
  femto::Core::SQtensor Tamp("T2", ET_indices, femto::t2_symm());

  vector<femto::Core::SQterm> result; result.reserve(femto::Nterms());
  //vector<femto::Core::SQterm> CoreTerms; CoreTerms.reserve(femto::Core::Nterms());

  { // contra1

  // One-body part ....
  for(vector<femto::Core::SQindex>::iterator p = Ps.begin();p != Ps.end();++p){
    for(vector<femto::Core::SQindex>::iterator q = Qs.begin();q != Qs.end();++q){

      std::vector<femto::Core::SQindex*> EH_indices;
      EH_indices.push_back(&(*p));
      EH_indices.push_back(&(*q));

      femto::Core::SQtensor h("h", EH_indices, femto::h1_symm());
      femto::Core::sfGen EH(EH_indices);

      // < Psi | EL h EH T2 ER | Psi >
      vector<femto::Core::SQtensor> ten1;
      ten1.push_back(EL1);
      ten1.push_back(h);
      ten1.push_back(EH);
      ten1.push_back(Tamp);
      ten1.push_back(ER);

      femto::Core::SQterm term1(1.0, coeff1, ten1);
      vector<femto::Core::SQterm> batch; batch.push_back(term1);
      femto::Core::normalOrder(&batch); 

      for(vector<femto::Core::SQterm>::iterator b_it = batch.begin();b_it != batch.end();++b_it){ 
        b_it->contractkDeltas(); // Burst Kronecker's deltas 
        b_it->decomposeRDMVirt();    // Decompose RDM, or sfGen 
      } 
      vector<femto::Core::SQterm> temp;
      femto::Core::screenTerms(batch, &temp); // Screen terms with negligible factor

      for(vector<femto::Core::SQterm>::iterator b_it = temp.begin();b_it != temp.end();++b_it) 
        b_it->transform2RDM(); // Transform sfGen to RDM (only if isInCanonical is true)

      result.insert(result.end(), temp.begin(), temp.end()); //*SLOW* 

    } // End q
  } // End p

  // Two-body part ....
  int c = 0;
  for(vector<femto::Core::SQindex>::iterator p = Ps.begin();p != Ps.end();++p){
    for(vector<femto::Core::SQindex>::iterator q = Qs.begin();q != Qs.end();++q){
      for(vector<femto::Core::SQindex>::iterator r = Rs.begin();r != Rs.end();++r){
        for(vector<femto::Core::SQindex>::iterator s = Ss.begin();s != Ss.end();++s){

          std::vector<femto::Core::SQindex*> EH_indices;
          EH_indices.push_back(&(*p));
          EH_indices.push_back(&(*q));
          EH_indices.push_back(&(*r));
          EH_indices.push_back(&(*s));

          femto::Core::SQtensor V("V2", EH_indices, femto::h2_symm());
          femto::Core::sfGen EH(EH_indices);

          // < Psi | EL V EH T2 ER | Psi >
          vector<femto::Core::SQtensor> ten1;
          ten1.push_back(EL1);
          ten1.push_back(V);
          ten1.push_back(EH);
          ten1.push_back(Tamp);
          ten1.push_back(ER);
                
          femto::Core::SQterm term1(0.5, coeff1, ten1);

          vector<femto::Core::SQterm> batch; batch.push_back(term1);
          femto::Core::normalOrder(&batch); 

          for(vector<femto::Core::SQterm>::iterator b_it = batch.begin();b_it != batch.end();++b_it){ 
            b_it->contractkDeltas(); // Burst Kronecker's deltas 
            b_it->decomposeRDMVirt();    // Decompose RDM, or sfGen 
          }

          vector<femto::Core::SQterm> temp;
          femto::Core::screenTerms(batch, &temp); // Screen terms with negligible factor

     	  for(vector<femto::Core::SQterm>::iterator b_it = temp.begin();b_it != temp.end();++b_it) 
	    { cout << *b_it << endl; b_it->transform2RDM(); }// Transform sfGen to RDM (only if isInCanonical is true)
          result.insert(result.end(), temp.begin(), temp.end());

        } // End s
      } // End r
    } // End q
  } // End p
  } // End contra1

  { // contra2

  // One-body part ....
  for(vector<femto::Core::SQindex>::iterator p = Ps.begin();p != Ps.end();++p){
    for(vector<femto::Core::SQindex>::iterator q = Qs.begin();q != Qs.end();++q){

      std::vector<femto::Core::SQindex*> EH_indices;
      EH_indices.push_back(&(*p));
      EH_indices.push_back(&(*q));

      femto::Core::SQtensor h("h", EH_indices, femto::h1_symm());
      femto::Core::sfGen EH(EH_indices);

      // < Psi | EL h EH T2 ER | Psi >
      vector<femto::Core::SQtensor> ten1;
      ten1.push_back(EL2);
      ten1.push_back(h);
      ten1.push_back(EH);
      ten1.push_back(Tamp);
      ten1.push_back(ER);

      femto::Core::SQterm term1(0.5, coeff1, ten1);
      vector<femto::Core::SQterm> batch; batch.push_back(term1);
      femto::Core::normalOrder(&batch); 

      for(vector<femto::Core::SQterm>::iterator b_it = batch.begin();b_it != batch.end();++b_it){ 
        b_it->contractkDeltas(); // Burst Kronecker's deltas 
        b_it->decomposeRDMVirt();    // Decompose RDM, or sfGen 
      } 
      vector<femto::Core::SQterm> temp;
      femto::Core::screenTerms(batch, &temp); // Screen terms with negligible factor

      for(vector<femto::Core::SQterm>::iterator b_it = temp.begin();b_it != temp.end();++b_it) 
        b_it->transform2RDM(); // Transform sfGen to RDM (only if isInCanonical is true)

      result.insert(result.end(), temp.begin(), temp.end()); //*SLOW* 

    } // End q
  } // End p

  // Two-body part ....
  int c = 0;
  for(vector<femto::Core::SQindex>::iterator p = Ps.begin();p != Ps.end();++p){
    for(vector<femto::Core::SQindex>::iterator q = Qs.begin();q != Qs.end();++q){
      for(vector<femto::Core::SQindex>::iterator r = Rs.begin();r != Rs.end();++r){
        for(vector<femto::Core::SQindex>::iterator s = Ss.begin();s != Ss.end();++s){

          std::vector<femto::Core::SQindex*> EH_indices;
          EH_indices.push_back(&(*p));
          EH_indices.push_back(&(*q));
          EH_indices.push_back(&(*r));
          EH_indices.push_back(&(*s));

          femto::Core::SQtensor V("V2", EH_indices, femto::h2_symm());
          femto::Core::sfGen EH(EH_indices);

          // < Psi | EL V EH T2 ER | Psi >
          vector<femto::Core::SQtensor> ten1;
          ten1.push_back(EL2);
          ten1.push_back(V);
          ten1.push_back(EH);
          ten1.push_back(Tamp);
          ten1.push_back(ER);
                
          femto::Core::SQterm term1(0.25, coeff1, ten1);

          vector<femto::Core::SQterm> batch; batch.push_back(term1);
          femto::Core::normalOrder(&batch); 

          for(vector<femto::Core::SQterm>::iterator b_it = batch.begin();b_it != batch.end();++b_it){ 
            b_it->contractkDeltas(); // Burst Kronecker's deltas 
            b_it->decomposeRDMVirt();    // Decompose RDM, or sfGen 
          }

          vector<femto::Core::SQterm> temp;
          femto::Core::screenTerms(batch, &temp); // Screen terms with negligible factor

     	  for(vector<femto::Core::SQterm>::iterator b_it = temp.begin();b_it != temp.end();++b_it) 
	    { cout << *b_it << endl; b_it->transform2RDM(); }// Transform sfGen to RDM (only if isInCanonical is true)
          result.insert(result.end(), temp.begin(), temp.end());

        } // End s
      } // End r
    } // End q
  } // End p
  } // End contra1

//*OVERLAP*   // Calculation of overlap 
//*OVERLAP*   {
//*OVERLAP*     vector<femto::Core::SQtensor> ten;
//*OVERLAP*     ten.push_back(EL);
//*OVERLAP*     ten.push_back(Tamp);
//*OVERLAP*     ten.push_back(ER);
//*OVERLAP*   
//*OVERLAP*     vector<string> Ecas;
//*OVERLAP*   
//*OVERLAP*     femto::Core::SQterm term(1.0, Ecas, ten);
//*OVERLAP*     vector<femto::Core::SQterm> batch; batch.push_back(term);
//*OVERLAP*     femto::Core::normalOrder(&batch);
//*OVERLAP*   
//*OVERLAP*     for(vector<femto::Core::SQterm>::iterator b_it = batch.begin();b_it != batch.end();++b_it){ 
//*OVERLAP*       b_it->contractkDeltas(); // Burst Kronecker's deltas 
//*OVERLAP*       b_it->decomposeRDMVirt();    // Decompose RDM, or sfGen 
//*OVERLAP*     }
//*OVERLAP*   
//*OVERLAP*     vector<femto::Core::SQterm> temp;
//*OVERLAP*     femto::Core::screenTerms(batch, &temp); // Screen terms with negligible factor
//*OVERLAP*   
//*OVERLAP*     for(vector<femto::Core::SQterm>::iterator b_it = temp.begin();b_it != temp.end();++b_it) 
//*OVERLAP*       { cout << *b_it << endl; b_it->transform2RDM(); }// Transform sfGen to RDM (only if isInCanonical is true)
//*OVERLAP*     result.insert(result.end(), temp.begin(), temp.end());
//*OVERLAP*   }

  vector<femto::Core::SQterm> combined_result;
  //  femto::Core::combineTerms(result, &combined_result);
  combined_result = result;

//*NO_CORE*   int count = 0;
//*NO_CORE*   cout << endl;
//*NO_CORE*   cout << "< RESULT >" << endl;
//*NO_CORE*   for(vector<femto::Core::SQterm>::iterator t = combined_result.begin();t != combined_result.end();++t){
//*NO_CORE*     //cout << count << " : " << *t << endl;
//*NO_CORE*     cout << boost::format("%5d : ") % count << *t << endl;
//*NO_CORE*     ++count;
//*NO_CORE*   }
//*NO_CORE* 
//*NO_CORE*   cout << endl;
//*NO_CORE*   int cnt = 0;
//*NO_CORE*   for(vector<femto::Core::SQterm>::iterator t = combined_result.begin();t != combined_result.end();++t, ++cnt){
//*NO_CORE*     t->convert2LaTeX();
//*NO_CORE*     if((cnt+1)%5==0) cout << endl;
//*NO_CORE*   } // End t
//*NO_CORE*   cout << endl;
//*NO_CORE* 
//*NO_CORE*   // ********* Construct Sigma *************
//*NO_CORE*   // Convert Dirac 2 Mulliken, do not do twice!
//*NO_CORE*   int count2 = 0;
//*NO_CORE*   cout << endl;
//*NO_CORE*   cout << "< RESULT in Mulliken notation >" << endl;
//*NO_CORE*   for(vector<femto::Core::SQterm>::iterator t = combined_result.begin();t != combined_result.end();++t){
//*NO_CORE*     for(size_t i = 0;i < t->get_tensors().size();++i){
//*NO_CORE*       if(t->get_tensors()[i].get_name() == "V2") t->get_tensors_ptr()[i]->convertD2M();
//*NO_CORE*       if(femto::Core::is_RDM(t->get_tensors()[i].get_name())) t->get_tensors_ptr()[i]->convertD2M();
//*NO_CORE*     } // End i
//*NO_CORE*     cout << boost::format("%5d : ") % count2 << *t << endl;
//*NO_CORE*     ++count2;
//*NO_CORE*   }
//*NO_CORE* 
//*NO_CORE*   femto::Core::SQtensor Sig("S2", EL_indices, femto::Core::u4_symm());  
//*NO_CORE*   //Sig.convert2LaTeX();
//*NO_CORE*   femto::Core::factorize(Sig, combined_result, "sigma_oovv_oovv", true, "V2", "T2");
//*NO_CORE*   // ********* Construct Sigma *************
  
  cout << endl << "Decompose RDMs ..... " << endl;
  int cnt = 0;
  vector<femto::Core::SQterm> result2;
  for(vector<femto::Core::SQterm>::iterator t = combined_result.begin();t != combined_result.end();++t){
    vector<femto::Core::SQterm> batch;
    //cout << boost::format("%5d : ") % cnt++ << *t << endl;
    femto::Core::decomposeRDMCore(*t, &batch);
    vector<femto::Core::SQterm>::iterator new_t = batch.begin();
    for(;new_t != batch.end();++new_t){
      new_t->contractkDeltas();
      new_t->transform2RDM();
    } // End new_t
    result2.insert(result2.end(), batch.begin(), batch.end());
  } // End t

  vector<femto::Core::SQterm> combined_result2; //cout << "Size " << result2.size() << endl;
  femto::Core::combineTerms(result2, &combined_result2);
  //femto::Core::combineTerms(result2); combined_result2 = result2; //*TEST*
  result2 = combined_result2;
  combined_result2.clear();
  femto::Core::combineTerms(result2, &combined_result2);

  cout << endl;
  cnt = 0;
  
  // TeX conversion
  ofstream formulas("formulas.tex");
  femto::Core::SQtensor Sig("S2", EL_indices1, femto::u4_symm());
  formulas << Sig.convert2LaTeX() << "&+=& ";  
  for(vector<femto::Core::SQterm>::iterator t = combined_result2.begin();t != combined_result2.end();++t, ++cnt){
    formulas << t->convert2LaTeX();
    if((cnt+1)%6==0) formulas << "\\\\" << endl << "& &";
  } // End t
  formulas << endl;

  int count2 = 0;
  cout << endl;
  cout << "< RESULT >" << endl;
  for(vector<femto::Core::SQterm>::iterator t = combined_result2.begin();t != combined_result2.end();++t){
    cout << boost::format("%5d : ") % count2 << *t << endl;
    ++count2;
  }
  cout << endl;

  // ****** Conversion of tensor into the code starts ********
  count2 = 0;
  cout << endl;
  cout << "< RESULT2 >" << endl;
  // Convert Dirac 2 Mulliken, do not do twice!
  for(vector<femto::Core::SQterm>::iterator t = combined_result2.begin();t != combined_result2.end();++t){
    for(size_t i = 0;i < t->get_tensors().size();++i){
      if(t->get_tensors()[i].get_name() == "V2") t->get_tensors_ptr()[i]->convertD2M();
      if(femto::is_RDM(t->get_tensors()[i].get_name())) t->get_tensors_ptr()[i]->convertD2M();
    } // End i
    cout << boost::format("%5d : ") % count2 << *t << endl;
    ++count2;
  }

  count2 = 0;
  cout << endl;
  cout << "< RESULT > SR-case" << endl;
  for(vector<femto::Core::SQterm>::iterator t = combined_result2.begin();t != combined_result2.end();++t){
    t->set_numConst(t->get_numConst() / 6);
    cout << boost::format("%5d : ") % count2 << *t << endl;
    ++count2;
  }
  cout << endl;

  // ********* Construct Sigma *************
  string thisname("sigma_cid");
  femto::Reaktor::SQreaktor gen(Sig, combined_result2, thisname, true, "V2", "T2");
  gen.set_use_cumulant(true);
  //gen.set_V_D(false);
  gen.generate(femto::Reaktor::Orz, femto::Reaktor::Factorize);
  vector<femto::Core::SQterm> terms(gen.get_inTerms());
  // ********* Construct Sigma *************

  cout << "< inTerms >" << endl;
  int count3(0);
  for(vector<femto::Core::SQterm>::iterator t = terms.begin();t != terms.end();++t){
    cout << boost::format("%5d : ") % count3 << *t << endl;
    ++count3;    
  }

}

