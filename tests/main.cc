
#include <iostream>
#include <vector>
#include <boost/format.hpp>
#include <femto.hpp>
#include <SQterm.hpp>
#include <SQindex.hpp>
#include <SQtensor.hpp>
#include <SQbinary.hpp>
#include <SQreaktor.hpp>

using namespace std;

int main(){

  std::vector<femto::Core::SQindex*> EL_indices;
  std::vector<femto::Core::SQindex*> ER_indices;

  femto::Core::SQindex A(femto::Core::SQindex("a", femto::virt));
  femto::Core::SQindex B(femto::Core::SQindex("b", femto::virt, true));
  femto::Core::SQindex C(femto::Core::SQindex("c", femto::virt));
  femto::Core::SQindex D(femto::Core::SQindex("d", femto::virt, true));

  femto::Core::SQindex I(femto::Core::SQindex("i", femto::act));
  femto::Core::SQindex J(femto::Core::SQindex("j", femto::act, true));
  femto::Core::SQindex K(femto::Core::SQindex("k", femto::act));
  femto::Core::SQindex L(femto::Core::SQindex("l", femto::act, true));

  vector<femto::Core::SQindex> Ps;
  vector<femto::Core::SQindex> Qs;
  vector<femto::Core::SQindex> Rs;
  vector<femto::Core::SQindex> Ss;
  for(size_t space = 0;space < 3;++space){
    Ps.push_back(femto::Core::SQindex("p", (femto::char_state)space, true));
    Qs.push_back(femto::Core::SQindex("q", (femto::char_state)space, true));
    Rs.push_back(femto::Core::SQindex("r", (femto::char_state)space, true));
    Ss.push_back(femto::Core::SQindex("s", (femto::char_state)space, true));
  }

  EL_indices.push_back(&I);
  EL_indices.push_back(&K);  
  EL_indices.push_back(&A);
  EL_indices.push_back(&C);  

  ER_indices.push_back(&B);
  ER_indices.push_back(&D);
  ER_indices.push_back(&J);
  ER_indices.push_back(&L);

  femto::Core::sfGen EL(EL_indices);
  femto::Core::sfGen ER(ER_indices);

  // T2 amplitude
  femto::Core::SQtensor Tamp("T2", ER_indices, femto::u4_symm());

  // Unit coefficient
  vector<string> coeff1;
  //coeff1.push_back("");

  vector<femto::Core::SQterm> result; result.reserve(femto::Nterms());
  //vector<femto::SQterm> CoreTerms; CoreTerms.reserve(femto::Nterms());

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
      ten1.push_back(EL);
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
          ten1.push_back(EL);
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

  vector<femto::Core::SQterm> combined_result;
  femto::Core::combineTerms(result, &combined_result);

  int count = 0;
  cout << endl;
  cout << "< RESULT >" << endl;
  for(vector<femto::Core::SQterm>::iterator t = combined_result.begin();t != combined_result.end();++t){
    //cout << count << " : " << *t << endl;
    cout << boost::format("%5d : ") % count << *t << endl;
    ++count;
  }

//*TEST*   cout << endl;
//*TEST*   cout << "Perm test .... " << endl;
//*TEST*   vector<int> ti, tj;
//*TEST*   for(size_t i = 0;i < 6;++i){
//*TEST*     ti.push_back(i);
//*TEST*     tj.push_back(i);
//*TEST*   }
//*TEST*   tj.erase(tj.begin()+3);
//*TEST*   tj.push_back(3);
//*TEST*   tj.erase(tj.begin()+2);
//*TEST*   tj.push_back(2);
//*TEST* 
//*TEST*   for(size_t i = 0;i < 6;++i){
//*TEST*     cout << ti[i] << ", " << tj[i] << endl;;
//*TEST*   }
//*TEST* 
//*TEST*   cout << endl << femto::get_num_perms(ti,tj) << endl;

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

  vector<femto::Core::SQterm> combined_result2;
  femto::Core::combineTerms(result2, &combined_result2);

  int count2 = 0;
  cout << endl;
  cout << "< RESULT >" << endl;
  for(vector<femto::Core::SQterm>::iterator t = combined_result2.begin();t != combined_result2.end();++t){
    //cout << count << " : " << *t << endl;
    cout << boost::format("%5d : ") % count2 << *t << endl;
    ++count2;
  }

  count2 = 0;
  cout << endl;
  cout << "< RESULT2 >" << endl;
  for(vector<femto::Core::SQterm>::iterator t = combined_result2.begin();t != combined_result2.end();++t){
    cout << boost::format("%5d : ") % count2 << *t << endl;
    ++count2;
  }

  femto::Core::makeFock(combined_result2);

  femto::Core::SQtensor Sig("S2", EL_indices, femto::u4_symm());
  string thisname("sigma_oovv_oovv_test");
  femto::Reaktor::SQreaktor gen(Sig, combined_result2, thisname, true, "V2", "T2");
  gen.generate(femto::Reaktor::Orz, femto::Reaktor::Factorize);

}
