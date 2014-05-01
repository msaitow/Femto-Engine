
#include <iostream>
#include <vector>
#include <fstream>
#include <boost/format.hpp>
#include <femto.hpp>
#include <SQterm.hpp>
#include <SQindex.hpp>
#include <SQtensor.hpp>
#include <Orz.hpp>

using namespace std;

// MM""""""""`M          dP oo                            
// MM  mmmmmmmM          88                               
// M`      MMMM .d8888b. 88 dP 88d888b. .d8888b. .d8888b. 
// MM  MMMMMMMM 88'  `"" 88 88 88'  `88 Y8ooooo. 88ooood8 
// MM  MMMMMMMM 88.  ... 88 88 88.  .88       88 88.  ... 
// MM        .M `88888P' dP dP 88Y888P' `88888P' `88888P' 
// MMMMMMMMMMMM                88                         

// Code to generate all the tensorial contraction codes for fic-MRCI

int main(){

  vector<vector<femto::SQindex*> > EL_indices;
  vector<vector<femto::SQindex*> > ER_indices;
  vector<vector<femto::SQindex*> > ET_indices;
  vector<femto::SQindex*> EL;
  vector<femto::SQindex*> ER;
  vector<femto::SQindex*> ET;

  femto::SQindex A = femto::SQindex("a", femto::virt);
  femto::SQindex B = femto::SQindex("b", femto::virt, true);
  femto::SQindex C = femto::SQindex("c", femto::virt);
  femto::SQindex D = femto::SQindex("d", femto::virt, true);

  femto::SQindex I = femto::SQindex("i", femto::act);
  femto::SQindex J = femto::SQindex("j", femto::act, true);
  femto::SQindex K = femto::SQindex("k", femto::act);
  femto::SQindex L = femto::SQindex("l", femto::act, true);
  femto::SQindex M = femto::SQindex("m", femto::act);
  femto::SQindex N = femto::SQindex("n", femto::act, true);

  femto::SQindex W = femto::SQindex("w", femto::core);
  femto::SQindex X = femto::SQindex("x", femto::core, true);
  femto::SQindex Y = femto::SQindex("y", femto::core);
  femto::SQindex Z = femto::SQindex("z", femto::core, true);

  vector<femto::SQindex> Ps;
  vector<femto::SQindex> Qs;
  vector<femto::SQindex> Rs;
  vector<femto::SQindex> Ss;
  for(size_t space = 0;space < 3;++space){
    Ps.push_back(femto::SQindex("p", (femto::char_state)space, true));
    Qs.push_back(femto::SQindex("q", (femto::char_state)space, true));
    Rs.push_back(femto::SQindex("r", (femto::char_state)space, true));
    Ss.push_back(femto::SQindex("s", (femto::char_state)space, true));
  }

  // External operator (left-side)
  EL.push_back(&I);
  EL.push_back(&K);  
  EL.push_back(&A);
  EL.push_back(&C);  
  EL_indices.push_back(EL);
  EL.clear();
  
  // Semiinternal operator (left-side)
  EL.push_back(&I);
  EL.push_back(&K);  
  EL.push_back(&M);
  EL.push_back(&A);  
  EL_indices.push_back(EL);
  EL.clear();

  // (c,c)->(v,v) operator (left-side)
  EL.push_back(&W);
  EL.push_back(&Y);  
  EL.push_back(&A);
  EL.push_back(&C);  
  EL_indices.push_back(EL);
  EL.clear();

  // External operator (right-side)
  ER.push_back(&B);
  ER.push_back(&D);
  ER.push_back(&J);
  ER.push_back(&L);
  ER_indices.push_back(ER);
  ER.clear();

  ET.push_back(&J);
  ET.push_back(&L);
  ET.push_back(&B);
  ET.push_back(&D);
  ET_indices.push_back(ET);
  ET.clear();

  // Semiinternal operator (right-side)
  ER.push_back(&N);
  ER.push_back(&B);
  ER.push_back(&L);
  ER.push_back(&J);
  ER_indices.push_back(ER);
  ER.clear();

  ET.push_back(&L);
  ET.push_back(&J);
  ET.push_back(&N);
  ET.push_back(&B);
  ET_indices.push_back(ET);
  ET.clear();

  // (c,c)->(v,v) operator (right-side)
  ER.push_back(&B);
  ER.push_back(&D);
  ER.push_back(&X);
  ER.push_back(&Z);
  ER_indices.push_back(ER);
  ER.clear();

  ET.push_back(&X);
  ET.push_back(&Z);
  ET.push_back(&B);
  ET.push_back(&D);
  ET_indices.push_back(ET);
  ET.clear();

  vector<string> names;
  names.push_back("oovv");
  names.push_back("ooov");
  names.push_back("ccvv");

  // File to write the sigma equations
  ofstream formulas("formulas.tex");

  // Generation of the diagonal elements for test
  cout << "* Calculation of the overlap vector .... " << endl << endl;
  for(size_t num = 0;num < EL_indices.size();++num){

    femto::sfGen EL(EL_indices[num]);
    femto::sfGen ER(ER_indices[num]);
  
    cout << "* " << num << " <" + names[num] + "/" << names[num] << ">" << endl << endl;
  
    // T2 amplitude
    femto::SQtensor Tamp("T2", ET_indices[num], femto::u4_symm());
    vector<femto::SQterm> result; result.reserve(femto::Nterms());
    //vector<femto::SQterm> CoreTerms; CoreTerms.reserve(femto::Nterms());
    
    vector<femto::SQtensor> ten;
    ten.push_back(EL);
    ten.push_back(Tamp);
    ten.push_back(ER);
  
    vector<string> coeff;
  
    femto::SQterm term(1.0, coeff, ten);
    vector<femto::SQterm> batch; batch.push_back(term);
    femto::normalOrder(&batch);
  
    for(vector<femto::SQterm>::iterator b_it = batch.begin();b_it != batch.end();++b_it){ 
      b_it->contractkDeltas(); // Burst Kronecker's deltas 
      b_it->decomposeRDMVirt();    // Decompose RDM, or sfGen 
      b_it->masquerade();
    }
  
    vector<femto::SQterm> temp;
    femto::screenTerms(batch, &temp); // Screen terms with negligible factor
  
    for(vector<femto::SQterm>::iterator b_it = temp.begin();b_it != temp.end();++b_it) 
      { cout << *b_it << endl; b_it->transform2RDM(); }// Transform sfGen to RDM (only if isInCanonical is true)
    result.insert(result.end(), temp.begin(), temp.end());

    cout << endl << "Decompose RDMs ..... " << endl;
    int cnt = 0;
    vector<femto::SQterm> result2;
    for(vector<femto::SQterm>::iterator t = result.begin();t != result.end();++t){
      vector<femto::SQterm> batch;
      //cout << boost::format("%5d : ") % cnt++ << *t << endl;
      femto::decomposeRDMCore(*t, &batch);
      vector<femto::SQterm>::iterator new_t = batch.begin();
      for(;new_t != batch.end();++new_t){
        new_t->contractkDeltas();
        new_t->transform2RDM();
      } // End new_t
      result2.insert(result2.end(), batch.begin(), batch.end());
    } // End t
  
    vector<femto::SQterm> combined_result2;
    femto::combineTerms(result2, &combined_result2);
  
    // ****** Conversion of tensor into the code starts ********
    int count2 = 0;
    cout << endl;
    cout << "< RESULT2 >" << endl;
    // Convert Dirac 2 Mulliken, do not do twice!
    for(vector<femto::SQterm>::iterator t = combined_result2.begin();t != combined_result2.end();++t){
      for(size_t i = 0;i < t->get_tensors().size();++i){
        if(t->get_tensors()[i].get_name() == "V2") t->get_tensors_ptr()[i]->convertD2M();
        if(femto::is_RDM(t->get_tensors()[i].get_name())) t->get_tensors_ptr()[i]->convertD2M();
      } // End i
      cout << boost::format("%5d : ") % count2 << *t << endl;
      ++count2;
    }
      
    // ********* Construct Sigma *************
    string thisname("overlap_" + names[num]);
    femto::SQtensor O2("O2", EL_indices[num], femto::u4_symm());
    femto::simpleloops(O2, combined_result2, thisname, true, "V2", "T2");
    // ********* Construct Sigma *************
  } // End num

  // Construction of the diagonal elements of the Hamiltonian
  int num_ham = 0;
  cout << "* Calculation of the diagonal elements starts ..... " << endl << endl;
  for(vector<vector<femto::SQindex*> >::iterator Li = EL_indices.begin();Li != EL_indices.end();++Li){
    vector<femto::SQindex*> Lid;
    Lid.push_back(Li->at(2));
    Lid.push_back(Li->at(3));
    Lid.push_back(Li->at(0));
    Lid.push_back(Li->at(1));

    femto::sfGen EL(*Li);
    femto::sfGen ER(Lid);

    // T2 amplitude
    size_t num_l = (size_t)(Li-EL_indices.begin());
    size_t num_r = (size_t)(Li-EL_indices.begin());
    femto::SQtensor Tamp("T2", ET_indices[num_r], femto::u4_symm());
  
    cout << "* " << num_ham++ << " L:" + names[num_l] + "/R:" + names[num_r] << endl << endl;
  
    // Unit coefficient
    vector<string> coeff1;
  
    vector<femto::SQterm> result; result.reserve(femto::Nterms());
    //vector<femto::SQterm> CoreTerms; CoreTerms.reserve(femto::Nterms());
  
    // One-body part ....
    for(vector<femto::SQindex>::iterator p = Ps.begin();p != Ps.end();++p){
      for(vector<femto::SQindex>::iterator q = Qs.begin();q != Qs.end();++q){
  
        std::vector<femto::SQindex*> EH_indices;
        EH_indices.push_back(&(*p));
        EH_indices.push_back(&(*q));
  
        femto::SQtensor h("h", EH_indices, femto::h1_symm());
        femto::sfGen EH(EH_indices);
  
        // < Psi | EL h EH T2 ER | Psi >
        vector<femto::SQtensor> ten1;
        ten1.push_back(EL);
        ten1.push_back(h);
        ten1.push_back(EH);
        ten1.push_back(ER);
  
        femto::SQterm term1(1.0, coeff1, ten1);
        vector<femto::SQterm> batch; batch.push_back(term1);
        femto::normalOrder(&batch); 
  
        for(vector<femto::SQterm>::iterator b_it = batch.begin();b_it != batch.end();++b_it){ 
          b_it->contractkDeltas(); // Burst Kronecker's deltas 
          b_it->decomposeRDMVirt();    // Decompose RDM, or sfGen 
        } 
        vector<femto::SQterm> temp;
        femto::screenTerms(batch, &temp); // Screen terms with negligible factor
  
        for(vector<femto::SQterm>::iterator b_it = temp.begin();b_it != temp.end();++b_it) 
          b_it->transform2RDM(); // Transform sfGen to RDM (only if isInCanonical is true)
  
        result.insert(result.end(), temp.begin(), temp.end()); //*SLOW* 

        if(names[num_l] == "ooov"){
          // < Psi | EL ER EH | Psi >
          vector<femto::SQtensor> ten2;
          ten2.push_back(EL);
          ten2.push_back(ER);
          ten2.push_back(h);
          ten2.push_back(EH);
  
          femto::SQterm term2(-1.0, coeff1, ten2);
          vector<femto::SQterm> batch2; batch2.push_back(term2);
          femto::normalOrder(&batch2); 
  
          for(vector<femto::SQterm>::iterator b_it = batch2.begin();b_it != batch2.end();++b_it){ 
            b_it->contractkDeltas(); // Burst Kronecker's deltas 
            b_it->decomposeRDMVirt();    // Decompose RDM, or sfGen 
          } 
          vector<femto::SQterm> temp2;
          femto::screenTerms(batch2, &temp2); // Screen terms with negligible factor
  
          for(vector<femto::SQterm>::iterator b_it = temp2.begin();b_it != temp2.end();++b_it) 
            b_it->transform2RDM(); // Transform sfGen to RDM (only if isInCanonical is true)
  
          result.insert(result.end(), temp2.begin(), temp2.end()); //*SLOW* 
        } // End if

      } // End q
    } // End p
  
    // Two-body part ....
    int c = 0;
    for(vector<femto::SQindex>::iterator p = Ps.begin();p != Ps.end();++p){
      for(vector<femto::SQindex>::iterator q = Qs.begin();q != Qs.end();++q){
        for(vector<femto::SQindex>::iterator r = Rs.begin();r != Rs.end();++r){
          for(vector<femto::SQindex>::iterator s = Ss.begin();s != Ss.end();++s){
  
            std::vector<femto::SQindex*> EH_indices;
            EH_indices.push_back(&(*p));
            EH_indices.push_back(&(*q));
            EH_indices.push_back(&(*r));
            EH_indices.push_back(&(*s));
  
            femto::SQtensor V("V2", EH_indices, femto::h2_symm());
            femto::sfGen EH(EH_indices);
  
            // < Psi | EL V EH ER | Psi >
            vector<femto::SQtensor> ten1;
            ten1.push_back(EL);
            ten1.push_back(V);
            ten1.push_back(EH);
            ten1.push_back(ER);
                  
            femto::SQterm term1(0.5, coeff1, ten1);
  
            vector<femto::SQterm> batch; batch.push_back(term1);
            femto::normalOrder(&batch); 
  
            for(vector<femto::SQterm>::iterator b_it = batch.begin();b_it != batch.end();++b_it){ 
              b_it->contractkDeltas(); // Burst Kronecker's deltas 
              b_it->decomposeRDMVirt();    // Decompose RDM, or sfGen 
            }
  
            vector<femto::SQterm> temp;
            femto::screenTerms(batch, &temp); // Screen terms with negligible factor
  
            for(vector<femto::SQterm>::iterator b_it = temp.begin();b_it != temp.end();++b_it) 
  	        { b_it->masquerade(); cout << *b_it << endl; b_it->transform2RDM(); }// Transform sfGen to RDM (only if isInCanonical is true)
            result.insert(result.end(), temp.begin(), temp.end());

            if(names[num_l] == "ooov"){
              // < Psi | EL ER V EH | Psi >
              vector<femto::SQtensor> ten2;
              ten2.push_back(EL);
              ten2.push_back(ER);
              ten2.push_back(V);
              ten2.push_back(EH);
                  
              femto::SQterm term2(-0.5, coeff1, ten2);
  	        //cout << "this " << term2 << endl;
              vector<femto::SQterm> batch2; batch2.push_back(term2);
              femto::normalOrder(&batch2); 
  
              for(vector<femto::SQterm>::iterator b_it = batch2.begin();b_it != batch2.end();++b_it){ 
                b_it->contractkDeltas(); // Burst Kronecker's deltas 
                b_it->decomposeRDMVirt();    // Decompose RDM, or sfGen 
              }
  
              vector<femto::SQterm> temp2;
              femto::screenTerms(batch2, &temp2); // Screen terms with negligible factor
  
       	      for(vector<femto::SQterm>::iterator b_it = temp2.begin();b_it != temp2.end();++b_it)
              // Transform sfGen to RDM (only if isInCanonical is true) 
  	        { b_it->masquerade(); cout << *b_it << endl; b_it->transform2RDM(); }
  
              result.insert(result.end(), temp2.begin(), temp2.end());
  	    } // End if

          } // End s
        } // End r
      } // End q
    } // End p

    if(names[num_l] == "ooov"){
    // Zeroth-body part
    // Ecas < Psi | EL T2 ER | Psi >
      vector<femto::SQtensor> ten;
      ten.push_back(EL);
      ten.push_back(ER);

      vector<string> Ecas;
      Ecas.push_back("Ecas");

      femto::SQterm term(1.0, Ecas, ten);
      vector<femto::SQterm> batch; batch.push_back(term);
      femto::normalOrder(&batch);

      for(vector<femto::SQterm>::iterator b_it = batch.begin();b_it != batch.end();++b_it){ 
        b_it->contractkDeltas(); // Burst Kronecker's deltas 
        b_it->decomposeRDMVirt();    // Decompose RDM, or sfGen 
        b_it->masquerade();
      }

      vector<femto::SQterm> temp;
      femto::screenTerms(batch, &temp); // Screen terms with negligible factor

      for(vector<femto::SQterm>::iterator b_it = temp.begin();b_it != temp.end();++b_it) 
        { cout << *b_it << endl; b_it->transform2RDM(); }// Transform sfGen to RDM (only if isInCanonical is true)
      result.insert(result.end(), temp.begin(), temp.end());

    } // End if

    vector<femto::SQterm> combined_result;
    femto::combineTerms(result, &combined_result);
    
    cout << endl << "Decompose RDMs ..... " << endl;
    int cnt = 0;
    vector<femto::SQterm> result2;
    for(vector<femto::SQterm>::iterator t = combined_result.begin();t != combined_result.end();++t){
      vector<femto::SQterm> batch;
      //cout << boost::format("%5d : ") % cnt++ << *t << endl;
      femto::decomposeRDMCore(*t, &batch);
      vector<femto::SQterm>::iterator new_t = batch.begin();
      for(;new_t != batch.end();++new_t){
        new_t->contractkDeltas();
        new_t->transform2RDM();
      } // End new_t
      result2.insert(result2.end(), batch.begin(), batch.end());
    } // End t
  
    vector<femto::SQterm> combined_result2;
    femto::combineTerms(result2, &combined_result2);
  
    cout << endl;
    cnt = 0;
        
    // TeX conversion
    femto::SQtensor Hdiag("Hdiag", *Li, femto::u4_symm());
    formulas << "<" + names[num_l] + "\\" + names[num_r] + ">" << endl; 
    formulas << Hdiag.convert2LaTeX() << "&=& ";  
    for(vector<femto::SQterm>::iterator t = combined_result2.begin();t != combined_result2.end();++t, ++cnt){
      formulas << t->convert2LaTeX();
      if((cnt+1)%6==0) formulas << "\\\\" << endl << "& &";
    } // End t
    formulas << endl << endl;
      
    int count2 = 0;
    cout << endl;
    cout << "< RESULT >" << endl;
    for(vector<femto::SQterm>::iterator t = combined_result2.begin();t != combined_result2.end();++t){
      cout << boost::format("%5d : ") % count2 << *t << endl;
      ++count2;
    }
    cout << endl;
  
    // ****** Conversion of tensor into the code starts ********
    count2 = 0;
    cout << endl;
    cout << "< RESULT2 >" << endl;
    // Convert Dirac 2 Mulliken, do not do twice!
    for(vector<femto::SQterm>::iterator t = combined_result2.begin();t != combined_result2.end();++t){
      for(size_t i = 0;i < t->get_tensors().size();++i){
        if(t->get_tensors()[i].get_name() == "V2") t->get_tensors_ptr()[i]->convertD2M();
        if(femto::is_RDM(t->get_tensors()[i].get_name())) t->get_tensors_ptr()[i]->convertD2M();
      } // End i
      cout << boost::format("%5d : ") % count2 << *t << endl;
      ++count2;
    }
      
    count2 = 0;
    cout << endl;
    vector<femto::SQterm> res;
    femto::find_JK_pairs(combined_result2, res);
    cout << "< RESULT3 >" << endl;
    for(vector<femto::SQterm>::iterator t = res.begin();t != res.end();++t){
      cout << boost::format("%5d : ") % count2 << *t << endl;
      ++count2;
    }
  
    int count = 0;
    cout << endl;
    cout << "< RESULT >" << endl;
    for(vector<femto::SQterm>::iterator t = res.begin();t != res.end();++t){
      //cout << count << " : " << *t << endl;
      cout << boost::format("%5d : ") % count << *t << endl;
      ++count;
    }
  
    // ********* Construct Sigma *************
    string thisname("diag_" + names[num_l]);
    femto::simpleloops(Hdiag, res, thisname, true, "V2", "T2");
    // ********* Construct Sigma *************
  } // End Li  

  num_ham = 0;
  cout << "* Calculation of the Sigma_{aa'}^{ee'} += <Psi0|EL H T0 |Psi0> ..... " << endl << endl;
  for(vector<vector<femto::SQindex*> >::iterator Li = EL_indices.begin();Li != EL_indices.end();++Li){
    femto::sfGen EL(*Li);

    // T2 amplitude
    size_t num_l = (size_t)(Li-EL_indices.begin());
    size_t num_r = (size_t)(Li-EL_indices.begin());
  
    cout << "* " << num_ham++ << " <" + names[num_l] + "/g>" << endl << endl;
  
    // T0 amplitude
    vector<string> coeffT0;
    coeffT0.push_back("T0");

    vector<femto::SQterm> result; result.reserve(femto::Nterms());
    //vector<femto::SQterm> CoreTerms; CoreTerms.reserve(femto::Nterms());
  
    // One-body part ....
    for(vector<femto::SQindex>::iterator p = Ps.begin();p != Ps.end();++p){
      for(vector<femto::SQindex>::iterator q = Qs.begin();q != Qs.end();++q){
  
        std::vector<femto::SQindex*> EH_indices;
        EH_indices.push_back(&(*p));
        EH_indices.push_back(&(*q));
  
        femto::SQtensor h("h", EH_indices, femto::h1_symm());
        femto::sfGen EH(EH_indices);
  
        // < Psi | EL h EH T0 | Psi >
        vector<femto::SQtensor> ten1;
        ten1.push_back(EL);
        ten1.push_back(h);
        ten1.push_back(EH);
  
        femto::SQterm term1(1.0, coeffT0, ten1);
        vector<femto::SQterm> batch; batch.push_back(term1);
        femto::normalOrder(&batch); 
  
        for(vector<femto::SQterm>::iterator b_it = batch.begin();b_it != batch.end();++b_it){ 
          b_it->contractkDeltas(); // Burst Kronecker's deltas 
          b_it->decomposeRDMVirt();    // Decompose RDM, or sfGen 
        } 
        vector<femto::SQterm> temp;
        femto::screenTerms(batch, &temp); // Screen terms with negligible factor
  
        for(vector<femto::SQterm>::iterator b_it = temp.begin();b_it != temp.end();++b_it) 
          b_it->transform2RDM(); // Transform sfGen to RDM (only if isInCanonical is true)
  
        result.insert(result.end(), temp.begin(), temp.end()); //*SLOW* 

      } // End q
    } // End p
  
    // Two-body part ....
    int c = 0;
    for(vector<femto::SQindex>::iterator p = Ps.begin();p != Ps.end();++p){
      for(vector<femto::SQindex>::iterator q = Qs.begin();q != Qs.end();++q){
        for(vector<femto::SQindex>::iterator r = Rs.begin();r != Rs.end();++r){
          for(vector<femto::SQindex>::iterator s = Ss.begin();s != Ss.end();++s){
  
            std::vector<femto::SQindex*> EH_indices;
            EH_indices.push_back(&(*p));
            EH_indices.push_back(&(*q));
            EH_indices.push_back(&(*r));
            EH_indices.push_back(&(*s));
  
            femto::SQtensor V("V2", EH_indices, femto::h2_symm());
            femto::sfGen EH(EH_indices);
  
            // < Psi | EL V EH T0 | Psi >
            vector<femto::SQtensor> ten1;
            ten1.push_back(EL);
            ten1.push_back(V);
            ten1.push_back(EH);
                  
            femto::SQterm term1(0.5, coeffT0, ten1);
  
            vector<femto::SQterm> batch; batch.push_back(term1);
            femto::normalOrder(&batch); 
  
            for(vector<femto::SQterm>::iterator b_it = batch.begin();b_it != batch.end();++b_it){ 
              b_it->contractkDeltas(); // Burst Kronecker's deltas 
              b_it->decomposeRDMVirt();    // Decompose RDM, or sfGen 
            }
  
            vector<femto::SQterm> temp;
            femto::screenTerms(batch, &temp); // Screen terms with negligible factor
  
            for(vector<femto::SQterm>::iterator b_it = temp.begin();b_it != temp.end();++b_it) 
  	      { b_it->masquerade(); cout << *b_it << endl; b_it->transform2RDM(); }// Transform sfGen to RDM (only if isInCanonical is true)
            result.insert(result.end(), temp.begin(), temp.end());

          } // End s
        } // End r
      } // End q
    } // End p
  
    vector<femto::SQterm> combined_result;
    femto::combineTerms(result, &combined_result);
    
    cout << endl << "Decompose RDMs ..... " << endl;
    int cnt = 0;
    vector<femto::SQterm> result2;
    for(vector<femto::SQterm>::iterator t = combined_result.begin();t != combined_result.end();++t){
      vector<femto::SQterm> batch;
      //cout << boost::format("%5d : ") % cnt++ << *t << endl;
      femto::decomposeRDMCore(*t, &batch);
      vector<femto::SQterm>::iterator new_t = batch.begin();
      for(;new_t != batch.end();++new_t){
        new_t->contractkDeltas();
        new_t->transform2RDM();
      } // End new_t
      result2.insert(result2.end(), batch.begin(), batch.end());
    } // End t
  
    vector<femto::SQterm> combined_result2;
    femto::combineTerms(result2, &combined_result2);
  
    cout << endl;
    cnt = 0;
        
    // TeX conversion
    femto::SQtensor Sig("S2", *Li, femto::u4_symm());
    formulas << "<" + names[num_l] + "\\g>" << endl; 
    formulas << Sig.convert2LaTeX() << "&=& ";  
    for(vector<femto::SQterm>::iterator t = combined_result2.begin();t != combined_result2.end();++t, ++cnt){
      formulas << t->convert2LaTeX();
      if((cnt+1)%6==0) formulas << "\\\\" << endl << "& &";
    } // End t
    formulas << endl << endl;
      
    int count2 = 0;
    cout << endl;
    cout << "< RESULT >" << endl;
    for(vector<femto::SQterm>::iterator t = combined_result2.begin();t != combined_result2.end();++t){
      cout << boost::format("%5d : ") % count2 << *t << endl;
      ++count2;
    }
    cout << endl;
  
    // ****** Conversion of tensor into the code starts ********
    count2 = 0;
    cout << endl;
    cout << "< RESULT2 >" << endl;
    // Convert Dirac 2 Mulliken, do not do twice!
    for(vector<femto::SQterm>::iterator t = combined_result2.begin();t != combined_result2.end();++t){
      for(size_t i = 0;i < t->get_tensors().size();++i){
        if(t->get_tensors()[i].get_name() == "V2") t->get_tensors_ptr()[i]->convertD2M();
        if(femto::is_RDM(t->get_tensors()[i].get_name())) t->get_tensors_ptr()[i]->convertD2M();
      } // End i
      cout << boost::format("%5d : ") % count2 << *t << endl;
      ++count2;
    }
      
    count2 = 0;
    cout << endl;
    vector<femto::SQterm> res;
    femto::find_JK_pairs(combined_result2, res);
    cout << "< RESULT3 >" << endl;
    for(vector<femto::SQterm>::iterator t = res.begin();t != res.end();++t){
      cout << boost::format("%5d : ") % count2 << *t << endl;
      ++count2;
    }
  
    int count = 0;
    cout << endl;
    cout << "< RESULT >" << endl;
    for(vector<femto::SQterm>::iterator t = res.begin();t != res.end();++t){
      //cout << count << " : " << *t << endl;
      cout << boost::format("%5d : ") % count << *t << endl;
      ++count;
    }
  
    // ********* Construct Sigma *************
    string thisname("sigma_" + names[num_l] + "_g");
    femto::factorize(Sig, res, thisname, true, "V2", "T2");
    // ********* Construct Sigma *************
  } // End Li  

  num_ham = 0;
  cout << "* Calculation of the Sigma_{0} += <Psi0| H T2 ER |Psi0> ..... " << endl << endl;
  for(vector<vector<femto::SQindex*> >::iterator Ri = ER_indices.begin();Ri != ER_indices.end();++Ri){
    femto::sfGen ER(*Ri);

    // T2 amplitude
    size_t num_r = (size_t)(Ri-ER_indices.begin());
    femto::SQtensor Tamp("T2", ET_indices[num_r], femto::u4_symm());
  
    cout << "* " << num_ham++ << " <g/" + names[num_r] + ">" << endl << endl;
  
    vector<string> coeff1;

    vector<femto::SQterm> result; result.reserve(femto::Nterms());
    //vector<femto::SQterm> CoreTerms; CoreTerms.reserve(femto::Nterms());
  
    // One-body part ....
    for(vector<femto::SQindex>::iterator p = Ps.begin();p != Ps.end();++p){
      for(vector<femto::SQindex>::iterator q = Qs.begin();q != Qs.end();++q){
  
        std::vector<femto::SQindex*> EH_indices;
        EH_indices.push_back(&(*p));
        EH_indices.push_back(&(*q));
  
        femto::SQtensor h("h", EH_indices, femto::h1_symm());
        femto::sfGen EH(EH_indices);
  
        // < Psi | h EH T2 ER | Psi >
        vector<femto::SQtensor> ten1;
        ten1.push_back(h);
        ten1.push_back(EH);
        ten1.push_back(Tamp);
        ten1.push_back(ER);
  
        femto::SQterm term1(1.0, coeff1, ten1);
        vector<femto::SQterm> batch; batch.push_back(term1);
        femto::normalOrder(&batch); 
  
        for(vector<femto::SQterm>::iterator b_it = batch.begin();b_it != batch.end();++b_it){ 
          b_it->contractkDeltas(); // Burst Kronecker's deltas 
          b_it->decomposeRDMVirt();    // Decompose RDM, or sfGen 
        } 
        vector<femto::SQterm> temp;
        femto::screenTerms(batch, &temp); // Screen terms with negligible factor
  
        for(vector<femto::SQterm>::iterator b_it = temp.begin();b_it != temp.end();++b_it) 
          b_it->transform2RDM(); // Transform sfGen to RDM (only if isInCanonical is true)
  
        result.insert(result.end(), temp.begin(), temp.end()); //*SLOW* 

      } // End q
    } // End p
  
    // Two-body part ....
    int c = 0;
    for(vector<femto::SQindex>::iterator p = Ps.begin();p != Ps.end();++p){
      for(vector<femto::SQindex>::iterator q = Qs.begin();q != Qs.end();++q){
        for(vector<femto::SQindex>::iterator r = Rs.begin();r != Rs.end();++r){
          for(vector<femto::SQindex>::iterator s = Ss.begin();s != Ss.end();++s){
  
            std::vector<femto::SQindex*> EH_indices;
            EH_indices.push_back(&(*p));
            EH_indices.push_back(&(*q));
            EH_indices.push_back(&(*r));
            EH_indices.push_back(&(*s));
  
            femto::SQtensor V("V2", EH_indices, femto::h2_symm());
            femto::sfGen EH(EH_indices);
  
            // < Psi | V EH T2 ER | Psi >
            vector<femto::SQtensor> ten1;
            ten1.push_back(V);
            ten1.push_back(EH);
            ten1.push_back(Tamp);
            ten1.push_back(ER);
                  
            femto::SQterm term1(0.5, coeff1, ten1);
  
            vector<femto::SQterm> batch; batch.push_back(term1);
            femto::normalOrder(&batch); 
  
            for(vector<femto::SQterm>::iterator b_it = batch.begin();b_it != batch.end();++b_it){ 
              b_it->contractkDeltas(); // Burst Kronecker's deltas 
              b_it->decomposeRDMVirt();    // Decompose RDM, or sfGen 
            }
  
            vector<femto::SQterm> temp;
            femto::screenTerms(batch, &temp); // Screen terms with negligible factor
  
            for(vector<femto::SQterm>::iterator b_it = temp.begin();b_it != temp.end();++b_it) 
  	        { b_it->masquerade(); cout << *b_it << endl; b_it->transform2RDM(); }// Transform sfGen to RDM (only if isInCanonical is true)
            result.insert(result.end(), temp.begin(), temp.end());

          } // End s
        } // End r
      } // End q
    } // End p
  
    vector<femto::SQterm> combined_result;
    femto::combineTerms(result, &combined_result);
    
    cout << endl << "Decompose RDMs ..... " << endl;
    int cnt = 0;
    vector<femto::SQterm> result2;
    for(vector<femto::SQterm>::iterator t = combined_result.begin();t != combined_result.end();++t){
      vector<femto::SQterm> batch;
      //cout << boost::format("%5d : ") % cnt++ << *t << endl;
      femto::decomposeRDMCore(*t, &batch);
      vector<femto::SQterm>::iterator new_t = batch.begin();
      for(;new_t != batch.end();++new_t){
        new_t->contractkDeltas();
        new_t->transform2RDM();
      } // End new_t
      result2.insert(result2.end(), batch.begin(), batch.end());
    } // End t
  
    vector<femto::SQterm> combined_result2;
    femto::combineTerms(result2, &combined_result2);
  
    cout << endl;
    cnt = 0;
        
    // TeX conversion
    vector<femto::SQindex*> null;
    femto::Symmetry null_symm;
    femto::SQtensor Sig("S0", null, null_symm);
    formulas << "<g\\" + names[num_r] + ">" << endl; 
    formulas << Sig.convert2LaTeX() << "&=& ";  
    for(vector<femto::SQterm>::iterator t = combined_result2.begin();t != combined_result2.end();++t, ++cnt){
      formulas << t->convert2LaTeX();
      if((cnt+1)%6==0) formulas << "\\\\" << endl << "& &";
    } // End t
    formulas << endl << endl;
      
    int count2 = 0;
    cout << endl;
    cout << "< RESULT >" << endl;
    for(vector<femto::SQterm>::iterator t = combined_result2.begin();t != combined_result2.end();++t){
      cout << boost::format("%5d : ") % count2 << *t << endl;
      ++count2;
    }
    cout << endl;
  
    // ****** Conversion of tensor into the code starts ********
    count2 = 0;
    cout << endl;
    cout << "< RESULT2 >" << endl;
    // Convert Dirac 2 Mulliken, do not do twice!
    for(vector<femto::SQterm>::iterator t = combined_result2.begin();t != combined_result2.end();++t){
      for(size_t i = 0;i < t->get_tensors().size();++i){
        if(t->get_tensors()[i].get_name() == "V2") t->get_tensors_ptr()[i]->convertD2M();
        if(femto::is_RDM(t->get_tensors()[i].get_name())) t->get_tensors_ptr()[i]->convertD2M();
      } // End i
      cout << boost::format("%5d : ") % count2 << *t << endl;
      ++count2;
    }
      
    count2 = 0;
    cout << endl;
    vector<femto::SQterm> res;
    femto::find_JK_pairs(combined_result2, res);
    cout << "< RESULT3 >" << endl;
    for(vector<femto::SQterm>::iterator t = res.begin();t != res.end();++t){
      cout << boost::format("%5d : ") % count2 << *t << endl;
      ++count2;
    }
  
    int count = 0;
    cout << endl;
    cout << "< RESULT >" << endl;
    for(vector<femto::SQterm>::iterator t = res.begin();t != res.end();++t){
      //cout << count << " : " << *t << endl;
      cout << boost::format("%5d : ") % count << *t << endl;
      ++count;
    }
  
    // ********* Construct Sigma *************
    string thisname("sigma_g_" + names[num_r]);
    femto::factorize(Sig, res, thisname, false, "V2", "T2");
    // ********* Construct Sigma *************
  } // End Ri  

  // Construction of Sigma_{aa'}^{ee'} += <Psi0|EL H TR ER |Psi0>
  int num_block = 0;
  cout << "* Construction of Sigma_{aa'}^{ee'} += <Psi0|EL H TR ER |Psi0> ..... " << endl << endl;
  for(vector<vector<femto::SQindex*> >::iterator Li = EL_indices.begin();Li != EL_indices.end();++Li){
    for(vector<vector<femto::SQindex*> >::iterator Ri = ER_indices.begin();Ri != ER_indices.end();++Ri){
      femto::sfGen EL(*Li);
      femto::sfGen ER(*Ri);

      // T2 amplitude
      size_t num_l = (size_t)(Li-EL_indices.begin());
      size_t num_r = (size_t)(Ri-ER_indices.begin());
      femto::SQtensor Tamp("T2", ET_indices[num_r], femto::u4_symm());
    
      cout << "* " << num_block++ << " <" + names[num_l] + "/" + names[num_r] + ">" << endl << endl;
    
      // Unit coefficient
      vector<string> coeff1;
    
      vector<femto::SQterm> result; result.reserve(femto::Nterms());
      //vector<femto::SQterm> CoreTerms; CoreTerms.reserve(femto::Nterms());
    
      // One-body part ....
      for(vector<femto::SQindex>::iterator p = Ps.begin();p != Ps.end();++p){
        for(vector<femto::SQindex>::iterator q = Qs.begin();q != Qs.end();++q){
    
          std::vector<femto::SQindex*> EH_indices;
          EH_indices.push_back(&(*p));
          EH_indices.push_back(&(*q));
    
          femto::SQtensor h("h", EH_indices, femto::h1_symm());
          femto::sfGen EH(EH_indices);
    
          // < Psi | EL h EH T2 ER | Psi >
          vector<femto::SQtensor> ten1;
          ten1.push_back(EL);
          ten1.push_back(h);
          ten1.push_back(EH);
          ten1.push_back(Tamp);
          ten1.push_back(ER);
    
          femto::SQterm term1(1.0, coeff1, ten1);
          vector<femto::SQterm> batch; batch.push_back(term1);
          femto::normalOrder(&batch); 
    
          for(vector<femto::SQterm>::iterator b_it = batch.begin();b_it != batch.end();++b_it){ 
            b_it->contractkDeltas(); // Burst Kronecker's deltas 
            b_it->decomposeRDMVirt();    // Decompose RDM, or sfGen 
          } 
          vector<femto::SQterm> temp;
          femto::screenTerms(batch, &temp); // Screen terms with negligible factor
    
          for(vector<femto::SQterm>::iterator b_it = temp.begin();b_it != temp.end();++b_it) 
            b_it->transform2RDM(); // Transform sfGen to RDM (only if isInCanonical is true)
    
          result.insert(result.end(), temp.begin(), temp.end()); //*SLOW* 
    
          if(names[num_l] == names[num_r] && names[num_r] == "ooov"){
            // < Psi | EL T2 ER EH | Psi >
            vector<femto::SQtensor> ten2;
            ten2.push_back(EL);
            ten2.push_back(Tamp);
            ten2.push_back(ER);
            ten2.push_back(h);
            ten2.push_back(EH);
    
            femto::SQterm term2(-1.0, coeff1, ten2);
            vector<femto::SQterm> batch2; batch2.push_back(term2);
            femto::normalOrder(&batch2); 
    
            for(vector<femto::SQterm>::iterator b_it = batch2.begin();b_it != batch2.end();++b_it){ 
              b_it->contractkDeltas(); // Burst Kronecker's deltas 
              b_it->decomposeRDMVirt();    // Decompose RDM, or sfGen 
            } 
            vector<femto::SQterm> temp2;
            femto::screenTerms(batch2, &temp2); // Screen terms with negligible factor
    
            for(vector<femto::SQterm>::iterator b_it = temp2.begin();b_it != temp2.end();++b_it) 
              b_it->transform2RDM(); // Transform sfGen to RDM (only if isInCanonical is true)
    
            result.insert(result.end(), temp2.begin(), temp2.end()); //*SLOW* 
          } // End if
    
        } // End q
      } // End p
    
      // Two-body part ....
      int c = 0;
      for(vector<femto::SQindex>::iterator p = Ps.begin();p != Ps.end();++p){
        for(vector<femto::SQindex>::iterator q = Qs.begin();q != Qs.end();++q){
          for(vector<femto::SQindex>::iterator r = Rs.begin();r != Rs.end();++r){
            for(vector<femto::SQindex>::iterator s = Ss.begin();s != Ss.end();++s){
    
              std::vector<femto::SQindex*> EH_indices;
              EH_indices.push_back(&(*p));
              EH_indices.push_back(&(*q));
              EH_indices.push_back(&(*r));
              EH_indices.push_back(&(*s));
    
              femto::SQtensor V("V2", EH_indices, femto::h2_symm());
              femto::sfGen EH(EH_indices);
    
              // < Psi | EL V EH T2 ER | Psi >
              vector<femto::SQtensor> ten1;
              ten1.push_back(EL);
              ten1.push_back(V);
              ten1.push_back(EH);
              ten1.push_back(Tamp);
              ten1.push_back(ER);
                    
              femto::SQterm term1(0.5, coeff1, ten1);
    
              vector<femto::SQterm> batch; batch.push_back(term1);
              femto::normalOrder(&batch); 
    
              for(vector<femto::SQterm>::iterator b_it = batch.begin();b_it != batch.end();++b_it){ 
                b_it->contractkDeltas(); // Burst Kronecker's deltas 
                b_it->decomposeRDMVirt();    // Decompose RDM, or sfGen 
              }
    
              vector<femto::SQterm> temp;
              femto::screenTerms(batch, &temp); // Screen terms with negligible factor
    
              for(vector<femto::SQterm>::iterator b_it = temp.begin();b_it != temp.end();++b_it) 
    	        { b_it->masquerade(); cout << *b_it << endl; b_it->transform2RDM(); }// Transform sfGen to RDM (only if isInCanonical is true)
              result.insert(result.end(), temp.begin(), temp.end());
    
              if(names[num_l] == names[num_r] && names[num_r] == "ooov"){
                // < Psi | EL T2 ER V EH | Psi >
                vector<femto::SQtensor> ten2;
                ten2.push_back(EL);
                ten2.push_back(Tamp);
                ten2.push_back(ER);
                ten2.push_back(V);
                ten2.push_back(EH);
                    
                femto::SQterm term2(-0.5, coeff1, ten2);
    	        //cout << "this " << term2 << endl;
                vector<femto::SQterm> batch2; batch2.push_back(term2);
                femto::normalOrder(&batch2); 
    
                for(vector<femto::SQterm>::iterator b_it = batch2.begin();b_it != batch2.end();++b_it){ 
                  b_it->contractkDeltas(); // Burst Kronecker's deltas 
                  b_it->decomposeRDMVirt();    // Decompose RDM, or sfGen 
                }
    
                vector<femto::SQterm> temp2;
                femto::screenTerms(batch2, &temp2); // Screen terms with negligible factor
    
         	for(vector<femto::SQterm>::iterator b_it = temp2.begin();b_it != temp2.end();++b_it)
                // Transform sfGen to RDM (only if isInCanonical is true) 
    	          { b_it->masquerade(); cout << *b_it << endl; b_it->transform2RDM(); }
    
                result.insert(result.end(), temp2.begin(), temp2.end());
    	      } // End if
    
            } // End s
          } // End r
        } // End q
      } // End p
    
      if(names[num_l] == names[num_r] && names[num_l] == "ooov"){
      // Zeroth-body part
      // Ecas < Psi | EL T2 ER | Psi >
        vector<femto::SQtensor> ten;
        ten.push_back(EL);
        ten.push_back(Tamp);
        ten.push_back(ER);
  
        vector<string> Ecas;
        Ecas.push_back("Ecas");
  
        femto::SQterm term(1.0, Ecas, ten);
        vector<femto::SQterm> batch; batch.push_back(term);
        femto::normalOrder(&batch);
  
        for(vector<femto::SQterm>::iterator b_it = batch.begin();b_it != batch.end();++b_it){ 
          b_it->contractkDeltas(); // Burst Kronecker's deltas 
          b_it->decomposeRDMVirt();    // Decompose RDM, or sfGen 
          b_it->masquerade();
        }
  
        vector<femto::SQterm> temp;
        femto::screenTerms(batch, &temp); // Screen terms with negligible factor
  
        for(vector<femto::SQterm>::iterator b_it = temp.begin();b_it != temp.end();++b_it) 
          { cout << *b_it << endl; b_it->transform2RDM(); }// Transform sfGen to RDM (only if isInCanonical is true)
        result.insert(result.end(), temp.begin(), temp.end());
      } // End if

      vector<femto::SQterm> combined_result;
      femto::combineTerms(result, &combined_result);
      
      cout << endl << "Decompose RDMs ..... " << endl;
      int cnt = 0;
      vector<femto::SQterm> result2;
      for(vector<femto::SQterm>::iterator t = combined_result.begin();t != combined_result.end();++t){
        vector<femto::SQterm> batch;
        //cout << boost::format("%5d : ") % cnt++ << *t << endl;
        femto::decomposeRDMCore(*t, &batch);
        vector<femto::SQterm>::iterator new_t = batch.begin();
        for(;new_t != batch.end();++new_t){
          new_t->contractkDeltas();
          new_t->transform2RDM();
        } // End new_t
        result2.insert(result2.end(), batch.begin(), batch.end());
      } // End t
    
      vector<femto::SQterm> combined_result2;
      femto::combineTerms(result2, &combined_result2);
    
      cout << endl;
      cnt = 0;
          
      // TeX conversion
      femto::SQtensor Sig("S2", *Li, femto::u4_symm());
      formulas << "<" + names[num_l] + "\\" + names[num_r] + ">" << endl; 
      formulas << Sig.convert2LaTeX() << "&=& ";  
      for(vector<femto::SQterm>::iterator t = combined_result2.begin();t != combined_result2.end();++t, ++cnt){
        formulas << t->convert2LaTeX();
        if((cnt+1)%6==0) formulas << "\\\\" << endl << "& &";
      } // End t
      formulas << endl << endl;
        
      int count2 = 0;
      cout << endl;
      cout << "< RESULT >" << endl;
      for(vector<femto::SQterm>::iterator t = combined_result2.begin();t != combined_result2.end();++t){
        cout << boost::format("%5d : ") % count2 << *t << endl;
        ++count2;
      }
      cout << endl;
    
      // ****** Conversion of tensor into the code starts ********
      count2 = 0;
      cout << endl;
      cout << "< RESULT2 >" << endl;
      // Convert Dirac 2 Mulliken, do not do twice!
      for(vector<femto::SQterm>::iterator t = combined_result2.begin();t != combined_result2.end();++t){
        for(size_t i = 0;i < t->get_tensors().size();++i){
          if(t->get_tensors()[i].get_name() == "V2") t->get_tensors_ptr()[i]->convertD2M();
          if(femto::is_RDM(t->get_tensors()[i].get_name())) t->get_tensors_ptr()[i]->convertD2M();
        } // End i
        cout << boost::format("%5d : ") % count2 << *t << endl;
        ++count2;
      }
        
      count2 = 0;
      cout << endl;
      vector<femto::SQterm> res;
      femto::find_JK_pairs(combined_result2, res);
      cout << "< RESULT3 >" << endl;
      for(vector<femto::SQterm>::iterator t = res.begin();t != res.end();++t){
        cout << boost::format("%5d : ") % count2 << *t << endl;
        ++count2;
      }
    
      int count = 0;
      cout << endl;
      cout << "< RESULT >" << endl;
      for(vector<femto::SQterm>::iterator t = res.begin();t != res.end();++t){
        //cout << count << " : " << *t << endl;
        cout << boost::format("%5d : ") % count << *t << endl;
        ++count;
      }
    
      // ********* Construct Sigma *************
      string thisname("sigma_" + names[num_l] + "_" + names[num_r]);
      femto::factorize(Sig, res, thisname, true, "V2", "T2");
      // ********* Construct Sigma *************
    } // End Li
  } // End Ri

}

