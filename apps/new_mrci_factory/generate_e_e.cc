
#include "Eclipse.hpp"

// MM""""""""`M          dP oo                            
// MM  mmmmmmmM          88                               
// M`      MMMM .d8888b. 88 dP 88d888b. .d8888b. .d8888b. 
// MM  MMMMMMMM 88'  `"" 88 88 88'  `88 Y8ooooo. 88ooood8 
// MM  MMMMMMMM 88.  ... 88 88 88.  .88       88 88.  ... 
// MM        .M `88888P' dP dP 88Y888P' `88888P' `88888P' 
// MMMMMMMMMMMM                88                         

using namespace std;

//#define _ONLY_OV
//#define _DIAG_ONLY

void
Eclipse::generate_e_e()
{

  // MS: I forgot why I added the statement below, but maybe it's just nothing (2014/02/07)
  //cout << "!!!!!!!!!!!! Now debugging mode !!!!!!!!!!!!" << endl;

  // Construction of Sigma_{aa'}^{ee'} += <Psi0|EL H TR ER |Psi0>
  int num_block = 0;
  string name("< Construction of Sigma_{aa'}^{ee'} += <Psi0|EL H TR ER |Psi0> > \n\n");
  cout << name;
  status += name;

  print_eqns << " ++ Sigma vector of fic-MRCI (e-e) ::" << endl << endl;

  for(auto Li = EL_indices.begin();Li != EL_indices.end();++Li){
    for(auto Ri = ER_indices.begin();Ri != ER_indices.end();++Ri){

      Femto::Core::sfGen EL(*Li);
      Femto::Core::sfGen ER(*Ri);

      // T2 amplitude
      size_t num_l = (size_t)(Li-EL_indices.begin());
      size_t num_r = (size_t)(Ri-ER_indices.begin());
      //Femto::Symmetry symm(Femto::u4_symm());
      Femto::Symmetry symm(Femto::t2_symm());
//*--      {
//*-- 	// If the the amplitude is of T_{XX'}^{YY'} type excitation, allow permutational symmetry
//*--	if(ET_indices[num_r].at(0)->get_char() == ET_indices[num_r].at(1)->get_char() && 
//*--	   ET_indices[num_r].at(2)->get_char() == ET_indices[num_r].at(3)->get_char()) symm = Femto::t2_symm();
//*--	//if(names[num_l] == "ccvv" && names[num_r] == "covv")  symm = Femto::t2_symm();
//*--	//if(names[num_l] == "covv" && names[num_r] == "covv")  symm = Femto::t2_symm();
//*--      }       
      Femto::Core::SQtensor Tamp("T2", ET_indices[num_r], symm);
    
      ostringstream stm_l;
      ostringstream stm_r;
      ostringstream stm1;
      stm_l << num_l;
      stm_r << num_r;
      stm1  << num_block++;
      string title("* " + stm1.str() + " <" + names[num_l] + "/" + names[num_r] + ">  ");
      cout << "* " << num_block << " <" + names[num_l] + "/" + names[num_r] + ">" << endl << endl;
      status += title;

#ifdef _DIAG_ONLY        
      if(num_l != num_r /*|| names[num_l] != "cooo" /*&& names[num_l] != "ooov"*/) { cout << "DIAG ONLY" << endl; continue; } 
#endif

#ifdef _ONLY_OV
      if((names[num_l] != "oovv" && names[num_l] != "ooov") || (names[num_r] != "oovv" && names[num_r] != "ooov")) { cout << "ONLY OV" << endl; continue; }
#endif

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
          ten1.push_back(Tamp);
          ten1.push_back(ER);

          double numcoeff = 1.0;
          if(names[num_l] == names[num_r] && (names[num_r] == "ooov" || names[num_l] == "cooo"))
            numcoeff = 0.5;           
          Femto::Core::SQterm term1(numcoeff, coeff1, ten1);
          vector<Femto::Core::SQterm> batch; batch.push_back(term1);

#ifndef _NEW_CON
	  Femto::Core::normalOrder(&batch); 
#else
          if(names[num_l] == names[num_r] && (names[num_r] == "ooov" || names[num_l] == "cooo")){
	    Femto::Core::normalOrderComm(&batch); 
	  } else{
	    Femto::Core::normalOrder(&batch); 
	  }
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

          if(names[num_l] == names[num_r] && (names[num_r] == "ooov" || names[num_l] == "cooo")){

#ifndef _NEW_CON    
            // < Psi | EL T2 ER h EH | Psi >
            vector<Femto::Core::SQtensor> ten2;
            ten2.push_back(EL);
            ten2.push_back(Tamp);
            ten2.push_back(ER);
            ten2.push_back(h);
            ten2.push_back(EH);
    
            Femto::Core::SQterm term2(-0.5, coeff1, ten2);
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
#endif

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

#ifndef _NEW_CON    
            // -1/2 <Psi| EL_tilde T ER_tilde h Eh |Psi>
            vector<Femto::Core::SQtensor> ten3;
            ten3.push_back(EL_tilde);
            ten3.push_back(Tamp);
            ten3.push_back(ER_tilde);
            ten3.push_back(h);
            ten3.push_back(EH);

            Femto::Core::SQterm term3(-0.5, coeff1, ten3);
            vector<Femto::Core::SQterm> batch3; batch3.push_back(term3);
            Femto::Core::normalOrder(&batch3); 
    
            for(vector<Femto::Core::SQterm>::iterator b_it = batch3.begin();b_it != batch3.end();++b_it){ 
              b_it->contractkDeltas(); // Burst Kronecker's deltas 
              b_it->decomposeRDMVirt();    // Decompose RDM, or sfGen 
            } 
            vector<Femto::Core::SQterm> temp3;
            Femto::Core::screenTerms(batch3, &temp3); // Screen terms with negligible factor
    
            for(vector<Femto::Core::SQterm>::iterator b_it = temp3.begin();b_it != temp3.end();++b_it) 
              b_it->transform2RDM(); // Transform sfGen to RDM (only if isInCanonical is true)
    
            result.insert(result.end(), temp3.begin(), temp3.end()); //*SLOW*
#endif   
 
            // +1/2 <Psi| EL_tilde h Eh T ER_tilde |Psi>
            vector<Femto::Core::SQtensor> ten4;
            ten4.push_back(EL_tilde);
            ten4.push_back(h);
            ten4.push_back(EH);
            ten4.push_back(Tamp);
            ten4.push_back(ER_tilde);
    
            Femto::Core::SQterm term4(+0.5, coeff1, ten4);
            vector<Femto::Core::SQterm> batch4; batch4.push_back(term4);

#ifndef _NEW_CON    
            Femto::Core::normalOrder(&batch4); 
#else 
            Femto::Core::normalOrderComm(&batch4); 
#endif

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
    
              // +1/2 (or +1/4) < Psi | EL V EH T2 ER | Psi >
              vector<Femto::Core::SQtensor> ten1;
              ten1.push_back(EL);
              ten1.push_back(V);
              ten1.push_back(EH);
              ten1.push_back(Tamp);
              ten1.push_back(ER);

              double numcoeff(0.5);
              if(names[num_l] == names[num_r] && (names[num_r] == "ooov" || names[num_l] == "cooo"))
                numcoeff = 0.25;           
	      
              Femto::Core::SQterm term1(numcoeff, coeff1, ten1);
    
              vector<Femto::Core::SQterm> batch; batch.push_back(term1);

#ifdef _NEW_CON
              if(names[num_l] == names[num_r] && (names[num_r] == "ooov" || names[num_l] == "cooo"))
		Femto::Core::normalOrderComm(&batch);
	      else
		Femto::Core::normalOrder(&batch);
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
    	        { b_it->masquerade(); cout << *b_it << endl; b_it->transform2RDM(); }// Transform sfGen to RDM (only if isInCanonical is true)
              result.insert(result.end(), temp.begin(), temp.end());
    
              if(names[num_l] == names[num_r] && (names[num_r] == "ooov"  || names[num_l] == "cooo")){

#ifndef _NEW_CON
                // -1/4 < Psi | EL T2 ER V EH | Psi >
                vector<Femto::Core::SQtensor> ten2;
                ten2.push_back(EL);
                ten2.push_back(Tamp);
                ten2.push_back(ER);
                ten2.push_back(V);
                ten2.push_back(EH);
                    
                Femto::Core::SQterm term2(-0.25, coeff1, ten2);
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
#endif

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

#ifndef _NEW_CON
                // -1/4 <Psi| EL_tilde T ER_tilde V EH |Psi>
                vector<Femto::Core::SQtensor> ten3;
                ten3.push_back(EL_tilde);
                ten3.push_back(Tamp);
                ten3.push_back(ER_tilde);
                ten3.push_back(V);
                ten3.push_back(EH);
                    
                Femto::Core::SQterm term3(-0.25, coeff1, ten3);
    	        //cout << "this " << term2 << endl;
                vector<Femto::Core::SQterm> batch3; batch3.push_back(term3);
                Femto::Core::normalOrder(&batch3); 
    
                for(vector<Femto::Core::SQterm>::iterator b_it = batch3.begin();b_it != batch3.end();++b_it){ 
                  b_it->contractkDeltas(); // Burst Kronecker's deltas 
                  b_it->decomposeRDMVirt();    // Decompose RDM, or sfGen 
                }
    
                vector<Femto::Core::SQterm> temp3;
                Femto::Core::screenTerms(batch3, &temp3); // Screen terms with negligible factor
    
         	for(vector<Femto::Core::SQterm>::iterator b_it = temp3.begin();b_it != temp3.end();++b_it)
                // Transform sfGen to RDM (only if isInCanonical is true) 
    	          { b_it->masquerade(); cout << *b_it << endl; b_it->transform2RDM(); }
    
                result.insert(result.end(), temp3.begin(), temp3.end());
#endif

                // +1/4 <Psi| EL_tilde V EH T ER_tilde |Psi>
                vector<Femto::Core::SQtensor> ten4;
                ten4.push_back(EL_tilde);
                ten4.push_back(V);
                ten4.push_back(EH);
                ten4.push_back(Tamp);
                ten4.push_back(ER_tilde);
                    
                Femto::Core::SQterm term4(+0.25, coeff1, ten4);
    	        //cout << "this " << term2 << endl;
                vector<Femto::Core::SQterm> batch4; batch4.push_back(term4);

#ifdef _NEW_CON
                Femto::Core::normalOrderComm(&batch4); 
#else
                Femto::Core::normalOrder(&batch4); 
#endif
    
                for(vector<Femto::Core::SQterm>::iterator b_it = batch4.begin();b_it != batch4.end();++b_it){ 
                  b_it->contractkDeltas(); // Burst Kronecker's deltas 
                  b_it->decomposeRDMVirt();    // Decompose RDM, or sfGen 
                }
    
                vector<Femto::Core::SQterm> temp4;
                Femto::Core::screenTerms(batch4, &temp4); // Screen terms with negligible factor
    
         	for(vector<Femto::Core::SQterm>::iterator b_it = temp4.begin();b_it != temp4.end();++b_it)
                // Transform sfGen to RDM (only if isInCanonical is true) 
    	          { b_it->masquerade(); cout << *b_it << endl; b_it->transform2RDM(); }
    
                result.insert(result.end(), temp4.begin(), temp4.end());

    	      } // End if
    
            } // End s
          } // End r
        } // End q
      } // End p
    
      if(names[num_l] == names[num_r] && (names[num_l] == "ooov" || names[num_l] == "cooo")){
      // Zeroth-body part
      // Ecas < Psi | EL T2 ER | Psi >
        vector<Femto::Core::SQtensor> ten;
        ten.push_back(EL);
        ten.push_back(Tamp);
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
      //Femto::Core::combineTerms(result, &combined_result);
      
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
      symm = Femto::u4_symm();
      {
	pair<Femto::char_state, Femto::char_state> from_to;
	from_to.first  = Li->at(0)->get_char();
	from_to.second = Li->at(2)->get_char();

	// If the the amplitude and sigma are of T_{XX'}^{YY'} type excitation, allow permutational symmetry
	if(from_to.first  == Li->at(0)->get_char()               &&
	   from_to.first  == Li->at(1)->get_char()               &&
	   from_to.first  == ET_indices[num_r].at(0)->get_char() &&
	   from_to.first  == ET_indices[num_r].at(1)->get_char() &&
	   from_to.second == Li->at(2)->get_char()               &&
	   from_to.second == Li->at(3)->get_char()               &&
	   from_to.second == ET_indices[num_r].at(2)->get_char() &&
	   from_to.second == ET_indices[num_r].at(3)->get_char())
	  symm = Femto::t2_symm();

      }       
      Femto::Core::SQtensor Sig("S2", *Li, symm);
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

      // Count the number of the 4-RDM
      const std::string name_d4(Femto::RDM_name()+"4");
      std::vector<Femto::Core::SQterm> d4s;
      for(std::vector<Femto::Core::SQterm>::iterator t = combined_result2.begin();t != combined_result2.end();++t){
        bool has_d4(false);
	for(size_t num_t = 0;num_t < t->get_tensors().size();++num_t){
	  if(t->get_tensors()[num_t].get_name() == name_d4) {
	    ++nums_rdm4.second;
            has_d4 = true;
	  } // End if
	} // End num_t
        if(has_d4) d4s.push_back(*t);
      } // End t

      // Take terms with 4-RDM 
      if(d4s.size()){
        d4_terms.insert(std::map<std::string, std::vector<Femto::Core::SQterm> >::value_type(names[num_l] + "/" + names[num_r], d4s));
	  //d4_terms[names[num_l] + "/" + names[num_r]] = d4s;
      } // End if

      ostringstream stm2;
      //*OLD* total_num_terms += count;
      //*OLD* stm2 << count;
      //*OLD* status += stm2.str() + " terms are generated .... \n\n";

#ifdef _CODE_GENE    
      // ********* Construct Sigma *************
      string thisname("sigma_" + names[num_l] + "_" + names[num_r]);
      Femto::Reaktor::SQreaktor gen(Sig, combined_result2, thisname, true, "V2", "T2");
      if(names[num_l] == "ccoo" && names[num_r] == "cooo") gen.set_oldgemm(true);
      gen.set_use_cumulant(true);
      //*OLD* if(names[num_l] != names[num_r]) //*This is OK, not the test*
      //*OLD*   gen.set_V_D(false);
      gen.set_V_D(true); // <-- NEW 2013/01/08
      gen.set_D4C(true);
      if(names[num_l] == "ooov" && names[num_r] == "cooo") gen.set_T_D(true);
      if(names[num_l].find("c") != string::npos || names[num_r].find("c") != string::npos) gen.set_guard_core(true);
      if(names[num_l].find("o") != string::npos || names[num_r].find("o") != string::npos) gen.set_guard_act(true);
      gen.init_header("mrci_header");
      gen.generate(Femto::Reaktor::Orz, Femto::Reaktor::Factorize);

      // Pick up terms with 4-RDMs, which are not transformed to Cn tensor (n = 4,5,6) 
      {
	std::vector<Femto::Core::SQterm> term_d4(gen.get_terms_d4());
	// Take terms with 4-RDM 
	if(term_d4.size()){
	  d4_terms_truely.insert(std::map<std::string, std::vector<Femto::Core::SQterm> >::value_type(names[num_l] + "/" + names[num_r], term_d4));
	} // End if        
      } // End scope

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
#endif
    } // End Li
  } // End Ri

}
