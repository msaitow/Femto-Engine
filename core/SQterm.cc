//
//  SQterm.cc
//  
//
//  Created by Masaaki Saitow on 12/06/30.
//  Copyright (c) 2012 Masaaki Saitow. All rights reserved.
//

#include <cmath>
#include <vector>
#include <algorithm>
#include <map>
#include <boost/format.hpp>
#include <Femto.hpp>
#include <SQterm.hpp>
#include <SQtensor.hpp>
#include <SQindex.hpp>

//#define _DEBUG_TERM

//#define _SORT
//#define _NEW_FACT

using namespace std;

namespace Femto { namespace Core { 

  // *********************************************************
  // 
  // *********************************************************
  SQterm::SQterm()
  {}

  // *********************************************************
  // 
  // *********************************************************
  SQterm::~SQterm()
  {}

  // *********************************************************
  // 
  // *********************************************************
  SQterm::SQterm(const double numConst, const vector<string> Consts, 
                 const vector<SQtensor> Tensors, const bool isInCanonical)
    : numConst_(numConst),
      Consts_(Consts),
      isInCanonical_(isInCanonical)
  {
    int num_kDeltas = 0;
    vector<SQtensor> Tnon_commutes;
    for(size_t I = 0;I < Tensors.size();++I){
      if(Tensors[I].isCommutable()){
        // Kronecker's delta with same indices aren't needed (2012/10/22)
	if(Tensors[I].get_name() == kDelta_name()){
	  if(find(Tensors_.begin(), Tensors_.end(), Tensors[I]) == Tensors_.end()) 
	    Tensors_.push_back(Tensors[I]);
          else ++num_kDeltas;
	}
	else Tensors_.push_back(Tensors[I]);
      }
      else   Tnon_commutes.push_back(Tensors[I]);
    }
    if(Tensors_.size()+Tnon_commutes.size()+num_kDeltas != Tensors.size()){
      cout << "Algorithmic Error .... " << endl;
      abort();
    }
    sort(Consts_.begin(), Consts_.end());
    sort(Tensors_.begin(), Tensors_.end());

    Tensors_.insert(Tensors_.end(), Tnon_commutes.begin(), Tnon_commutes.end());

    set_summedBody();

    // Count number of sfGen
    int count = 0;
    for(size_t I = 0;I < Tensors_.size();++I){
      if(is_sfGen(Tensors_[I].get_name())) ++count;
    } // End I
    // Already normal ordered
    if(count == 0 || count == 1) set_isInCanonical(true);

#ifdef _DEBUG_TERM
    cout << "Content of the Tensors_ .... " << endl;               //*TEST* 
    for(size_t i = 0;i < Tensors_.size();++i) cout << Tensors_[i]; //*TEST*
    cout << endl << endl;                                          //*TEST*
#endif

  }

  // *********************************************************
  // 
  // *********************************************************
  SQterm::SQterm(const SQterm &obj)
    : isInCanonical_(obj.isInCanonical_),
      numConst_(obj.numConst_),
      Consts_(obj.Consts_),
      Tensors_(obj.Tensors_)
  { set_summedBody(); }

  // *********************************************************
  // 
  // *********************************************************
  void SQterm::set_summedBody()
  {
    // Firstly, assign all the &summedIndices_[:] to &tempIndices[:]
    vector<SQindex> tempIndices(summedIndices_);
    for(size_t i = 0;i < tempIndices.size();++i){
      for(size_t I = 0;I < Tensors_.size();++I){
        vector<SQindex*> temp1(Tensors_[I].get_indices());
        for(size_t j = 0;j < temp1.size();++j){
          if(*temp1[j]==tempIndices[i]) Tensors_[I].put_indices(j, &tempIndices[i]);
	} // End j
      } // End I
    } // End i
    summedIndices_.clear();

    // Search for the dummy indices .... 
    for(size_t I = 0;I < Tensors_.size();++I){
      vector<SQindex*> temp1(Tensors_[I].get_indices());
      vector<SQindex> indices;
      for(size_t i = 0;i < temp1.size();++i) indices.push_back(*temp1[i]);
      for(size_t j = 0;j < indices.size();++j){
        if(find(summedIndices_.begin(), summedIndices_.end(), indices[j])==summedIndices_.end()) {
          summedIndices_.push_back(*temp1[j]);
        } // End if
      } // End j
    } // End I
    
    // Replace all the dummy indices, which are shared among all the terms, with those
    // copied in the summedIndices
    for(size_t i = 0;i < summedIndices_.size();++i){
      SQindex temp_i(summedIndices_[i]);
      for(size_t I = 0;I < Tensors_.size();++I){
        vector<SQindex*> temp1(Tensors_[I].get_indices());
        vector<SQindex> indices;        
        for(size_t j = 0;j < temp1.size();++j) indices.push_back(*temp1[j]);
        for(size_t j = 0;j < indices.size();++j){
          if(indices[j] == temp_i) Tensors_[I].put_indices(j, &summedIndices_[i]);
	}
      } // End I
    } // End i

//*OLD*      // Set the dummy names for all the dummy indices
//*OLD*      int Ccount = 0;
//*OLD*      int Ocount = 0;
//*OLD*      int Vcount = 0;
//*OLD*      for(size_t i = 0;i < summedIndices_.size();++i){
//*OLD*        // In case of core
//*OLD*        if(summedIndices_[i].get_char() == 0 && summedIndices_[i].get_isSummed()){
//*OLD*          ++Ccount;
//*OLD*          ostringstream stm;
//*OLD*          stm << Ccount;
//*OLD*          summedIndices_[i].put_index("c" + stm.str());
//*OLD*        }
//*OLD*        else
//*OLD*        // In case of active
//*OLD*  	if(summedIndices_[i].get_char() == 1 && summedIndices_[i].get_isSummed()){
//*OLD*          ++Ocount;
//*OLD*          ostringstream stm;
//*OLD*          stm << Ocount;
//*OLD*          summedIndices_[i].put_index("o" + stm.str());
//*OLD*        }
//*OLD*        else
//*OLD*        // In case of virtual
//*OLD*  	if(summedIndices_[i].get_char() == 2 && summedIndices_[i].get_isSummed()){
//*OLD*          ++Vcount;
//*OLD*          ostringstream stm;
//*OLD*          stm << Vcount;
//*OLD*          summedIndices_[i].put_index("v" + stm.str());
//*OLD*        }
//*OLD*        
//*OLD*      } // End i    

  }


  // *********************************************************
  // 
  // *********************************************************
  void SQterm::masquerade()
  {
     // Set the dummy names for all the dummy indices
     int Ccount(0);
     int Ocount(0);
     int Vcount(0);
     for(size_t i = 0;i < summedIndices_.size();++i){
       // In case of core
       if(summedIndices_[i].get_char() == 0 && summedIndices_[i].get_isSummed()){
         ++Ccount;
         ostringstream stm;
         stm << Ccount;
         summedIndices_[i].put_index("c" + stm.str());
       }
       else
       // In case of active
 	if(summedIndices_[i].get_char() == 1 && summedIndices_[i].get_isSummed()){
         ++Ocount;
         ostringstream stm;
         stm << Ocount;
         summedIndices_[i].put_index("o" + stm.str());
       }
       else
       // In case of virtual
 	if(summedIndices_[i].get_char() == 2 && summedIndices_[i].get_isSummed()){
         ++Vcount;
         ostringstream stm;
         stm << Vcount;
         summedIndices_[i].put_index("v" + stm.str());
       }     
     } // End i
     for(auto t = Tensors_.begin();t != Tensors_.end();++t){
       t->sortIndices();
     } // End t    

  }


//*OLD*   // *********************************************************
//*OLD*   // Mostly perfect, but with Factorize .... 
//*OLD*   // *********************************************************
//*OLD*   void SQterm::set_summedBody()
//*OLD*   {
//*OLD*     // Firstly, assign all the &summedIndices_[:] to &tempIndices[:]
//*OLD*     vector<SQindex> tempIndices(summedIndices_);
//*OLD*     for(size_t i = 0;i < tempIndices.size();++i){
//*OLD*       for(size_t I = 0;I < Tensors_.size();++I){
//*OLD*         vector<SQindex*> temp1 = Tensors_[I].get_indices();
//*OLD*         for(size_t j = 0;j < temp1.size();++j){
//*OLD*           if(*temp1[j]==tempIndices[i]) Tensors_[I].put_indices(j, &tempIndices[i]);
//*OLD* 	} // End j
//*OLD*       } // End I
//*OLD*     } // End i
//*OLD*     summedIndices_.clear();
//*OLD* 
//*OLD*     // Search for the dummy indices .... 
//*OLD*     for(size_t I = 0;I < Tensors_.size();++I){
//*OLD*       vector<SQindex*> temp1 = Tensors_[I].get_indices();
//*OLD*       vector<SQindex> indices;
//*OLD*       for(size_t i = 0;i < temp1.size();++i) indices.push_back(*temp1[i]);
//*OLD*       for(size_t j = 0;j < indices.size();++j){
//*OLD*         if(find(summedIndices_.begin(), summedIndices_.end(), indices[j])==summedIndices_.end()) {
//*OLD*           summedIndices_.push_back(*temp1[j]);
//*OLD*         } // End if
//*OLD*       } // End j
//*OLD*     } // End I
//*OLD*     
//*OLD*     // Replace all the dummy indices, which are shared among all the terms, with those
//*OLD*     // copied in the summedIndices
//*OLD*     for(size_t i = 0;i < summedIndices_.size();++i){
//*OLD*       SQindex temp_i = summedIndices_[i];
//*OLD*       for(size_t I = 0;I < Tensors_.size();++I){
//*OLD*         vector<SQindex*> temp1 = Tensors_[I].get_indices();
//*OLD*         vector<SQindex> indices;        
//*OLD*         for(size_t j = 0;j < temp1.size();++j) indices.push_back(*temp1[j]);
//*OLD*         for(size_t j = 0;j < indices.size();++j){
//*OLD*           if(indices[j] == temp_i) Tensors_[I].put_indices(j, &summedIndices_[i]);
//*OLD* 	}
//*OLD*       } // End I
//*OLD*     } // End i
//*OLD* 
//*OLD*      // Set the dummy names for all the dummy indices
//*OLD*      int Ccount = 0;
//*OLD*      int Ocount = 0;
//*OLD*      int Vcount = 0;
//*OLD*      for(size_t i = 0;i < summedIndices_.size();++i){
//*OLD*        // In case of core
//*OLD*        if(summedIndices_[i].get_char() == 0 && summedIndices_[i].get_isSummed()){
//*OLD*          ++Ccount;
//*OLD*          ostringstream stm;
//*OLD*          stm << Ccount;
//*OLD*          summedIndices_[i].put_index("c" + stm.str());
//*OLD*        }
//*OLD*        else
//*OLD*        // In case of active
//*OLD*  	if(summedIndices_[i].get_char() == 1 && summedIndices_[i].get_isSummed()){
//*OLD*          ++Ocount;
//*OLD*          ostringstream stm;
//*OLD*          stm << Ocount;
//*OLD*          summedIndices_[i].put_index("o" + stm.str());
//*OLD*        }
//*OLD*        else
//*OLD*        // In case of virtual
//*OLD*  	if(summedIndices_[i].get_char() == 2 && summedIndices_[i].get_isSummed()){
//*OLD*          ++Vcount;
//*OLD*          ostringstream stm;
//*OLD*          stm << Vcount;
//*OLD*          summedIndices_[i].put_index("v" + stm.str());
//*OLD*        }
//*OLD*        
//*OLD*      } // End i    
//*OLD* 
//*OLD*   }


//*OLD*   // *********************************************************
//*OLD*   // *THE SECOND OLDEST* Wick expansion works with this .....
//*OLD*   // *********************************************************
//*OLD*   void SQterm::set_summedBody()
//*OLD*   {
//*OLD* 
//*OLD*     // Firstly, assign all the &summedIndices_[:] to &tempIndices[:]
//*OLD*     vector<SQindex> tempIndices(summedIndices_);
//*OLD*     for(size_t i = 0;i < tempIndices.size();++i){
//*OLD*       for(size_t I = 0;I < Tensors_.size();++I){
//*OLD*         vector<SQindex*> temp1 = Tensors_[I].get_indices();
//*OLD*         for(size_t j = 0;j < temp1.size();++j){
//*OLD*           if(*temp1[j]==tempIndices[i]) Tensors_[I].put_indices(j, &tempIndices[i]);
//*OLD* 	} // End j
//*OLD*       } // End I
//*OLD*     } // End i
//*OLD* 
//*OLD*     summedIndices_.clear();
//*OLD* 
//*OLD*     if(DEBUG_){    
//*OLD*     // Print *TEST* pointers of the associated indices in all the tensors 
//*OLD*     cout << endl;
//*OLD*     cout << "Before .... " << endl;
//*OLD*     for(size_t I = 0;I < Tensors_.size();++I){
//*OLD*       cout << "Current Tensor ... " << Tensors_[I] << endl;
//*OLD*       vector<SQindex*> temp1 = Tensors_[I].get_indices();
//*OLD*       for(size_t i = 0;i < temp1.size();++i){
//*OLD*         cout << i << " " << temp1[i] << ", " << *(temp1[i]) << endl;
//*OLD*       } // End i
//*OLD*       cout << endl;
//*OLD*     } // End I
//*OLD*     cout << endl;
//*OLD*     }
//*OLD* 
//*OLD*     // Search for the dummy indices .... 
//*OLD*     for(size_t I = 0;I < Tensors_.size();++I){
//*OLD*       vector<SQindex*> temp1 = Tensors_[I].get_indices();
//*OLD*       vector<SQindex> indices;
//*OLD*       for(size_t i = 0;i < temp1.size();++i) indices.push_back(*temp1[i]);
//*OLD*       for(size_t j = 0;j < indices.size();++j){
//*OLD*         if(indices[j].get_isSummed() && find(summedIndices_.begin(), summedIndices_.end(), indices[j])==summedIndices_.end()) {
//*OLD*           summedIndices_.push_back(*temp1[j]);
//*OLD*         } // End if
//*OLD*       } // End j
//*OLD*     } // End I
//*OLD*     
//*OLD*     // Replace all the dummy indices, which are shared among all the terms, with those
//*OLD*     // copied in the summedIndices
//*OLD*     for(size_t i = 0;i < summedIndices_.size();++i){
//*OLD*       SQindex temp_i = summedIndices_[i];
//*OLD*       for(size_t I = 0;I < Tensors_.size();++I){
//*OLD*         vector<SQindex*> temp1 = Tensors_[I].get_indices();
//*OLD*         vector<SQindex> indices;        
//*OLD*         for(size_t j = 0;j < temp1.size();++j) indices.push_back(*temp1[j]);
//*OLD*         for(size_t j = 0;j < indices.size();++j){
//*OLD*           if(indices[j] == temp_i) Tensors_[I].put_indices(j, &summedIndices_[i]);
//*OLD* 	}
//*OLD*       } // End I
//*OLD*     } // End i
//*OLD* 
//*OLD* //*TEST*     cout << endl;
//*OLD* //*TEST*     for(size_t i = 0;i < summedIndices_.size();++i) cout << summedIndices_[i] << endl;
//*OLD* 
//*OLD*     if(DEBUG_){
//*OLD*     // Print *TEST* pointers of the associated indices in all the tensors 
//*OLD*     cout << "After .... " << endl;
//*OLD*     for(size_t I = 0;I < Tensors_.size();++I){
//*OLD*       cout << "Current Tensor ... " << Tensors_[I] << endl;
//*OLD*       vector<SQindex*> temp1 = Tensors_[I].get_indices();
//*OLD*       for(size_t i = 0;i < temp1.size();++i){
//*OLD*         cout << i << " " << temp1[i] << ", " << *(temp1[i]) << endl;
//*OLD*       } // End i
//*OLD*       cout << endl;
//*OLD*     } // End I
//*OLD*     cout << endl;
//*OLD* 
//*OLD*     // Print *TEST* pointers of indices in summedIndices_
//*OLD*     cout << endl;
//*OLD*     cout << "summedIndices_ .... " << endl;
//*OLD*     for(size_t I = 0;I < summedIndices_.size();++I) 
//*OLD*       cout << I << " " << &(summedIndices_[I]) << ", " << summedIndices_[I] << endl;
//*OLD*     cout << endl;
//*OLD*     }
//*OLD* 
//*OLD*     // Set the dummy names for all the dummy indices
//*OLD*     int Ccount = 0;
//*OLD*     int Ocount = 0;
//*OLD*     int Vcount = 0;
//*OLD*     for(size_t i = 0;i < summedIndices_.size();++i){
//*OLD*       // In case of core
//*OLD*       if(summedIndices_[i].get_char() == 0){
//*OLD*         ++Ccount;
//*OLD*         ostringstream stm;
//*OLD*         stm << Ccount;
//*OLD*         summedIndices_[i].put_index("c" + stm.str());
//*OLD*       }
//*OLD*       else
//*OLD*       // In case of active
//*OLD*       if(summedIndices_[i].get_char() == 1){
//*OLD*         ++Ocount;
//*OLD*         ostringstream stm;
//*OLD*         stm << Ocount;
//*OLD*         summedIndices_[i].put_index("o" + stm.str());
//*OLD*       }
//*OLD*       else
//*OLD*       // In case of virtual
//*OLD*       if(summedIndices_[i].get_char() == 2){
//*OLD*         ++Vcount;
//*OLD*         ostringstream stm;
//*OLD*         stm << Vcount;
//*OLD*         summedIndices_[i].put_index("v" + stm.str());
//*OLD*       }
//*OLD*       
//*OLD*     } // End i    
//*OLD* 
//*OLD*   }


//*OLD*
//*OLD*   // *********************************************************
//*OLD*   // *THE OLDEST VERSION* May be bugged ..... 
//*OLD*   // *********************************************************
//*OLD*   void SQterm::set_summedBody()
//*OLD*   {
//*OLD*     summedIndices_.clear();
//*OLD* 
//*OLD*     if(DEBUG_){    
//*OLD*     // Print *TEST* pointers of the associated indices in all the tensors 
//*OLD*     cout << endl;
//*OLD*     cout << "Before .... " << endl;
//*OLD*     for(size_t I = 0;I < Tensors_.size();++I){
//*OLD*       cout << "Current Tensor ... " << Tensors_[I] << endl;
//*OLD*       vector<SQindex*> temp1 = Tensors_[I].get_indices();
//*OLD*       for(size_t i = 0;i < temp1.size();++i){
//*OLD*         cout << i << " " << temp1[i] << ", " << *(temp1[i]) << endl;
//*OLD*       } // End i
//*OLD*       cout << endl;
//*OLD*     } // End I
//*OLD*     cout << endl;
//*OLD*     }
//*OLD* 
//*OLD*     // Search for the dummy indices .... 
//*OLD*     for(size_t I = 0;I < Tensors_.size();++I){
//*OLD*       vector<SQindex*> temp1 = Tensors_[I].get_indices();
//*OLD*       vector<SQindex> indices;
//*OLD*       for(size_t i = 0;i < temp1.size();++i) indices.push_back(*temp1[i]);
//*OLD*       for(size_t j = 0;j < indices.size();++j){
//*OLD*         if(indices[j].get_isSummed() && find(summedIndices_.begin(), summedIndices_.end(), indices[j])==summedIndices_.end()) {
//*OLD*           summedIndices_.push_back(*temp1[j]);
//*OLD*         } // End if
//*OLD*       } // End j
//*OLD*     } // End I
//*OLD*     
//*OLD*     // Replace all the dummy indices, which are shared among all the terms, with those
//*OLD*     // copied in the summedIndices
//*OLD*     for(size_t I = 0;I < Tensors_.size();++I){
//*OLD*       vector<SQindex*> temp1 = Tensors_[I].get_indices();
//*OLD*       vector<SQindex> indices;
//*OLD*       for(size_t i = 0;i < temp1.size();++i) indices.push_back(*temp1[i]);
//*OLD*       for(size_t j = 0;j < indices.size();++j){
//*OLD*         if(indices[j].get_isSummed()) {
//*OLD*           
//*OLD* 	} // End if        
//*OLD*       } // End j
//*OLD*     } //End I
//*OLD* 
//*OLD*     for(size_t i = 0;i < summedIndices_.size();++i){
//*OLD*       SQindex temp_i = summedIndices_[i];
//*OLD*       for(size_t I = 0;I < Tensors_.size();++I){
//*OLD*         vector<SQindex*> temp1 = Tensors_[I].get_indices();
//*OLD*         vector<SQindex> indices;        
//*OLD*         for(size_t j = 0;j < temp1.size();++j) indices.push_back(*temp1[j]);
//*OLD*         for(size_t j = 0;j < indices.size();++j){
//*OLD*           if(indices[j] == temp_i) Tensors_[I].put_indices(j, &summedIndices_[i]);
//*OLD* 	}
//*OLD*       } // End I
//*OLD*     } // End i
//*OLD* 
//*OLD* //*TEST*     cout << endl;
//*OLD* //*TEST*     for(size_t i = 0;i < summedIndices_.size();++i) cout << summedIndices_[i] << endl;
//*OLD* 
//*OLD*     if(DEBUG_){
//*OLD*     // Print *TEST* pointers of the associated indices in all the tensors 
//*OLD*     cout << "After .... " << endl;
//*OLD*     for(size_t I = 0;I < Tensors_.size();++I){
//*OLD*       cout << "Current Tensor ... " << Tensors_[I] << endl;
//*OLD*       vector<SQindex*> temp1 = Tensors_[I].get_indices();
//*OLD*       for(size_t i = 0;i < temp1.size();++i){
//*OLD*         cout << i << " " << temp1[i] << ", " << *(temp1[i]) << endl;
//*OLD*       } // End i
//*OLD*       cout << endl;
//*OLD*     } // End I
//*OLD*     cout << endl;
//*OLD* 
//*OLD*     // Print *TEST* pointers of indices in summedIndices_
//*OLD*     cout << endl;
//*OLD*     cout << "summedIndices_ .... " << endl;
//*OLD*     for(size_t I = 0;I < summedIndices_.size();++I) 
//*OLD*       cout << I << " " << &(summedIndices_[I]) << ", " << summedIndices_[I] << endl;
//*OLD*     cout << endl;
//*OLD*     }
//*OLD* 
//*OLD*     // Set the dummy names for all the dummy indices
//*OLD*     int Ccount = 0;
//*OLD*     int Ocount = 0;
//*OLD*     int Vcount = 0;
//*OLD*     for(size_t i = 0;i < summedIndices_.size();++i){
//*OLD*       // In case of core
//*OLD*       if(summedIndices_[i].get_char() == 0){
//*OLD*         ++Ccount;
//*OLD*         ostringstream stm;
//*OLD*         stm << Ccount;
//*OLD*         summedIndices_[i].put_index("c" + stm.str());
//*OLD*       }
//*OLD*       else
//*OLD*       // In case of active
//*OLD*       if(summedIndices_[i].get_char() == 1){
//*OLD*         ++Ocount;
//*OLD*         ostringstream stm;
//*OLD*         stm << Ocount;
//*OLD*         summedIndices_[i].put_index("o" + stm.str());
//*OLD*       }
//*OLD*       else
//*OLD*       // In case of virtual
//*OLD*       if(summedIndices_[i].get_char() == 2){
//*OLD*         ++Vcount;
//*OLD*         ostringstream stm;
//*OLD*         stm << Vcount;
//*OLD*         summedIndices_[i].put_index("v" + stm.str());
//*OLD*       }
//*OLD*       
//*OLD*     } // End i    
//*OLD* 
//*OLD*   }

//*INCOMPLETE*   // *********************************************************
//*INCOMPLETE*   // 
//*INCOMPLETE*   // *********************************************************
//*INCOMPLETE*   void SQterm::erase_duplication()
//*INCOMPLETE*   {
//*INCOMPLETE*     vector<int> count;   count.resize((int)summedIndices_.size());
//*INCOMPLETE*     vector<int> location;location.resize((int)summedIndices_.size());
//*INCOMPLETE*     for(size_t i = 0;i < summedIndices_.size()-1;++i){
//*INCOMPLETE*       for(size_t j = i+1;j < summedIndices_.size();++j){
//*INCOMPLETE*         if(summedBody_[i]==summedBody_[j]) { ++count[i]; location[i] = j; }
//*INCOMPLETE*       } // End j    
//*INCOMPLETE*     } // End i
//*INCOMPLETE*     for(size_t i = 0;i < summedIndices.size();++i){
//*INCOMPLETE*       if(count[i]){
//*INCOMPLETE*         for(vector<SQtensor>::iterator t = Tensors_.begin();t != Tensors.end();++t){
//*INCOMPLETE* 	  if(*(t->get_indices()) == summedBody_[i]) t->put_indices(i, &summedBody[i]);
//*INCOMPLETE* 	} // End t
//*INCOMPLETE*         summedIndices
//*INCOMPLETE*       } // End if
//*INCOMPLETE*     } // End for  
//*INCOMPLETE*   }

  // *********************************************************
  // 
  // *********************************************************
  vector<SQtensor> SQterm::get_tensors() const
  { return Tensors_; }

  // *********************************************************
  // 
  // *********************************************************
  vector<SQtensor*> SQterm::get_tensors_ptr()
  {
    vector<SQtensor*> retval;
    for(size_t i = 0;i < Tensors_.size();++i) retval.push_back(&(Tensors_[i])); 
    return retval; 
  }

  // *********************************************************
  // 
  // *********************************************************
  double SQterm::get_numConst() const
  { return numConst_; }

  // *********************************************************
  // 
  // *********************************************************
  vector<string> SQterm::get_Consts() const
  { return Consts_; }

  // *********************************************************
  // 
  // *********************************************************
  bool SQterm::get_isInCanonical() const
  { return isInCanonical_; }

  // *********************************************************
  // 
  // *********************************************************
  void SQterm::set_isInCanonical(const bool isInCanonical)
  { isInCanonical_ = isInCanonical; }

  // *********************************************************
  // 
  // *********************************************************
  void SQterm::set_tensors(const vector<SQtensor> &Tensors)
  {
    vector<SQtensor> Tcommutes, Tnon_commutes;
    Tcommutes.reserve(Tensors.size());
    Tnon_commutes.reserve(Tensors.size());
    for(size_t I = 0;I < Tensors.size();++I){
      if(Tensors[I].isCommutable()) Tcommutes.push_back(Tensors[I]);
      else                          Tnon_commutes.push_back(Tensors[I]);
    }
    if(Tcommutes.size()+Tnon_commutes.size() != Tensors.size()){
      cout << "Algorithmic Error .... " << endl;
      abort();
    }
    sort(Tcommutes.begin(), Tcommutes.end());

    Tensors_.clear();
    Tensors_.reserve(Tensors.size());

    // Modified 2012/10/23
    Tensors_.insert(Tensors_.begin(), Tcommutes.begin(),     Tcommutes.end());
    Tensors_.insert(Tensors_.end(),   Tnon_commutes.begin(), Tnon_commutes.end());
    
//*OLD*     for(size_t i = 0;i < Tcommutes.size();++i) Tensors_.push_back(Tcommutes[i]);
//*OLD*     for(size_t i = 0;i < Tnon_commutes.size();++i) Tensors_.push_back(Tnon_commutes[i]);

    set_summedBody();
    
  }

  // *********************************************************
  // 
  // *********************************************************
  SQterm SQterm::operator=(const SQterm &obj){
    isInCanonical_ = obj.isInCanonical_;
    numConst_      = obj.numConst_;     
    Consts_        = obj.Consts_;       
    Tensors_       = obj.Tensors_;      

    set_summedBody();

    return (*this);
  }

//*SLEEP_FOR_AWHILE*   // *********************************************************
//*SLEEP_FOR_AWHILE*   // 
//*SLEEP_FOR_AWHILE*   // *********************************************************
//*SLEEP_FOR_AWHILE*   SQterm SQterm::operator+(const SQterm &obj){
//*SLEEP_FOR_AWHILE*     if(!isAdditive(*this, obj)){
//*SLEEP_FOR_AWHILE*       cout << "Terms are not additive." << endl;
//*SLEEP_FOR_AWHILE*       abort();
//*SLEEP_FOR_AWHILE*     }
//*SLEEP_FOR_AWHILE* 
//*SLEEP_FOR_AWHILE*     //    SQterm retval(*this);
//*SLEEP_FOR_AWHILE*     numConst_ += obj.numConst_;
//*SLEEP_FOR_AWHILE*     //    cout << this->numConst_ + obj.numConst_ << endl; //*TEST*
//*SLEEP_FOR_AWHILE* //*OLD*     retval.set_numConst(this->numConst_ + obj.numConst_);
//*SLEEP_FOR_AWHILE* //*OLD*     retval.set_tensors(this->Tensors_);
//*SLEEP_FOR_AWHILE* //*OLD*     retval.set_isInCanonical(this->isInCanonical_);
//*SLEEP_FOR_AWHILE* //*OLD*     retval.set_Consts(this->Consts_);
//*SLEEP_FOR_AWHILE* //*OLD*     retval.set_summedBody();
//*SLEEP_FOR_AWHILE* 
//*SLEEP_FOR_AWHILE*     return *this;
//*SLEEP_FOR_AWHILE*   }

//*OLD*   // *********************************************************
//*OLD*   // 
//*OLD*   // *********************************************************
//*OLD*   SQterm SQterm::operator+(const SQterm &obj){
//*OLD*     if(!isAdditive(*this, obj)){
//*OLD*       cout << "Terms are not additive." << endl;
//*OLD*       abort();
//*OLD*     }
//*OLD* 
//*OLD*     SQterm retval(*this);
//*OLD*     retval.set_numConst(this->numConst_ + obj.numConst_);
//*OLD*     //    cout << this->numConst_ + obj.numConst_ << endl; //*TEST*
//*OLD* //*OLD*     retval.set_numConst(this->numConst_ + obj.numConst_);
//*OLD* //*OLD*     retval.set_tensors(this->Tensors_);
//*OLD* //*OLD*     retval.set_isInCanonical(this->isInCanonical_);
//*OLD* //*OLD*     retval.set_Consts(this->Consts_);
//*OLD* //*OLD*     retval.set_summedBody();
//*OLD* 
//*OLD*     return retval;
//*OLD*   }

  // *********************************************************
  // 
  // *********************************************************
  bool SQterm::operator==(const SQterm &obj)
  {
    if(numConst_ != obj.numConst_) return false;
    if(Consts_   != obj.Consts_)   return false;
    if(Tensors_  != obj.Tensors_)  return false;

    return true;
  }

  // *********************************************************
  // 
  // *********************************************************
  void SQterm::set_numConst(const double num)
  { numConst_ = num; }

  // *********************************************************
  // 
  // *********************************************************
  void SQterm::set_Consts(const vector<string> consts)
  { Consts_ = consts; }

  // *********************************************************
  // 
  // *********************************************************
  ostream& operator <<(std::ostream &os, const SQterm &t)
  {
    os << boost::format("(%14.8f)") % t.get_numConst() << " ";

    vector<string> strs(t.get_Consts());
    for(size_t i = 0;i < strs.size();++i) os << strs[i] << " ";

    vector<SQtensor> tensors(t.get_tensors());    
    for(size_t i = 0;i < tensors.size();++i) os << tensors[i];
    return os;
  }

  // *********************************************************
  // 
  // *********************************************************
  vector<SQindex*> SQterm::get_summedBody()
  {
    vector<SQindex*> retval;
    auto ind = summedIndices_.begin();
    for(;ind != summedIndices_.end();++ind) retval.push_back(&(*ind));
    return retval;
  }

  // *********************************************************
  // If the term is in canonical order, transform sfGen to RDM 
  // *********************************************************
  void SQterm::transform2RDM()
  {
    //set_summedBody(); //*TEST* 
    auto ten = Tensors_.begin();
    if(!get_isInCanonical()) return;
    bool set_flag = false; 
    for(;ten != Tensors_.end();){
      if(is_sfGen(ten->get_name())){
        RDM temp(ten->get_indices());
        Tensors_.erase(ten);
        Tensors_.push_back(temp);
        set_flag = true;        
      } // End if
      else ++ten;
    } // End ten
    if(set_flag) {
      vector<SQtensor> commT;
      vector<SQtensor> ncommT;
      for(size_t i = 0;i < Tensors_.size();++i)
        if(Tensors_[i].isCommutable())  commT.push_back(Tensors_[i]);
	else                           ncommT.push_back(Tensors_[i]);
 
      sort(commT.begin(), commT.end());
      commT.insert(commT.end(), ncommT.begin(), ncommT.end());
      Tensors_ = commT;
      //set_summedBody(); //*TEST*
    }
  }

  // *********************************************************
  // Contract Kronecker deltas
  // *********************************************************
  void SQterm::contractkDeltas()
  {
    auto t = Tensors_.begin();
    for(; t != Tensors_.end();){
      // In case of delta_{p}^{p} ....
      if((t->get_name()==kDelta_name()) && (*(t->get_indices()[0]) == *(t->get_indices()[1])))
        t = Tensors_.erase(t);
      // In case of both indices of kDelta are dummy ....
      else if((t->get_name()==kDelta_name()) && !(t->get_indices()[0]->get_isSummed()) && !(t->get_indices()[1]->get_isSummed())){
        if(t->get_indices()[0]->get_char() != t->get_indices()[1]->get_char()){
	  numConst_ = 0; // This term shoud be zero
	  return;
	} // End if
	++t;
      } // End else if
      // In case of either of indices in kDelta is dummy .... 
      else if((t->get_name()==kDelta_name()) && (t->get_indices()[0]->get_isSummed() || t->get_indices()[1]->get_isSummed())){
        if(t->get_indices()[0]->get_char() != t->get_indices()[1]->get_char()){
          numConst_ = 0;
	  return;
	} // End if
//*TEST* 	++t; //*TEST*
	else{
	  //           set_summedBody(); //*TEST*
          SQindex *killed;
          SQindex *killer;
          if(t->get_indices()[1]->get_isSummed()) { killed = t->get_indices()[1]; killer = t->get_indices()[0];}
          else                                    { killed = t->get_indices()[0]; killer = t->get_indices()[1];}
 	  // Kill the index to be killed 
          auto killed_ptr = find(summedIndices_.begin(), summedIndices_.end(), *killed);
          if(killed_ptr == summedIndices_.end()){
            cout << "Algorithmic Error ..... " << endl;
            abort();
	  }
	  else{
            *killed_ptr = *killer; 
	  }
 	  t = Tensors_.erase(t);
//*TEST*           cout << endl;                         //*TEST*
//*TEST*           cout << "Killed : " << *killed << " Kill : " << *kill << endl; //*TEST*
//*TEST*           print_summedBody(); cout << endl;//*TEST*
//*TEST*           ++t; //*TEST* 
	} // End else
      } // End else if
      else ++t;

    } // End for
    //set_summedBody(); //*SLOW* ???
    // This should be improved like  
  }

  // *********************************************************
  // Decompose sfGen, or RDM (only if isInCanonical = true)
  // *********************************************************
  void SQterm::decomposeRDMVirt()
  {
    if(!isInCanonical_) return;
    auto t = Tensors_.begin();
    for(;t != Tensors_.end();++t){
      if(is_sfGen(t->get_name()) || is_RDM(t->get_name())){
        vector<SQindex*> indices(t->get_indices());
        auto i = indices.begin();
        for(;i != indices.end();++i){
          // In case of core .... 
          if((*i)->get_char() == (char_state)0) {
	    /* Do for core indices */
	  } // End if
          // In case of virtual ....
          else if((*i)->get_char() == (char_state)2) {
            numConst_ = 0;
	    return;
	  } // End else if
	} // End for
      } // End if
    } // End for
  }

  // *********************************************************
  // 
  // *********************************************************
  void SQterm::print_summedBody()
  {
    cout << ">> summedBody <<" << endl;
    auto i = summedIndices_.begin();
    int count = 0;
    for(;i != summedIndices_.end();++i)
      cout << (boost::format("[%10d] ") % count++) << *i << "(" << &(*i) << ") " <<  (i->get_isExt() ? "E" : "I") << endl; 
  }

  // *********************************************************
  // 
  // *********************************************************
  string SQterm::convert2LaTeX() const
  {
    string term("");
    if     (numConst_ == 0.0) return "";

    ostringstream stm;
    stm << fabs(numConst_);
    term += (numConst_ > 0.0 ? "+ " : "- ");
    if(fabs(numConst_) != 1.0) term += stm.str() + " ";
    for(size_t i = 0;i < Consts_.size();++i)
      if(Consts_[i] != "") term += Consts_[i] + " ";

    for(size_t i = 0;i < Tensors_.size();++i)
      term += Tensors_[i].convert2LaTeX();
    return term;
  }

  // *********************************************************
  // Return whether two SQterms are additive or not
  // *********************************************************
  // *NOTE* :: If there are common tensor type in terms, the current algorithm 
  //           may not work perfectly (2012/10/22)
  bool isAdditive(SQterm &a, SQterm &b)
  {
    // Compare constants
    vector<string> Const_a(a.get_Consts());
    vector<string> Const_b(b.get_Consts());
    if(Const_a != Const_b) return false;

#ifndef _NEW_FACT
    // Modified 2012/10/24
    return isFactorizable(a, b);
#else
    return isFactorizable2(a, b);
#endif

//* --- *     // Compare names of tensors and associated indices
//* --- *     vector<SQtensor> ten_a(a.get_tensors());
//* --- *     vector<SQtensor> ten_b(b.get_tensors());
//* --- *     if(ten_a.size() != ten_b.size()) return false;
//* --- *     for(size_t i = 0;i < ten_a.size();++i){
//* --- *       if(ten_a[i].get_name()           != ten_b[i].get_name())           
//* --- *         return false;
//* --- *       if(ten_a[i].get_indices().size() != ten_b[i].get_indices().size())
//* --- *         return false;
//* --- *     } // End i
//* --- * 
//* --- *     // Count all the dummies in each orbital group
//* --- *     vector<SQindex*> a_indices(a.get_summedBody());
//* --- *     vector<SQindex*> c_ptr;
//* --- *     vector<SQindex*> o_ptr;
//* --- *     vector<SQindex*> v_ptr;
//* --- *     for(vector<SQindex*>::iterator i = a_indices.begin();i != a_indices.end();++i){
//* --- *       if     ((*i)->get_char()==(char_state)0 && (*i)->get_isSummed()) c_ptr.push_back(*i); 
//* --- *       else if((*i)->get_char()==(char_state)1 && (*i)->get_isSummed()) o_ptr.push_back(*i); 
//* --- *       else if((*i)->get_char()==(char_state)2 && (*i)->get_isSummed()) v_ptr.push_back(*i);
//* --- *     } // End i
//* --- * 
//* --- *     // Before comparing the tensors, masquerade the other one!
//* --- *     b.masquerade();
//* --- * 
//* --- *     // Permute pairs of indices in each orbital group
//* --- *     IIvector c_perms = makePermutations((int)c_ptr.size());
//* --- *     IIvector o_perms = makePermutations((int)o_ptr.size());
//* --- *     IIvector v_perms = makePermutations((int)v_ptr.size());
//* --- *     int c_max = (c_perms.size() ? c_perms.size() : 1);
//* --- *     int o_max = (o_perms.size() ? o_perms.size() : 1);
//* --- *     int v_max = (v_perms.size() ? v_perms.size() : 1);
//* --- *     for(int ic = 0;ic < c_max;++ic){
//* --- *       for(int io = 0;io < o_max;++io){
//* --- *         for(int iv = 0;iv < v_max;++iv){
//* --- *           // Replace each name of index
//* --- *           if(c_perms.size())
//* --- *             //cout << c_perms.size() << endl;
//* --- *             for(size_t i = 0;i < c_ptr.size();++i){
//* --- *               stringstream num; num << c_perms[ic][i]+1;
//* --- *               c_ptr[i]->put_index("c"+num.str());
//* --- *             }
//* --- * 	  if(o_perms.size())
//* --- *             for(size_t i = 0;i < o_ptr.size();++i){
//* --- *               stringstream num; num << o_perms[io][i]+1;
//* --- *               o_ptr[i]->put_index("o"+num.str());
//* --- *             }
//* --- *           if(v_perms.size())
//* --- *             for(size_t i = 0;i < v_ptr.size();++i){
//* --- *               stringstream num; num << v_perms[iv][i]+1;
//* --- *               v_ptr[i]->put_index("v"+num.str()); 
//* --- *             }
//* --- * 	    if(a.get_tensors() == b.get_tensors()) return true;
//* --- *         } // End iv
//* --- *       } // End io
//* --- *     } // End ic
//* --- * 
//* --- *     return false;
  }


  // *Created (2012/10/24)*
  // *********************************************************
  // Return whether two SQterms are factorizable or not
  // *********************************************************
  // *NOTE* :: If there are common tensor type in terms, the current algorithm 
  //           may not work perfectly (2012/10/22)
  bool isFactorizable(SQterm &a, SQterm &b)
  {

    // Compare names of tensors and associated indices
    vector<SQtensor> ten_a(a.get_tensors());
    vector<SQtensor> ten_b(b.get_tensors());
    if(ten_a.size() != ten_b.size()) return false;
    for(size_t i = 0;i < ten_a.size();++i){
      if(ten_a[i].get_name()           != ten_b[i].get_name())           
        return false;
      if(ten_a[i].get_indices().size() != ten_b[i].get_indices().size())
        return false;
    } // End i

//*UNCOMPLETE*     // Consider the which tensor holds the non-dummy indices
//*UNCOMPLETE*     map<SQindex*, vector<string> > a_ind2ten, b_ind2ten;
//*UNCOMPLETE*     for(size_t a = 0;a < ten_a.size();++a){
//*UNCOMPLETE*       for(size_t i = 0;i < ten_a[a].get_indices().size();++i){
//*UNCOMPLETE*         if()
//*UNCOMPLETE*       } // End i 
//*UNCOMPLETE*     } // End a

    // Count all the dummies in each orbital group
    vector<SQindex*> a_indices(a.get_summedBody());
    vector<SQindex*> c_ptr;
    vector<SQindex*> o_ptr;
    vector<SQindex*> v_ptr;
    for(auto i = a_indices.begin();i != a_indices.end();++i){
      if     ((*i)->get_char()==Femto::core && (*i)->get_isSummed()) c_ptr.push_back(*i); 
      else if((*i)->get_char()==Femto::act  && (*i)->get_isSummed()) o_ptr.push_back(*i); 
      else if((*i)->get_char()==Femto::virt && (*i)->get_isSummed()) v_ptr.push_back(*i);
    } // End i

    // Before comparing the tensors, masquerade the other one!
    b.masquerade();

    // Permute pairs of indices in each orbital group
    IIvector c_perms = makePermutations((int)c_ptr.size());
    IIvector o_perms = makePermutations((int)o_ptr.size());
    IIvector v_perms = makePermutations((int)v_ptr.size());
    int c_max = (c_perms.size() ? c_perms.size() : 1);
    int o_max = (o_perms.size() ? o_perms.size() : 1);
    int v_max = (v_perms.size() ? v_perms.size() : 1);
    for(int ic = 0;ic < c_max;++ic){
      for(int io = 0;io < o_max;++io){
        for(int iv = 0;iv < v_max;++iv){
          // Replace each name of index
          if(c_perms.size())
            //cout << c_perms.size() << endl;
            for(size_t i = 0;i < c_ptr.size();++i){
              stringstream num; num << c_perms[ic][i]+1;
              c_ptr[i]->put_index("c"+num.str());
            }
	  if(o_perms.size())
            for(size_t i = 0;i < o_ptr.size();++i){
              stringstream num; num << o_perms[io][i]+1;
              o_ptr[i]->put_index("o"+num.str());
            }
          if(v_perms.size())
            for(size_t i = 0;i < v_ptr.size();++i){
              stringstream num; num << v_perms[iv][i]+1;
              v_ptr[i]->put_index("v"+num.str()); 
            }
#ifdef _SORT
	    vector<SQtensor> a_tensors(a.get_tensors()); sort(a_tensors.begin(), a_tensors.end()); //*TEST* 
	    //vector<SQtensor> b_tensors(b.get_tensors()); sort(b_tensors.begin(), b_tensors.end()); //*TEST*
	    if(a_tensors == b.get_tensors()) return true; //*TEST
#else
	    if(a.get_tensors() == b.get_tensors()) return true;
#endif
        } // End iv
      } // End io
    } // End ic

    return false;
  }


  // *Created (2012/11/14)*
  // *********************************************************
  // Return whether two SQterms are factorizable or not
  // *********************************************************
  // *NOTE* :: This will work perfectly, even if there are common tensors
  //           But it seems to be somewhat heavier
  bool isFactorizable2(SQterm &a, SQterm &b)
  {
    // Compare names of tensors and associated indices
    vector<SQtensor> ten_a(a.get_tensors());
    vector<SQtensor> ten_b(b.get_tensors());
    if(ten_a.size() != ten_b.size()) return false;

    //Divide the tensors by the names   
    typedef map<string, vector<SQtensor*> > ten_ptr; // name of tensor and the positions 
    ten_ptr a_ten, b_ten;
    for(size_t i = 0;i < a.get_tensors().size();++i) 
      a_ten[a.get_tensors()[i].get_name()].push_back(a.get_tensors_ptr()[i]);
    for(size_t i = 0;i < b.get_tensors().size();++i) 
      b_ten[b.get_tensors()[i].get_name()].push_back(b.get_tensors_ptr()[i]);
    for(auto i = a_ten.begin();i != a_ten.end();++i){
      string a_key(i->first);
      if(b_ten.find(a_key) == b_ten.end())        return false;
      if(b_ten[a_key].size() != i->second.size()) return false;
    } // End i

    int count(0);
    vector<IIvector> i_perms(a_ten.size());
    for(auto i = a_ten.begin();i != a_ten.end();++i,++count) i_perms[count] = makePermutations((int)i->second.size());

    // Count all the dummies in each orbital group
    vector<SQindex*> a_indices(a.get_summedBody());
    vector<SQindex*> c_ptr;
    vector<SQindex*> o_ptr;
    vector<SQindex*> v_ptr;
    for(auto i = a_indices.begin();i != a_indices.end();++i){
      if     ((*i)->get_char()==(char_state)0 && (*i)->get_isSummed()) c_ptr.push_back(*i); 
      else if((*i)->get_char()==(char_state)1 && (*i)->get_isSummed()) o_ptr.push_back(*i); 
      else if((*i)->get_char()==(char_state)2 && (*i)->get_isSummed()) v_ptr.push_back(*i);
    } // End i

    // Before comparing the tensors, masquerade the other one!
    b.masquerade();    

    // Permute pairs of indices in each orbital group
    IIvector c_perms(makePermutations((int)c_ptr.size()));
    IIvector o_perms(makePermutations((int)o_ptr.size()));
    IIvector v_perms(makePermutations((int)v_ptr.size()));
    int c_max = (c_perms.size() ? c_perms.size() : 1);
    int o_max = (o_perms.size() ? o_perms.size() : 1);
    int v_max = (v_perms.size() ? v_perms.size() : 1);
    for(int ic = 0;ic < c_max;++ic){
      for(int io = 0;io < o_max;++io){
        for(int iv = 0;iv < v_max;++iv){
          // Replace each name of index
          if(c_perms.size())
            for(size_t i = 0;i < c_ptr.size();++i){
              stringstream num; num << c_perms[ic][i]+1;
              c_ptr[i]->put_index("c"+num.str());
            }
	  if(o_perms.size())
            for(size_t i = 0;i < o_ptr.size();++i){
              stringstream num; num << o_perms[io][i]+1;
              o_ptr[i]->put_index("o"+num.str());
            }
          if(v_perms.size())
            for(size_t i = 0;i < v_ptr.size();++i){
              stringstream num; num << v_perms[iv][i]+1;
              v_ptr[i]->put_index("v"+num.str()); 
            }
	  // Compare all the members of the corresponding subgroups
          // in every possible orders
	  {  
            bool all_ok(true);
            int count(0);
            for(auto i = a_ten.begin();i != a_ten.end();++i,++count){
              bool ok;
              //IIvector i_perms(makePermutations((int)i->second.size()));
	      for(int num_p = 0;num_p < i_perms[count].size();++num_p){
                ok = true;
		for(size_t num = 0;num < i->second.size();++num){
		  if(!(*(b_ten[i->first].at(num)) == *(i->second.at(i_perms[count][num_p][num])))) ok = false;
		} // End num
                if(ok) break;
	      } // End p
              if(!ok) { all_ok = false; break; }
	    } // End i
            if(all_ok) return true;
	  } // End scope
        } // End iv
      } // End io
    } // End ic

    return false;    
  }


//*OLD*   // *********************************************************
//*OLD*   // Return whether two SQterms are additive or not
//*OLD*   // *********************************************************
//*OLD*   bool isAdditive(SQterm &a, const SQterm &b)
//*OLD*   {
//*OLD*     // Compare constants
//*OLD*     vector<string> Const_a = a.get_Consts();
//*OLD*     vector<string> Const_b = b.get_Consts();
//*OLD*     if(Const_a != Const_b) return false;
//*OLD* 
//*OLD*     // Compare names of tensors and associated indices
//*OLD*     vector<SQtensor> ten_a = a.get_tensors();
//*OLD*     vector<SQtensor> ten_b = b.get_tensors();
//*OLD*     if(ten_a.size() != ten_b.size()) return false;
//*OLD*     for(size_t i = 0;i < ten_a.size();++i){
//*OLD*       if(ten_a[i].get_name()           != ten_b[i].get_name())           
//*OLD*         return false;
//*OLD*       if(ten_a[i].get_indices().size() != ten_b[i].get_indices().size())
//*OLD*         return false;
//*OLD*     } // End i
//*OLD* 
//*OLD*     // Count all the dummies in each orbital group
//*OLD*     vector<SQindex*> a_indices = a.get_summedBody();
//*OLD*     vector<SQindex*> c_ptr;
//*OLD*     vector<SQindex*> o_ptr;
//*OLD*     vector<SQindex*> v_ptr;
//*OLD*     for(vector<SQindex*>::iterator i = a_indices.begin();i != a_indices.end();++i){
//*OLD*       if     ((*i)->get_char()==(char_state)0 && (*i)->get_isSummed()) c_ptr.push_back(*i); 
//*OLD*       else if((*i)->get_char()==(char_state)1 && (*i)->get_isSummed()) o_ptr.push_back(*i); 
//*OLD*       else if((*i)->get_char()==(char_state)2 && (*i)->get_isSummed()) v_ptr.push_back(*i);
//*OLD*     } // End i
//*OLD* 
//*OLD*     // Permute pairs of indices in each orbital group
//*OLD*     IIvector c_perms = makePermutations((int)c_ptr.size());
//*OLD*     IIvector o_perms = makePermutations((int)o_ptr.size());
//*OLD*     IIvector v_perms = makePermutations((int)v_ptr.size());
//*OLD*     int c_max = (c_perms.size() ? c_perms.size() : 1);
//*OLD*     int o_max = (o_perms.size() ? o_perms.size() : 1);
//*OLD*     int v_max = (v_perms.size() ? v_perms.size() : 1);
//*OLD*     for(int ic = 0;ic < c_max;++ic){
//*OLD*       for(int io = 0;io < o_max;++io){
//*OLD*         for(int iv = 0;iv < v_max;++iv){
//*OLD*           // Replace each name of index
//*OLD*           if(c_perms.size())
//*OLD*             //cout << c_perms.size() << endl;
//*OLD*             for(size_t i = 0;i < c_ptr.size();++i){
//*OLD*               stringstream num; num << c_perms[ic][i]+1;
//*OLD*               c_ptr[i]->put_index("c"+num.str());
//*OLD*             }
//*OLD* 	  if(o_perms.size())
//*OLD*             for(size_t i = 0;i < o_ptr.size();++i){
//*OLD*               stringstream num; num << o_perms[io][i]+1;
//*OLD*               o_ptr[i]->put_index("o"+num.str());
//*OLD*             }
//*OLD*           if(v_perms.size())
//*OLD*             for(size_t i = 0;i < v_ptr.size();++i){
//*OLD*               stringstream num; num << v_perms[iv][i]+1;
//*OLD*               v_ptr[i]->put_index("v"+num.str()); 
//*OLD*             }
//*OLD* 	    if(a.get_tensors() == b.get_tensors()) return true;
//*OLD*         } // End iv
//*OLD*       } // End io
//*OLD*     } // End ic
//*OLD* 
//*OLD*     return false;
//*OLD*   }


//*OLD*   // *********************************************************
//*OLD*   // Return whether two SQterms are additive or not
//*OLD*   // *********************************************************
//*OLD*   bool isAdditive(const SQterm &a, const SQterm &b)
//*OLD*   {
//*OLD*     vector<SQtensor> ten_a = a.get_tensors();
//*OLD*     vector<SQtensor> ten_b = b.get_tensors();
//*OLD*     if(ten_a != ten_b) return false;
//*OLD* 
//*OLD*     vector<string> Const_a = a.get_Consts();
//*OLD*     vector<string> Const_b = b.get_Consts();
//*OLD*     if(Const_a != Const_b) return false;
//*OLD* 
//*OLD*     return true;
//*OLD*   }

  // *********************************************************
  // Delete all the terms in inTerms with zero factors
  // *********************************************************
  void screenTerms(const vector<SQterm> &inTerms, vector<SQterm> *outTerms, const double crit)
  {
    outTerms->reserve(inTerms.size());
    for(size_t i = 0;i < inTerms.size();++i){
      if(fabs(inTerms[i].get_numConst()) > crit) outTerms->push_back(inTerms[i]);
    } // End i
  }

  // *********************************************************
  // Delete all the terms in inTerms with zero factors
  // *********************************************************
  vector<SQterm> screenTerms(const vector<SQterm> &inTerms, const double crit)
  {
    vector<SQterm> retval;
    for(size_t i = 0;i < inTerms.size();++i){
      if(fabs(inTerms[i].get_numConst()) > crit) retval.push_back(inTerms[i]);
    } // End for 
    return retval;
  }

//*SLOW*   vector<SQterm> termChop(vector<SQterm> &inTerms, const double crit)
//*SLOW*   {
//*SLOW*     vector<SQterm>::iterator t = inTerms.begin();
//*SLOW*     for(;t != inTerms.end();){
//*SLOW*       if(fabs(t->get_numConst()) <= crit) t = inTerms.erase(t);
//*SLOW*       else ++t; 
//*SLOW*     } // End for
//*SLOW*     return inTerms;
//*SLOW*   }

  // *********************************************************
  // Combine terms with isAdditive()==true
  // *********************************************************
  vector<SQterm> combineTerms(vector<SQterm> &inTerms)
  {
    vector<SQterm> retval;
    int count(0);
    for(auto t1 = inTerms.begin();t1 != inTerms.end();++t1, ++count){
      //retval.push_back(*t1);
      //vector<SQterm>::iterator t_ret = retval.end(); --t_ret;
      for(auto t2 = inTerms.begin();t2 != inTerms.end();){
        if(isAdditive(*t1, *t2) && t1!=t2) {
          t1->set_numConst(t1->get_numConst()+t2->get_numConst()); 
          t2 = inTerms.erase(t2);
          //t1->masquerade();
	} // End if
        else ++t2;
      } // End for
      if(fabs(t1->get_numConst())) retval.push_back(*t1);
    } // End for
    if((int)inTerms.size() != count){
      cout << "Something is wrong in combining terms .... " << endl;
      abort();
    } // End if

    return retval;
  }

  // *********************************************************
  // Combine terms with isAdditive()==true
  // *********************************************************
  void combineTerms(vector<SQterm> &inTerms, vector<SQterm> *outTerms)
  {
    int count = 0;
    outTerms->reserve(inTerms.size());
    for(auto t1 = inTerms.begin();t1 != inTerms.end();++t1, ++count){
      for(auto t2 = inTerms.begin();t2 != inTerms.end();){
        if(isAdditive(*t1, *t2) && t1!=t2) { 
          t1->set_numConst(t1->get_numConst()+t2->get_numConst()); 
          t2 = inTerms.erase(t2);
          //t1->masquerade();
        } // Enf if
        else ++t2;
      } // End t2
      if(fabs(t1->get_numConst())) outTerms->push_back(*t1);
    } // End for
    //cout << " : " << inTerms.size() << ", " << count << endl;
    if((int)inTerms.size() != count){
      cout << "Something is wrong in combining terms .... " << endl;
      abort();
    } // End if
  }

  // *********************************************************
  // Combine terms with isAdditive()==true
  // *********************************************************
  void combineTerms2(vector<SQterm> &inTerms, vector<SQterm> *outTerms)
  {
    int count = 0;
    outTerms->reserve(inTerms.size());
    vector<SQterm*> p_inTerms; p_inTerms.reserve(inTerms.size());
    for(auto t1 = inTerms.begin();t1 != inTerms.end();++t1) 
      p_inTerms.push_back(&(*t1));    
    for(auto t1 = p_inTerms.begin();t1 != p_inTerms.end();++t1, ++count){
      for(auto t2 = p_inTerms.begin();t2 != p_inTerms.end();){
        if(isAdditive(**t1, **t2) && t1!=t2) { 
          (*t1)->set_numConst((*t1)->get_numConst()+(*t2)->get_numConst()); 
          t2 = p_inTerms.erase(t2);
          //t1->masquerade();
        } // Enf if
        else ++t2;
      } // End t2
      if(fabs((*t1)->get_numConst())) outTerms->push_back(**t1);
    } // End for
    //cout << " : " << inTerms.size() << ", " << count << endl;
    if((int)p_inTerms.size() != count){
      cout << "Something is wrong in combining terms .... " << endl;
      abort();
    } // End if
  }

  // *********************************************************
  // Combine terms with isAdditive()==true
  // This implementation may be the most elegant one, but no warrentry!
  // *********************************************************
  void combineTerms_new(vector<SQterm> &inTerms)
  {
    for(auto t1 = inTerms.begin();t1 != inTerms.end();){
      for(auto t2 = inTerms.begin();t2 != inTerms.end();){
        if(isAdditive(*t1, *t2) && t1!=t2) { 
          t1->set_numConst(t1->get_numConst()+t2->get_numConst()); 
          t2 = inTerms.erase(t2);
          //t1->masquerade();
        } // Enf if
        else ++t2;
      } // End t2
      if(!fabs(t1->get_numConst())) t1 = inTerms.erase(t1);
      else ++t1;
    } // End for
  }

//*BUG*   // *********************************************************
//*BUG*   // Combine terms with isAdditive()==true
//*BUG*   // *********************************************************
//*BUG*   void combineTerms(const vector<SQterm> &inTerms, vector<SQterm> *outTerms)
//*BUG*   {
//*BUG*     outTerms->reserve(inTerms.size());
//*BUG*     vector<SQterm> temp(inTerms);
//*BUG*     for(size_t i = 0;i < temp.size();++i){
//*BUG*       outTerms->push_back(temp[i]);
//*BUG*       for(size_t j = i+1;j < inTerms.size();++j){
//*BUG*         if(isAdditive(inTerms[j], temp[i])){
//*BUG*           vector<SQterm>::iterator t = outTerms->end(); --t;
//*BUG*           *t = *t + inTerms[j]; 
//*BUG*         } // End if
//*BUG*       } // End j    
//*BUG*     } // End i
//*BUG*   }

  // *********************************************************
  // Decompose the RDM to contract core indices
  // *********************************************************
  void decomposeRDMCore(SQterm &inTerm, vector<SQterm> *outTerms)
  {
    if(!inTerm.get_isInCanonical()) return;
    vector<SQtensor> sfGen_list;
    vector<SQtensor> other_list;

    vector<SQtensor> tensors(inTerm.get_tensors());
    //cout << inTerm << endl;
    for(size_t t = 0;t < tensors.size();++t){
      if(is_sfGen(tensors[t].get_name()) || is_RDM(tensors[t].get_name())) 
           sfGen_list.push_back(tensors[t]);
      else other_list.push_back(tensors[t]);
    } // End t
    if(!sfGen_list.size()) return;

    SQtensor e1(*sfGen_list.begin());
    int o1 = (int)(e1.get_indices().size())/2;

    vector<int> CindCre;
    vector<int> CindDes;
    for(size_t i = 0;i < (size_t)o1;++i){
      if(e1.get_indices()[i   ]->get_char()==(char_state)0) CindCre.push_back(i);
      if(e1.get_indices()[i+o1]->get_char()==(char_state)0) CindDes.push_back(i+o1);      
    } // End i

    int nc = (int)CindCre.size();
    outTerms->reserve(Nterms());
    if(nc != CindDes.size()) return;
    if(nc == 0) return outTerms->push_back(inTerm);

    // New rank of sfGen
    int newOrder(o1 - nc);
    IIvector perms(makePermutations(nc));
    for(auto perm = perms.begin();perm != perms.end();++perm){
      // Compute all the pairs of contraction ....
      // Example: (conPairs[0][p], conPairs[1][p]) is the p-th contraction pair
      IIvector conPairs;
      Ivector p;
      conPairs.push_back(p);
      conPairs.push_back(p);
      for(size_t i = 0;i < nc;++i){
        conPairs[0].push_back(CindCre[(*perm)[i]]);
        conPairs[1].push_back(CindDes[i]);
      } // End i

      // Evaluate constant prefactor
      int prefactor = 1;
      vector<int> ind;
      for(int i = 0;i < nc;++i){
        auto i_index(find(ind.begin(), ind.end(), i));

        if(i_index == ind.end()){
          ind.push_back(i);
          int i0(conPairs[0][i]);
          int i1(conPairs[1][i] - o1);
          if(i1 == i0){
            prefactor = prefactor * 2;
            continue;
	  } // End if
          // While i1 in conPair[0] ////////////
          auto i_index(find(conPairs[0].begin(), conPairs[0].end(), i1));
          for(;i_index!=conPairs[0].end();){
            ind.push_back((int)(i_index-conPairs[0].begin()));
            i1 = conPairs[1][(size_t)(i_index-conPairs[0].begin())] - o1;
            if(i1==i0){ prefactor = prefactor * 2; break; }
            i_index = find(conPairs[0].begin(), conPairs[0].end(), i1);
	  } // End i_index
	} // End if
      } // End i

      // Initialize the term's tensor list
      vector<SQtensor> tensorList(other_list);
      vector<SQindex*> indexList;
      SQindex *ptr = NULL;
      for(int i = 0;i < 2*newOrder;++i) indexList.push_back(ptr);

      // Polulate the index list for the new excitation operator
      // Also, create a kronecker delta for each contraction pair
      int count = 0;
      SQtensor* kDelta_ptr(NULL);
      for(int i1 = 0;i1 < o1;++i1){
	auto i_index(find(conPairs[0].begin(), conPairs[0].end(), i1));  
        if(i_index != conPairs[0].end()){
          int i2 = conPairs[1][(size_t)(i_index-conPairs[0].begin())];
          SQindex* ind1(e1.get_indices()[(size_t)i1]);
          SQindex* ind2(e1.get_indices()[(size_t)i2]);
          vector<SQindex*> ind_delta;
          ind_delta.push_back(ind1);
          ind_delta.push_back(ind2);
          tensorList.push_back(kDelta(ind_delta));
          kDelta_ptr = &(tensorList[tensorList.size()-1]);
	} // End if
        else{
          indexList[count] = e1.get_indices()[i1];
          int i2(i1 + o1);
          auto i_index(find(conPairs[1].begin(), conPairs[1].end(), i2));
          for(;i_index!=conPairs[1].end();) {
            i2 = conPairs[0][(size_t)(i_index-conPairs[1].begin())] + o1;
            i_index = find(conPairs[1].begin(), conPairs[1].end(), i2);
	  } // End i_index
          indexList[count+newOrder] = e1.get_indices()[i2];
          ++count;
	} // End else
      } // End i1

      // Ensure that all slots in the index list have been filled
      for(auto ind = indexList.begin();ind != indexList.end();++ind)
        if(*ind == NULL){
          cout << "There is at least one unassigned index in the new spin-free unitary group generator" << endl;
          abort(); 
	} // End if

      // If kDelta has two indices belong to different orbital group, eliminate this term.
      if(kDelta_ptr != NULL){
        if(kDelta_ptr->get_indices()[0]->get_char() != kDelta_ptr->get_indices()[1]->get_char())
          continue;
      } // End if

      // Add the new excitation operator to the tensor list
      if(indexList.size()!=0) tensorList.push_back(sfGen(indexList));
      
      // Determine the sign 
      IIvector indPairs;
      Ivector q;
      indPairs.push_back(q);      
      indPairs.push_back(q);      
      for(int i = 0;i < o1;++i){
        indPairs[0].push_back(i);
	auto i_index(find(conPairs[0].begin(), conPairs[0].end(), i));
        if(i_index != conPairs[0].end())
          indPairs[1].push_back(conPairs[1][(size_t)(i_index-conPairs[0].begin())]-o1);
        else{
          SQindex* a(e1.get_indices()[i]);
          auto a_index(find(indexList.begin(), indexList.end(), a));
          SQindex* b(indexList[(size_t)(a_index-indexList.begin())+newOrder]);
          vector<SQindex*> temp;
          for(size_t j = o1;j < 2*o1;++j) temp.push_back(e1.get_indices()[j]);
          auto b_index(find(temp.begin(), temp.end(), b));
          indPairs[1].push_back((int)(b_index-temp.begin()));
	} // End else
      } // End i
      //cout << "************+++++++++ " << endl;     
      int nperm(get_num_perms(indPairs[0], indPairs[1]));
      double sign(nperm%2 ? -1 : 1);

      outTerms->push_back(SQterm(sign*prefactor*inTerm.get_numConst(), inTerm.get_Consts(), tensorList)); 
      //cout << outTerms->back() << endl;
    } // End perm

  }


//*OLD*   // *********************************************************
//*OLD*   // Decompose the RDM to contract core indices
//*OLD*   // *********************************************************
//*OLD*   void decomposeRDMCore(SQterm &inTerm, vector<SQterm> *outTerms)
//*OLD*   {
//*OLD*     if(!inTerm.get_isInCanonical()) return;
//*OLD*     vector<SQtensor> sfGen_list;
//*OLD*     vector<SQtensor> other_list;
//*OLD* 
//*OLD*     vector<SQtensor> tensors(inTerm.get_tensors());
//*OLD*     //cout << inTerm << endl;
//*OLD*     for(size_t t = 0;t < tensors.size();++t){
//*OLD*       if(is_sfGen(tensors[t].get_name()) || is_RDM(tensors[t].get_name())) 
//*OLD*            sfGen_list.push_back(tensors[t]);
//*OLD*       else other_list.push_back(tensors[t]);
//*OLD*     } // End t
//*OLD*     if(!sfGen_list.size()) return;
//*OLD* 
//*OLD*     SQtensor e1(*sfGen_list.begin());
//*OLD*     int o1 = (int)(e1.get_indices().size())/2;
//*OLD* 
//*OLD*     vector<int> CindCre;
//*OLD*     vector<int> CindDes;
//*OLD*     for(size_t i = 0;i < (size_t)o1;++i){
//*OLD*       if(e1.get_indices()[i   ]->get_char()==(char_state)0) CindCre.push_back(i);
//*OLD*       if(e1.get_indices()[i+o1]->get_char()==(char_state)0) CindDes.push_back(i+o1);      
//*OLD*     } // End i
//*OLD* 
//*OLD*     int nc = (int)CindCre.size();
//*OLD*     outTerms->reserve(Nterms());
//*OLD*     if(nc != CindDes.size()) return;
//*OLD*     if(nc == 0) return outTerms->push_back(inTerm);
//*OLD* 
//*OLD*     // New rank of sfGen
//*OLD*     int newOrder = o1 - nc;
//*OLD*     IIvector perms = makePermutations(nc);
//*OLD*     for(IIvector::iterator perm = perms.begin();perm != perms.end();++perm){
//*OLD*       // Compute all the pairs of contraction ....
//*OLD*       // Example: (conPairs[0][p], conPairs[1][p]) is the p-th contraction pair
//*OLD*       IIvector conPairs;
//*OLD*       Ivector p;
//*OLD*       conPairs.push_back(p);
//*OLD*       conPairs.push_back(p);
//*OLD*       for(size_t i = 0;i < nc;++i){
//*OLD*         conPairs[0].push_back(CindCre[(*perm)[i]]);
//*OLD*         conPairs[1].push_back(CindDes[i]);
//*OLD*       } // End i
//*OLD* 
//*OLD*       // Evaluate constant prefactor
//*OLD*       int prefactor = 1;
//*OLD*       vector<int> ind;
//*OLD*       for(int i = 0;i < nc;++i){
//*OLD*         vector<int>::iterator i_index = find(ind.begin(), ind.end(), i);
//*OLD* 
//*OLD*         if(i_index == ind.end()){
//*OLD*           ind.push_back(i);
//*OLD*           int i0 = conPairs[0][i];
//*OLD*           int i1 = conPairs[1][i] - o1;
//*OLD*           if(i1 == i0){
//*OLD*             prefactor = prefactor * 2;
//*OLD*             continue;
//*OLD* 	  } // End if
//*OLD*           // While i1 in conPair[0] ////////////
//*OLD*           Ivector::iterator i_index = find(conPairs[0].begin(), conPairs[0].end(), i1);
//*OLD*           for(;i_index!=conPairs[0].end();){
//*OLD*             ind.push_back((int)(i_index-conPairs[0].begin()));
//*OLD*             i1 = conPairs[1][(size_t)(i_index-conPairs[0].begin())] - o1;
//*OLD*             if(i1==i0){ prefactor = prefactor * 2;break; }
//*OLD*             i_index = find(conPairs[0].begin(), conPairs[0].end(), i1);
//*OLD* 	  } // End i_index
//*OLD* 	} // End if
//*OLD*       } // End i
//*OLD* 
//*OLD*       // Initialize the term's tensor list
//*OLD*       vector<SQtensor> tensorList(other_list);
//*OLD*       vector<SQindex*> indexList;
//*OLD*       SQindex *ptr = NULL;
//*OLD*       for(int i = 0;i < 2*newOrder;++i) indexList.push_back(ptr);
//*OLD* 
//*OLD*       // Polulate the index list for the new excitation operator
//*OLD*       // Also, create a kronecker delta for each contraction pair
//*OLD*       int count = 0;
//*OLD*       for(int i1 = 0;i1 < o1;++i1){
//*OLD* 	Ivector::iterator i_index = find(conPairs[0].begin(), conPairs[0].end(), i1);  
//*OLD*         if(i_index != conPairs[0].end()){
//*OLD*           int i2 = conPairs[1][(size_t)(i_index-conPairs[0].begin())];
//*OLD*           SQindex* ind1 = e1.get_indices()[(size_t)i1];
//*OLD*           SQindex* ind2 = e1.get_indices()[(size_t)i2];
//*OLD*           vector<SQindex*> ind_delta;
//*OLD*           ind_delta.push_back(ind1);
//*OLD*           ind_delta.push_back(ind2);
//*OLD*           tensorList.push_back(kDelta(ind_delta));
//*OLD* 	} // End if
//*OLD*         else{
//*OLD*           indexList[count] = e1.get_indices()[i1];
//*OLD*           int i2 = i1 + o1;
//*OLD*           Ivector::iterator i_index = find(conPairs[1].begin(), conPairs[1].end(), i2);
//*OLD*           for(;i_index!=conPairs[1].end();) {
//*OLD*             i2 = conPairs[0][(size_t)(i_index-conPairs[1].begin())] + o1;
//*OLD*             i_index = find(conPairs[1].begin(), conPairs[1].end(), i2);
//*OLD* 	  } // End i_index
//*OLD*           indexList[count+newOrder] = e1.get_indices()[i2];
//*OLD*           ++count;
//*OLD* 	} // End else
//*OLD*       } // End i1
//*OLD* 
//*OLD*       // Ensure that all slots in the index list have been filled
//*OLD*       for(vector<SQindex*>::iterator ind = indexList.begin();ind != indexList.end();++ind)
//*OLD*         if(*ind == NULL){
//*OLD*           cout << "There is at least one unassigned index in the new spin-free unitary group generator" << endl;
//*OLD*           abort(); 
//*OLD* 	} // End if
//*OLD* 
//*OLD*       // Add the new excitation operator to the tensor list
//*OLD*       if(indexList.size()!=0) tensorList.push_back(sfGen(indexList));
//*OLD*       
//*OLD*       // Determine the sign 
//*OLD*       IIvector indPairs;
//*OLD*       Ivector q;
//*OLD*       indPairs.push_back(q);      
//*OLD*       indPairs.push_back(q);      
//*OLD*       for(int i = 0;i < o1;++i){
//*OLD*         indPairs[0].push_back(i);
//*OLD* 	Ivector::iterator i_index = find(conPairs[0].begin(), conPairs[0].end(), i);
//*OLD*         if(i_index != conPairs[0].end())
//*OLD*           indPairs[1].push_back(conPairs[1][(size_t)(i_index-conPairs[0].begin())]-o1);
//*OLD*         else{
//*OLD*           SQindex* a = e1.get_indices()[i];
//*OLD*           vector<SQindex*>::iterator a_index = find(indexList.begin(), indexList.end(), a);
//*OLD*           SQindex* b = indexList[(size_t)(a_index-indexList.begin())+newOrder];
//*OLD*           vector<SQindex*> temp;
//*OLD*           for(size_t j = o1;j < 2*o1;++j) temp.push_back(e1.get_indices()[j]);
//*OLD*           vector<SQindex*>::iterator b_index = find(temp.begin(), temp.end(), b);
//*OLD*           indPairs[1].push_back((int)(b_index-temp.begin()));
//*OLD* 	} // End else
//*OLD*       } // End i
//*OLD*       //cout << "************+++++++++ " << endl;     
//*OLD*       int nperm = get_num_perms(indPairs[0], indPairs[1]);
//*OLD*       double sign = (nperm%2 ? -1 : 1);
//*OLD* 
//*OLD*       outTerms->push_back(SQterm(sign*prefactor*inTerm.get_numConst(), inTerm.get_Consts(), tensorList)); 
//*OLD*       //cout << outTerms->back() << endl;
//*OLD*     } // End perm
//*OLD* 
//*OLD*   }

  // *********************************************************
  // Created 2012/11/26 (not tested yet)
  // *********************************************************
  //SQterm SQterm::operator*(const SQterm &obj)
  SQterm times(SQterm &a, SQterm &b)
  {
    SQterm retval;
    // Merge Tensors_
    vector<SQtensor> ten_a(a.get_tensors());
    vector<SQtensor> ten_b(b.get_tensors());
    ten_a.insert(ten_a.end(), ten_b.begin(), ten_b.end());
    
    // Rename the dummy indices in b to work with those of a
    vector<int> a_dummies(3); // core=0, active=1, virtual=2
    for(size_t num = 0;num < a.get_summedBody().size();++num)
      if(a.get_summedBody()[num]->get_isSummed()){
	if     (a.get_summedBody()[num]->get_char() == core) ++a_dummies[0];
	else if(a.get_summedBody()[num]->get_char() == act ) ++a_dummies[1];
	else if(a.get_summedBody()[num]->get_char() == virt) ++a_dummies[2];
      } // End if

    vector<int> b_dummies(3);  // core=0, active=1, virtual=2
    for(size_t num = 0;num < b.get_summedBody().size();++num)    
      if(b.get_summedBody()[num]->get_isSummed()){
        ostringstream stm;
	if     (b.get_summedBody()[num]->get_char() == core) { 
          stm << a_dummies[0] + (b_dummies[0]++);
          b.get_summedBody()[num]->put_index("c"+stm.str());
	} // End if
	else if(b.get_summedBody()[num]->get_char() == act ) {
          stm << a_dummies[1] + (b_dummies[1]++);
          b.get_summedBody()[num]->put_index("o"+stm.str());
	} // End if
	else if(b.get_summedBody()[num]->get_char() == virt) {
          stm << a_dummies[2] + (b_dummies[2]++);
          b.get_summedBody()[num]->put_index("v"+stm.str());
	} // End if
      } // End if

    retval.set_tensors(ten_a);
    b.masquerade(); // Revert all the names of dummy indices
    
    // Merge Consts_
    vector<string> Consts_a(a.get_Consts());
    vector<string> Consts_b(b.get_Consts());
    Consts_a.insert(Consts_a.end(), Consts_b.begin(), Consts_b.end());
    retval.set_Consts(Consts_a);

    // Merge numConst_
    retval.set_numConst((a.get_numConst())*(b.get_numConst()));

    return retval;
  }

}} // Femto::core

