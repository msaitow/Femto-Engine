//
//  SQindex.cc
//  
//
//  Created by Masaaki Saitow on 12/03/26.
//  Copyright (c) 2012 Masaaki Saitow. All rights reserved.
//

#include <iostream>
#include <SQindex.hpp>

using namespace std;

namespace Femto { namespace Core {
  
  // *********************************************************
  // 
  // *********************************************************
  SQindex::SQindex(const SQindex &obj)
    : index_(obj.index_),
      charactar_(obj.charactar_),
      isSummed_(obj.isSummed_),
      isExt_(obj.isExt_)
  { 
    //if((int)charactar_!=0 && (int)charactar_!=1 && (int)charactar_!=2 && (int)charactar_!=3){
    if((int)charactar_ < 0 || (int)charactar_ > 10 ){
      cout << "Incorrect orbital type is specified ... " << endl;
      abort();
    } 
  }

  // *********************************************************
  // 
  // *********************************************************
  SQindex::SQindex(const string name, const enum char_state nature, 
                   const bool summedFlag, const bool extFlag)
    : index_(name), 
      charactar_(nature),
      isSummed_(summedFlag),
      isExt_(extFlag)
  {
    if((int)charactar_!=0 && (int)charactar_!=1 && (int)charactar_!=2 && /* (int)charactar_!=3 &&*/
       (int)charactar_!=4 && (int)charactar_!=5 && (int)charactar_!=6 &&
       (int)charactar_!=7 && (int)charactar_!=8 && (int)charactar_!=9 && (int)charactar_!=10){
      cout << "Incorrect orbital type is specified ... " << endl;
      cout << "If you are interested in use of the auxiliary index for DF/RI method, call SQdf_aux instead of SQindex" << endl;
      abort();
    } 
  }

  // *********************************************************
  // 
  // *********************************************************
  string 
  SQindex::get_index() const{ return index_; }

  // *********************************************************
  // 
  // *********************************************************  
  enum char_state
  SQindex::get_char() const{ return charactar_; }

  // *********************************************************
  // 
  // *********************************************************
  bool 
  SQindex::get_isSummed() const{ return isSummed_; }

  // *********************************************************
  // 
  // *********************************************************
  bool 
  SQindex::get_isExt() const{ return isExt_; }

//*ALTERNATIVE*   // *********************************************************
//*ALTERNATIVE*   // 
//*ALTERNATIVE*   // *********************************************************
//*ALTERNATIVE*   bool
//*ALTERNATIVE*   SQindex::operator<(const SQindex &other) const
//*ALTERNATIVE*   {
//*ALTERNATIVE*     bool retval;
//*ALTERNATIVE*     if      (charactar_ >  other.charactar_) retval = true;
//*ALTERNATIVE*     else if (charactar_ <  other.charactar_) retval = false;
//*ALTERNATIVE*     else if (charactar_ == other.charactar_){
//*ALTERNATIVE*       if      (index_ >  other.index_) retval = true;
//*ALTERNATIVE*       else if (index_ <  other.index_) retval = false;
//*ALTERNATIVE*       else if (index_ == other.index_){
//*ALTERNATIVE*         if      (isSummed_ >  other.isSummed_) retval = true;
//*ALTERNATIVE*         else if (isSummed_ <  other.isSummed_) retval = false;
//*ALTERNATIVE*         else if (isSummed_ == other.isSummed_){
//*ALTERNATIVE*           if      (isExt_ >  other.isExt_) retval = true;
//*ALTERNATIVE*           else if (isExt_ <  other.isExt_) retval = false;
//*ALTERNATIVE*           else if (isExt_ == other.isExt_) retval = false;
//*ALTERNATIVE* 	}
//*ALTERNATIVE*       }
//*ALTERNATIVE*     }
//*ALTERNATIVE*     return retval;
//*ALTERNATIVE*   }

  // *********************************************************
  // 
  // *********************************************************
  bool
  SQindex::operator<(const SQindex &other) const
  {
    bool retval;
    if      (charactar_ <  other.charactar_) retval = true;
    else if (charactar_ >  other.charactar_) retval = false;
    else if (charactar_ == other.charactar_){
      if      (index_ <  other.index_) retval = true;
      else if (index_ >  other.index_) retval = false;
      else if (index_ == other.index_){
        if      (isSummed_ <  other.isSummed_) retval = true;
        else if (isSummed_ >  other.isSummed_) retval = false;
        else if (isSummed_ == other.isSummed_){
          if      (isExt_ <  other.isExt_) retval = true;
          else if (isExt_ >  other.isExt_) retval = false;
          else if (isExt_ == other.isExt_) retval = false;
	}
      }
    }
    return retval;
  }

  // *********************************************************
  // 
  // *********************************************************
  bool
  SQindex::operator>(const SQindex &other) const
  {
    bool retval = true;
    if((*this)==other) retval = false;
    else if((*this)<other) retval = false;    
    return retval;
  }
 
  // *********************************************************
  // 
  // *********************************************************
  bool
  SQindex::operator==(const SQindex &other) const
  {
    //bool retval = true;
    if (index_     != other.index_)     return false;
    if (charactar_ != other.charactar_) return false;
    if (isSummed_  != other.isSummed_)  return false;
    return true;
  }

  // *********************************************************
  // 
  // *********************************************************
  bool
  SQindex::operator!=(const SQindex &other) const
  {
    if(*this == other) return false;
    return true;
  }

  // *********************************************************
  // 
  // *********************************************************
  bool
  SQindex::operator%=(const SQindex &other) const
  {
    bool retval = true;
    if      (index_     != other.index_)     return false;
    else if (charactar_ != other.charactar_) return false;
    else if (isSummed_  != other.isSummed_)  return false;
    else if (isExt_     != other.isExt_)     return false;
    else return true;
  }

  // *********************************************************
  // 
  // *********************************************************
  void
  SQindex::put_index(const string name) 
  { 
    if(!isSummed_){ cout << "Name of non-dummy index is unchangeable" << endl; abort(); }
    index_ = name; 
  }

  // *********************************************************
  // 
  // *********************************************************
  void
  SQindex::switch_isExt(const bool extFlag) { isExt_ = extFlag; }
 
  // *********************************************************
  // 
  // *********************************************************
  SQindex&
  SQindex::operator=(const SQindex &other){
    if(this == &other) return *this;
    index_     = other.index_;
    charactar_ = other.charactar_;
    isSummed_  = other.isSummed_;
    isExt_     = other.isExt_;

    //if((int)charactar_!=0 && (int)charactar_!=1 && (int)charactar_!=2 && (int)charactar_!=3){
    if((int)charactar_ < 0 || (int)charactar_ > 10 ){
      cout << "Incorrect orbital type is specified ... " << endl;
      abort();
    } 

    return *this;
  }

  // *********************************************************
  // Added 2012/11/27  :: Basically it's same to operator =, 
  // but not always recommendable to use
  // *********************************************************
  void
  SQindex::copy(const SQindex &other){
    index_     = other.index_;
    charactar_ = other.charactar_;
    isSummed_  = other.isSummed_;
    isExt_     = other.isExt_;
  }

  // *********************************************************
  // 
  // *********************************************************
  ostream&
  operator<<(std::ostream &os, const SQindex &i)
  {
    std::string IndType;
    //*-- Spin-free index --*//
    if      (i.charactar_ == core  ) IndType = "\"core\"";
    else if (i.charactar_ == act   ) IndType = "\"active\"";
    else if (i.charactar_ == virt  ) IndType = "\"virtual\""; 
    else if (i.charactar_ == aux   ) IndType = "\"auxiliary\"";
    //*-- Alpha index --*// 
    else if (i.charactar_ == core_a) IndType = "\"core (alpha)\"";
    else if (i.charactar_ == act_a ) IndType = "\"active (alpha)\"";
    else if (i.charactar_ == virt_a) IndType = "\"virtual (alpha)\"";
    //*-- Beta index --*// 
    else if (i.charactar_ == core_b) IndType = "\"core (beta)\"";
    else if (i.charactar_ == act_b ) IndType = "\"active (bta)\"";
    else if (i.charactar_ == virt_b) IndType = "\"virtual (beta)\""; 
    //*-- AO basis index --*//
    else if (i.charactar_ == ao    ) IndType = "\"ao basis\""; 
    os << "@[" << i.index_ << ", " << IndType << "]";
    return os;
  }

}} //Femto::core
