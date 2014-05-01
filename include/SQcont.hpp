//
//  SQcont.hpp
//  
//  Small template that hold a stream of objects like SQtensor and so on.
//
//  Created by Masaaki Saitow on 14/04/10.
//  Copyright (c) 2012 Masaaki Saitow. All rights reserved.
//

#pragma once

#include <vector>
#include <iostream>
#include <algorithm>
#include <boost/format.hpp>
#include <boost/assign.hpp>

// :::::::::: Declaration ::::::::::
namespace Femto {

  ////////////////////////////////////////
  // Abstraction for std::for_each
  ////////////////////////////////////////
  template<typename Tcont, typename Tfunc>
  Tfunc for_each(Tcont &cont, Tfunc func)
  { return std::for_each(cont.begin(), cont.end(), func); }

  //////////////////////////////////////////////////
  // Standard container template used in libFemto
  //////////////////////////////////////////////////
  template<class T>
  class SQcont{

  ///////////////////////////////////////////////////////////
  ////// Future example /////////////////////////////////////
  ///////////////////////////////////////////////////////////
  private:
    std::vector<T> p_;
  public:
    const std::vector<T> &p() const { return p_; };
    std::vector<T> &p() { return p_; };
  ///////////////////////////////////////////////////////////

  public:
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructors
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    SQcont() {}
    SQcont(std::vector<T> const & ts) { p_ = ts; }
    SQcont(SQcont<T> const & ts)      { p_ = ts.p(); }
    SQcont(T const & t1)
    { 
      using namespace boost::assign;
      p_ += t1; 
    }

    SQcont(T const & t1, T const & t2)
    { 
      using namespace boost::assign;
      p_ += t1, t2;
    }    

    SQcont(T const & t1, T const & t2, T const & t3)    
    { 
      using namespace boost::assign;
      p_ += t1, t2, t3;
    }    

    SQcont(T const & t1, T const & t2, T const & t3, T const & t4)    
    { 
      using namespace boost::assign;
      p_ += t1, t2, t3, t4; 
    }    

    SQcont(T const & t1, T const & t2, T const & t3, T const & t4, T const & t5)    
    { 
      using namespace boost::assign;
      p_ += t1, t2, t3, t4, t5; 
    }    

    SQcont(T const & t1, T const & t2, T const & t3, T const & t4, T const & t5, T const & t6)    
    { 
      using namespace boost::assign;
      p_ += t1, t2, t3, t4, t5, t6; 
    }    

    SQcont(T const & t1, T const & t2, T const & t3, T const & t4, T const & t5, T const & t6, T const & t7)    
    { 
      using namespace boost::assign;
      p_ += t1, t2, t3, t4, t5, t6, t7; 
    }    

    SQcont(T const & t1, T const & t2, T const & t3, T const & t4, T const & t5, T const & t6, T const & t7, T const & t8)
    { 
      using namespace boost::assign;
      p_ += t1, t2, t3, t4, t5, t6, t7, t8; 
    }    

    SQcont(T const & t1, T const & t2, T const & t3, T const & t4, T const & t5, T const & t6, T const & t7, T const & t8, T const & t9)    
    { 
      using namespace boost::assign;
      p_ += t1, t2, t3, t4, t5, t6, t7, t8, t9; 
    }    

    SQcont(T const & t1, T const & t2, T const & t3, T const & t4, T const & t5, T const & t6, T const & t7, T const & t8, T const & t9, T const & t10)
    { 
      using namespace boost::assign;
      p_ += t1, t2, t3, t4, t5, t6, t7, t8, t9, t10; 
    }    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //////////////// Interfaces ////////////////////////////////////////
    typename std::vector<T>::iterator         begin()  { return p_.begin();  }
    typename std::vector<T>::iterator         end()    { return p_.end();    }
    typename std::vector<T>::reverse_iterator rbegin() { return p_.rbegin(); }
    typename std::vector<T>::reverse_iterator rend()   { return p_.rend();   }

    typename std::vector<T>::const_iterator cbegin()  const { return p_.cbegin();  }
    typename std::vector<T>::const_iterator cend()    const { return p_.cend();    }
    typename std::vector<T>::const_reverse_iterator rbegin() const { return p_.rbegin(); }
    typename std::vector<T>::const_reverse_iterator rend()   const { return p_.rend();   }

    typename std::vector<T>::iterator erase(typename std::vector<T>::iterator t) { return p_.erase(t); }

    void clear() { p_.clear(); }
    void reserve(const size_t num) { p_.reserve(num); }
    void p_resize(const size_t num) { p_.resize(num); }
    void pop_back() { p_.pop_back(); }
    void push_back(T const & t) { p_.push_back(t); }
    void insert(typename std::vector<T>::iterator str, typename std::vector<T>::iterator bgn, typename std::vector<T>::iterator end)
    { p_.insert(str, bgn, end); }
    void insert(typename std::vector<T>::iterator pos, T const & t)
    { p_.insert(pos, t); }
    size_t size() const { return p_.size(); }

    bool operator == (SQcont<T> const & t) { return (this->p_ == t.p_ ? true  : false); }
    bool operator != (SQcont<T> const & t) { return (this->p_ != t.p_ ? true  : false); }
    bool operator  > (SQcont<T> const & t) { return (this->p_  > t.p_ ? true  : false); }
    bool operator  < (SQcont<T> const & t) { return (this->p_  < t.p_ ? true  : false); }

    //bool &operator >= (SQcont<T> const & t) { return (*this) >= t; }
    //bool &operator <= (SQcont<T> const & t) { return (*this) <= t; }

    T & back() { return p_.back(); }
    const T & back() const { return p_.back(); }
    T at(size_t num_t) { return p_[num_t]; }
    T &operator [] (size_t num_t) { return p_[num_t]; }
    T const &operator [] (size_t num_t) const { return p_[num_t]; }
    SQcont<T> & operator = (SQcont<T> const & ts) { p_ = ts.p_; return *this; }
    ////////////////////////////////////////////////////////////////////

    int count(T const & t);
    SQcont<T> &operator <= (T const & t); 
    SQcont<T> &operator += (SQcont<T> & ts); 
    SQcont<T> &operator += (SQcont<T>   ts); 
    SQcont<T> &operator +  (SQcont<T> & t) { SQcont<T> out(p_); out += t; return out;}

    SQcont<T> diff_set(SQcont<T> & t);

  };
  
  template <class T>
  std::ostream & operator <<  (std::ostream &os, SQcont<T> const &t);
  
  template <class T>
  std::ostream & operator <<= (std::ostream &os, SQcont<T> const &t);
  
} // Femto 

// :::::::::: Implementation ::::::::::
namespace Femto {

  // *********************************************************
  // 
  // *********************************************************
  template <class T>
  int SQcont<T>::count (T const &t)
  {
    int num(0);
    for(auto i = p_.cbegin();i != p_.cend();++i) if(*i == t) ++num;
    return num;
  }

  // *********************************************************
  // 
  // *********************************************************
  template <class T>
  SQcont<T> & SQcont<T>::operator += (SQcont<T> & ts)
  {
    insert(p_.end(), ts.p().begin(), ts.p().end());
    return *this;
  }

  // *********************************************************
  // 
  // *********************************************************
  template <class T>
  SQcont<T> & SQcont<T>::operator += (SQcont<T> ts)
  {
    insert(p_.end(), ts.p().begin(), ts.p().end());
    return *this;
  }

  // *********************************************************
  // 
  // *********************************************************
  template <class T>
  SQcont<T> & SQcont<T>::operator <= (T const &t)
  {
    p_.push_back(t);
    return *this;
  }

  // *********************************************************
  // Make difference set of t1 (this) to t2
  // *********************************************************
  template<class T>
  SQcont<T> SQcont<T>::diff_set(SQcont<T> & t2)
  {
    SQcont<T> outs;
    for(auto s = p_.begin();s != p_.end();++s) 
      if(!t2.count(*s) && !outs.count(*s)) outs <= *s;
    return outs;
  }

//  // *********************************************************
//  // 
//  // *********************************************************
//  template <class T>
//  T & SQcont<T>::operator [] (size_t num_t)
//  { return p_[num_t]; }

  // *********************************************************
  // 
  // *********************************************************
  template <class T>
  std::ostream & operator <<= (std::ostream &os, SQcont<T> const &t)
  {
    for(auto i = t.p().cbegin();i != t.p().cend();++i)
      os << *i;

    return os;
  }

  // *********************************************************
  // 
  // *********************************************************
  template <class T>
  std::ostream & operator << (std::ostream &os, SQcont<T> const &t)
  {
    int num(0);
    for(auto i = t.p().begin();i != t.p().end();++i){
      os << (boost::format("  [%4d ] ") % num++) << *i << std::endl;
    } // End i

    return os;
  }

  /////////////////////////////////////////
  // Utilities
  /////////////////////////////////////////
  template<class T>
  typename std::vector<T>::iterator search(SQcont<T> const & ts, T const &t)
  { return std::find(ts.begin(), ts.end(), t); }

  //////////////////////////////////////////////////
  // Make union of two SQconts, t1 and t2
  //////////////////////////////////////////////////
  template<typename T>
  SQcont<T> make_union(SQcont<T> const &t1, SQcont<T> const &t2)
  {
    SQcont<T> outs;
    for(auto t = t1.cbegin();t != t1.cend();++t) if(!outs.count(*t)) outs <= *t;
    for(auto t = t2.cbegin();t != t2.cend();++t) if(!outs.count(*t)) outs <= *t;
    return outs;
  }

  //////////////////////////////////////////////////
  // Make intersection of two SQconts, t1 and t2
  //////////////////////////////////////////////////
  template<typename T>
  SQcont<T> make_intersection(SQcont<T> &t1, SQcont<T> &t2)
  {
    SQcont<T> outs;
    for(auto t = t1.cbegin();t != t1.cend();++t) if(t2.count(*t) && !outs.count(*t)) outs <= *t;
    for(auto t = t2.cbegin();t != t2.cend();++t) if(t1.count(*t) && !outs.count(*t)) outs <= *t;
    return outs;
  }

  //////////////////////////////////////////////////////////
  // Make symmetric difference of two SQconts, t1 and t2
  //////////////////////////////////////////////////////////
  template<typename T>
  SQcont<T> make_symmdiff(SQcont<T> &t1, SQcont<T> &t2)
  {
    SQcont<T> outs;
    for(auto t = t1.cbegin();t != t1.cend();++t) if(!t2.count(*t) && !outs.count(*t)) outs <= *t;
    for(auto t = t2.cbegin();t != t2.cend();++t) if(!t1.count(*t) && !outs.count(*t)) outs <= *t;
    return outs;
  }

  ////////////////////////////////////////////////////////////
  // Eliminate the duplicated member in the given SQcont<T>
  ////////////////////////////////////////////////////////////
  template<typename T>
  SQcont<T> eliminate_duplicate(SQcont<T> &t)
  {
    SQcont<T> outs;
    for(auto i = t.cbegin();i != t.cend();++i) if(!outs.count(*i)) outs <= *i;
    return outs;
  }

} // Femto
 

