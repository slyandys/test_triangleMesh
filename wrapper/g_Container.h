// ==========================================================================
// $Id: g_Container.h 3389 2012-10-19 08:23:25Z jlang $
// Wrapper code to interface 
// ==========================================================================
// (C)opyright:
//
//   Jochen Lang
//   SITE, University of Ottawa
//   800 King Edward Ave.
//   Ottawa, On., K1N 6N5
//   Canada.
//   http://www.site.uottawa.ca
//
// Creator: Jochen Lang
// Email:   jlang@site.uottawa.ca
// ==========================================================================
// $Rev: 3389 $
// $LastChangedBy: jlang $
// $LastChangedDate: 2012-10-19 04:23:25 -0400 (Fri, 19 Oct 2012) $
#ifndef WRAPPER_G_CONTAINER_H
#define WRAPPER_G_CONTAINER_H

#include <vector>

// Should use a namespace

template <class T>
class g_Container : public std::vector<T> {
 public:
  // T& operator[] (int index);
  // const T& operator[] (int index) const;
  /* void insert(T& item) {
    push_back( item );
  }
  */

  void insert(const T& item) {
    push_back( item );
  }

  void append(const g_Container<T>& oContainer) {
    std::vector<T>::insert(this->end(), oContainer.begin(), oContainer.end()); 
  }
  // int id() const;
  int numberOfItems() const {
    return this->size();
  }
  // void clear();
};

  
#endif
