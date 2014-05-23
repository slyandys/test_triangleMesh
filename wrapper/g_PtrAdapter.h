// ==========================================================================
// $Id: g_PtrAdapter.h 3389 2012-10-19 08:23:25Z jlang $
// Wrapper code to interface g_plane
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
// ==========================================================================
#ifndef WRAPPER_G_PTRADAPTER_H
#define WRAPPER_G_PTRADAPTER_H

// Should use a namespace
#include <set>


#include "g_Element.h"
#include "g_PEdge.h"


template <class T> class g_PtrAdapter {
  T& d_container;
 public:
  g_PtrAdapter(T& _container) : d_container( _container ) {} 
  void removeDuplicates() {
    // Put entries in a set and get them back
    std::set<typename T::value_type> tSet;
    for ( typename T::iterator cIt=d_container.begin();
	cIt != d_container.end(); ++cIt ) {
      tSet.insert( *cIt );
    }
    d_container.clear(); // delete old pointers
    for ( typename std::set<typename T::valuetype>::iterator uniqueIt = tSet.begin();
	  uniqueIt != tSet.end(); ++uniqueIt ) { 
      d_container.insert( *uniqueIt );
    }
    return;
  }
};

// should these work by ID and if what do we need to do if two elements are the same
template <> class g_PtrAdapter<g_Element*> {
  g_ElementContainer& d_container;
 public:
  g_PtrAdapter(g_ElementContainer& _container) : d_container( _container ) {};
  void removeDuplicates() {
    // Put entries in a set and get them back
    std::set<g_Element*> eSet;
    for ( g_ElementContainer::iterator eIt=d_container.begin();
	eIt != d_container.end(); ++eIt ) {
      eSet.insert( *eIt );
    }
    d_container.clear(); // delete old pointers
    for ( std::set<g_Element*>::iterator uniqueIt = eSet.begin();
	  uniqueIt != eSet.end(); ++uniqueIt ) { 
      d_container.insert( *uniqueIt );
    }
    return;
  }
};

template <> class g_PtrAdapter<g_PEdge*> {
  g_PEdgeContainer& d_container;
 public:
  g_PtrAdapter(g_PEdgeContainer& _container) : d_container( _container ) {};
  void removeDuplicates() {
    // Put entries in a set and get them back
    std::set<g_PEdge*> eSet;
    for ( g_PEdgeContainer::iterator eIt=d_container.begin();
	eIt != d_container.end(); ++eIt ) {
      eSet.insert( *eIt );
    }
    d_container.clear(); // delete old pointers
    for ( std::set<g_PEdge*>::iterator uniqueIt = eSet.begin();
	  uniqueIt != eSet.end(); ++uniqueIt ) { 
      d_container.insert( *uniqueIt );
    }
    return;
  }
};

template <> class g_PtrAdapter<g_Node*> {
  g_NodeContainer& d_container;
 public:
  g_PtrAdapter(g_NodeContainer& _container) : d_container( _container ) {};
  void removeDuplicates() {
    // Put entries in a set and get them back
    std::set<g_Node*> nSet;
    for ( g_NodeContainer::iterator nIt=d_container.begin();
	nIt != d_container.end(); ++nIt ) {
      nSet.insert( *nIt );
    }
    d_container.clear(); // delete old pointers
    for ( std::set<g_Node*>::iterator uniqueIt = nSet.begin();
	  uniqueIt != nSet.end(); ++uniqueIt ) { 
      d_container.insert( *uniqueIt );
    }
    return;
  }
};


#endif
