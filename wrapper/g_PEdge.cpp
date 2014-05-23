// ==========================================================================
// $Id: g_PEdge.cpp 3390 2012-10-19 17:09:57Z jlang $
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
// $Rev: 3390 $
// $LastChangedBy: jlang $
// $LastChangedDate: 2012-10-19 13:09:57 -0400 (Fri, 19 Oct 2012) $
// ==========================================================================
#ifndef WRAPPER_G_PEDGE_CPP
#define WRAPPER_G_PEDGE_CPP

// Should use a namespace

#include "g_PEdge.h"

using std::cerr;
using std::endl;


g_PEdge::g_PEdge() {}


g_PEdge::g_PEdge( const g_Node& _first, const g_Node& _last) : d_first(_first), d_last(_last) {
}

void g_PEdge::nodes( g_Node& _first, g_Node& _last ) {
  d_first = _first; d_last = _last;
} 

g_Node g_PEdge::firstNode() {
  return d_first;
}

g_Node g_PEdge::lastNode() {
  return d_last;
}

const g_ElementContainer& g_PEdge::elements() {
  // get the elements for first and last 
  g_ElementContainer fElements = d_first.elements();
  g_ElementContainer lElements = d_last.elements();
  d_faces.clear();
  for ( g_ElementContainer::const_iterator ffIt=fElements.begin(); 
	ffIt != fElements.end(); ++ffIt ) {
    for ( g_ElementContainer::const_iterator lfIt=lElements.begin(); 
	  lfIt != lElements.end(); ++lfIt ) {
      if ( (*ffIt)->id() == (*lfIt)->id() ) d_faces.push_back( *ffIt );
    }
  }
  return d_faces;
}


double g_PEdge::length() {
  return d_first.d_coordinate.DistanceTo( d_last.d_coordinate );
}
  


bool g_PEdge::operator<(const g_PEdge &oEdge) const {
  int p0 = (d_first.id()<d_last.id()) ? d_first.id() : d_last.id();
  int p1 = (d_first.id()<d_last.id()) ? d_last.id() : d_first.id();
  int oP0 = (oEdge.d_first.id()<oEdge.d_last.id()) ? oEdge.d_first.id() : oEdge.d_last.id();
  int oP1 = (oEdge.d_first.id()<oEdge.d_last.id()) ? oEdge.d_last.id() : oEdge.d_first.id();
  if ( p0 < oP0 ||
       ( p0 == oP0 && p1 < oP1 )) {
    return true;
  } 
  return false;
}


#endif
