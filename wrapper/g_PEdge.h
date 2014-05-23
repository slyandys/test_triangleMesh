// ==========================================================================
// $Id: g_PEdge.h 3829 2013-04-09 08:38:58Z jlang $
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
// $Rev: 3829 $
// $LastChangedBy: jlang $
// $LastChangedDate: 2013-04-09 04:38:58 -0400 (Tue, 09 Apr 2013) $
// ==========================================================================
#ifndef WRAPPER_G_PEDGE_H
#define WRAPPER_G_PEDGE_H

// Should use a namespace

#include "g_Node.h"
#include "g_Element.h"
#include "g_Container.h"


class g_PEdge {
  g_Node d_first, d_last;
  g_ElementContainer d_faces;
  
 public: 
  g_PEdge();
  g_PEdge( const g_Node&, const g_Node& );
  g_Node firstNode();
  g_Node lastNode();
  void nodes(g_Node&, g_Node&);
  const g_ElementContainer& elements();
  double length();

  bool operator<(const g_PEdge &oEdge) const;
};

typedef g_Container<g_PEdge *> g_PEdgeContainer;

#endif
