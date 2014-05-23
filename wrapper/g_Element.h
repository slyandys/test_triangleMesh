// ==========================================================================
// Wrapper code to interface g_Element
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
// $Rev: 3834 $
// $LastChangedBy: jlang $
// $LastChangedDate: 2013-04-11 16:03:39 -0400 (Thu, 11 Apr 2013) $
// ==========================================================================
#ifndef WRAPPER_G_ELEMENT_H
#define WRAPPER_G_ELEMENT_H

#include "g_Container.h"
class g_Node;
class g_PEdge;
class g_Part;

typedef g_Container<g_PEdge *> g_PEdgeContainer;
typedef g_Container<g_Node *> g_NodeContainer;


// Should use a namespace

class g_Element {
 protected:
  g_Part* d_part;
  int d_id; // -1
  g_NodeContainer d_nodes;
  friend class g_Part;
  friend class g_Node;
 public:
  g_Element() : d_id(-1), d_part(0) {}
  void node( g_Node* );
  //  g_Node* nodes();
  int id() const;
  g_NodeContainer nodes() const;
  g_PEdgeContainer pEdges() const;
  void emptyNodeList();
  void replaceNodeAt( int&, g_Node*& );
};

typedef g_Container<g_Element *> g_ElementContainer;


#endif
