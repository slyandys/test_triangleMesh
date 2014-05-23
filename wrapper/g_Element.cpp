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
#ifndef WRAPPER_G_ELEMENT_CPP
#define WRAPPER_G_ELEMENT_CPP

#include "g_Element.h"
#include "g_Node.h"
#include "g_PEdge.h"

using std::cerr;
using std::endl;

void g_Element::node( g_Node* n ) {
  d_nodes.insert(n);
  n->element(*this);
}

int g_Element::id() const {
  return d_id;
}

g_NodeContainer g_Element::nodes() const {
  return d_nodes;
}


// Make an edge list for current face
g_PEdgeContainer g_Element::pEdges() const {
  g_PEdgeContainer edges;
  g_NodeContainer::const_iterator prev = d_nodes.begin(); 
  for ( g_NodeContainer::const_iterator nIt = ++d_nodes.begin(); 
	nIt != d_nodes.end(); ++nIt ) {
    // cerr << "At: " << (*nIt)->id() << endl;
    edges.insert(new g_PEdge(**prev,**nIt));
    prev = nIt;
  }
  edges.insert(new g_PEdge(**(--d_nodes.end()),**d_nodes.begin()));
  return edges;
}

void g_Element::emptyNodeList() {
#if 0
  // May leave pointers to this node around
  for ( g_NodeContainer::iterator nIt = ++d_nodes.begin(); 
	nIt != d_nodes.end(); ++nIt ) {
    (*nIt)->removeElement( this );
  }
#endif
  d_nodes.clear();
  return;
}

void g_Element::replaceNodeAt( int& index, g_Node*& n) {
  d_nodes[index] = n;
  return;
}




#endif
