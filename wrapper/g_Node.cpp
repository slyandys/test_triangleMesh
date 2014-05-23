// ==========================================================================
// Wrapper code to interface g_Node
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
// $Rev: 3407 $
// $LastChangedBy: jlang $
// $LastChangedDate: 2012-11-02 14:45:16 -0400 (Fri, 02 Nov 2012) $
// ==========================================================================
#ifndef WRAPPER_G_NODE_CPP
#define WRAPPER_G_NODE_CPP

#include "g_Node.h"
#include "g_Element.h"
#include "g_Part.h"

g_Node::g_Node() : d_id(-1), d_coordinate(0, 0, 0)
{ }

g_Node::g_Node( g_Vector _coord ) : d_id(-1), d_coordinate(_coord) 
{ }


void g_Node::coordinate(const g_Vector& _coord) {
  d_coordinate = _coord;
}

g_Vector& g_Node::coordinate() {
  return d_coordinate; 
}


g_Vector g_Node::coordinate() const {
  return d_coordinate; 
}

int g_Node::id() const {
  return d_id;
}

void g_Node::id( int _id ) {
  d_id = _id;
  return;
}


g_ElementContainer g_Node::elements() {
  return d_elements;
}

void g_Node::print(ostream & os) {
  os << "Node " << d_id << " : ";
  for ( g_ElementContainer::iterator fIt=d_elements.begin(); fIt != d_elements.end(); ++fIt ) { 
    os << (*fIt)->id() << " ";
  }
  os << endl;
}

void g_Node::removeElement(g_Element* _ele) {
#if 0
  // Removing an element from a node seems to leave the node in the element?
  // removes this g_Node from _ele
  for ( g_NodeContainer::iterator nIt=_ele->d_nodes.begin(); nIt != _ele->d_nodes.end(); ++nIt ) { 
    if ( d_id == (*nIt)->id() ) {
      cerr << "Delete node " << (*nIt)->id() << " from " << _ele->id() << endl;
      _ele->d_nodes.erase( nIt );
      break;
    }
  }
#endif
  // remove _ele from our container
  for ( g_ElementContainer::iterator fIt=d_elements.begin(); fIt != d_elements.end(); ++fIt ) { 
    if ((*fIt)->id() == _ele->id()) {
#ifdef DEBUG_ELEMENT
      cerr << "Delete element " << (*fIt)->id() << " from " << d_id << endl;
#endif
      d_elements.erase(fIt);
      break;
    }
  }
  return;
}

void g_Node::removeElement(g_Element& _ele) {
  removeElement( &_ele );
  return;
}

void g_Node::element(g_Element& _ele) {
  // Add element
  d_elements.push_back(&_ele);
  return;
}

bool g_Node::operator==(const g_Node& _oNode) const {
  return (this == &_oNode);
}

bool g_Node::operator!=(const g_Node& _oNode) const {
  return !(this == &_oNode);
}

#endif
