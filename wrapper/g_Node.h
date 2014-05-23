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
#ifndef WRAPPER_G_NODE_H
#define WRAPPER_G_NODE_H

#include "g_Container.h"
#include "g_Vector.h"

using std::ostream;

class g_Element;
class g_Part;

typedef g_Container<g_Element *> g_ElementContainer;

// Should use a namespace

class g_Node {
 protected:
  g_Part* d_part;
  int d_id; // -1
  g_Vector d_coordinate;
  g_ElementContainer d_elements;
  friend class g_PEdge;
  friend class g_Part;
 public:
  g_Node();
  g_Node( g_Vector );
  void coordinate(const g_Vector& );
  // g_Vector* coordinate();
  g_Vector& coordinate();
  g_Vector coordinate() const;
  int id() const;
  void id( int _id );
  g_ElementContainer elements();
  void removeElement(g_Element*);
  void removeElement(g_Element&);
  void element(g_Element&);
  // void element();
  // void pEdges();
  void print(ostream & );
  bool operator==(const g_Node& ) const;
  bool operator!=(const g_Node& ) const;
};

typedef g_Container<g_Node *> g_NodeContainer;




#endif
