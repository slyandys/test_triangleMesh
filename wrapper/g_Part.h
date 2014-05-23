// ==========================================================================
// $Id: g_Part.h 3406 2012-10-31 16:23:40Z jlang $
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
// $Rev: 3406 $
// $LastChangedBy: jlang $
// $LastChangedDate: 2012-10-31 12:23:40 -0400 (Wed, 31 Oct 2012) $
// ==========================================================================
#ifndef WRAPPER_G_PART_H
#define WRAPPER_G_PART_H

// Should use a namespace

#include "g_Node.h"
#include "g_Element.h"
#include "g_PEdge.h"

#include <set>

class g_Part {
  // Make sure we can access g_Part from related types
  friend class g_Node;
  friend class g_Element;
  friend class g_PEdge;

  // copy all nodes and elements
  void deep_copy(const g_Part& );

 protected:
  // Store vertices
  int d_nnodes;
  g_NodeContainer d_nodes;
  // Store faces
  int d_nelements;
  g_ElementContainer d_elements;
  // Edges will be generated on the fly if needed

  // This will delete all pointers of the nodes and elements
  void reset();

 public:
  g_Part();
  g_Part(const g_Part& _oPart);
  virtual ~g_Part(); 

  g_Part& operator=( const g_Part& _oPart );

  g_PEdgeContainer uniqueEdges();
  g_NodeContainer nodes();
  g_ElementContainer elements();

  void element( g_Element* ); // adds a new face to the mesh; element will have handles to vertices
  void node( g_Node* ); // adds a vertex to the mesh

  inline int getNumberOfNodes() const {
    return d_nnodes;
  }
  // Destructor which frees as side effect all container elements
  // ~g_Part();
};
  

#endif
