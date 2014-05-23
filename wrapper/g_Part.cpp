// ==========================================================================
// $Id: g_Part.cpp 3829 2013-04-09 08:38:58Z jlang $
// Wrapper code to interface g_part
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
#ifndef WRAPPER_G_PART_CPP
#define WRAPPER_G_PART_CPP

// Should use a namespace

#include "g_Node.h"
#include "g_Element.h"
#include "g_PEdge.h"

#include "g_Part.h"

#include <set>

using std::cerr;
using std::endl;

struct PEdgeComp {
  bool operator() ( g_PEdge* const & lhs, g_PEdge* const & rhs ) {
    return *lhs < *rhs;     
  }
};


// Make a deep copy
void g_Part::deep_copy(const g_Part& _oPart) {
  for ( g_NodeContainer::const_iterator nIt = _oPart.d_nodes.begin(); 
	nIt != _oPart.d_nodes.end();
	++nIt ) {
    this->node(new g_Node((*nIt)->coordinate()));
  }
  for ( g_ElementContainer::const_iterator eIt = _oPart.d_elements.begin(); 
	eIt != _oPart.d_elements.end();
	++eIt ) {
    g_Element* element = new g_Element( );
    element->node(d_nodes[(*eIt)->d_nodes[0]->id()-1]);
    element->node(d_nodes[(*eIt)->d_nodes[1]->id()-1]);
    element->node(d_nodes[(*eIt)->d_nodes[2]->id()-1]);
    this->element(element);
  }
}

g_Part::g_Part(const g_Part& _oPart) :   d_nnodes(0),
					 d_nelements(0) {
  deep_copy(_oPart);
}

  
g_Part::g_Part() : d_nnodes(0), d_nelements(0) {
} 


g_Part::~g_Part() {
  // cerr << "Dtor g_Part" << endl;
  for ( g_NodeContainer::iterator nIt = d_nodes.begin(); 
	nIt != d_nodes.end();
	++nIt ) {
    delete (*nIt);
  }
  for ( g_ElementContainer::iterator eIt = d_elements.begin(); 
	eIt != d_elements.end();
	++eIt ) {
    delete (*eIt);
  }
} 


g_Part& g_Part::operator=(const g_Part& _oPart) 
{
  if ( this != &_oPart ) {
    reset();
    deep_copy(_oPart);
  }
  return *this;
}




void g_Part::reset() {
  for ( g_NodeContainer::iterator nIt = d_nodes.begin(); 
	nIt != d_nodes.end();
	++nIt ) {
    delete (*nIt);
  }
  for ( g_ElementContainer::iterator eIt = d_elements.begin(); 
	eIt != d_elements.end();
	++eIt ) {
      delete (*eIt);
  }
  d_nodes.clear(); d_nnodes = 0; d_elements.clear(); d_nelements = 0;
}


/**
 * Build an edge list with unique edges and return
 */
g_PEdgeContainer g_Part::uniqueEdges() {
  // loop over all faces 
  std::set<g_PEdge*,PEdgeComp> edgeSet;
  for ( g_ElementContainer::iterator fIt=d_elements.begin();
	fIt != d_elements.end(); ++fIt ) {
    // loop over the edges and add
    // cerr << "Face: " << (*fIt)->id() << endl;
    g_PEdgeContainer f_edges = (*fIt)->pEdges();
    for ( g_PEdgeContainer::iterator eIt = f_edges.begin();
	  eIt != f_edges.end(); ++eIt ) {
      // cerr <<  (*eIt)->firstNode().id() << " - " << (*eIt)->lastNode().id() << endl;
      std::pair<std::set<g_PEdge*,PEdgeComp>::iterator,bool> added = edgeSet.insert( *eIt );
      // inserted element already existed?
      if ( !added.second ) {
	delete *eIt; // avoid memory leaks
      }
    }
  }
  g_PEdgeContainer allEdges;
  for ( std::set<g_PEdge*,PEdgeComp>::iterator eUnique_it = edgeSet.begin();
	eUnique_it != edgeSet.end(); ++eUnique_it ) { 
    allEdges.insert( *eUnique_it );
  }
  return allEdges;
}


/** 
 * Return a copy of the node container
 */
g_NodeContainer g_Part::nodes() {
  return d_nodes;
}


/** 
 * Return a copy of element container
 */
g_ElementContainer g_Part::elements() {
  return d_elements;
}


/**
 *  Adds a new face to the mesh; element will have handles to vertices
 */
void g_Part::element( g_Element* face) {
  // add a face to trimesh
  g_NodeContainer nodes = face->nodes();
  // Add g_Element if triangle
  if ( nodes.numberOfItems() == 3 &&
       nodes[0]->id() != -1 &&
       nodes[1]->id() != -1 &&
       nodes[2]->id() != -1 ) {
    // Update g_Element
    face->d_part = this;
    // id starts at 1
    face->d_id = ++d_nelements;
    d_elements.push_back(face);
  }
  return;
}


/**
 *  Adds a vertex to the mesh
 */
void g_Part::node( g_Node* v) {
  // Add g_Node
  v->d_part = this;
  // id starts at 1
  v->d_id = ++d_nnodes;
  d_nodes.push_back(v);
  return;
}

/*
g_Part::~g_Part() {
  for ( g_ElementContainer::iterator fIt=d_elements.begin();
	fIt != d_elements.end(); ++fIt ) {

  }
  for ( g_NodeContainer::iterator vIt=d_nodes.begin();
	vIt != d_elements.end(); ++vIt ) {
  }
}
*/

#endif
