// ==========================================================================
// $Id: TriangleMesh.h 3765 2013-03-15 19:25:00Z jlang $
// Wrapper code to interface TriangleMesh
// TriangleMesh implementation must support VRML reading, normal
// generation, setting of normals, vertex (node) access.
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
// $Rev: 3765 $
// $LastChangedBy: jlang $
// $LastChangedDate: 2013-03-15 15:25:00 -0400 (Fri, 15 Mar 2013) $
// ==========================================================================
#ifndef WRAPPER_TRIANGLE_MESH_H
#define WRAPPER_TRIANGLE_MESH_H

#include "TriMeshSimple.h"

#include <vector>
using std::vector;

#include <set>
using std::set;
using std::multiset;

#include <iostream>
using std::cout;
using std::endl;

#include <fstream>
using std::ifstream;
using std::ofstream;

#include <string>

#include "g_Part.h"
#include "g_Vector.h"

#include "geodesic_algorithm_exact.h"

// typedef g_Container<g_PEdge *> g_PEdgeContainer;
// typedef g_Container<g_Node *> g_NodeContainer;


template <class T>
struct svector:vector<T> {
  inline svector() : vector<T>()
  {} 
  inline svector(int n) : vector<T>(n) 
  {}    
  // already in vector
  // T operator[](int); 
};


// Should use a namespace
class Coordinate : public TriMeshSimple::Point {
  // point is really a vec3f which is a typedef from VectorT
  // geometry subpackage
 public:
  inline double x() {
    return (*this)[0];
  }
  inline double y() {
    return (*this)[1];
  }
  inline double z() {
    return (*this)[2];
  }
  inline double DistanceTo(const g_Vector& _oVec) const {
    return g_Vector(*this).DistanceTo(_oVec);
  }
  /* should just work
  inline Coordinate operator-(Coordinate&) {
  }
  */
};


// Should this be a g_Node?
/* class Node { */
/*  public: */
/*   Coordinate coordinate(); */
/* }; */
  
class Color : public TriMeshSimple::Color {
 public:
  inline Color() : TriMeshSimple::Color() {};
  inline Color( float _r, float _g, float _b ) 
    : TriMeshSimple::Color( _r, _g, _b ) {}
  inline Color( const TriMeshSimple::Color& _oCol ) : TriMeshSimple::Color( _oCol ) {}
};


class TriangleMesh : public g_Part {
 private:
  geodesic::Mesh d_geo_mesh;

  TriMeshSimple* d_tms;
  TriMeshSimple::VertexHandle* d_vhandleA;

  void calcTriMeshSimple(bool addColor=false, bool addNormals=false);

  svector<g_Vector> d_normals;
  svector<Color> d_color;

 public:
  TriangleMesh();
  TriangleMesh( const g_Part& );
  TriangleMesh( const TriangleMesh& );
  ~TriangleMesh();
  TriangleMesh& operator=( const TriangleMesh& );

  void init( ifstream&, const std::string& _ext = ".PLY" );
  void exportVRML( ofstream&, const std::string& _ext = ".PLY");

  void setColorToMesh(svector<Color>&);

  void calculateVertexNormals( bool _keep = true );

  void initGeodesicDistanceCalculation();
  double getGeodesicDistance(int, double&, svector<double>&);

  size_t getNumberOfTriangles();

  double getMeshResolution();

  void SetNormals(int _id, int _dim, double _val );
  void Normals_Resize();
  // get normal at node with id
  g_Vector getNormalAtNode(int _id);
};


#endif
