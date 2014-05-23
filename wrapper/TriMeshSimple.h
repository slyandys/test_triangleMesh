// ==========================================================================
// $Id: TriMeshSimple.h 3389 2012-10-19 08:23:25Z jlang $
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
#ifndef WRAPPER_TRI_MESH_SIMPLE_H
#define  WRAPPER_TRI_MESH_SIMPLE_H

#include "OpenMesh/Core/IO/Options.hh"
#include "OpenMesh/Core/IO/MeshIO.hh"
#include "OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh"
#include "OpenMesh/Core/Geometry/VectorT.hh"


// Define double traits
struct DoubleTraits : OpenMesh::DefaultTraits
{
  // Make Point and Normal a vector of doubles
  typedef OpenMesh::Vec3d Point;
  typedef OpenMesh::Vec3d Normal;
  // Make Color float 
  typedef OpenMesh::Vec3f Color;
};

// Define my mesh with the new traits!
typedef OpenMesh::TriMesh_ArrayKernelT<DoubleTraits> TriMeshSimple;


#endif
