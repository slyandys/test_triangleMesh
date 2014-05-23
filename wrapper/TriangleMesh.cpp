// ==========================================================================
// $Id: TriangleMesh.cpp 3829 2013-04-09 08:38:58Z jlang $
// Wrapper code to interface TriangleMesh
// TriangleMesh implementation does not support VRML IO -- use PLY instead.
// Supports normal generation, setting of normals, vertex (node) access.
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
#ifndef WRAPPER_TRIANGLE_MESH_CPP
#define WRAPPER_TRIANGLE_MESH_CPP

// #define _DEBUG_TRIANGLE_MESH_

#include "TriangleMesh.h"

using std::cerr;

TriangleMesh::TriangleMesh() : g_Part(), d_tms(0), d_vhandleA(0) {
}

TriangleMesh::TriangleMesh( const g_Part& _part ) : g_Part( _part ), d_tms(0), d_vhandleA(0) {
#ifdef DEBUG_TRIANGLE_MESH
  cerr << "TriangleMesh::TriangleMesh( const g_Part& )" << endl;
#endif
}

TriangleMesh::TriangleMesh( const TriangleMesh& _oMesh ) : g_Part( _oMesh ), d_tms(0), d_vhandleA(0),
							   d_normals( _oMesh.d_normals), d_color( _oMesh.d_color)
{
#ifdef DEBUG_TRIANGLE_MESH
  cerr << "TriangleMesh::TriangleMesh( const TriangleMesh& )" << endl;
#endif
  if (_oMesh.d_tms != 0 || _oMesh.d_vhandleA != 0 ) {
    calcTriMeshSimple( d_color.size() == d_nnodes, d_normals.size() == d_nnodes );
  }
}


TriangleMesh& TriangleMesh::operator=( const TriangleMesh& _oMesh ) {
  if ( this != &_oMesh ) {
    // Clean up this
    if (d_tms != 0) delete d_tms;
    if (d_vhandleA != 0) delete[] d_vhandleA;
    d_tms = 0;
    d_vhandleA = 0;
    // need to assign g_Part
    g_Part::operator=( _oMesh );
    d_normals = _oMesh.d_normals;
    d_color = _oMesh.d_color;
    if (_oMesh.d_tms != 0 || _oMesh.d_vhandleA != 0 ) {
      calcTriMeshSimple( d_color.size() == d_nnodes, d_normals.size() == d_nnodes );
    }
  }
  return *this;
}


TriangleMesh::~TriangleMesh() {
  if (d_tms != 0) delete d_tms;
  if (d_vhandleA != 0) delete[] d_vhandleA;
}

void TriangleMesh::init( ifstream& _ifs, const std::string& _ext ) {
  if ( d_tms != 0 ) { 
    delete d_tms;
    d_tms = 0;
  }
  // Keep tms local fo consistency
  TriMeshSimple tms; 
  tms.request_vertex_normals();
  tms.request_vertex_colors();
  // Read with the help of OpenMesh
  OpenMesh::IO::Options ropt =  OpenMesh::IO::Options::Default;
  ropt += OpenMesh::IO::Options::Binary;
  ropt += OpenMesh::IO::Options::VertexColor;
  ropt += OpenMesh::IO::Options::VertexNormal;
  if (!OpenMesh::IO::read_mesh<TriMeshSimple>(tms, _ifs, _ext, ropt)) {
    std::cerr << "ifstream read error\n";
    exit(1);
  }
  // Now initialize g_Part
  reset();
  // Create a map from vhandle to id
  std::map<TriMeshSimple::VertexHandle,int> vMap;
  // Loop over all vertices
  for (TriMeshSimple::VertexIter vIt = tms.vertices_begin();
       vIt != tms.vertices_end(); ++vIt) {
    this->node(new g_Node( tms.point(*vIt)));
#ifdef _DEBUG_TRIANGLE_MESH_    
    cerr << "Added Node: " << d_nnodes << endl; 
#endif
    vMap[vIt.handle()] = d_nnodes;
  }
  // Loop over all faces 
  for (TriMeshSimple::FaceIter fIt=tms.faces_begin();
       fIt!=tms.faces_end(); ++fIt) {
    g_Element* ele = new g_Element();
#ifdef _DEBUG_TRIANGLE_MESH_    
    cerr << "Add Face: ";
#endif
    // Loop over all vertices in the current face
    for (TriMeshSimple::FaceVertexIter fvIt=tms.fv_iter(fIt.handle());
	 fvIt; ++fvIt) {
#ifdef _DEBUG_TRIANGLE_MESH_    
      cerr << vMap[fvIt.handle()]-1 << " ";
#endif
      ele->node(d_nodes[vMap[fvIt.handle()]-1]);
    } 
    this->element(ele);
#ifdef _DEBUG_TRIANGLE_MESH_    
    cerr << " in face: " << d_nelements << endl; 
#endif
  }
  if (ropt.check(OpenMesh::IO::Options::VertexColor)) {
    cerr << "Colors read from file." << endl;
    d_color.clear();
    // Loop over all vertices and get vertex colors
    for (TriMeshSimple::VertexIter vIt = tms.vertices_begin(); 
	 vIt != tms.vertices_end(); ++vIt) {
      d_color.push_back(tms.color(*vIt));
    }
  }
  if (ropt.check(OpenMesh::IO::Options::VertexNormal)) {
    cerr << "Vertex normals read from file." << endl;
    d_normals.clear();
    // Loop over all vertices and get vertex normals
    for (TriMeshSimple::VertexIter vIt = tms.vertices_begin(); 
	 vIt != tms.vertices_end(); ++vIt) {
      d_normals.push_back(tms.normal(*vIt));
    }
  }
 return;
}


void TriangleMesh::calcTriMeshSimple( bool addColor, bool addNormals ) {
  // Convert to TriMeshSimple
  if (d_tms != 0) {
    delete d_tms;
    d_tms = 0;
  }
  d_tms = new TriMeshSimple();
  if (d_vhandleA != 0) delete[] d_vhandleA;
  d_vhandleA = new TriMeshSimple::VertexHandle[d_nnodes];
  int index = 0; 
  for ( g_NodeContainer::const_iterator nIt = d_nodes.begin();
	nIt != d_nodes.end(); ++nIt ) {
#ifdef _DEBUG_TRIANGLE_MESH_    
    cerr << (*nIt)->coordinate() << endl;
#endif
    d_vhandleA[index++] = d_tms->add_vertex(TriMeshSimple::Point((*nIt)->coordinate()));
  }
  std::vector<TriMeshSimple::VertexHandle> faceVhandles;
  for ( g_ElementContainer::const_iterator fIt=d_elements.begin(); fIt != d_elements.end(); ++fIt ) {
    g_NodeContainer faceNodes = (*fIt)->nodes(); 
    faceVhandles.clear();
#ifdef _DEBUG_TRIANGLE_MESH_    
    cerr << "Make face: ";
#endif
    for ( g_NodeContainer::const_iterator nIt = faceNodes.begin(); 
	  nIt != faceNodes.end(); ++nIt ) {
#ifdef _DEBUG_TRIANGLE_MESH_    
      cerr << (*nIt)->id() - 1 << " ";
#endif
      faceVhandles.push_back(d_vhandleA[(*nIt)->id() - 1]);
    }
#ifdef _DEBUG_TRIANGLE_MESH_    
    cerr << endl;
#endif
    d_tms->add_face(faceVhandles);
  }
  if ( addColor && d_color.size() == d_nnodes ) {
    d_tms->request_vertex_colors();
    svector<Color>::const_iterator cIt = d_color.begin();
    for (TriMeshSimple::VertexIter vIt = d_tms->vertices_begin();
	 vIt != d_tms->vertices_end(); ++vIt) {
      // cerr << "Setting color!" << endl;
      d_tms->set_color( vIt, *cIt );
      cIt++;
    }
  }
  if ( addNormals && d_normals.size() == d_nnodes ) {
    d_tms->request_vertex_normals();
    svector<g_Vector>::const_iterator nIt = d_normals.begin();
    for (TriMeshSimple::VertexIter vIt = d_tms->vertices_begin();
	 vIt != d_tms->vertices_end(); ++vIt) {
      // cerr << "Setting normal!" << endl;
      d_tms->set_normal( vIt, *nIt );
      nIt++;
    }
  }
}


void TriangleMesh::exportVRML( ofstream& _ofs,  const std::string& _ext ) {
  calcTriMeshSimple(d_color.size() == d_nnodes, d_normals.size() == d_nnodes);
  OpenMesh::IO::Options wopt = OpenMesh::IO::Options::Default;
  wopt += OpenMesh::IO::Options::Binary;
  if ( d_color.size() == d_nnodes ) {
    wopt += OpenMesh::IO::Options::VertexColor;
  }
  if ( d_normals.size() == d_nnodes ) {
    wopt += OpenMesh::IO::Options::VertexNormal;
  }
  if (!OpenMesh::IO::write_mesh(*d_tms, _ofs, _ext, wopt )) {
    std::cerr << "ofstream write error -- mesh not saved\n";
    return;
  }
  return;
}


void TriangleMesh::setColorToMesh(svector<Color>& _col) {
  d_color = _col;
  return;
}

void TriangleMesh::calculateVertexNormals( bool _keep ) {
  // could arguably be "and" because without faces no face and hence no
  // vertex normals
  if ( d_nelements > 0 || d_normals.size() != d_nnodes ) { 
    calcTriMeshSimple(d_color.size() == d_nnodes, false); // This will create a d_tms
    d_tms->request_vertex_normals();
    d_tms->request_face_normals();
    // assure we have vertex normals
    if (!d_tms->has_vertex_normals()) {
      std::cerr << "ERROR: Standard vertex property 'Normals' not available!\n";
      return;
    }
    if (!d_tms->has_face_normals()) {
      std::cerr << "ERROR: Standard face property 'Normals' not available!\n";
      return;
    }
    // OpenMesh to calculate vertex normals
    d_normals.clear();
    d_tms->update_face_normals();
    d_tms->update_vertex_normals();
    // Loop over all vertices and get vertex normals
#ifdef DEBUG_TRIANGLE_MESH
    int nCnt=0;
    cerr << "Normals: ";
    cerr << "(Nodes: " << d_nnodes << ")" << endl;
#endif
    for (TriMeshSimple::VertexIter vIt = d_tms->vertices_begin(); 
	 vIt != d_tms->vertices_end(); ++vIt) {
      d_normals.push_back(d_tms->normal(*vIt));
#ifdef DEBUG_TRIANGLE_MESH
      cerr << nCnt << ": " << d_normals[nCnt] << endl;
      ++nCnt;
#endif
    }
    // dispose the tms normals -- is this needed?
    // d_tms->release_face_normals();
    // d_tms->release_vertex_normals();
  }
  return;
}

void TriangleMesh::initGeodesicDistanceCalculation() {
  // xyzxyz... list of point coordinates
  std::vector<double> points;	
  // v0v1v2... 0-indexed face list
  std::vector<unsigned> faces;
  
  for ( g_NodeContainer::const_iterator nIt = d_nodes.begin();
	nIt != d_nodes.end(); ++nIt ) {
    g_Vector coord = (*nIt)->coordinate();
    points.push_back(coord.x());
    points.push_back(coord.y());
    points.push_back(coord.z());
  }
  for ( g_ElementContainer::const_iterator fIt=d_elements.begin(); 
	fIt != d_elements.end(); ++fIt ) {
    g_NodeContainer faceNodes = (*fIt)->nodes(); 
    for ( g_NodeContainer::const_iterator nIt = faceNodes.begin(); 
	  nIt != faceNodes.end(); ++nIt ) {
      faces.push_back((*nIt)->id()-1);
    }
  }
   
  d_geo_mesh.initialize_mesh_data(points, faces);
  return;
}


// Change to index 1 start 
// _radius for maximum calculation
// size of _geodesicDists should be the number of nodes in this mesh
double TriangleMesh::getGeodesicDistance(int _index, double& _radius, 
					 svector<double>& _geodesicDists) {

  geodesic::GeodesicAlgorithmExact algorithm(&d_geo_mesh);
  // Could consider the below approximation for memory and speed
  // geodesic::GeodesicAlgorithmDijkstra algorithm(&d_geo_mesh);		
  // geodesic::GeodesicAlgorithmSubdivision algorithm(&d_geo_mesh,2);

  // pack the source
  std::vector<geodesic::SurfacePoint> source;
  source.push_back(geodesic::SurfacePoint(&d_geo_mesh.vertices()[_index-1]));
  // pack targets -- need all vertice
  std::vector<geodesic::SurfacePoint> targets;		
  for ( int tIndex=0; tIndex < d_nnodes; ++tIndex ) {
    targets.push_back(geodesic::SurfacePoint(&d_geo_mesh.vertices()[tIndex])); 
  }
  // init from source
  algorithm.propagate(source, _radius); //cover the mesh up to radius
#ifdef DEBUG_TRIANGLE_MESH
  algorithm.print_statistics();
#endif
  // compute the paths to all targets
  std::vector<geodesic::SurfacePoint> path;
  for(int i=0; i<targets.size(); ++i) {
    algorithm.trace_back(targets[i], path);
    // find the path
    _geodesicDists[i] = geodesic::length( path );
#ifdef DEBUG_MESH_SIMPLIFICATION
    geodesic::print_info_about_path(path);
#endif
  }

  return 0.0;
}

size_t TriangleMesh::getNumberOfTriangles() {
  return d_nelements; 
};

double TriangleMesh::getMeshResolution() {
  double res = 0.0;
  int cnt=0;

  // Loop over all edges and calculate average length
  for ( g_ElementContainer::iterator fIt=d_elements.begin();
	fIt != d_elements.end(); ++fIt ) {
    // loop over the edges 
    g_PEdgeContainer f_edges = (*fIt)->pEdges();
    for ( g_PEdgeContainer::iterator eIt = f_edges.begin();
	  eIt != f_edges.end(); ++eIt ) {
      res *= (static_cast<double>(cnt)/(cnt+1));
      res += (*eIt)->length()/(++cnt);
#ifdef _DEBUG_TRIANGLE_MESH_    
      cerr << "At " << cnt+1 << " " << res << "+";
      cerr << (*eIt)->length();
      cerr << " = " << res << endl;
#endif
      delete *eIt;
      *eIt = NULL;
    }
    f_edges.clear();
  }
  return res;
}

// set normal at node with 0-index 
void TriangleMesh::SetNormals(int _id, int _dim, double _val ) {
  if ( d_normals.size() != d_nnodes ) {
    d_normals.resize( d_nnodes );
  }
  (d_normals[_id])[_dim] = _val;
  return;
}

void TriangleMesh::Normals_Resize() {
  for ( svector<g_Vector>::iterator nIt = d_normals.begin();
	nIt != d_normals.end(); ++nIt ) {
    nIt->Normalize();
  }
  return;
}

// get normal at node with index 1..N 
g_Vector TriangleMesh::getNormalAtNode(int _id) {
  assert(d_normals.size()>=_id);
  return d_normals[_id-1];
}

#endif
