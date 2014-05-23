/** Simplistic test */
#include "TriangleMesh.h"


#include <fstream>

using std::cout;
using std::endl;

using std::ofstream;
using std::ifstream;

int main(int argc, char *argv[] ) { 

  g_Part* mesh = new g_Part();

  // Create 4 Nodes
  // node 1 -- g_Part starts counting at 1!
  g_Vector coord( -1.0, -1.0, 0.0 );
  g_Node* node = new g_Node( coord );
  mesh->node( node );
  // node 2
  coord.Set( 1.0, -1.0, 0.0 );
  node = new g_Node( coord );
  mesh->node( node );
  // node 3
  coord.Set( 0.0, 0.0, 2.0 );
  node = new g_Node( coord );
  mesh->node( node );
  // node 4
  coord.Set( 0.0, 1.0, 0.0 );
  node = new g_Node( coord );
  mesh->node( node );
  cout << "4 Nodes created." << endl;

  // Make 4 faces for a tetrahedron
  // face 1-2-3
  g_Element* element = new g_Element();
  element->node(mesh->nodes()[0]);
  element->node(mesh->nodes()[1]);
  element->node(mesh->nodes()[2]);
  mesh->element(element);
  // face 2-4-3
  element = new g_Element();
  element->node(mesh->nodes()[1]);
  element->node(mesh->nodes()[3]);
  element->node(mesh->nodes()[2]);
  mesh->element(element);
  // face 3-4-1
  element = new g_Element();
  element->node(mesh->nodes()[2]);
  element->node(mesh->nodes()[3]);
  element->node(mesh->nodes()[0]);
  mesh->element(element);
  // face 1-4-2
  element = new g_Element();
  element->node(mesh->nodes()[0]);
  element->node(mesh->nodes()[3]);
  element->node(mesh->nodes()[1]);
  mesh->element(element);
  cout << "4 Faces created." << endl;
  // Get faces
  g_ElementContainer faces = mesh->elements();
  for ( g_ElementContainer::const_iterator fIt = faces.begin(); fIt != faces.end(); ++fIt ) {
    g_NodeContainer nodes = (*fIt)->nodes();
    cout << "Face: " << (*fIt)->id() << " Nodes: ";
    for ( g_NodeContainer::const_iterator nIt = nodes.begin(); nIt != nodes.end(); ++nIt ) {
      cout << (*nIt)->id() << " ";
    }
    cout << endl;
  }
  // Get edges
  g_PEdgeContainer edges = mesh->uniqueEdges();
  cout << "Got edges." << endl;
  for ( g_PEdgeContainer::const_iterator eIt = edges.begin(); eIt != edges.end(); ++eIt ) {
    int id0 = (*eIt)->firstNode().id();
    int id1 = (*eIt)->lastNode().id();
    g_Vector coord0 = mesh->nodes()[id0-1]->coordinate();
    g_Vector coord1 = mesh->nodes()[id1-1]->coordinate();
    double edgeLength = coord0.DistanceTo(coord1);
    cout << "Edge: " << id0 << " to " << id1 
	 << " Coord: [" <<  coord0 << "], [" << coord1 
	 << "] Length: " << edgeLength << endl; 
    g_ElementContainer neighbors = (*eIt)->firstNode().elements();
    cout << "Neighbors: ";
    for ( g_ElementContainer::const_iterator eleIt = neighbors.begin(); 
	  eleIt != neighbors.end(); ++eleIt ) {
      cout << (*eleIt)->id()-1 << " ";
    }
    cout << endl;
  }
  // Create normals for the mesh
  TriangleMesh trim(*mesh);
  trim.calculateVertexNormals();
  int nn = trim.getNumberOfNodes();
  for ( int nId=0; nId < nn; ++nId ) {
    cout << "Normal at " << nId+1 << " : " << trim.getNormalAtNode( nId +1) << endl;
  }
  // Get mesh resolution
  cout << "Mesh resolution is: " << trim.getMeshResolution() << endl;
  // Compute geodesic distances 
  trim.initGeodesicDistanceCalculation();
  svector<double> result(trim.getNumberOfNodes());
  double radius = 2.0;
  trim.getGeodesicDistance( 1, radius, result );
  // Set colors to mesh
  svector<Color> vColor;
  vColor.push_back( Color( 0, 0, 1 ));
  vColor.push_back( Color( 0, 1, 0 ));
  vColor.push_back( Color( 1, 0, 0 ));
  vColor.push_back( Color( 0.5, 0.5, 0.5 ));
  trim.setColorToMesh( vColor );

  // Save the mesh to file
  ofstream ofs("tetra.ply", std::fstream::binary );
  trim.exportVRML( ofs ); 
  ofs.close();
  cerr << "Stream closed" << endl;
  // delete the nodes and elements of the mesh -- not really
  delete mesh;
  cerr << "Mesh deleted" << endl;

  // read it back in
  ifstream ifs("tetra.ply", std::fstream::binary );
  trim.init(ifs);
  ifs.close();
  cerr << "Read back" << endl;


  if ( argc > 1 ) {
    // read a big mesh
    cerr << "Reading: " << argv[1] << endl;
    ifs.open(argv[1],std::fstream::binary);
    trim.init(ifs);
    ifs.close();
    trim.initGeodesicDistanceCalculation();
    result.clear();
    result.resize(trim.getNumberOfNodes());
    radius = 200.0;
    trim.getGeodesicDistance( 1, radius, result );
    cout << "Resolution: " << trim.getMeshResolution() << endl;
 }

  return 0;
}
