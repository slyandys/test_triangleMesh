#ifndef MESHSIMP
#define MESHSIMP

#include <vector>
#include <g_Part.h>
#include <g_PEdge.h>
#include <g_PtrAdapter.h>
#include <math.h>
#include <queue>
// #include <windows.h>
#include "MdLibrary.h"
#include "TriangleMesh.h"
#include "tetgen.h" 

#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_cost_function.h>
#include <vnl/algo/vnl_lbfgsb.h>

#define GRAD_TOL 0.0001
#define MDS_THRESHOLD_DIVIDE0 0.000000001

using namespace std;

//Structure containing the collapsing tree:
struct treeNode
{
	int id;
	//store ids (1-based) of parent and second children (half-edge collapse means first child is the node itself)
	int parent;
	vector<int> child2;
	//store the neighbor ids (1-based) before collapsing along with the weights. This is used when recomputing verex positions.
	vector <int> neighborIndices;
	vector <int> elementIndices;
	vector <double> omegas;
	vector <double> bs;
};

//Structure containing edges and weights:
struct priorityEdge
{
	g_Node * firstNode;
	g_Node * lastNode;
	double weight;
};

//Class does half-edge collapses to only keep the vertices given in indices.
class MeshSimplification
{
public:
	MeshSimplification();
	MeshSimplification(g_Part * mesh);
	~MeshSimplification();

	void initNewMesh(g_Part * mesh);

	//the container contains the vertices corresponding to the indices to keep:
	g_Part * reconstructMesh(g_NodeContainer cont, vector<int> & indicesToKeep, bool minimizeKilianEnergy = true);
	//For multi-level approach, we may want to move more vertices during minimization:
	g_Part * reconstructMesh(g_NodeContainer cont, vector<int> & indicesToKeep, vector<int> & indicesForMinimization, bool minKilian);

	//collapse until only numNodes nodes are left:
	g_Part * getSimplifiedMesh(int numNodes, bool checkIntersection = false);
	//collapse until only numNodes nodes are left. Give some vertices that should be collapsed (they get higher priority)
	g_Part * getSimplifiedMesh(int numNodes, vector<bool> importantToCollapse, bool checkIntersection = false);
	//collapse until numNodes are left AND keep indices:
	g_Part * getSimplifiedMesh(int numNodes, vector<int> & indicesToKeep, bool checkIntersection = false);
	//only keep indices:
	g_Part * getSimplifiedMesh(vector<int> & indicesToKeep, bool checkIntersection = false,
		bool keepNonCollapsables = false);
	vector <int> getAdditionalIndices(){return additionalIndicesKept;}
	vector <int> getAllIndices(){return allIndicesKept;}

	//Method to get the mean value coordinate of a vertex in the deformed mesh (NEEDS TO HAVE SAME STRUCTURE AS ORIGINAL MESH)
	g_Vector getMeanValueCoordinate(g_Part * deformedMesh, int index);

	g_Part * getOriginalMesh() {return originalMesh;}

	//minimize energy for a given mesh:
	bool minimizeEnergy(g_Part *& newMesh, vector<int> & indicesToKeep);

	//Define what cost to use. Default: geomCost.
	void useVolumeCost(){geomCost = false;};
	void setGeomCost(){geomCost = true;};

private:
	g_PEdge * findBestHalfEdge(priority_queue<priorityEdge, vector<priorityEdge>, 
		greater<priorityEdge> > & edgesToCollapse, bool checkIntersection);
	bool isLegalMove(g_PEdge * halfEdge, bool checkIntersection);

	//Cost by volume
	double computeHalfEdgeCostVolume(g_PEdge * halfEdge, bool checkIntersection);
	//Cost according to Garland and Heckbert:
	double computeHalfEdgeCostGeom(g_PEdge * halfEdge, bool checkIntersection);

	void halfEdgeCollapse(int vertexId1, int vertexId2, vector<int> & newNeighbors);
	void createFinalMesh(vector<int> & indicesToKeep, bool allDeleted, bool keepNonCollapsables);
	bool isBoundaryNode(int nodeIdx);
	bool leadsSelfIntersection(g_PEdge * halfEdge);
	//Check required for mean value encoding:
	bool leadsProjectionProblems(g_PEdge * halfEdge);
	vector<treeNode> treeStructure;

	vector<treeNode> originalMeshWeights;
	void storeOriginalMeshWeights();
	void initOriginal();

	//necessary to reconstruct vertex positions of collapsed nodes after motion:
	void storeMeanValueWeights(int index, vector<treeNode> & structure);
	//only call this function when the parent of this node is the root)
	g_Vector computeMeanValueVertexPos(int index, vector<treeNode> & structure);
	//after all vertices are put back, minimize an energy:
	bool minimizeEnergy(vector<int> & indicesToKeep);
	bool minimizeEnergyKilian(vector<int> & indicesToKeep);

	//Get the (numeric) partial derivative of the edge (id1, id2) w.r.t. id1 at the coordinate index
	double getPartialDerivative(int id1, int id2, int index);

	//test if the mesh is valid: just a test method for now
	bool isMeshValid(g_Part * mesh);

	g_Part * originalMesh, * workingMesh, * finalMesh;
	vector<int> additionalIndicesKept;
	vector<int> allIndicesKept;

	double perturbConst, resolution;
	bool geomCost;

	class MeanValueCostFunction : public vnl_cost_function {
		public:
			MeanValueCostFunction(MeshSimplification * simplify, int * indexToKeep, int dimension) : vnl_cost_function(dimension) 
			{
				this->simplify=simplify;
				this->indexToKeep = indexToKeep;
			}
		    virtual void compute (vnl_vector< double > const &x, double *f, vnl_vector< double > *g);
			
		private:
		    MeshSimplification * simplify;
			int * indexToKeep;
	};
	friend void MeanValueCostFunction::compute(vnl_vector< double > const &x, double *f, vnl_vector< double > *g);

	class KilianCostFunction : public vnl_cost_function {
		public:
			KilianCostFunction(MeshSimplification * simplify, int * indexToKeep, int dimension) : vnl_cost_function(dimension) 
			{
				this->simplify=simplify;
				this->indexToKeep = indexToKeep;
			}
		    virtual void compute (vnl_vector< double > const &x, double *f, vnl_vector< double > *g);
			
		private:
		    MeshSimplification * simplify;
			int * indexToKeep;
	};
	friend void KilianCostFunction::compute(vnl_vector< double > const &x, double *f, vnl_vector< double > *g);
};

#endif
