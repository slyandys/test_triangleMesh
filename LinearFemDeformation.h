#ifndef LINEAR_FEM_DEFORMATION
#define LINEAR_FEM_DEFORMATION

#include "TriangleMesh.h"
#include "tetgen.h" 
#include "MeshSimplification.h"
#include "CorrespondSurfaceAndVolumetricTemplate.h"
#include "RBFDeform.h"

#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_cost_function.h>
#include <vnl/algo/vnl_lbfgsb.h>

#include <ANN/ANN.h>

#include "clapack_defines.h"

#define THRESHOLD_DIV_ZERO 0.000000001
#define SIMPLIFIED_NUM 500 // 2000 set to 500 for hand JL

#define MAX_NUM_VERT_IN_TET 20000 
#define ARG "pqa" 
#define ARG_Q "pqaQ" 
#define ARG_O "pqaO" 

#define NUM_ITERATIONS 3
#define MAX_OFFSET 2

#define ACCURACY 1 //0 for lower accuracy and 1 for medium accuracy
#define OUTPUT_OPT_INFO 1 // set to 1 for extra output from lbfgsb and to 0 otherwise

//#define _FORCE_DIRECTION_NORMAL_TO_TRI //rotate the force vector such that it's normal to contact triangle


/*
This class is used to optimize a linear FEM model with respect to a known deformation, a known contact point, and a known force direction.
Note that the classes TriangleMesh, g_Part, g_Node, and g_Element can maybe be replaced using the open source trimesh library.
The class MeshSimplification performs edge collapses to get a multi-resolution hierarchy and can probably be replaced by some open source library.
For the collapse, you may also choose an option that uses the command line of MeshLab. Note that for this option, a lot of things are hard coded, so
you may need to adapt things.
*/

class LinearFemDeformation
{
public:
	LinearFemDeformation();
	LinearFemDeformation(TriangleMesh * inMesh, TriangleMesh * tempMesh = NULL, char * exportTet = NULL, int mode = 0, g_Node * contactNode = NULL);
	~LinearFemDeformation();

	void initialize();

	//Methods to set information: need to be called first:
	void setParameters(double youngsModulus, double poissonRatio);
	void setContactPtAndForce(double * contactPointAndForceIn, int numContactPoints = 1);
	void setNewPositions(g_NodeContainer newPositions, vector<int> validNodes);

	//Test if the mesh is ok: returns false if no intersection and true otherwise
	//Mode = 0: uses MeshSimplification to simplify the mesh, mode = 1: uses MeshLab to simplify the mesh, and mode = 2 expects tempMesh to be given and
	//merely computes the indices.
	bool intersectionTetgen(TriangleMesh * inMesh, TriangleMesh *& tempMesh, int mode = 0);

	//Method computes an FEM simulation: set Parameters needs to be called first!
	bool computeSimulation(bool updateSurfaceNodes, bool useNormalOffsets = false, double * forceVector = NULL);

	//Method optimizes material parameters and then computes FEM simulation: setContactPtAndForce and setNewPositions needs to be called first
	bool optimizeAndSimulate(bool useNormalOffsets = false);

	//Only call this method after computeSimulation or optimizeAndSimulate
	vector<int> getPredictedIndices();

	static vector<int> findClosestIndex(g_Vector vec, g_NodeContainer closeNodes, double contactRadius = 0);

	//Only call this method after computeSimulation or optimizeAndSimulate
	void updateCoordinates(g_NodeContainer * deformedPositions, char * exportName = NULL);
	void setCoordinates(char * filename);

	//Only call these methods after computeSimulation or optimizeAndSimulate
	g_NodeContainer getUpdatedNodes()
	{
		return updatedPositions;
	}

	//Methods to check for validity of data
	static int fem_isnan(double x)
	{
	   return x != x;
	}

	static int fem_isinf(double x)
	{
	   if ((x == x) && ((x - x) != 0.0)) return (x < 0.0 ? -1 : 1);
	   else return 0;
	}

private:
	//Method sets force vector, offet vector, modifiedK matrix, and modifiedF vector based on the newPositions and the validNodes 
	//(this can then be used to compute the simulation)
	//If a forceVector is given, use this vector as input forces. Otherwise, it is assumed that forces have been set before.
	bool computeSimulationInfo(double *& modifiedK, double *& modifiedF, bool useNormalOffsets, double * forceVector = NULL);

	//Method computes offsets based on an RBF
	void computeOffsetsRBF(bool useNormalOffsets);

	//Method sets result of simulation to offset vector of previously invalid nodes
	void updateOffsets(double *& modifiedF);

	//Methods set the force / offset vectors according to the contact points / observed nodes. Returns a boolean vector indicating if the force / offset at point i is known.
	//For setOffsets, there is an option to set all offsets at surface nodes (even at unobserved ones!)
	vector<bool> setForces();
	vector<bool> setOffsets(bool useNormalOffsets, bool useOnlyObservedNodes = true);

	//Method changes newPositions and validNodes
	bool computeInverseOfMapping(bool useNormalOffsets = false);

	//If derivative = -1: computes the matrices K or localK used for FEM. 
	//If derivative = 0: computes the partial derivative of matrices K or localK w.r.t. poissonRatio. 
	//If derivative = 1: computes the partial derivative of matrices K or localK w.r.t. youngsModulus. 
	bool computeMatrixK(int derivative = -1);
	void computeLocalMatrixK(double *& localK, int tetrahedronId, int derivative = -1);

	//Only works for matrixDimension <= 4
	double computeDeterminant(double * matrix, int matrixDimension);

	// Internal class used in the minimization process
	// Needed to be declared here because of the "friending"
	class FemCostFunction : public vnl_cost_function {
		public:
			FemCostFunction(LinearFemDeformation * fem, vector<bool> validForces) : vnl_cost_function(2) 
			{
				this->fem = fem;
				this->validForces = validForces;
			};
		    virtual void compute (vnl_vector< double > const &x, double *f, vnl_vector< double > *g);

			//For optimization w.r.t. Poisson ratio and Young's modulus
			double getEnergyAndGradient(double youngsModulusTest, double poissonRatioTest, double *& gradient);

		private:
		    LinearFemDeformation * fem;
			vector<bool> validForces;
	};
	friend void FemCostFunction::compute(vnl_vector< double > const &x, double *f, vnl_vector< double > *g);

	//Create a mesh to use for tetgen
	void createTetgenMesh(tetgenio *& in, TriangleMesh * tempMesh);

	TriangleMesh * mesh, * lowResMesh;
	g_NodeContainer allNodes;
	vector< vector<int> > allTetrahedra;
	
	vector<int> mappingIndices, invMappingIndices, contactNodes;
	vector<double> mappingWeights, invMappingWeights, mappingOffsets, invMappingOffsets;
	
	int numNodes, numTetrahedra, numContactPoints;
	double * matrixK, * offsets, * forces;
	vector<double> lengthRBFOffsets;

	double youngsModulus, poissonRatio;

	//Stored information because of how lbfgsb code is called:
	g_NodeContainer updatedPositions;
	vector<int> validNodes;
	double * contactPointAndForce;
};

#endif