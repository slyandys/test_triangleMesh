#ifndef FINE_FITTING_TRACKING_LOCAL
#define FINE_FITTING_TRACKING_LOCAL

#include "util_wrapper.h"

#include "MeshSimplification.h"
#include "LinearFemDeformation.h"

#include <vector>
#include <cmath>
using namespace std; 

#define MAX_NUM_ITER 1000
#define ENERGY_TOLERANCE 0.0001
#define MAX_ANGLE 60 

#define MAX_DIST 5.0 

#define PI 3.14159265

#define MD_THRESHOLD_DIVIDE 10.0

/*
This class is used to deform a template mesh given using the class TriangleMesh (an extension of g_Part) to a frame (points + normals). 
Note that the classes TriangleMesh, g_Part, g_Node, and g_Element can maybe be replaced using the open source trimesh library.
The class MeshSimplification performs edge collapses to get a multi-resolution hierarchy and can probably be replaced by some open source library.

NOTE: THIS CLASS USES QUATERNIONS TO EXPRESS THE TRANSFORMATIONS
*/

class FineFittingTrackingLocalFrame
{
public:
	
	FineFittingTrackingLocalFrame();
	~FineFittingTrackingLocalFrame();

	void setTemplate(TriangleMesh * mesh);
	void resetTemplateCoo(g_NodeContainer nodes, char * filename = NULL);
	void setFrame(TriangleMesh * mesh, 
		      g_Vector probePos=g_Vector(), g_Vector probeOrient=g_Vector(), double depth = 0.0, double radius=-1);

	// set the position and orientation for the probe in the first frame
	// implies that the probe does not slip and remains attached to the location on the mesh
	bool setProbeConstraint( const g_Vector& probePos, const g_Vector& probeOrient );

	TriangleMesh * computeDeformedTemplate(char * filename = NULL, bool rigidAligment = false, bool * selfIntersectReturn = NULL, bool onlyUseMarkers = false);

	double * getTransformationParameters(){return transformationParameters;}

	//Use this method to adjust the geometry of the mesh to the shape predicted by the FEM simulation: 
	//changes the nodes to the transformed coordinate
	//ONLY CALL THIS AFTER computeDeformedTemplate was called
	void adjustTransformations(g_Part * targetMesh, vector<int> indicesModifiedWithFEM, bool * selfIntersectReturn = NULL);

	//Export the transformations: only call this AFTER adjustTransformations
	void exportTransformationParameters(char * filename);

	//This method can be used to define a set of manually clicked markers to achieve a better alignment to the first frame
	void setMarkers(char * markersToFirstFrame, int numMarkers, double scaleFactor);
	void unsetMarkers();

	vector<int> getNearestNeighbors(bool excludeBoundary);

	//Size of neighborhood used for access energy
	double SIZE_NBHD;

private:

	void postProcess(g_Part * currentMesh);

	void smoothTransformations(int index, char * filename = NULL);

	void alignRigidlyToFrame(bool onlyUseMarkers);
	void computeTranslationMat(double tx, double ty, double tz, double *& mat);
	void computeRotationMat(double tx, double ty, double tz, double angle, double *& mat);
	void computeScalingMat(double s, double *& mat);
	void computeTranslationGrad(double tx, double ty, double tz, double *& mat, int index);
	void computeRotationGrad(double tx, double ty, double tz, double angle, double *& mat, int index);
	void computeScalingGrad(double s, double *& mat);

	void fitToData(g_Part * mesh, vector<g_Vector> normalVectors, vector< vector<int> > clusterIds, vector<int> * isFeature = NULL, bool onlyUseMarkers = false);
	double solveOptimization(g_Part * mesh, vector<g_Vector> normalVectors, vector< vector<int> > clusterIds, vector<int> * isFeature = NULL, bool onlyUseMarkers = false);

	void computeHierarchy(int numberLevels);
	
	void computeNearestNeighbors(g_Part * mesh, vector<g_Vector> normalVectors, vector<int> * isFeature = NULL, bool onlyUseMarkers = false,
		bool excludeBoundary = false);
	double nearestNeighborEnergy(g_Part * mesh, bool useFeatures = false);
	void nearestNeighborGradient(g_Part * mesh, double *& g);
	
	void computeAccessNeighbors(g_Part * mesh, vector<int> * indices = NULL);
	double accessEnergy(g_Part * mesh);
	void accessGradient(g_Part * mesh, double *& g);
	
	void computeMdNeighbors(g_Part * mesh, vector<int> * indices = NULL);
	double mdEnergy(g_Part * mesh);
	void mdGradient(g_Part * mesh, double *& g);
	
	void computeMdNeighborsFixedBarycenters(g_Part * mesh, vector<int> * indices = NULL);
	double mdEnergyFixedBarycenters(g_Part * mesh);
	void mdGradientFixedBarycenters(g_Part * mesh, double *& g);
	
	double distanceResolutionEnergy(g_Part * mesh, g_Part * reconstructedMesh, vector<int> indices);
	void distanceResolutionGradient(g_Part * mesh, g_Part * reconstructedMesh, vector<int> indices, double *& g);

	g_Vector computeTransformedPoint(double * X, g_Vector point, int index = -1);
	
	//Member Variables:
	g_PEdgeContainer allEdges;
	int smallestMeshNum;
	double * transformationParameters;
	double nnWeight, accessWeight, mdWeight;
	double largestDist;
	TriangleMesh * templateMesh, * frameMesh, * resultMesh;
	g_Part * smallestMesh;
	vector<MeshSimplification *> hierarchy;
	//Store the normal vectors of the hierarchy
	vector< vector<g_Vector> > meshNormals;
	vector<int> nearestNeighbors;
	vector<g_Vector> approxNearestPointPlaneNeighbors;
	vector< vector<int> > isFeatureInHierarchy;
	vector< vector<int> > indicesForIntersection;
	vector< vector<int> > indicesForSimpleRegularization;
	vector< vector<int> > indicesForMd;
	vector<int> numUsedNeighbors;
	vector<int> triangleIdsForMd;
	vector<double> barycentersForMd;

	//kd-tree for frameMesh:
	bool kdInitialized;
	ANNpointArray		dataPts;
	ANNkd_tree*			kdTree;

	// Points for probe locations are stored after this index in frameMesh
	int probeNodeBegin;
	// use of constraints for probe location
	bool useProbeConstraint;
	std::vector< vector<bool> > constraintNodesHierarchy;
	
	// Internal classes used in the minimization process
	// Needed to be declared here because of the "friending"
	class RigidAlignmentCostFunction : public vnl_cost_function {
		public:
			RigidAlignmentCostFunction(FineFittingTrackingLocalFrame * tracker) : vnl_cost_function(8) {this->tracker=tracker;};
		    virtual void compute (vnl_vector< double > const &x, double *f, vnl_vector< double > *g);
			bool recomputeNeighbors;
			bool featuresOnly;
		private:
		    FineFittingTrackingLocalFrame * tracker;
			vector<int> fixedNN;
	};
	friend void RigidAlignmentCostFunction::compute(vnl_vector< double > const &x, double *f, vnl_vector< double > *g);

	class SolveOptimizationCostFunction : public vnl_cost_function {
		public:
			SolveOptimizationCostFunction(FineFittingTrackingLocalFrame * tracker, g_Part * mesh, vector<g_Vector> normalVectors, vector< vector<int> > clusterIds, vector<int> * isFeature,
				int dimension) : vnl_cost_function(dimension) 
			{
				this->tracker = tracker;
				this->mesh = mesh;
				this->normalVectors = normalVectors;
				this->clusterIds = clusterIds;
				this->isFeature = isFeature;
			}
		    virtual void compute (vnl_vector< double > const &x, double *f, vnl_vector< double > *g);
			
		private:
		    FineFittingTrackingLocalFrame * tracker;
			g_Part * mesh; 
			vector<g_Vector> normalVectors;
			vector< vector<int> > clusterIds; 
			vector<int> * isFeature;
	};
	friend void SolveOptimizationCostFunction::compute(vnl_vector< double > const &x, double *f, vnl_vector< double > *g);

	class AdjustTransformationCostFunction : public vnl_cost_function {
		public:
			AdjustTransformationCostFunction(FineFittingTrackingLocalFrame * tracker, g_Part * mesh, g_Part * targetMesh, vector<int> indices, vector<int> indicesLowRes,
				int dimension) : vnl_cost_function(dimension) 
			{
				this->tracker=tracker;
				this->mesh = mesh;
				this->targetMesh = targetMesh;
				this->indices = indices;
				this->indicesLowRes = indicesLowRes;
			}
		    virtual void compute (vnl_vector< double > const &x, double *f, vnl_vector< double > *g);
			
		private:
		    FineFittingTrackingLocalFrame * tracker;
			g_Part * mesh, * targetMesh; 
			vector<int> indices, indicesLowRes;
	};
	friend void AdjustTransformationCostFunction::compute(vnl_vector< double > const &x, double *f, vnl_vector< double > *g);

	class SmoothTransformationCostFunction : public vnl_cost_function {
		public:
			SmoothTransformationCostFunction(FineFittingTrackingLocalFrame * tracker, g_Part * mesh, g_Part * reconstructedMesh, vector<int> newIndices,
				int dimension) : vnl_cost_function(dimension) 
			{
				this->tracker=tracker;
				this->mesh = mesh;
				this->reconstructedMesh = reconstructedMesh;
				this->newIndices = newIndices;
			}
		    virtual void compute (vnl_vector< double > const &x, double *f, vnl_vector< double > *g);
			
		private:
		    FineFittingTrackingLocalFrame * tracker;
			g_Part * mesh, * reconstructedMesh; 
			vector<int> newIndices;
	};
	friend void SmoothTransformationCostFunction::compute(vnl_vector< double > const &x, double *f, vnl_vector< double > *g);
};

#endif
