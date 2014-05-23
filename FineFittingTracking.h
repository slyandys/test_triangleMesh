#ifndef FINE_FITTING_TRACKING
#define FINE_FITTING_TRACKING

#include "util_wrapper.h"

#include "MeshSimplification.h"
#include "LinearFemDeformation.h"

#include <vector>
#include <math.h>
using namespace std; 

#define MAX_NUM_ITER 1000
#define ENERGY_TOLERANCE 0.0001
#define MAX_ANGLE 30 // was 60 for dino model and buste

//Size of neighborhood used for access energy
#define SIZE_NBHD 10 // was 3 for dino model // was 10 for buste
#define MAX_DIST 1.0 // was 5.0 for dino model

/*
This class is used to deform a template mesh given using the class TriangleMesh (an extension of g_Part) to a frame (points + normals). 
Note that the classes TriangleMesh, g_Part, g_Node, and g_Element can maybe be replaced using the open source trimesh library.
The class MeshSimplification performs edge collapses to get a multi-resolution hierarchy and can probably be replaced by some open source library.
*/

class FineFittingTracking
{
public:
	
	FineFittingTracking();
	~FineFittingTracking();

	void setTemplate(TriangleMesh * mesh);
	void resetTemplateCoo(g_NodeContainer nodes, char * filename = NULL);
	void setFrame(TriangleMesh * mesh);

	TriangleMesh * computeDeformedTemplate(char * filename = NULL, bool rigidAligment = false, bool * selfIntersectReturn = NULL, bool onlyUseMarkers = false);

	double * getTransformationParameters(){return transformationParameters;}
	void setInitialTransformationParameters(double * param);

	//Use this method to adjust the geometry of the mesh to the shape predicted by the FEM simulation: 
	//changes the nodes to the transformed coordinate
	//ONLY CALL THIS AFTER computeDeformedTemplate was called
	void adjustTransformations(g_Part * targetMesh, vector<int> indicesModifiedWithFEM, bool * selfIntersectReturn = NULL);

	//Set information about parts from a file: this is not used
	void setClustersFromFile(char * filename);

	//Export the transformations: only call this AFTER adjustTransformations
	void exportTransformationParameters(char * filename);

	//Set information about feature flow fields from a file: this is not used
	void setFeatures(char * featureFile);
	vector<int> getFeatures(int index = 0);

	//This method can be used to define a set of manually clicked markers to achieve a better alignment to the first frame
	void setMarkers(char * markersToFirstFrame, int numMarkers);
	void unsetMarkers();

	vector<int> getNearestNeighbors(bool symmetric);

private:

	void postProcessGeometry(vector<int> * isFeature = NULL);
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
		bool symmetric = false);
	double isometricEnergy(g_Part * mesh);
	void isometricGradient(g_Part * mesh, double *& g);
	double nearestNeighborEnergy(g_Part * mesh);
	void nearestNeighborGradient(g_Part * mesh, double *& g);
	double regularizationEnergy(g_Part * mesh, vector<int> * indices = NULL);
	void regularizationGradient(g_Part * mesh, double *& g);
	void computeAccessNeighbors(g_Part * mesh, vector<int> * indices = NULL);
	double accessEnergy(g_Part * mesh);
	void accessGradient(g_Part * mesh, double *& g);
	void computeMdNeighbors(g_Part * mesh, vector<int> * indices = NULL);
	double mdEnergy(g_Part * mesh);
	void mdGradient(g_Part * mesh, double *& g);
	double rigidEnergy(g_Part * mesh, vector<int> * indices = NULL);
	void rigidGradient(g_Part * mesh, double *& g);
	double distanceResolutionEnergy(g_Part * mesh, g_Part * reconstructedMesh, vector<int> indices);
	void distanceResolutionGradient(g_Part * mesh, g_Part * reconstructedMesh, vector<int> indices, double *& g);
	double clusterEnergy(g_Part * mesh, vector< vector<int> > clusterIds);
	void clusterGradient(g_Part * mesh, double *& g, vector< vector<int> > clusterIds);

	g_Vector computeTransformedPoint(g_Part * mesh, int index);
	
	//Member Variables:
	g_PEdgeContainer allEdges;
	int smallestMeshNum;
	double * transformationParameters, * transformationParameterInitialization;
	double isoWeight, nnWeight, regWeight, rigidWeight, accessWeight, clusterWeight, mdWeight;
	double largestDist;
	TriangleMesh * templateMesh, * frameMesh, * resultMesh;
	g_Part * smallestMesh;
	vector<MeshSimplification *> hierarchy;
	//Store the normal vectors of the hierarchy
	vector< vector<g_Vector> > meshNormals;
	vector<int> nearestNeighbors;
	vector< vector<int> > isFeatureInHierarchy;
	vector< vector<int> > indicesForIntersection;
	vector< vector< vector<int> > > clusterIndicesInHierarchy;
	vector<int> triangleIdsForMd;
	vector<double> barycentersForMd;

	//kd-tree for frameMesh:
	bool kdInitialized;
	ANNpointArray		dataPts;
	ANNkd_tree*			kdTree;

	// Internal classes used in the minimization process
	// Needed to be declared here because of the "friending"
	class RigidAlignmentCostFunction : public vnl_cost_function {
		public:
			RigidAlignmentCostFunction(FineFittingTracking * tracker) : vnl_cost_function(8) {this->tracker=tracker;};
		    virtual void compute (vnl_vector< double > const &x, double *f, vnl_vector< double > *g);
			bool recomputeNeighbors;
			bool featuresOnly;
		private:
		    FineFittingTracking * tracker;
			vector<int> fixedNN;
	};
	friend void RigidAlignmentCostFunction::compute(vnl_vector< double > const &x, double *f, vnl_vector< double > *g);

	class SolveOptimizationCostFunction : public vnl_cost_function {
		public:
			SolveOptimizationCostFunction(FineFittingTracking * tracker, g_Part * mesh, vector<g_Vector> normalVectors, vector< vector<int> > clusterIds, vector<int> * isFeature,
				int dimension) : vnl_cost_function(dimension) 
			{
				this->tracker=tracker;
				this->mesh = mesh;
				this->normalVectors = normalVectors;
				this->clusterIds = clusterIds;
				this->isFeature = isFeature;
			}
		    virtual void compute (vnl_vector< double > const &x, double *f, vnl_vector< double > *g);
			
		private:
		    FineFittingTracking * tracker;
			g_Part * mesh; 
			vector<g_Vector> normalVectors;
			vector< vector<int> > clusterIds; 
			vector<int> * isFeature;
	};
	friend void SolveOptimizationCostFunction::compute(vnl_vector< double > const &x, double *f, vnl_vector< double > *g);

	class AdjustTransformationCostFunction : public vnl_cost_function {
		public:
			AdjustTransformationCostFunction(FineFittingTracking * tracker, g_Part * mesh, g_Part * targetMesh, vector<int> indices, vector<int> indicesLowRes,
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
		    FineFittingTracking * tracker;
			g_Part * mesh, * targetMesh; 
			vector<int> indices, indicesLowRes;
	};
	friend void AdjustTransformationCostFunction::compute(vnl_vector< double > const &x, double *f, vnl_vector< double > *g);

	class PostProcessCostFunction : public vnl_cost_function {
		public:
			PostProcessCostFunction(FineFittingTracking * tracker, g_Part * mesh, g_Part * meanValueCooMesh, vector<int> invisibleIndices,
				int dimension) : vnl_cost_function(dimension) 
			{
				this->tracker=tracker;
				this->mesh = mesh;
				this->meanValueCooMesh = meanValueCooMesh;
				this->invisibleIndices = invisibleIndices;
			}
		    virtual void compute (vnl_vector< double > const &x, double *f, vnl_vector< double > *g);
			
		private:
		    FineFittingTracking * tracker;
			g_Part * mesh, * meanValueCooMesh; 
			vector<int> invisibleIndices;
	};
	friend void PostProcessCostFunction::compute(vnl_vector< double > const &x, double *f, vnl_vector< double > *g);

	class SmoothTransformationCostFunction : public vnl_cost_function {
		public:
			SmoothTransformationCostFunction(FineFittingTracking * tracker, g_Part * mesh, g_Part * reconstructedMesh, vector<int> newIndices,
				int dimension) : vnl_cost_function(dimension) 
			{
				this->tracker=tracker;
				this->mesh = mesh;
				this->reconstructedMesh = reconstructedMesh;
				this->newIndices = newIndices;
			}
		    virtual void compute (vnl_vector< double > const &x, double *f, vnl_vector< double > *g);
			
		private:
		    FineFittingTracking * tracker;
			g_Part * mesh, * reconstructedMesh; 
			vector<int> newIndices;
	};
	friend void SmoothTransformationCostFunction::compute(vnl_vector< double > const &x, double *f, vnl_vector< double > *g);
};

#endif
