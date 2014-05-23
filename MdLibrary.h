#ifndef MD_LIB
#define MD_LIB

#define MD_DEFAULT_EPSILON 0.000001

#include <algorithm>
#include <omp.h>
#include "g_Vector.h"
#include "g_Node.h"
#include "g_Container.h"
#include "g_PEdge.h"
#include "g_PtrAdapter.h"
#include "g_Plane.h"
#include <vector>
#include <time.h>

#define USE_ANN 1
#define MEASURE_TIME 0

#include <ANN/ANN.h>

using namespace std;

class MinimumDistance{
public:
	// If barycenters is not NULL, the barycentric coordinates of the closest points on the segments (p1, p2 or p1, u) and (p3, p4 of p3, v) are returned.
	static double minimumDistance(g_Vector& p1, g_Vector& p2, g_Vector& p3, g_Vector& p4, bool bSquared = true, double ** barycenters = NULL);
	static double minimumDistance(g_Vector& p1, g_Vector& u, g_Vector& p3, g_Vector& v, double d2121, bool bSquared = true, double ** barycenters = NULL);

	// Computes the shortest distance between t1 and t2. If returnSmallest is not NULL, the case that led to the smallest
	// distance is returned. If barycenters is not NULL, the barycentric coordinates of the points leading to the smallest
	// distance on t1 and t2 are returned (barycenters on t1 followed by barycenters on t2).
	static double minimumDistance(g_Element & t1, g_Element & t2, int * returnSmallest = NULL, double ** barycenters = NULL);

	static double boundingBoxDist(g_Vector& p1, g_Vector& p2, g_Vector& p3, g_Vector& p4);
	static void mdGrad_p1(g_Vector& p1, g_Vector& p2, g_Vector& p3, g_Vector& p4, g_Vector& grad);
	static bool intersect(g_Element & t1, g_Element & t2, double threshold = 0.00001);
	static bool intersectInPlane(g_Element & t, g_Vector & p1, g_Vector & p2, double threshold = 0.00001);

	//Non-static members:
	MinimumDistance()
	{
		edges = NULL;
		nodes = NULL;
		meshElements = NULL;
	}
	MinimumDistance(g_PEdgeContainer * edges, g_NodeContainer * nodes)
	{
		this->edges = edges;
		this->nodes = nodes;
	}
	MinimumDistance(const g_ElementContainer * elements, const g_NodeContainer * nodes)
	{
		meshElements = elements;
		this->nodes = nodes;
	}
	~MinimumDistance()
	{
		edges = NULL;
		nodes = NULL;
		meshElements = NULL;
	}

	const g_NodeContainer * getNodes(){return nodes;}

	//Segment-segment distances:
	double evaluateSegments();
	void gradientSegments(double*& gradient);
	//Only take one point into consideration (for projection):
	double evaluateSegments(int id);
	void gradientSegments(int id, double*& gradient);
	double getMaxStepSizeSegments();
	double getMaxStepSizeSegments(int id);

	//Triangle-triangle distances:
	double evaluateTriangles();
	void gradientTriangles(double*& gradient);
	double evaluateTriangles(int id);
	void gradientTriangles(int id, double*& gradient);
	double getMaxStepSizeTriangles();
	double getMaxStepSizeTriangles(int id);
	bool selfIntersection(double threshold = 0.00001);
	//Methods if only ONE nearest neighor per triangle is used for evergy and gradient:
	void updateNN();
	void updateNNKdSpherical();
	void updateNNKd();
	double evaluateTrianglesOneNN();
	void gradientTrianglesOneNN(double *& gradient);

	void updateNNKd(int id);
	double evaluateTrianglesOneNN(int id);
	void gradientTrianglesOneNN(int id, double *& gradient);

private:
	static void initializeIdentity(double *& matrix);
	static void matrixMultiply(g_Vector vec1, g_Vector vec2, double *& result);
	static g_Vector matrixMultiply(double *& matrix, g_Vector vec);
	static void transposeMatrix(double *& matrix);
	static void scalarMultiply(double scalarMultiplier, double *& result);
	static void matrixAdd(double *& mat1, double *& mat2, double *& result);
	static void addSegmentGrad(g_Node * p1, g_Node * p2, g_Node * p3, g_Node * p4, double *& gradient, int id = -1);
	static void addPointPlaneGrad(g_Element * t, g_Node * p, double *& gradient, int id = -1);
	static bool projectionInTriangle(g_Element t, g_Vector p, double ** barycenters = NULL);

	const g_NodeContainer * nodes;
	g_PEdgeContainer * edges;
	const g_ElementContainer * meshElements;

	//make it possible to keep the nearest neighbor and its distance:
	vector<int> nearestNeighbor;
	vector<double> distNN;
};

#endif
