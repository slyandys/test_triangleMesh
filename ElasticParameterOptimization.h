#ifndef ELASTIC_PARAMETER_OPTIMIZATION
#define ELASTIC_PARAMETER_OPTIMIZATION

#include "TriangleMesh.h"
#include "tetgen.h" 
#include "MeshSimplification.h"
#include "CorrespondSurfaceAndVolumetricTemplate.h"

#include "clapack_defines.h"

#define THRESHOLD_DIV_ZERO 0.000000001
#define SIMPLIFIED_NUM 500 // 2000 set to 500 for hand JL

#define MAX_NUM_VERT_IN_TET 20000 
#define ARG "pqa" 
#define ARG_Q "pqaQ" 
#define ARG_O "pqaO" 


// 1e10/0/250000/80000//80685 //initial values
#define INITIAL_LAMBDA 80685 // 1100 // 138600 // 2.5e5
// 1e10/10350/5200/10000//8965
#define INITIAL_MU 8965  // 1e4 // 9900 // 1e-5 // 5200
#define INITIAL_RESIDUAL 1

#define SUPPORT_SURF_THR 0.001 //distance from support surface in meters (1mm for foamblock and ear , 1cm for redball)
#define MAX_DIST 1.0E6 //maximum number used for error

//#define  _FIX_DISPLACEMENT_DIRECTION_ //for redball data
#define _FIX_FORCE_DIRECTION_
#define WITH_FORCE_ERROR


/* This class is used to optimize the linear elasticity parameters using a two step approach.
Having known both the nodal displacement and force for contact nodes, force for interior nodes,
and just the displacement for non-contact surface nodes, we can calculate the displacement for 
interior nodes. Since we have the nodal displacement for the whole tetmesh, we can make use of
contact nodes' forces to estimate the elasticity parameters. */



class ElasticParameterOptimization
{

public:
	ElasticParameterOptimization();
	//mode = 0 : tetrahedralization by tetgen, mode = 1 : tetrahedralization read from input vtk file
	ElasticParameterOptimization(TriangleMesh * inMesh, g_Node * contactPoint, double cRadius, TriangleMesh * lresMesh = NULL, int mode = 0, char * supportSurfFileName = NULL, char * solidmeshFileName = NULL);
	~ElasticParameterOptimization();

	//computes lambda and mu using least square
	void computeElasticityParameters(int numFrame, int maxIter=20);
	//calculates the nodal displacement from initial frame mesh to current frame
	void computeDisplacement(TriangleMesh * frMesh);
	void readForce(char * forceFileName);

	class NLMinimization{
	public:
		NLMinimization(ElasticParameterOptimization * const EPO);
		double operator()(double * x);

	private:
		ElasticParameterOptimization * epo;
		int iter;
	};

	
	
	


private:
	//void computeDisplacement(double* Qcalc);
	double computeDisplacementError(double* Qcalc, int iter=-1);
	void computeFreeSurfDisplacement(double* Qcalc, double Lambda, double Mu);
	double * computeReducedYoungsModulus(double poisson = 0.45);
	void computeContactArea();
	void computeCalculatedForce(double* Fcalc, double Lambda, double Mu);
	double computeForceError(double* Fcalc);
	
	
	void InitMeshfromFile( char * solidmeshFileName = NULL);
	void printParameters(char * filename = NULL);
	void print(double * matrix, int m, int n, char * filename);
	void exportSolidMesh(char * filename, double* Qout=NULL );

	bool intersectionTetgen(TriangleMesh * inMesh, TriangleMesh *& tempMesh);
	void createTetgenMesh(tetgenio *& in, TriangleMesh * tempMesh);

	void computeMatrixK(double youngsModulus, double poissonRatio, int derivative);
	void computeLocalMatrixK(double *& localK, int tetrahedronId, double youngsModulus, double poissonRatio, int derivative);
	double computeDeterminant(double * matrix, int matrixDimension);

	void setForce();
	void setSupportSurfaceNodes(char * supportSurfFileName);
	int ISContactNode(int Node);
	void optimizationSetup(int maxIter);

	
	
	TriangleMesh * lowResMesh, * templateMesh, * frameMesh;
	g_NodeContainer allNodes;
	vector< vector<int> > allTetrahedra;
	vector<int> interNodes, contactNodes, supportSurfNodes, freeSurfNodes;
	vector<double> contactWeights, contactOffset;
	vector<int> mapStartIndices;
	vector<double> mapStartWeights, mapStartOffsets;
	vector< set<int> > trianglesOnSurface;


	
	int numNodes, numTetrahedra, numContactPoints, numInterNodes, numSurfNodes, numFrame, numSuppNodes, numFreeSNodes;
	double contactTriArea, contactRadius;
	double * matrixK, * Q, * F, * nodeError, * normalError, * lameResult;

	double d_Fmax, d_Umax;


	vector<double> frameForce, lambda, mu, residualsqrt,  overallError, overallErrorDisp, overallErrorForce;

	



	class LinearLameOpt{
	public:
		LinearLameOpt(ElasticParameterOptimization * const EPO, int maxIt);
		~LinearLameOpt();
		void computeLameParam();

	private:
		//compute decomposed matrix K for least square optimization
		void computeLocalMatrix_H_J(double *& localH, double *& localJ, int tetrahedronId);
		void computeMatrix_H_J(double *& matrixH, double *& matrixJ);

		ElasticParameterOptimization * epo;
		long int dimension, dim1, dim2, variabledim, nrhs, info, lwork;
		double alpha, beta;
		char trans;
		int numUsedEqn, maxIter;
		double * work, * matrixH, * matrixJ, * HQ, * JQ, * A, * b;
	};
	LinearLameOpt * linLameOpt;


	class NonLinearLameOpt{
	public:
		NonLinearLameOpt(ElasticParameterOptimization * const EPO, int maxIt);
		~NonLinearLameOpt();
		void computeLameParam();

		double minlambdaRange, minmuRange, minpoisson;

	private:
		ElasticParameterOptimization * epo;
		//NLMinimization  nlm;
		int n, konvge, kcount, icount, numres, ifault;
		double ynewlo, reqmin;
		double * start, * xmin, * step;

	};
	NonLinearLameOpt * nonlinLameOpt;

};


#endif
