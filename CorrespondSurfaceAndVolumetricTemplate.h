#ifndef CORRES_TEMPLATES
#define CORRES_TEMPLATES

#include "util_wrapper.h"

#include "TriangleMesh.h"
#include "g_Plane.h"

#include "clapack_defines.h"

/*Class to compute the correspondences between a surface mesh and a tetrahedral mesh */
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

class CorrespondTemplates
{
public:
	CorrespondTemplates();
	CorrespondTemplates(char * surfTempOriginal, char * surfTempRigidAligned, char * volumetricNodes, char * volumetricElements, int numNodes, int numElements);
	~CorrespondTemplates();

	void setAlignedSurface(TriangleMesh * surf);
	void setTetModel(g_NodeContainer volumetricMeshNodes, vector< vector<int> > allTetrahedra);

	void deformVolumetricMeshRigidly(char * outfile = NULL);
	void computeMappings(char * outfileMappingTet = NULL, char * outfileMappingSurf = NULL, g_Node * contactPoint = NULL);
	void importMappings(char * infileMappingTet, char * infileMappingSurf, int numNodesTet, int numNodesSurf);

	vector<set<int> > getTriOnSurf(){return trianglesOnSurface;}
	vector<int> getContactPoints(){return contactTriangleId;}
	vector<double> getContactWeights(){return contactWeights;}
	vector<double> getContactOffset(){return contactoffset;}
	vector<int> getMapTet2SurfID(){return maptet2surfID;}
	vector<int> getMappingIndicesTet(){return mappingIndicesTet;}
	vector<double> getMappingWeightsTet(){return mappingWeightsTet;}
	vector<double> getMappingOffsetsTet(){return offsetsTet;}
	vector<int> getMappingIndicesSurf(){return mappingIndicesSurf;}
	vector<double> getMappingWeightsSurf(){return mappingWeightsSurf;}
	vector<double> getMappingOffsetsSurf(){return offsetsSurf;}

	g_NodeContainer getTetNodes(){return volumetricMeshNodes;}
	TriangleMesh * getOriginal(){return originalSurfaceTemplate;}
	TriangleMesh * getAligned(){return rigidlyAlignedSurfaceTemplate;}

	//Export the surface only:
	TriangleMesh * computeSurface(char * exportFilename = NULL, bool useMapping = false );

	TriangleMesh * correspondSurface(TriangleMesh * meshToDeform, TriangleMesh * frame, char * exportFilename);

private:
	void findSurfaceNodes();
	bool computeRigidAlignment(TriangleMesh * originalMesh, TriangleMesh * alignedMesh, g_NodeContainer contToBeAligned);
	// g_Vector computeClosestPoint(g_Node * vertex, g_Element * triangle);
	// vector<double> computeBarycentrics(g_Vector coordinate, g_Vector coordinateOnSurface, g_Element * triangle);
	vector<double> computeBarycentrics( g_Vector& coordinate, 
					    g_Vector& coordinateOnSurface, 
					    g_Element * const triangle);
	void saveTetMesh(const char *fName, bool useMapping=false );

	TriangleMesh * originalSurfaceTemplate, * rigidlyAlignedSurfaceTemplate;
	g_NodeContainer volumetricMeshNodes;
	vector< vector<int> > allTetrahedra;
	vector< set<int> > trianglesOnSurface;

	vector<int> surfaceNodes, maptet2surfID;
	bool meshAligned;

	//Mapping from tetrahedral mesh to surface mesh
	vector<int> mappingIndicesTet;
	vector<double> mappingWeightsTet;
	vector<double> offsetsTet;

	//Mapping from surface mesh to tetrahedral mesh
	vector<int> mappingIndicesSurf, contactTriangleId;
	vector<double> mappingWeightsSurf, contactWeights;
	vector<double> offsetsSurf, contactoffset;
};

#endif
