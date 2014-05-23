#include "CorrespondSurfaceAndVolumetricTemplate.h"

#include "distPoint3Triangle3.h"

using std::min;
using std::max;

using Wm5::computeClosestPoint;


CorrespondTemplates::CorrespondTemplates()
{
	meshAligned = false;
	originalSurfaceTemplate = NULL;
	rigidlyAlignedSurfaceTemplate = NULL;
}

CorrespondTemplates::CorrespondTemplates(char * surfTempOriginal, char * surfTempRigidAligned, char * volumetricNodes, char * volumetricElements, 
										 int numNodes, int numElements)
{
        originalSurfaceTemplate = new TriangleMesh();
	if (!importMeshWrapper(surfTempOriginal,originalSurfaceTemplate)) {
		cout<<"Problem opening "<<surfTempOriginal<<endl;
		delete originalSurfaceTemplate;
		originalSurfaceTemplate = 0;
		return;
	}

        rigidlyAlignedSurfaceTemplate = new TriangleMesh();
	if (!importMeshWrapper(surfTempRigidAligned,rigidlyAlignedSurfaceTemplate)) {
		cout<<"Problem opening "<<surfTempOriginal<<endl;
		delete rigidlyAlignedSurfaceTemplate;
		rigidlyAlignedSurfaceTemplate = 0;
		return; // will leak
	}

	int i, intRead1, intRead2, intRead3, intRead4;
	float floatRead1, floatRead2, floatRead3;

	FILE * fp = fopen(volumetricNodes, "r");
	if(fp == NULL)
	{
		cout<<"File "<<volumetricNodes<<" not found"<<endl;
		return;
	}
	for(i = 0; i < numNodes; i++)
	{
		fscanf(fp, "%*d %*d %f %f %f ", &floatRead1, &floatRead2, &floatRead3);
		g_Node * node = new g_Node();
		node->coordinate(g_Vector(floatRead1, floatRead2, floatRead3));
		volumetricMeshNodes.insert(node);
	}
	fclose(fp);

	fp = fopen(volumetricElements, "r");
	if(fp == NULL)
	{
		cout<<"File "<<volumetricElements<<" not found"<<endl;
		return;
	}
	for(i = 0; i < numElements; i++)
	{
		fscanf(fp, "%*d %*d %*d %d %d %d %d ", &intRead1, &intRead2, &intRead3, &intRead4);
		vector<int> tetrahedron;
		tetrahedron.push_back(intRead1-1);
		tetrahedron.push_back(intRead2-1);
		tetrahedron.push_back(intRead3-1);
		tetrahedron.push_back(intRead4-1);
		allTetrahedra.push_back(tetrahedron);
	}
	fclose(fp);

	meshAligned = false;
}

CorrespondTemplates::~CorrespondTemplates()
{
	if(originalSurfaceTemplate != NULL)
	{
		delete originalSurfaceTemplate;
		originalSurfaceTemplate = NULL;
	}
	if(rigidlyAlignedSurfaceTemplate != NULL)
	{
		delete rigidlyAlignedSurfaceTemplate;
		rigidlyAlignedSurfaceTemplate = NULL;
	}
}

void CorrespondTemplates::setAlignedSurface(TriangleMesh * surf)
{
	originalSurfaceTemplate = new TriangleMesh(*surf);
	rigidlyAlignedSurfaceTemplate = new TriangleMesh(*surf);
	meshAligned = true;
}

void CorrespondTemplates::setTetModel(g_NodeContainer meshNodes, vector< vector<int> > allTetrahedra)
{
	volumetricMeshNodes.clear();
	for(int i = 0; i < meshNodes.numberOfItems(); i++)
		volumetricMeshNodes.insert(meshNodes[i]);
	this->allTetrahedra = allTetrahedra;
}

void CorrespondTemplates::deformVolumetricMeshRigidly(char * outfile)
{
	if(!meshAligned)
	{
		//Align the volumetric mesh using the known alignment of the template mesh
		computeRigidAlignment(originalSurfaceTemplate, rigidlyAlignedSurfaceTemplate, volumetricMeshNodes);

		meshAligned = true;
	}

	//Export
	if(outfile != NULL)
	{
		FILE * fp = fopen(outfile, "w");
		if(fp == NULL)
		{
			cout<<"File "<<outfile<<" could not be created"<<endl;
			return;
		}
		for(int i = 0; i < volumetricMeshNodes.numberOfItems(); i++)
			fprintf(fp, "%d -1 %f %f %f\n", i, volumetricMeshNodes[i]->coordinate().x(), volumetricMeshNodes[i]->coordinate().y(), volumetricMeshNodes[i]->coordinate().z());
		fclose(fp);
	}
}

bool CorrespondTemplates::computeRigidAlignment(TriangleMesh * originalMesh, TriangleMesh * alignedMesh, g_NodeContainer contToBeAligned)
{
	long int i, numberOfNodes, dimension, dimensionH, lwork, info, determinant;
	double alpha, beta;
	char job, trans;
	bool returnval;
	numberOfNodes = originalMesh->getNumberOfNodes();
	dimension = 3;
	dimensionH = 4;
	job = 'A';
	trans = 'N';
	alpha = 1.0;
	beta = 0.0;
	lwork = 2*numberOfNodes;
	double * work = new double[lwork];
	double * A = new double[numberOfNodes * dimensionH];
	double * b = new double[numberOfNodes * dimension];
	double * S = new double[dimension];
	double * U = new double[dimension*dimension];
	double * VT = new double[dimension*dimension];
	double * R = new double[dimensionH*dimension];

	returnval = true;

	//populate the arrays A and b:
	for(i = 0; i < numberOfNodes; i++)
	{
		A[i] = originalMesh->nodes()[i]->coordinate().x();
		A[numberOfNodes + i] = originalMesh->nodes()[i]->coordinate().y();
		A[2*numberOfNodes + i] = originalMesh->nodes()[i]->coordinate().z();
		A[3*numberOfNodes + i] = 1.0;
		b[i] = alignedMesh->nodes()[i]->coordinate().x();
		b[numberOfNodes + i] = alignedMesh->nodes()[i]->coordinate().y();
		b[2*numberOfNodes + i] = alignedMesh->nodes()[i]->coordinate().z();
	}
	//compute an estimate of R:
	clapack::dgels_(&trans, &numberOfNodes, &dimensionH, &dimension, A, &numberOfNodes, b, &numberOfNodes, work, &lwork, &info);
	if(info != 0)
	{
		cout<<"Problem with estimating R "<<info<<endl;
		goto align_EXIT;
	}
	//make sure that R is a valid matrix (orthonormal): b contains R
	delete [] work;
	lwork = max(3 * min(dimension, dimensionH) + max(dimension, dimensionH), 5*min(dimension, dimensionH));
	work = new double[lwork];
	//only copy 3 by 3 submatrix to R:
	for(i = 0; i < dimension; i++)
	{
		R[i] = b[i];
		R[dimension + i] = b[numberOfNodes + i];
		R[2*dimension + i] = b[2*numberOfNodes + i];
	}
	//copy the modified version back to b (as help storage) and then to R in the correct order
	R[10] = R[8]; R[9] = R[7]; R[8] = R[6];
	R[6] = R[5]; R[5] = R[4]; R[4] = R[3];
	R[3] = b[3]; R[7] = b[numberOfNodes+3]; R[11] = b[2*numberOfNodes+3];

	//transform the coordinates:
	numberOfNodes = contToBeAligned.numberOfItems();
	delete [] A;
	A = new double[numberOfNodes * dimensionH];
	//populate the array A:
	for(i = 0; i < numberOfNodes; i++)
	{
		A[i] = contToBeAligned[i]->coordinate().x();
		A[numberOfNodes + i] = contToBeAligned[i]->coordinate().y();
		A[2*numberOfNodes + i] = contToBeAligned[i]->coordinate().z();
		A[3*numberOfNodes + i] = 1.0;
	}
	//do the rigid transformation and set it to part:
	clapack::dgemm_(&trans, &trans, &numberOfNodes, &dimension, &dimensionH, &alpha, A, &numberOfNodes, R, &dimensionH, &beta,
		b, &numberOfNodes);
	for(i = 0; i < numberOfNodes; i++)
	{
		contToBeAligned[i]->coordinate(g_Vector(b[i], b[numberOfNodes + i], b[2*numberOfNodes + i]));
	}

align_EXIT:
	delete [] work;
	delete [] A;
	delete [] b;
	delete [] S;
	delete [] U;
	delete [] VT;
	delete [] R;

	return returnval;
}

TriangleMesh * CorrespondTemplates::correspondSurface(TriangleMesh * meshToDeform, TriangleMesh * frame, char * exportFilename)
{
	int i, j;
	double dist, minDist, threshold;
	g_Vector globallyClosest, closestPoint;

	threshold = 20.0 * frame->getMeshResolution();

	for(i = 0; i < meshToDeform->getNumberOfNodes(); i++)
	{
		for(j = 0; j < frame->getNumberOfTriangles(); j++)
		{
		  // closestPoint = computeClosestPoint(meshToDeform->nodes()[i], frame->elements()[j]);
		  // dist = closestPoint.DistanceTo(meshToDeform->nodes()[i]->coordinate());
		        g_Vector closestPoint;
		        dist = computeClosestPoint(*meshToDeform->nodes()[i], *frame->elements()[j], closestPoint );
			dist = sqrt(dist);
			if((j == 0) || (dist < minDist))
			{
				minDist = dist;
				globallyClosest = closestPoint;
			}
		}

		if(minDist < threshold) meshToDeform->nodes()[i]->coordinate(globallyClosest);
	}

	exportMeshWrapper(exportFilename, meshToDeform);

	return meshToDeform;
}


void CorrespondTemplates::computeMappings(char * outfileMappingTet, char * outfileMappingSurf, g_Node * contactPoint)
{
	//Compute the mappings using barycentric coordinates
	findSurfaceNodes();
	maptet2surfID = surfaceNodes;
	// output tetmesh surface
	computeSurface( "tet_mesh_surf.wrl" );
	deformVolumetricMeshRigidly();

	int i, j, k, numNodesTet, numElemsTet, numNodesSurf, numElemsSurf;
	double dist, minDist, d, mind, distContactPoint, minDistContactPoint;
	numNodesTet = volumetricMeshNodes.numberOfItems();
	numElemsTet = trianglesOnSurface.size();
	numNodesSurf = rigidlyAlignedSurfaceTemplate->getNumberOfNodes();
	numElemsSurf = rigidlyAlignedSurfaceTemplate->getNumberOfTriangles();

	//First compute the mapping from the tetrahedral mesh to the surface mesh
	mappingIndicesTet.resize(3*numNodesTet);
	mappingWeightsTet.resize(3*numNodesTet);
	offsetsTet.resize(numNodesTet);

	for(i = 0; i < numNodesTet; i++)
	{
		if(surfaceNodes[i] == -1)
		{
			mappingIndicesTet[3*i] = mappingIndicesTet[3*i+1] = mappingIndicesTet[3*i+2] = -1;
			mappingWeightsTet[3*i] = mappingWeightsTet[3*i+1] = mappingWeightsTet[3*i+2] = -1;
		}
		else
		{
			for(j = 0; j < numElemsSurf; j++)
			{
			  // g_Vector closestPoint = computeClosestPoint(volumetricMeshNodes[i], rigidlyAlignedSurfaceTemplate->elements()[j]);
			  // dist = closestPoint.DistanceTo(volumetricMeshNodes[i]->coordinate());
		          g_Vector closestPoint;
			  dist = computeClosestPoint(*volumetricMeshNodes[i], *rigidlyAlignedSurfaceTemplate->elements()[j], closestPoint );
			  dist = sqrt(dist);

				if((j == 0) || (dist < minDist))
				{
					minDist = dist;
					
					mappingIndicesTet[3*i] = rigidlyAlignedSurfaceTemplate->elements()[j]->nodes()[0]->id()-1;
					mappingIndicesTet[3*i+1] = rigidlyAlignedSurfaceTemplate->elements()[j]->nodes()[1]->id()-1;
					mappingIndicesTet[3*i+2] = rigidlyAlignedSurfaceTemplate->elements()[j]->nodes()[2]->id()-1;

					//----------------------
					for(int n = 0; n < 3; n++){
						d = volumetricMeshNodes[i]->coordinate().DistanceTo(rigidlyAlignedSurfaceTemplate->elements()[j]->nodes()[n]->coordinate());
						if((n == 0) || (d < mind)){
							mind = d;
							maptet2surfID[i] = rigidlyAlignedSurfaceTemplate->elements()[j]->nodes()[n]->id()-1;
						}
					}
					//------------------------

					vector<double> barycentrics = computeBarycentrics(volumetricMeshNodes[i]->coordinate(), closestPoint, rigidlyAlignedSurfaceTemplate->elements()[j]);
					mappingWeightsTet[3*i] = barycentrics[0];
					mappingWeightsTet[3*i+1] = barycentrics[1];
					mappingWeightsTet[3*i+2] = barycentrics[2];
					offsetsTet[i] = barycentrics[3];
				}
			}
		}
	}

	//Export
	if(outfileMappingTet != NULL)
	{
#if 0
		FILE * fp = fopen(outfileMappingTet, "w");
		if(fp == NULL)
		{
			cout<<"Problem writing "<<outfileMappingTet<<endl;
			return;
		}
		for(i = 0; i < numNodesTet; i++)
			fprintf(fp, "%d %d %d %f %f %f %f\n", mappingIndicesTet[3*i], mappingIndicesTet[3*i+1], mappingIndicesTet[3*i+2], 
				mappingWeightsTet[3*i], mappingWeightsTet[3*i+1], mappingWeightsTet[3*i+2], offsetsTet[i]);
		fclose(fp);
#else
		// Make a surface mesh using the tetmesh surface plus the coordinates from the surface mesh
		computeSurface( outfileMappingTet, true );
#endif
	}

	//Second compute the mapping from the surface mesh to the tetrahedral mesh
	mappingIndicesSurf.resize(3*numNodesSurf);
	mappingWeightsSurf.resize(3*numNodesSurf);
	offsetsSurf.resize(numNodesSurf);
	contactTriangleId.resize(3);
	contactWeights.resize(3);
	contactoffset.resize(1);

	for(i = 0; i < numNodesSurf; i++)
	{
		for(j = 0; j < numElemsTet; j++)
		{
			g_Element elem;
			set<int>::iterator triangleIt;
			for(triangleIt = trianglesOnSurface[j].begin(); triangleIt != trianglesOnSurface[j].end(); triangleIt++)
				elem.node(volumetricMeshNodes[*triangleIt]);

			// g_Vector closestPoint = computeClosestPoint(rigidlyAlignedSurfaceTemplate->nodes()[i], &elem);
			// dist = closestPoint.DistanceTo(rigidlyAlignedSurfaceTemplate->nodes()[i]->coordinate());
		        g_Vector closestPoint;
		        dist = computeClosestPoint(*rigidlyAlignedSurfaceTemplate->nodes()[i], elem, closestPoint );
			dist = sqrt(dist);

			if((j == 0) || (dist < minDist))
			{
				minDist = dist;

				k = 0;
				for(triangleIt = trianglesOnSurface[j].begin(); triangleIt != trianglesOnSurface[j].end(); triangleIt++, k++)
					mappingIndicesSurf[3*i+k] = *triangleIt;

				vector<double> barycentrics = computeBarycentrics(rigidlyAlignedSurfaceTemplate->nodes()[i]->coordinate(), closestPoint, &elem);
				mappingWeightsSurf[3*i] = barycentrics[0];
				mappingWeightsSurf[3*i+1] = barycentrics[1];
				mappingWeightsSurf[3*i+2] = barycentrics[2];
				offsetsSurf[i] = barycentrics[3];
			}


			//---------------------------
			if(i == 0 && contactPoint != NULL){
			        g_Vector closestContactPoint;
				distContactPoint = computeClosestPoint(*contactPoint, elem, closestContactPoint );
				distContactPoint = sqrt(distContactPoint);
				// g_Vector closestContactPoint = computeClosestPoint(contactPoint, &elem);
				// distContactPoint = closestContactPoint.DistanceTo(contactPoint->coordinate());
				if((j == 0) || (distContactPoint < minDistContactPoint))
				{
					minDistContactPoint = distContactPoint;

					k = 0;
					for(triangleIt = trianglesOnSurface[j].begin(); triangleIt != trianglesOnSurface[j].end(); triangleIt++, k++)
						contactTriangleId[3*i+k] = *triangleIt;

					vector<double> barycentricsContactPoint = computeBarycentrics(contactPoint->coordinate(), closestContactPoint, &elem);
					contactWeights[3*i] = barycentricsContactPoint[0];
					contactWeights[3*i+1] = barycentricsContactPoint[1];
					contactWeights[3*i+2] = barycentricsContactPoint[2];
					contactoffset[i] = barycentricsContactPoint[3];
				}
			}
			//----------------------------
		}
	}

	//Export
	if(outfileMappingSurf != NULL)
	{
#if 0 
		FILE * fp = fopen(outfileMappingSurf, "w");
		if(fp == NULL)
		{
			cout<<"Problem writing "<<outfileMappingSurf<<endl;
			return;
		}
		for(i = 0; i < numNodesSurf; i++)
			fprintf(fp, "%d %d %d %f %f %f %f\n", mappingIndicesSurf[3*i], mappingIndicesSurf[3*i+1], mappingIndicesSurf[3*i+2], 
				mappingWeightsSurf[3*i], mappingWeightsSurf[3*i+1], mappingWeightsSurf[3*i+2], offsetsSurf[i]);
		fclose(fp);
#else
		// Make a surface mesh using the tetmesh surface plus the coordinates from the surface mesh
		TriangleMesh surfTet;
		for (int i=0; i<numNodesSurf; ++i ) {
		  // std::cerr << i << " " << std::endl;
		  g_Vector coord =  mappingWeightsSurf[3*i] * volumetricMeshNodes[mappingIndicesSurf[3*i]]->coordinate();
		  coord += mappingWeightsSurf[3*i+1] * volumetricMeshNodes[mappingIndicesSurf[3*i+1]]->coordinate(); 
		  coord += mappingWeightsSurf[3*i+2] * volumetricMeshNodes[mappingIndicesSurf[3*i+2]]->coordinate();
		  surfTet.node(new g_Node( coord ));
		}
		// std::cerr << std::endl;
		g_ElementContainer eleSurf = rigidlyAlignedSurfaceTemplate->elements();
		g_NodeContainer newNodes = surfTet.nodes();
		for ( g_ElementContainer::const_iterator fIt=eleSurf.begin(); fIt != eleSurf.end(); ++fIt ) {
		  g_Element *ele = new g_Element(); 
		  g_NodeContainer faceNodes = (*fIt)->nodes(); 
		  for ( g_NodeContainer::const_iterator nIt = faceNodes.begin(); 
			nIt != faceNodes.end(); ++nIt ) {
		    ele->node(newNodes[(*nIt)->id()-1]);
		  }
		  surfTet.element(ele);
		}
		exportMeshWrapper( outfileMappingSurf, &surfTet ); 
#endif


	}
}

void CorrespondTemplates::importMappings(char * infileMappingTet, char * infileMappingSurf, int numNodesTet, int numNodesSurf)
{
	int i, intRead1, intRead2, intRead3;
	float floatRead1, floatRead2, floatRead3, floatRead4;

	mappingIndicesTet.clear();
	mappingWeightsTet.clear();
	offsetsTet.clear();

	FILE * fp = fopen(infileMappingTet, "r");
	if(fp == NULL)
	{
		cout<<"Problem reading "<<infileMappingTet<<endl;
		return;
	}
	for(i = 0; i < numNodesTet; i++)
	{
		fscanf(fp, "%d %d %d %f %f %f %f ", &intRead1, &intRead2, &intRead3, &floatRead1, &floatRead2, &floatRead3, &floatRead4);
		mappingIndicesTet.push_back(intRead1);
		mappingIndicesTet.push_back(intRead2);
		mappingIndicesTet.push_back(intRead3);
		mappingWeightsTet.push_back(floatRead1);
		mappingWeightsTet.push_back(floatRead2);
		mappingWeightsTet.push_back(floatRead3);
		offsetsTet.push_back(floatRead4);
	}
	fclose(fp);

	mappingIndicesSurf.clear();
	mappingWeightsSurf.clear();
	offsetsSurf.clear();

	fp = fopen(infileMappingSurf, "r");
	if(fp == NULL)
	{
		cout<<"Problem reading "<<infileMappingSurf<<endl;
		return;
	}
	for(i = 0; i < numNodesSurf; i++)
	{
		fscanf(fp, "%d %d %d %f %f %f %f ", &intRead1, &intRead2, &intRead3, &floatRead1, &floatRead2, &floatRead3, &floatRead4);
		mappingIndicesSurf.push_back(intRead1);
		mappingIndicesSurf.push_back(intRead2);
		mappingIndicesSurf.push_back(intRead3);
		mappingWeightsSurf.push_back(floatRead1);
		mappingWeightsSurf.push_back(floatRead2);
		mappingWeightsSurf.push_back(floatRead3);
		offsetsSurf.push_back(floatRead4);
	}
	fclose(fp);
}

void CorrespondTemplates::findSurfaceNodes() 
{
	int i, numTetNodes;

	numTetNodes = volumetricMeshNodes.numberOfItems();
	surfaceNodes.resize(numTetNodes);
	for(i = 0; i < numTetNodes; i++) surfaceNodes[i] = -1;

	trianglesOnSurface.clear();

	multiset< set<int> > trianglesInTetrahedralMesh;
	for(i = 0; i < allTetrahedra.size(); i++)
	{
		set<int> triangle1; triangle1.insert(allTetrahedra[i][0]); triangle1.insert(allTetrahedra[i][1]); triangle1.insert(allTetrahedra[i][2]);
		set<int> triangle2; triangle2.insert(allTetrahedra[i][0]); triangle2.insert(allTetrahedra[i][1]); triangle2.insert(allTetrahedra[i][3]);
		set<int> triangle3; triangle3.insert(allTetrahedra[i][0]); triangle3.insert(allTetrahedra[i][2]); triangle3.insert(allTetrahedra[i][3]);
		set<int> triangle4; triangle4.insert(allTetrahedra[i][1]); triangle4.insert(allTetrahedra[i][2]); triangle4.insert(allTetrahedra[i][3]);

		trianglesInTetrahedralMesh.insert(triangle1);
		trianglesInTetrahedralMesh.insert(triangle2);
		trianglesInTetrahedralMesh.insert(triangle3);
		trianglesInTetrahedralMesh.insert(triangle4);
	}

	multiset< set<int> >::iterator it;
	for(it = trianglesInTetrahedralMesh.begin(); it != trianglesInTetrahedralMesh.end(); it++)
	{
		//If a triangle occurs only once, set all of its nodes to be surface nodes
		if(trianglesInTetrahedralMesh.count(*it) == 1)
		{
			set<int>::iterator triangleIt;
			for(triangleIt = (*it).begin(); triangleIt != (*it).end(); triangleIt++)
				surfaceNodes[*triangleIt] = 1;
			trianglesOnSurface.push_back(*it);
		}
	}
}

TriangleMesh * CorrespondTemplates::computeSurface(char * exportFilename, bool useMapping)
{
	int i, j, counter, numTetNodes;
	TriangleMesh * resultMesh = new TriangleMesh();
	vector<int> surfaceIds;

	//Find the surface:
	findSurfaceNodes();

	numTetNodes = volumetricMeshNodes.numberOfItems();
	surfaceIds.resize(numTetNodes);

	g_NodeContainer nodesSurf;
	if ( useMapping ) 
	   nodesSurf = rigidlyAlignedSurfaceTemplate->nodes();

	//Add the nodes:
	counter = 0;
	for(i = 0; i < numTetNodes; i++)
	{
		if(surfaceNodes[i] == -1) 
			surfaceIds[i] = -1;
		else 
		{
			surfaceIds[i] = counter;
			counter++;

			g_Vector coord;
			if ( useMapping ) {
			  coord =  mappingWeightsTet[3*i] * nodesSurf[mappingIndicesTet[3*i]]->coordinate(); 
			  coord += mappingWeightsTet[3*i+1] * nodesSurf[mappingIndicesTet[3*i+1]]->coordinate(); 
			  coord += mappingWeightsTet[3*i+2] * nodesSurf[mappingIndicesTet[3*i+2]]->coordinate();
			} else {
			  coord = volumetricMeshNodes[i]->coordinate();
			}
			g_Node * newNode = new g_Node(coord); 
			resultMesh->node(newNode);
		}
	}

	//Add the triangles (find the correct orientation):
	for(i = 0; i < trianglesOnSurface.size(); i++)
	{
		g_Element * elem = new g_Element();
	
		for(j = 0; j < allTetrahedra.size(); j++)
		{
			set<int> triangle; triangle.insert(allTetrahedra[j][0]); triangle.insert(allTetrahedra[j][1]); triangle.insert(allTetrahedra[j][2]);
			if(trianglesOnSurface[i] == triangle)
			{
				elem->node(resultMesh->nodes()[surfaceIds[allTetrahedra[j][0]]]);
				elem->node(resultMesh->nodes()[surfaceIds[allTetrahedra[j][2]]]);
				elem->node(resultMesh->nodes()[surfaceIds[allTetrahedra[j][1]]]);
				break;
			}
			triangle.clear(); triangle.insert(allTetrahedra[j][0]); triangle.insert(allTetrahedra[j][1]); triangle.insert(allTetrahedra[j][3]);
			if(trianglesOnSurface[i] == triangle)
			{
				elem->node(resultMesh->nodes()[surfaceIds[allTetrahedra[j][0]]]);
				elem->node(resultMesh->nodes()[surfaceIds[allTetrahedra[j][1]]]);
				elem->node(resultMesh->nodes()[surfaceIds[allTetrahedra[j][3]]]);
				break;
			}
			triangle.clear(); triangle.insert(allTetrahedra[j][0]); triangle.insert(allTetrahedra[j][2]); triangle.insert(allTetrahedra[j][3]);
			if(trianglesOnSurface[i] == triangle)
			{
				elem->node(resultMesh->nodes()[surfaceIds[allTetrahedra[j][0]]]);
				elem->node(resultMesh->nodes()[surfaceIds[allTetrahedra[j][3]]]);
				elem->node(resultMesh->nodes()[surfaceIds[allTetrahedra[j][2]]]);
				break;
			}
			triangle.clear(); triangle.insert(allTetrahedra[j][1]); triangle.insert(allTetrahedra[j][2]); triangle.insert(allTetrahedra[j][3]);
			if(trianglesOnSurface[i] == triangle)
			{
				elem->node(resultMesh->nodes()[surfaceIds[allTetrahedra[j][1]]]);
				elem->node(resultMesh->nodes()[surfaceIds[allTetrahedra[j][2]]]);
				elem->node(resultMesh->nodes()[surfaceIds[allTetrahedra[j][3]]]);
				break;
			}
		}

		resultMesh->element(elem);
	}

	//Export:
	if(exportFilename != NULL)
	{
		exportMeshWrapper(exportFilename, resultMesh);
	}

	return resultMesh;
}




/*
g_Vector CorrespondTemplates::computeClosestPoint(g_Node * vertex, g_Element * triangle)
{
	//Check the projection if it is in the triangle:
	g_Plane plane(triangle->nodes()[0]->coordinate(), triangle->nodes()[1]->coordinate(), triangle->nodes()[2]->coordinate());
	g_Vector projection = plane.project(vertex->coordinate());
	g_Vector pos = vertex->coordinate();
	vector<double> bCoord = computeBarycentrics(pos,projection,triangle);
	g_Vector closestPoint;

	bCoord[0] = (sgn(bCoord[0]) == 0)?1e-39:bCoord[0];
	bCoord[1] = (sgn(bCoord[1]) == 0)?1e-39:bCoord[1];
	bCoord[2] = (sgn(bCoord[2]) == 0)?1e-39:bCoord[2];

	if ( bCoord[0] > 0 && bCoord[1] > 0 && bCoord[2] > 0 ) {
	  // Point projects inside
	  closestPoint = projection;
	} else {
	  if ( bCoord[0] * bCoord[1] * bCoord[2] < 0 ) {
	    // Closest to an edge -- or vertex; closestPoint calculation likely not correct
	    if ( bCoord[0] < 0 ) {
	      float scale = 1.0 / (bCoord[1] + bCoord[2]);
	      closestPoint = scale * bCoord[1] *  triangle->nodes()[1]->coordinate() +
		scale * bCoord[2] *  triangle->nodes()[2]->coordinate();
	    } else {
	      if ( bCoord[1] < 0 ) {
		float scale = 1.0 / (bCoord[0] + bCoord[2]);
		closestPoint = scale * bCoord[0] *  triangle->nodes()[0]->coordinate() +
		  scale * bCoord[2] *  triangle->nodes()[2]->coordinate();
	      } else {
		if ( bCoord[2] < 0 ) {
		  float scale = 1.0 / (bCoord[1] + bCoord[0]);
		  closestPoint = scale * bCoord[1] *  triangle->nodes()[1]->coordinate() +
		    scale * bCoord[0] *  triangle->nodes()[0]->coordinate();
		}
	      }
	    }
	  } else {
	    // We are closest to a vertex
	    if ( bCoord[0] > 0 ) {
	      closestPoint = triangle->nodes()[0]->coordinate();
	    } else {
	      if ( bCoord[1] > 0 ) {
		closestPoint = triangle->nodes()[1]->coordinate();
	      } else {
		if ( bCoord[2] > 0 ) 
		  closestPoint = triangle->nodes()[0]->coordinate();
	      }
	    }
	  }
	}
	return closestPoint;
}
*/

/*
g_Vector CorrespondTemplates::computeClosestPoint(g_Node * vertex, g_Element * triangle)
{
	int i, turn1, turn2, turn3;
	double dot, dist, minDist;
	g_Vector closestPoint, projection;

	//Check the three vertices
	for(i = 0; i < 3; i++)
	{
		dist = vertex->coordinate().DistanceTo(triangle->nodes()[i]->coordinate());
		if((i == 0) || (dist < minDist))
		{
			minDist = dist;
			closestPoint = triangle->nodes()[i]->coordinate();
		}
	}	

	//Check the projection if it is in the triangle:
	g_Plane plane(triangle->nodes()[0]->coordinate(), triangle->nodes()[1]->coordinate(), triangle->nodes()[2]->coordinate());
	projection = plane.project(vertex->coordinate());

	dot = (triangle->nodes()[1]->coordinate()-triangle->nodes()[0]->coordinate()).Cross(projection-triangle->nodes()[1]->coordinate()).Dot(plane.normal());
	if(dot < 0) turn1 = -1;
	else turn1 = 1;
	dot = (triangle->nodes()[2]->coordinate()-triangle->nodes()[1]->coordinate()).Cross(projection-triangle->nodes()[2]->coordinate()).Dot(plane.normal());
	if(dot < 0) turn2 = -1;
	else turn2 = 1;
	dot = (triangle->nodes()[0]->coordinate()-triangle->nodes()[2]->coordinate()).Cross(projection-triangle->nodes()[0]->coordinate()).Dot(plane.normal());
	if(dot < 0) turn3 = -1;
	else turn3 = 1;
	if(turn1 == turn2 && turn2 == turn3) 
	{
		dist = vertex->coordinate().DistanceTo(projection);
		if(dist < minDist)
		{
			minDist = dist;
			closestPoint = projection;
		}
	} 

	return closestPoint;
}
*/


vector<double> CorrespondTemplates::computeBarycentrics( g_Vector& coordinate, 
							 g_Vector& coordinateOnSurface, 
							 g_Element *  triangle)
{
	double areaABC, areaPBC, areaPCA;
	vector<double> barycentrics;

	g_Vector normal = (triangle->nodes()[1]->coordinate()-triangle->nodes()[0]->coordinate()).Cross(triangle->nodes()[2]->coordinate()-triangle->nodes()[0]->coordinate());
	normal.Normalize();

	// Compute twice area of triangle ABC
	areaABC = normal.Dot((triangle->nodes()[1]->coordinate()-triangle->nodes()[0]->coordinate()).Cross(triangle->nodes()[2]->coordinate()-triangle->nodes()[0]->coordinate()));

	// Compute barycentrics
	areaPBC = normal.Dot((triangle->nodes()[1]->coordinate()-coordinateOnSurface).Cross(triangle->nodes()[2]->coordinate()-coordinateOnSurface));
	barycentrics.push_back(areaPBC / areaABC);

	areaPCA = normal.Dot((triangle->nodes()[2]->coordinate()-coordinateOnSurface).Cross(triangle->nodes()[0]->coordinate()-coordinateOnSurface));
	barycentrics.push_back(areaPCA / areaABC);

	barycentrics.push_back(1.0 - barycentrics[0] - barycentrics[1]);

	// Compute the normal offset
	g_Vector projectedPoint = barycentrics[0]*triangle->nodes()[0]->coordinate() + barycentrics[1]*triangle->nodes()[1]->coordinate() + barycentrics[2]*triangle->nodes()[2]->coordinate();
	projectedPoint = projectedPoint - coordinate;
	barycentrics.push_back(normal.Dot(projectedPoint));

	return barycentrics;
}
