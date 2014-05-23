#include "FineFittingTracking.h"

FineFittingTracking::FineFittingTracking()
{
	smallestMeshNum = 1000;
	templateMesh = NULL;
	frameMesh = NULL;
	smallestMesh = NULL;
	hierarchy.clear();
	transformationParameters = NULL;
	kdInitialized = false;
	allEdges.clear();
	transformationParameterInitialization = NULL;
	resultMesh = NULL;
	clusterIndicesInHierarchy.clear();
	isFeatureInHierarchy.clear();
}

FineFittingTracking::~FineFittingTracking()
{
	int i;
	if(smallestMesh != NULL) 
	{
		delete smallestMesh;
		smallestMesh = NULL;
	}
	if(hierarchy.size() > 0)
	{
		for(i = 0; i < (int)hierarchy.size(); i++)
		{
			delete hierarchy[i];
		}
		hierarchy.clear();
	}
	if(transformationParameters != NULL)
	{
		delete [] transformationParameters;
		transformationParameters = NULL;
	}
	if(kdInitialized)
	{
		delete kdTree;
		annDeallocPts(dataPts);
		annClose();	
		kdInitialized = false;
	}
	for(i = 0; i < (int)allEdges.numberOfItems(); i++) delete allEdges[i];
	allEdges.clear();
	if(transformationParameterInitialization != NULL)
	{
		delete [] transformationParameterInitialization;
		transformationParameterInitialization = NULL;
	}
	if(resultMesh != NULL)
	{
		delete resultMesh;
		resultMesh = NULL;
	}
}

vector<int> FineFittingTracking::getNearestNeighbors(bool symmetric)
{
	vector<int> * isFeature;
	if(isFeatureInHierarchy.size() > 0) 
		isFeature = &(isFeatureInHierarchy[0]);
	else isFeature = NULL;

	computeNearestNeighbors(templateMesh, meshNormals[0], isFeature, false, symmetric); 

	return nearestNeighbors;
}

void FineFittingTracking::exportTransformationParameters(char * filename)
{
	if(templateMesh == NULL) return;

	FILE * fp = fopen(filename, "w");
	for (int i = 0; i < (int)templateMesh->getNumberOfNodes()* 12; i++)
		fprintf(fp, "%f ", transformationParameters[i]);
	fclose(fp);
}

void FineFittingTracking::setTemplate(TriangleMesh * mesh)
{
	templateMesh = mesh;
	largestDist = MAX_DIST * templateMesh->getMeshResolution(); 


/*/FOR FRAME AFTER A FEW FRAMES WERE COMPUTED ALREADY:
int numNodes = templateMesh->getNumberOfNodes();
int numLevels = (int) ceil(log((double)numNodes / (double)smallestMeshNum))+1;
if(hierarchy.size() == 0) computeHierarchy(numLevels);
numLevels = (int)hierarchy.size();
//END */
}

void FineFittingTracking::setInitialTransformationParameters(double * param)
{
	int i, dimension;

	if(templateMesh == NULL) return;

	dimension = templateMesh->nodes().numberOfItems() * 12;
	if(transformationParameterInitialization != NULL) delete [] transformationParameterInitialization;
	transformationParameterInitialization = new double[dimension];

	for(i = 0; i < dimension; i++) transformationParameterInitialization[i] = param[i];
}

void FineFittingTracking::resetTemplateCoo(g_NodeContainer nodes, char * filename)
{
	int i, j;

	g_Part * currentMesh, * previousMesh;

	//Set new coordinates to the template
	if(templateMesh == NULL) return;
	else if(templateMesh->getNumberOfNodes() > (int)nodes.numberOfItems()) return;
	else
	{
		for(i = 0; i < templateMesh->getNumberOfNodes(); i++)
			templateMesh->nodes()[i]->coordinate(nodes[i]->coordinate());
	}

	//Set new coordinates to the hierarchy
	if(hierarchy.size() > 0)
	{
		for(i = 0; i < (int)hierarchy.size(); i++)
		{
			currentMesh = hierarchy[i]->getOriginalMesh();

			if(i == 0)
			{
				for(j = 0; j < (int)currentMesh->nodes().numberOfItems(); j++)
					currentMesh->nodes()[j]->coordinate(nodes[j]->coordinate());
			}

			else
			{
				vector<int> indices = hierarchy[i-1]->getAdditionalIndices();
				for(j = 0; j < (int)currentMesh->nodes().numberOfItems(); j++)
					currentMesh->nodes()[j]->coordinate(previousMesh->nodes()[indices[j]]->coordinate());
			}

			if(filename != NULL)
			{
				char buffer[200];
				sprintf(buffer, "%s_level%d.wrl", filename, i);
				TriangleMesh mesh(*currentMesh);
				exportMeshWrapper( buffer, &mesh );
			}

			previousMesh = currentMesh;
		}
	}
}

void FineFittingTracking::setFrame(TriangleMesh * mesh)
{
	frameMesh = mesh;
	//Initialize the kd-tree:
	if(kdInitialized)
	{
		delete kdTree;
		annDeallocPts(dataPts);
		annClose();	
		kdInitialized = false;
	}

	dataPts = annAllocPts(frameMesh->getNumberOfNodes(), 3);
	for(int i = 0; i < frameMesh->getNumberOfNodes(); i++)
	{
		dataPts[i][0] = frameMesh->nodes()[i]->coordinate().x();
		dataPts[i][1] = frameMesh->nodes()[i]->coordinate().y();
		dataPts[i][2] = frameMesh->nodes()[i]->coordinate().z();
	}
	kdTree = new ANNkd_tree(dataPts, frameMesh->getNumberOfNodes(), 3);
	kdInitialized = true;
}

void FineFittingTracking::setClustersFromFile(char * filename) 
{
	if(templateMesh == NULL) return;

	int i, j, k, l, numParts, numPoints, index;
	bool found;
	vector< vector<int> > clusters;

	//Read the clusters:
	FILE * fp = fopen(filename, "r");
	if(fp == NULL) return;

	fscanf(fp, "Number parts %d ", &numParts);
	clusters.resize(numParts);
	
	for(i = 0; i < numParts; i++)
	{
		fscanf(fp, "Part %d: %d ", &index, &numPoints);
		clusters[i].clear();
		for(j = 0; j < numPoints; j++)
		{
			fscanf(fp, "%d ", &index);
			clusters[i].push_back(index);
		}
	}

	fclose(fp);

	//Compute the cluster neighbors for each point:
	clusterIndicesInHierarchy.clear();
	vector< vector<int> > clusterIndicesLevel0(templateMesh->getNumberOfNodes());
	
	for(i = 0; i < templateMesh->getNumberOfNodes(); i++)
	{
		set<int> setForPoint;
		for(j = 0; j < numParts; j++)
		{
			found  = false;
			for(k = 0; k < clusters[j].size(); k++)
			{
				if(i == clusters[j][k])
				{
					found = true;
					break;
				}
			}
			if(found)
			{
				for(k = 0; k < clusters[j].size(); k++) setForPoint.insert(clusters[j][k]);
			}
		}

		set<int>::iterator it;
		for(it = setForPoint.begin(); it != setForPoint.end(); it++)
			clusterIndicesLevel0[i].push_back(*it);
	}
	clusterIndicesInHierarchy.push_back(clusterIndicesLevel0);

	//Transfer to all hierarchy levels:
	for(i = 0; i < hierarchy.size(); i++)
	{
		vector<int> indices = hierarchy[i]->getAdditionalIndices();
		vector< vector<int> > clusterIndicesCurrentLevel(indices.size());

		for(j = 0; j < indices.size(); j++)
		{
			for(k = 0; k < clusterIndicesInHierarchy[i][indices[j]].size(); k++)
			{
				for(l = 0; l < indices.size(); l++) 
				{
					if(indices[l] == clusterIndicesInHierarchy[i][indices[j]][k])
						clusterIndicesCurrentLevel[j].push_back(l);
				}
			}
		}

		clusterIndicesInHierarchy.push_back(clusterIndicesCurrentLevel);
	}
}

void FineFittingTracking::setFeatures(char * featureFile)
{
	if((templateMesh == NULL) || (frameMesh == NULL)) return;

	int i, j, numNodes, numSamples, retVal, minId1, minId2, numFeatures, numAcceptedFeatures;
	double dist, minDist;
	g_Vector currentPoint;

	numNodes = templateMesh->getNumberOfNodes();
	numSamples = frameMesh->getNumberOfNodes();

	FILE * fp = fopen(featureFile, "rb");
	if(fp == NULL)
	{
		cout<<"Could not find file "<<featureFile<<endl;
		return;
	}

	isFeatureInHierarchy.clear();

	vector<int> isFeatureInTemplate (numNodes);
	for(i = 0; i < numNodes; i++) isFeatureInTemplate[i] = -1;
	
	//read the file
	float * floatBuffer = new float[6];
		
	retVal = 1;
	numFeatures = 0;
	numAcceptedFeatures = 0;
	largestDist = MAX_DIST*templateMesh->getMeshResolution(); 
	while(retVal != EOF)
	{
		if(fread(floatBuffer, sizeof(float), 6, fp) != 6) 
		{
			retVal = EOF;
			continue;
		}

		numFeatures++;

		//Find the point on the current template corresponding to the feature
		currentPoint.Set(floatBuffer[0], floatBuffer[1], floatBuffer[2]);
		for(i = 0; i < numNodes; i++)
		{
			dist = templateMesh->nodes()[i]->coordinate().DistanceTo(currentPoint);
			if((i == 0) || (dist < minDist))
			{
				minDist = dist;
				minId1 = i;
			}
		}

		//Only accept the feature if it is close to the template
		if(minDist < largestDist)
		{
			//Find the point on the data corresponding to the feature
			currentPoint.Set(floatBuffer[3], floatBuffer[4], floatBuffer[5]);
			for(i = 0; i < numSamples; i++)
			{
				dist = frameMesh->nodes()[i]->coordinate().DistanceTo(currentPoint);
				if((i == 0) || (dist < minDist))
				{
					minDist = dist;
					minId2 = i;
				}
			}
			if(minDist < largestDist)
			{
				isFeatureInTemplate[minId1] = minId2;
				numAcceptedFeatures++;
			}
		}
	}

	fclose(fp);

	cout<<"Number features "<<numFeatures<<" Number accepted "<<numAcceptedFeatures<<endl;

	isFeatureInHierarchy.push_back(isFeatureInTemplate);

	//Transfer to all hierarchy levels
	for(i = 0; i < hierarchy.size(); i++)
	{
		vector<int> indices = hierarchy[i]->getAllIndices();
		vector<int> isFeature(indices.size());

		for(j = 0; j < indices.size(); j++)
			isFeature[j] = isFeatureInHierarchy[i][indices[j]];

		isFeatureInHierarchy.push_back(isFeature);
	}

	delete [] floatBuffer;
}

vector<int> FineFittingTracking::getFeatures(int index)
{
	if(isFeatureInHierarchy.size() < index + 1)
	{
		vector<int> features;
		return features;
	}
	else return isFeatureInHierarchy[index];
}

void FineFittingTracking::setMarkers(char * markersToFirstFrame, int numMarkers)
{
	if((templateMesh == NULL) || (frameMesh == NULL)) return;

	int i, j, numNodes, numSamples, minId1, minId2;
	float floatRead [6];
	double dist, minDist;
	g_Vector currentPoint;

	numNodes = templateMesh->getNumberOfNodes();
	numSamples = frameMesh->getNumberOfNodes();

	FILE * fp = fopen(markersToFirstFrame, "r");
	if(fp == NULL)
	{
		cout<<"Could not find file "<<markersToFirstFrame<<endl;
		return;
	}

	isFeatureInHierarchy.clear();

	vector<int> isFeatureInTemplate (numNodes);
	for(i = 0; i < numNodes; i++) isFeatureInTemplate[i] = -1;
	
	//read the file		
	for(i = 0; i < numMarkers; i++)
	{
		fscanf(fp, "%f %f %f %f %f %f ", &(floatRead[0]), &(floatRead[1]), &(floatRead[2]), &(floatRead[3]), &(floatRead[4]), &(floatRead[5]));

		//Find the point on the current template corresponding to the feature
		currentPoint.Set(floatRead[0], floatRead[1], floatRead[2]);
		for(j = 0; j < numNodes; j++)
		{
			dist = templateMesh->nodes()[j]->coordinate().DistanceTo(currentPoint);
			if((j == 0) || (dist < minDist))
			{
				minDist = dist;
				minId1 = j;
			}
		}

		//Find the point on the data corresponding to the feature
		currentPoint.Set(floatRead[3], floatRead[4], floatRead[5]);
		for(j = 0; j < numSamples; j++)
		{
			dist = frameMesh->nodes()[j]->coordinate().DistanceTo(currentPoint);
			if((j == 0) || (dist < minDist))
			{
				minDist = dist;
				minId2 = j;
			}
		}
		isFeatureInTemplate[minId1] = minId2;
	}
	fclose(fp);

	isFeatureInHierarchy.push_back(isFeatureInTemplate);

	//Transfer to all hierarchy levels
	for(i = 0; i < hierarchy.size(); i++)
	{
		vector<int> indices = hierarchy[i]->getAllIndices();
		vector<int> isFeature(indices.size());

		for(j = 0; j < indices.size(); j++)
			isFeature[j] = isFeatureInHierarchy[i][indices[j]];

		isFeatureInHierarchy.push_back(isFeature);
	}
}

void FineFittingTracking::unsetMarkers()
{
	isFeatureInHierarchy.clear();
}

TriangleMesh * FineFittingTracking::computeDeformedTemplate(char * filename, bool rigidAligment, bool * selfIntersectReturn, bool onlyUseMarkers)
{
	int i, j, k, numNodes, numLevels, counter, resolutionIndex;
	g_Part * currentMesh;

	//First, align rigidly:
	if(rigidAligment) 
	{
		alignRigidlyToFrame(onlyUseMarkers);

		if(filename != NULL)
		{
			char buffer[200];
			sprintf(buffer, "%s_after_rigid_alignment.wrl", filename);
			exportMeshWrapper(buffer, templateMesh);
		}
	}

	//If not set, initialize to identity:
	if(transformationParameterInitialization == NULL)
	{
		transformationParameterInitialization = new double[templateMesh->nodes().numberOfItems() * 12];
		for(i = 0; i < (int)templateMesh->nodes().numberOfItems() * 12; i++) transformationParameterInitialization[i] = 0;
		for(i = 0; i < (int)templateMesh->nodes().numberOfItems(); i++)
		{
			transformationParameterInitialization[12*i] = 1.0;
			transformationParameterInitialization[12*i+4] = 1.0;
			transformationParameterInitialization[12*i+8] = 1.0;
		}
	}

	//Second, compute multi-resolution non-rigid alignment:
	numNodes = templateMesh->getNumberOfNodes();
	numLevels = (int) ceil(log((double)numNodes / (double)smallestMeshNum))+1;

	if(hierarchy.size() == 0) computeHierarchy(numLevels);
	numLevels = (int)hierarchy.size();

	currentMesh = smallestMesh;

	if(transformationParameters != NULL) delete [] transformationParameters;
	transformationParameters = new double[currentMesh->nodes().numberOfItems() * 12];
	for(i = 0; i < (int)currentMesh->nodes().numberOfItems(); i++)
	{
		resolutionIndex = i;
		for(j = numLevels-1; j >= 0; j--) resolutionIndex = hierarchy[j]->getAdditionalIndices()[resolutionIndex];
		for(j = 0; j < 12; j++) transformationParameters[12*i+j] = transformationParameterInitialization[12*resolutionIndex+j];
	}

	for(i = numLevels-1; i >= 0; i--)
	{
		TriangleMesh tempMesh(*currentMesh);
		largestDist = MAX_DIST * tempMesh.getMeshResolution();

		vector<int> * isFeature;
		if(isFeatureInHierarchy.size() > 0) 
			isFeature = &(isFeatureInHierarchy[i+1]);
		else isFeature = NULL;
		
		vector< vector<int> > clusterIds;
		if(clusterIndicesInHierarchy.size() > 0) clusterIds = clusterIndicesInHierarchy[i+1];
		fitToData(currentMesh, meshNormals[i+1], clusterIds, isFeature, onlyUseMarkers);

		if(filename != NULL)
		{
			char buffer[200];
			TriangleMesh test(*currentMesh);	
			for(j = 0; j < (int)test.nodes().numberOfItems(); j++)
				test.nodes()[j]->coordinate(computeTransformedPoint(&test, j));
			if(clusterIndicesInHierarchy.size() > 0)
			{
				srand((int)time(0));
				svector<Color> colors(test.getNumberOfNodes());
				for(j = 0; j < clusterIndicesInHierarchy[i+1].size(); j++)
				{
					Color col((float)rand()/(float)RAND_MAX, (float)rand()/(float)RAND_MAX, (float)rand()/(float)RAND_MAX); 
					for(k = 0; k < clusterIndicesInHierarchy[i+1][j].size(); k++)
						colors[clusterIndicesInHierarchy[i+1][j][k]] = col;
				}
				test.setColorToMesh(colors);
			}
			sprintf(buffer, "%s_after_deformation_level_%d.wrl", filename, i);
			exportMeshWrapper( buffer, &test );
		}

		//Update both the current mesh and the transformation parameters:
		currentMesh = hierarchy[i]->getOriginalMesh();

		double * newTrafoParams = new double[currentMesh->nodes().numberOfItems() * 12];
		vector<int> indices = hierarchy[i]->getAdditionalIndices();
		counter = 0;
		for(j = 0; j < (int)currentMesh->nodes().numberOfItems(); j++)
		{
			if(counter < indices.size())
			{
				if(indices[counter] == j)
				{
					for(k = 0; k < 12; k++) newTrafoParams[12*j+k] = transformationParameters[12*counter+k];
					counter++;
				}
				else
				{
					resolutionIndex = j;
					for(k = i-1; k >= 0; k--) resolutionIndex = hierarchy[k]->getAdditionalIndices()[resolutionIndex];
					for(k = 0; k < 12; k++) newTrafoParams[12*j+k] = transformationParameterInitialization[12*resolutionIndex+k];
				}
			}
			else
			{
				resolutionIndex = j;
				for(k = i-1; k >= 0; k--) resolutionIndex = hierarchy[k]->getAdditionalIndices()[resolutionIndex];
				for(k = 0; k < 12; k++) newTrafoParams[12*j+k] = transformationParameterInitialization[12*resolutionIndex+k];
			}
		}

		delete [] transformationParameters;
		transformationParameters = newTrafoParams;
		smoothTransformations(i, filename);

		if(filename != NULL)
		{
			char buffer[200];
			TriangleMesh test(*currentMesh);	
			for(j = 0; j < (int)test.nodes().numberOfItems(); j++)
				test.nodes()[j]->coordinate(computeTransformedPoint(&test, j));
			sprintf(buffer, "%s_after_expansion_level_%d.wrl", filename, i);
			exportMeshWrapper( buffer, &test );
		}

		cout<<"Finished level "<<i<<endl;
	}

	vector<int> * isFeature;
	if(isFeatureInHierarchy.size() > 0) 
		isFeature = &(isFeatureInHierarchy[0]);
	else isFeature = NULL;

	vector< vector<int> > clusterIds;
	if(clusterIndicesInHierarchy.size() > 0) clusterIds = clusterIndicesInHierarchy[0];
	
	fitToData(currentMesh, meshNormals[0], clusterIds, isFeature, onlyUseMarkers);

	if(filename != NULL)
	{
		char buffer[200];
		TriangleMesh test(*currentMesh);	
		for(j = 0; j < (int)test.nodes().numberOfItems(); j++)
			test.nodes()[j]->coordinate(computeTransformedPoint(&test, j));
		sprintf(buffer, "%s_before_post_processing.wrl", filename);
		exportMeshWrapper( buffer, &test );
	}

	//Post-process the mesh: not required if we use FEM
	//postProcessGeometry(isFeature);

	//NOW: Return the mesh to use with FEM
	if(selfIntersectReturn != NULL) 
	{
		//Test for self-intersections:
		TriangleMesh checkIntersectionMesh(*currentMesh);	
		for(i = 0; i < (int)checkIntersectionMesh.nodes().numberOfItems(); i++)
			checkIntersectionMesh.nodes()[i]->coordinate(computeTransformedPoint(&checkIntersectionMesh, i));
		MinimumDistance dist(&(checkIntersectionMesh.elements()), &(checkIntersectionMesh.nodes()));
		bool selfIntersect = dist.selfIntersection(0);

		*selfIntersectReturn = selfIntersect;
	}

	// Return the result:
	if(resultMesh != NULL) delete resultMesh;
	resultMesh = new TriangleMesh(*currentMesh);

	//Set the result (modify the mesh):
	for(i = 0; i < (int)resultMesh->nodes().numberOfItems(); i++)
		resultMesh->nodes()[i]->coordinate(computeTransformedPoint(currentMesh, i));

	return resultMesh;
}

void FineFittingTracking::adjustTransformations(g_Part * targetMesh, vector<int> indicesModifiedWithFEM, bool * selfIntersectReturn)
{
	int i, j;

	g_Part * mesh = templateMesh;

	vector<int> filteredIndices;
	for(i = 0; i < (int)indicesModifiedWithFEM.size(); i++)
	{
		if(mesh->nodes()[indicesModifiedWithFEM[i]]->coordinate().DistanceTo(targetMesh->nodes()[indicesModifiedWithFEM[i]]->coordinate()) < largestDist)
			filteredIndices.push_back(indicesModifiedWithFEM[i]);
	}

	int dimension = (int)indicesModifiedWithFEM.size() * 12;

	vnl_vector< double > x (dimension, 0.0);
	for(i = 0; i < (int)indicesModifiedWithFEM.size(); i++)
	{
		for(j = 0; j < 12; j++) x[12*i+j] = transformationParameters[12*indicesModifiedWithFEM[i]+j];
	}

	// Set the weights used here:
	nnWeight = 10.0;
	rigidWeight = 1.0;
	accessWeight = 10.0;
	mdWeight = 1.0;
		
	//Compute the access indices:
	if(accessWeight > 0) computeAccessNeighbors(mesh, &indicesModifiedWithFEM);
	if(mdWeight > 0) computeMdNeighbors(mesh, &indicesModifiedWithFEM);

	//perform the nonlinear optimization
	AdjustTransformationCostFunction costFct(this, mesh, targetMesh, indicesModifiedWithFEM, filteredIndices, dimension);
	vnl_lbfgsb * minimizer = new vnl_lbfgsb(costFct);
	minimizer->set_cost_function_convergence_factor ( 10000000 ); 
	minimizer->set_projected_gradient_tolerance ( 0.00001 ); 
	minimizer->set_max_function_evals(MAX_NUM_ITER);
	minimizer->minimize( x );

	//Copy the transformations:
	for(i = 0; i < (int)indicesModifiedWithFEM.size(); i++)
	{
		for(j = 0; j < 12; j++) transformationParameters[12*indicesModifiedWithFEM[i]+j] = x[12*i+j];
	}

	delete minimizer;

	//Change the coordinates of the input
	for(i = 0; i < (int)mesh->nodes().numberOfItems(); i++)
		targetMesh->nodes()[i]->coordinate(computeTransformedPoint(mesh, i));

	//Check for self-intersections
	if(selfIntersectReturn != NULL) 
	{
		//Test for self-intersections:
		MinimumDistance dist(&(targetMesh->elements()), &(targetMesh->nodes()));
		bool selfIntersect = dist.selfIntersection(0);

		*selfIntersectReturn = selfIntersect;
	}
}

void FineFittingTracking::AdjustTransformationCostFunction::compute(vnl_vector< double > const &x, double *f, vnl_vector< double > *g)
{
	int i, j, dimension;

	dimension = (int)indices.size() * 12;

	double * helpGrad1 = new double[mesh->nodes().numberOfItems() * 12];
	double * helpGrad2 = new double[mesh->nodes().numberOfItems() * 12];

	// Update the transformation parameters:
	for(i = 0; i < (int)indices.size(); i++)
	{
		for(j = 0; j < 12; j++) tracker->transformationParameters[12*indices[i]+j] = x[12*i+j];
	} 

	// Compute Energy:
	*f = 0;
	if(tracker->accessWeight > 0) *f += tracker->accessWeight * tracker->accessEnergy(mesh);
	if(tracker->mdWeight > 0) *f += tracker->mdWeight * tracker->mdEnergy(mesh);
	if(tracker->rigidWeight > 0) *f += tracker->rigidWeight * tracker->rigidEnergy(mesh, &indices);
	if(tracker->nnWeight > 0) *f += tracker->nnWeight * tracker->distanceResolutionEnergy(mesh, targetMesh, indicesLowRes);

	// Compute Gradient:
	for(i = 0; i < dimension; i++) (*g)[i] = 0;
	for(i = 0; i < mesh->nodes().numberOfItems() * 12; i++) helpGrad2[i] = 0;

	if(tracker->accessWeight > 0)
	{
		tracker->accessGradient(mesh, helpGrad1);
		for(i = 0; i < mesh->nodes().numberOfItems() * 12; i++) 
			helpGrad2[i] += tracker->accessWeight * helpGrad1[i];
	}
	if(tracker->mdWeight > 0)
	{
		tracker->mdGradient(mesh, helpGrad1);
		for(i = 0; i < mesh->nodes().numberOfItems() * 12; i++) 
			helpGrad2[i] += tracker->mdWeight * helpGrad1[i];
	}
	if(tracker->rigidWeight > 0)
	{
		tracker->rigidGradient(mesh, helpGrad1);
		for(i = 0; i < mesh->nodes().numberOfItems() * 12; i++)
			helpGrad2[i] += tracker->rigidWeight * helpGrad1[i];
	}
	if(tracker->nnWeight > 0) 
	{
		tracker->distanceResolutionGradient(mesh, targetMesh, indicesLowRes, helpGrad1);
		for(i = 0; i < (int)indicesLowRes.size(); i++)
		{
			for(j = 0; j < 12; j++) helpGrad2[12*indicesLowRes[i]+j] += tracker->nnWeight * helpGrad1[12*i+j];
		}
	}

	for(i = 0; i < (int)indices.size(); i++)
	{
		for(j = 0; j < 12; j++) (*g)[12*i+j] = helpGrad2[12*indices[i]+j];
	}

	delete [] helpGrad1;
	delete [] helpGrad2;
}

void FineFittingTracking::postProcessGeometry(vector<int> * isFeature)
{
	int i, j, numUpdateIterations;

	g_Part * mesh = hierarchy[0]->getOriginalMesh();
	vector<int> invisibleIndices;

	//Compute the points without nearest neighbors:
	computeNearestNeighbors(mesh, meshNormals[0], isFeature, false);
	for(i = 0; i < (int)mesh->nodes().numberOfItems(); i++)
	{
		if(nearestNeighbors[i] == -1) invisibleIndices.push_back(i);
	}

	int dimension = (int)(invisibleIndices.size()) * 12;

	vnl_vector<double> x(dimension, 0.0);
	for(i = 0; i < invisibleIndices.size(); i++)
	{
		for(j = 0; j < 12; j++) x[12*i+j] = transformationParameters[12*invisibleIndices[i]+j];
	}

	// Set the weights used here:
	accessWeight = 1.0;
	rigidWeight = 1.0;
	nnWeight = 100.0;

	// Compute the access indices
	computeAccessNeighbors(mesh, &invisibleIndices);

	// Do multiple energy minimizations using the mean value coordinates
	numUpdateIterations = 50;
	for(i = 0; i < numUpdateIterations; i++)
	{
		// Compute the deformed mesh:
		g_Part deformedMesh(*mesh);
		for(j = 0; j < (int)mesh->nodes().numberOfItems(); j++)
			deformedMesh.nodes()[j]->coordinate(computeTransformedPoint(mesh, j));

		// Compute the mean-value coordinates:
		g_Part meanValueCooMesh(deformedMesh);
		for(j = 0; j < (int)mesh->nodes().numberOfItems(); j++)
			meanValueCooMesh.nodes()[j]->coordinate(0.9*deformedMesh.nodes()[j]->coordinate() + 0.1*hierarchy[0]->getMeanValueCoordinate(&deformedMesh, j));

		//perform the nonlinear optimization
		PostProcessCostFunction costFct(this, mesh, &meanValueCooMesh, invisibleIndices, dimension);
		vnl_lbfgsb * minimizer = new vnl_lbfgsb(costFct);
		minimizer->set_cost_function_convergence_factor ( 10000000 );
		minimizer->set_projected_gradient_tolerance ( 0.00001 ); 
		minimizer->set_max_function_evals(MAX_NUM_ITER);
		minimizer->minimize( x );

		delete minimizer;
	}

	//Copy the transformations:
	for(i = 0; i < invisibleIndices.size(); i++)
	{
		for(j = 0; j < 12; j++) transformationParameters[12*invisibleIndices[i]+j] = x[12*i+j];
	}
}

void FineFittingTracking::PostProcessCostFunction::compute(vnl_vector< double > const &x, double *f, vnl_vector< double > *g)
{
	int j, k, dimension;

	dimension = (invisibleIndices.size()) * 12;

	double * helpGrad1 = new double[mesh->nodes().numberOfItems() * 12];
	double * helpGrad2 = new double[mesh->nodes().numberOfItems() * 12];

	// Update the transformation parameters:
	for(j = 0; j < invisibleIndices.size(); j++)
	{
		for(k = 0; k < 12; k++) tracker->transformationParameters[12*invisibleIndices[j]+k] = x[12*j+k];
	}

	// Compute Energy:
	*f = 0;
	*f += tracker->accessWeight * tracker->accessEnergy(mesh);
	*f += tracker->rigidWeight * tracker->rigidEnergy(mesh, &invisibleIndices);
	*f += tracker->nnWeight * tracker->distanceResolutionEnergy(mesh, meanValueCooMesh, invisibleIndices);

	// Compute Gradient:
	for(j = 0; j < dimension; j++) 
		(*g)[j] = 0;
	tracker->accessGradient(mesh, helpGrad1);
	for(j = 0; j < mesh->nodes().numberOfItems() * 12; j++) 
		helpGrad2[j] = tracker->accessWeight * helpGrad1[j];
	tracker->rigidGradient(mesh, helpGrad1);
	for(j = 0; j < mesh->nodes().numberOfItems() * 12; j++)
		helpGrad2[j] += tracker->rigidWeight * helpGrad1[j];
	tracker->distanceResolutionGradient(mesh, meanValueCooMesh, invisibleIndices, helpGrad1);
	for(j = 0; j < invisibleIndices.size(); j++)
	{
		for(k = 0; k < 12; k++) helpGrad2[12*invisibleIndices[j]+k] += tracker->nnWeight * helpGrad1[12*j+k];
	}

	for(j = 0; j < invisibleIndices.size(); j++)
	{
		for(k = 0; k < 12; k++) (*g)[12*j+k] = helpGrad2[12*invisibleIndices[j]+k];
	}

	delete [] helpGrad1;
	delete [] helpGrad2;
}

void FineFittingTracking::smoothTransformations(int index, char * filename)
{
	int i, j, counterId; 

	g_Part * mesh = hierarchy[index]->getOriginalMesh();
	vector<int> indices = hierarchy[index]->getAdditionalIndices();
	vector<int> newIndices;

	int dimension = (int)(mesh->nodes().numberOfItems() - indices.size()) * 12;

	counterId = 0;
	for(i = 0; i < (int)mesh->nodes().numberOfItems(); i++)
	{
		if(counterId < indices.size())
		{
			if(indices[counterId] == i) counterId++;
			else newIndices.push_back(i);
		}
		else newIndices.push_back(i);
	}

	vnl_vector<double> x( dimension, 0.0 );

	for(i = 0; i < newIndices.size(); i++)
	{
		for(j = 0; j < 12; j++) x[12*i+j] = transformationParameters[12*newIndices[i]+j];
	}

	// Set the weights used here:
	accessWeight = 1.0;
	mdWeight = 1.0;
	rigidWeight = 1.0;
	nnWeight = 100.0;

	// Compute the access indices
	if(accessWeight > 0) computeAccessNeighbors(mesh, &newIndices);
	if(mdWeight > 0) computeMdNeighbors(mesh, &newIndices);

	// Compute the mean-value coordinates of the new points:
	g_NodeContainer nodeCont;
	for(i = 0; i < indices.size(); i++)
	{
		g_Node * node = new g_Node(computeTransformedPoint(mesh, indices[i]));
		nodeCont.insert(node);
	}
	g_Part * reconstructedMesh = hierarchy[index]->reconstructMesh(nodeCont, indices, true);

	if(filename != NULL)
	{
		char buffer[200];
		sprintf(buffer, "%s_reconstructed_mesh_%d.wrl", filename, index);
		TriangleMesh exportMesh(*reconstructedMesh);
		exportMeshWrapper( buffer, &exportMesh );
	}

	for(i = 0; i < indices.size(); i++) delete nodeCont[i];

	//perform the nonlinear optimization
	SmoothTransformationCostFunction costFct(this, mesh, reconstructedMesh, newIndices, dimension);
	vnl_lbfgsb * minimizer = new vnl_lbfgsb(costFct);
	minimizer->set_cost_function_convergence_factor ( 10000000 );
	minimizer->set_projected_gradient_tolerance ( 0.00001 ); 
	minimizer->set_max_function_evals(MAX_NUM_ITER);
	minimizer->minimize( x );

	//Copy the transformations:
	for(i = 0; i < newIndices.size(); i++)
	{
		for(j = 0; j < 12; j++) transformationParameters[12*newIndices[i]+j] = x[12*i+j];
	}

	//delete stuff:
	delete minimizer;
}

void FineFittingTracking::SmoothTransformationCostFunction::compute(vnl_vector< double > const &x, double *f, vnl_vector< double > *g)
{
	int i, j, dimension;

	dimension = newIndices.size() * 12;

	double * helpGrad1 = new double[mesh->nodes().numberOfItems() * 12];
	double * helpGrad2 = new double[mesh->nodes().numberOfItems() * 12];

	// Update the transformation parameters:
	for(i = 0; i < newIndices.size(); i++)
	{
		for(j = 0; j < 12; j++) tracker->transformationParameters[12*newIndices[i]+j] = x[12*i+j];
	}

	// Compute Energy:
	*f = 0;
	if(tracker->accessWeight > 0) *f += tracker->accessWeight * tracker->accessEnergy(mesh);
	if(tracker->mdWeight > 0) *f += tracker->mdWeight * tracker->mdEnergy(mesh);
	if(tracker->rigidWeight > 0) *f += tracker->rigidWeight * tracker->rigidEnergy(mesh, &newIndices);
	if(tracker->nnWeight > 0) *f += tracker->nnWeight * tracker->distanceResolutionEnergy(mesh, reconstructedMesh, newIndices);

	// Compute Gradient:
	for(i = 0; i < dimension; i++) (*g)[i] = 0;
	for(i = 0; i < mesh->nodes().numberOfItems() * 12; i++) helpGrad2[i] = 0;

	if(tracker->accessWeight > 0) 
	{
		tracker->accessGradient(mesh, helpGrad1);
		for(i = 0; i < mesh->nodes().numberOfItems() * 12; i++) 
			helpGrad2[i] += tracker->accessWeight * helpGrad1[i];
	}
	if(tracker->mdWeight > 0) 
	{
		tracker->mdGradient(mesh, helpGrad1);
		for(i = 0; i < mesh->nodes().numberOfItems() * 12; i++) 
			helpGrad2[i] += tracker->mdWeight * helpGrad1[i];
	}
	if(tracker->rigidWeight > 0) 
	{
		tracker->rigidGradient(mesh, helpGrad1);
		for(i = 0; i < mesh->nodes().numberOfItems() * 12; i++)
			helpGrad2[i] += tracker->rigidWeight * helpGrad1[i];
	}
	if(tracker->nnWeight > 0) 
	{
		tracker->distanceResolutionGradient(mesh, reconstructedMesh, newIndices, helpGrad1);
		for(i = 0; i < newIndices.size(); i++)
		{
			for(j = 0; j < 12; j++) helpGrad2[12*newIndices[i]+j] += tracker->nnWeight * helpGrad1[12*i+j];
		}
	}

	for(i = 0; i < newIndices.size(); i++)
	{
		for(j = 0; j < 12; j++) (*g)[12*i+j] = helpGrad2[12*newIndices[i]+j];
	}

	delete [] helpGrad1;
	delete [] helpGrad2;
}

void FineFittingTracking::alignRigidlyToFrame(bool onlyUseMarkers)
{
	//Check that everything is initialized
	if(!kdInitialized) return;

	vnl_vector<double> x( 8, 0.0 );
	x[0] = x[4] = x[5] = x[6] = 1.0;
	x[1] = x[2] = x[3] = x[7] = 0.0;

	//initialize the parameters for the nonlinear optimization
	RigidAlignmentCostFunction costFct(this);
	
	vnl_lbfgsb * minimizer = new vnl_lbfgsb(costFct);
	minimizer->set_cost_function_convergence_factor ( 10000000 );
	minimizer->set_projected_gradient_tolerance ( 0.00001 ); 
	minimizer->set_max_function_evals(MAX_NUM_ITER);

	costFct.featuresOnly = true;
	minimizer->minimize( x );
	if(!onlyUseMarkers)
	{
		costFct.featuresOnly = false;
		for(int numNNSteps = 0; numNNSteps < 30; numNNSteps++)
		{
			costFct.recomputeNeighbors = true;
			minimizer->minimize( x );
		}
	}

	//Set the result (modify the templateMesh):
	double * transformation1 = new double[16];
	double * transformation2 = new double[16];
	double * transformation3 = new double[16];
	double * transformation4 = new double[16];
	double * partialResult1 = new double[16];
	double * partialResult2 = new double[16];
	double * trafo = new double[16];
	double * queryPt = new double[3];
	char trans = 'N';
	long int dimTrafo = 4;
	double alpha = 1.0;
	double beta = 0.0; 

	computeTranslationMat(-x[1], -x[2], -x[3], transformation1);
	computeRotationMat(x[4], x[5], x[6], x[7], transformation2);
	computeScalingMat(x[0], transformation3);
	computeTranslationMat(x[1], x[2], x[3], transformation4);
	clapack::dgemm_(&trans, &trans, &dimTrafo, &dimTrafo, &dimTrafo, &alpha, transformation4, &dimTrafo, transformation3, &dimTrafo, &beta, partialResult1, &dimTrafo);
	clapack::dgemm_(&trans, &trans, &dimTrafo, &dimTrafo, &dimTrafo, &alpha, transformation2, &dimTrafo, transformation1, &dimTrafo, &beta, partialResult2, &dimTrafo);
	clapack::dgemm_(&trans, &trans, &dimTrafo, &dimTrafo, &dimTrafo, &alpha, partialResult1, &dimTrafo, partialResult2, &dimTrafo, &beta, trafo, &dimTrafo);
	for(int i = 0; i < (int)templateMesh->getNumberOfNodes(); i++)
	{
		queryPt[0] = trafo[0]*templateMesh->nodes()[i]->coordinate().x()+trafo[4]*templateMesh->nodes()[i]->coordinate().y()+trafo[8]*templateMesh->nodes()[i]->coordinate().z()+trafo[12];
		queryPt[1] = trafo[1]*templateMesh->nodes()[i]->coordinate().x()+trafo[5]*templateMesh->nodes()[i]->coordinate().y()+trafo[9]*templateMesh->nodes()[i]->coordinate().z()+trafo[13];
		queryPt[2] = trafo[2]*templateMesh->nodes()[i]->coordinate().x()+trafo[6]*templateMesh->nodes()[i]->coordinate().y()+trafo[10]*templateMesh->nodes()[i]->coordinate().z()+trafo[14];
		templateMesh->nodes()[i]->coordinate(g_Vector(queryPt[0], queryPt[1], queryPt[2]));
	}

	//delete stuff:
	delete minimizer;
	delete [] transformation1;
	delete [] transformation2;
	delete [] transformation3;
	delete [] transformation4;
	delete [] partialResult1;
	delete [] partialResult2;
	delete [] trafo;
	delete [] queryPt;

	cout<<"Finished rigid alignment"<<endl;
}

void FineFittingTracking::RigidAlignmentCostFunction::compute(vnl_vector< double > const &x, double *f, vnl_vector< double > *g)
{
	int i, j; 

	double * transformation1 = new double[16];
	double * transformation2 = new double[16];
	double * transformation3 = new double[16];
	double * transformation4 = new double[16];
	double * partialResult1 = new double[16];
	double * partialResult2 = new double[16];
	double * trafo = new double[16];
	double * trafoGrad = new double[16*8];
	char trans = 'N';
	char transY = 'T';
	long int dimTrafo = 4;
	double alpha = 1.0;
	double beta = 0.0; 

	ANNpoint			queryPt;
	ANNidxArray			nnIdx;	
	ANNdistArray		dists;
	queryPt = annAllocPt(3);			
	nnIdx = new ANNidx[1];	
	dists = new ANNdist[1];

	tracker->templateMesh->calculateVertexNormals();
	tracker->frameMesh->calculateVertexNormals();

	// Compute energy and gradient:
	tracker->computeTranslationMat(-x[1], -x[2], -x[3], transformation1);
	tracker->computeRotationMat(x[4], x[5], x[6], x[7], transformation2);
	tracker->computeScalingMat(x[0], transformation3);
	tracker->computeTranslationMat(x[1], x[2], x[3], transformation4);
	clapack::dgemm_(&trans, &trans, &dimTrafo, &dimTrafo, &dimTrafo, &alpha, transformation4, &dimTrafo, transformation3, &dimTrafo, &beta, partialResult1, &dimTrafo);
	clapack::dgemm_(&trans, &trans, &dimTrafo, &dimTrafo, &dimTrafo, &alpha, transformation2, &dimTrafo, transformation1, &dimTrafo, &beta, partialResult2, &dimTrafo);
	clapack::dgemm_(&trans, &trans, &dimTrafo, &dimTrafo, &dimTrafo, &alpha, partialResult1, &dimTrafo, partialResult2, &dimTrafo, &beta, trafo, &dimTrafo);	
	for(i = 0; i < 8; i++)
	{
		if((i == 1) || (i == 2) || (i == 3)) tracker->computeTranslationGrad(-x[1], -x[2], -x[3], transformation1, i-1);
		else tracker->computeTranslationMat(-x[1], -x[2], -x[3], transformation1);
		if(i > 3) tracker->computeRotationGrad(x[4], x[5], x[6], x[7], transformation2, i-4);
		else tracker->computeRotationMat(x[4], x[5], x[6], x[7], transformation2);
		if(i == 0) tracker->computeScalingGrad(x[0], transformation3);
		else tracker->computeScalingMat(x[0], transformation3);
		if((i == 1) || (i == 2) || (i == 3)) tracker->computeTranslationGrad(x[1], x[2], x[3], transformation1, i-1);
		else tracker->computeTranslationMat(x[1], x[2], x[3], transformation4);
		clapack::dgemm_(&trans, &trans, &dimTrafo, &dimTrafo, &dimTrafo, &alpha, transformation4, &dimTrafo, transformation3, &dimTrafo, &beta, partialResult1, &dimTrafo);
		clapack::dgemm_(&trans, &trans, &dimTrafo, &dimTrafo, &dimTrafo, &alpha, transformation2, &dimTrafo, transformation1, &dimTrafo, &beta, partialResult2, &dimTrafo);
		clapack::dgemm_(&trans, &trans, &dimTrafo, &dimTrafo, &dimTrafo, &alpha, partialResult1, &dimTrafo, partialResult2, &dimTrafo, &beta, &(trafoGrad[16*i]), &dimTrafo);	
	}

	int numValid = 0;
	*f = 0;
	for(i = 0; i < 8; i++) (*g)[i] = 0;
	if(featuresOnly || recomputeNeighbors) fixedNN.clear();
	for(i = 0; i < (int)tracker->templateMesh->getNumberOfNodes(); i++)
	{
		//deform point and compute nearest neighbor:
		queryPt[0] = trafo[0]*tracker->templateMesh->nodes()[i]->coordinate().x()+trafo[4]*tracker->templateMesh->nodes()[i]->coordinate().y()+
			trafo[8]*tracker->templateMesh->nodes()[i]->coordinate().z()+trafo[12];
		queryPt[1] = trafo[1]*tracker->templateMesh->nodes()[i]->coordinate().x()+trafo[5]*tracker->templateMesh->nodes()[i]->coordinate().y()+
			trafo[9]*tracker->templateMesh->nodes()[i]->coordinate().z()+trafo[13];
		queryPt[2] = trafo[2]*tracker->templateMesh->nodes()[i]->coordinate().x()+trafo[6]*tracker->templateMesh->nodes()[i]->coordinate().y()+
			trafo[10]*tracker->templateMesh->nodes()[i]->coordinate().z()+trafo[14];

		if(featuresOnly)
		{
			if(tracker->isFeatureInHierarchy.size() > 0) fixedNN.push_back(tracker->isFeatureInHierarchy[0][i]);
			else fixedNN.push_back(-1);
		}

		else if(recomputeNeighbors)
		{
			//Set to the features in case it is not NULL. Otherwise, use kd-tree search.
			if(tracker->isFeatureInHierarchy.size() > 0) 
			{
				if(tracker->isFeatureInHierarchy[0][i] != -1) fixedNN.push_back(tracker->isFeatureInHierarchy[0][i]);
				else
				{
					tracker->kdTree->annkPriSearch(queryPt, 1, nnIdx, dists);
					if((tracker->templateMesh->nodes()[i]->coordinate().DistanceTo(tracker->frameMesh->nodes()[nnIdx[0]]->coordinate()) < tracker->largestDist) && 
						(tracker->templateMesh->getNormalAtNode(i+1).AngleBetween(tracker->frameMesh->getNormalAtNode(nnIdx[0]+1)) < MAX_ANGLE/180.0*3.14159))
						fixedNN.push_back(nnIdx[0]);
					else fixedNN.push_back(-1);
				}
			}
			else
			{
				tracker->kdTree->annkPriSearch(queryPt, 1, nnIdx, dists);
				if((tracker->templateMesh->nodes()[i]->coordinate().DistanceTo(tracker->frameMesh->nodes()[nnIdx[0]]->coordinate()) < tracker->largestDist) && 
					(tracker->templateMesh->getNormalAtNode(i+1).AngleBetween(tracker->frameMesh->getNormalAtNode(nnIdx[0]+1)) < MAX_ANGLE/180.0*3.14159))
					fixedNN.push_back(nnIdx[0]);
				else fixedNN.push_back(-1);
			}
		}

		if(fixedNN[i] != -1)
		{
			g_Vector diffVec(g_Vector(queryPt[0], queryPt[1], queryPt[2])-tracker->frameMesh->nodes()[fixedNN[i]]->coordinate());
			*f += diffVec.SquaredLength();

			//compute gradient:
			for(j = 0; j < 8; j++)
			{
				queryPt[0] = trafoGrad[16*j+0]*tracker->templateMesh->nodes()[i]->coordinate().x()+trafoGrad[16*j+4]*tracker->templateMesh->nodes()[i]->coordinate().y()+
					trafoGrad[16*j+8]*tracker->templateMesh->nodes()[i]->coordinate().z()+trafoGrad[16*j+12];
				queryPt[1] = trafoGrad[16*j+1]*tracker->templateMesh->nodes()[i]->coordinate().x()+trafoGrad[16*j+5]*tracker->templateMesh->nodes()[i]->coordinate().y()+
					trafoGrad[16*j+9]*tracker->templateMesh->nodes()[i]->coordinate().z()+trafoGrad[16*j+13];
				queryPt[2] = trafoGrad[16*j+2]*tracker->templateMesh->nodes()[i]->coordinate().x()+trafoGrad[16*j+6]*tracker->templateMesh->nodes()[i]->coordinate().y()+
					trafoGrad[16*j+10]*tracker->templateMesh->nodes()[i]->coordinate().z()+trafoGrad[16*j+14];
	
				(*g)[j] += 2.0*diffVec.Dot(g_Vector(queryPt[0], queryPt[1], queryPt[2]));
			}

			numValid ++;
		}	
	}
	// Normalize:
	if(numValid > 0)
	{
		*f = *f * (double)tracker->templateMesh->getNumberOfNodes() / (double)numValid;
		for(i = 0; i < 8; i++) (*g)[i] = (*g)[i] * (double)tracker->templateMesh->getNumberOfNodes() / (double)numValid;
	}

	//If neighbors were computed, set recomputation to false:
	if(recomputeNeighbors) recomputeNeighbors = false;

	//Free space
	delete [] transformation1;
	delete [] transformation2;
	delete [] transformation3;
	delete [] transformation4;
	delete [] partialResult1;
	delete [] partialResult2;
	delete [] trafo;
	delete [] trafoGrad;
	annDeallocPt(queryPt);
	delete [] nnIdx;
	delete [] dists;
}

void FineFittingTracking::computeTranslationMat(double tx, double ty, double tz, double *& mat)
{
	mat[0] = 1;		mat[4] = 0;		mat[8] = 0;			mat[12] = tx;
	mat[1] = 0;		mat[5] = 1;		mat[9] = 0;			mat[13] = ty;
	mat[2] = 0;		mat[6] = 0;		mat[10] = 1;		mat[14] = tz;
	mat[3] = 0;		mat[7] = 0;		mat[11] = 0;		mat[15] = 1;
}

void FineFittingTracking::computeRotationMat(double tx, double ty, double tz, double angle, double *& mat)
{
	double quaternion[4];
	g_Vector crossProd(tx, ty, tz);
	crossProd.Normalize();
	crossProd = crossProd * sin(angle/2.0);
	quaternion[0] = cos(angle/2.0);		quaternion[1] = crossProd.x();		quaternion[2] = crossProd.y();		quaternion[3] = crossProd.z();
		
	mat[0] = pow(quaternion[0],2.0) + pow(quaternion[1],2.0) - pow(quaternion[2],2.0) - pow(quaternion[3],2.0);		
	mat[4] = 2*quaternion[1]*quaternion[2] - 2*quaternion[0]*quaternion[3];			
	mat[8] = 2*quaternion[1]*quaternion[3] + 2*quaternion[0]*quaternion[2];				
	mat[12] = 0;
	mat[1] = 2*quaternion[1]*quaternion[2] + 2*quaternion[0]*quaternion[3];	
	mat[5] = pow(quaternion[0],2.0) - pow(quaternion[1],2.0) + pow(quaternion[2],2.0) - pow(quaternion[3],2.0);	
	mat[9] = 2*quaternion[2]*quaternion[3] - 2*quaternion[0]*quaternion[1];				
	mat[13] = 0;
	mat[2] = 2*quaternion[1]*quaternion[3] - 2*quaternion[0]*quaternion[2];	
	mat[6] = 2*quaternion[2]*quaternion[3] + 2*quaternion[0]*quaternion[1];	
	mat[10] = pow(quaternion[0],2.0) - pow(quaternion[1],2.0) - pow(quaternion[2],2.0) + pow(quaternion[3],2.0);	
	mat[14] = 0;
	mat[3] = 0;			mat[7] = 0;			mat[11] = 0;			mat[15] = 1;	
}

void FineFittingTracking::computeScalingMat(double scale, double *& mat)
{
	mat[0] = scale;		mat[4] = 0;			mat[8] = 0;				mat[12] = 0;
	mat[1] = 0;			mat[5] = scale;		mat[9] = 0;				mat[13] = 0;
	mat[2] = 0;			mat[6] = 0;			mat[10] = scale;		mat[14] = 0;
	mat[3] = 0;			mat[7] = 0;			mat[11] = 0;			mat[15] = 1;
}

void FineFittingTracking::computeTranslationGrad(double tx, double ty, double tz, double *& mat, int index)
{
	for(int i = 0; i < 16; i++) mat[i] = 0;
	if(index == 0) mat[12] = 1;
	else if(index == 1) mat[13] = 1;
	else mat[14] = 1;
}

void FineFittingTracking::computeRotationGrad(double tx, double ty, double tz, double alpha, double *& mat, int index)
{
	// Normalize:
	double axisX, axisY, axisZ;
	g_Vector axis(tx, ty, tz);
	axis.Normalize();
	axisX = axis.x();	axisY = axis.y(); axisZ = axis.z();

	// Compute things related to angles:
	double sinA = sin(alpha/2.0);
	double cosA = cos(alpha/2.0);

	if(index < 3)
	{
		double grad_axisX, grad_axisY, grad_axisZ;
		if(index == 0)
		{
			//Gradient w.r.t. tx
			grad_axisX = (ty*ty+tz*tz)/pow(tx*tx+ty*ty+tz*tz, 1.5);
			grad_axisY = (-tx*ty)/pow(tx*tx+ty*ty+tz*tz, 1.5);
			grad_axisZ = (-tx*tz)/pow(tx*tx+ty*ty+tz*tz, 1.5);
		}
		else if(index == 1)
		{
			//Gradient w.r.t. ty
			grad_axisX = (-tx*ty)/pow(tx*tx+ty*ty+tz*tz, 1.5);
			grad_axisY = (tx*tx+tz*tz)/pow(tx*tx+ty*ty+tz*tz, 1.5);
			grad_axisZ = (-ty*tz)/pow(tx*tx+ty*ty+tz*tz, 1.5);
		}
		else if(index == 2)
		{
			//Gradient w.r.t. tz
			grad_axisX = (-tx*tz)/pow(tx*tx+ty*ty+tz*tz, 1.5);
			grad_axisY = (-ty*tz)/pow(tx*tx+ty*ty+tz*tz, 1.5);
			grad_axisZ = (tx*tx+ty*ty)/pow(tx*tx+ty*ty+tz*tz, 1.5);
		}

		mat[0] = 2*axisX*sinA*sinA*grad_axisX - 2*axisY*sinA*sinA*grad_axisY - 2*axisZ*sinA*sinA*grad_axisZ;			
		mat[4] = 2*sinA*sinA*(grad_axisX*axisY+axisX*grad_axisY) - 2*cosA*sinA*grad_axisZ;			
		mat[8] = 2*sinA*sinA*(grad_axisX*axisZ+axisX*grad_axisZ) + 2*cosA*sinA*grad_axisY; 	
		mat[1] = 2*sinA*sinA*(grad_axisX*axisY+axisX*grad_axisY) + 2*cosA*sinA*grad_axisZ;			
		mat[5] = -2*axisX*sinA*sinA*grad_axisX + 2*axisY*sinA*sinA*grad_axisY - 2*axisZ*sinA*sinA*grad_axisZ;		
		mat[9] = 2*sinA*sinA*(grad_axisY*axisZ+axisY*grad_axisZ) - 2*cosA*sinA*grad_axisX;
		mat[2] = 2*sinA*sinA*(grad_axisX*axisZ+axisX*grad_axisZ) - 2*cosA*sinA*grad_axisY; 		
		mat[6] = 2*sinA*sinA*(grad_axisY*axisZ+axisY*grad_axisZ) + 2*cosA*sinA*grad_axisX;				
		mat[10] = -2*axisX*sinA*sinA*grad_axisX - 2*axisY*sinA*sinA*grad_axisY + 2*axisZ*sinA*sinA*grad_axisZ;
	}
	else 
	{
		double gradSinSin = 2*sinA*cosA*0.5;
		double gradSinCos = sinA*(-sinA)*0.5 + cosA*cosA*0.5;
		double gradCosCos = 2*cosA*(-sinA)*0.5;

		//Gradient w.r.t. angle:
		mat[0] = gradCosCos + axisX*axisX*gradSinSin - axisY*axisY*gradSinSin - axisZ*axisZ*gradSinSin;
		mat[4] = 2*axisX*axisY*gradSinSin - 2*axisZ*gradSinCos;
		mat[8] = 2*axisX*axisZ*gradSinSin + 2*axisY*gradSinCos;	
		mat[1] = 2*axisX*axisY*gradSinSin + 2*axisZ*gradSinCos;			
		mat[5] = gradCosCos - axisX*axisX*gradSinSin + axisY*axisY*gradSinSin - axisZ*axisZ*gradSinSin;
		mat[9] = 2*axisY*axisZ*gradSinSin - 2*axisX*gradSinCos;	
		mat[2] = 2*axisX*axisZ*gradSinSin - 2*axisY*gradSinCos;			
		mat[6] = 2*axisY*axisZ*gradSinSin + 2*axisX*gradSinCos;
		mat[10] = gradCosCos - axisX*axisX*gradSinSin - axisY*axisY*gradSinSin + axisZ*axisZ*gradSinSin;
	}

	mat[12] = 0;		mat[13] = 0;		mat[14] = 0;
	mat[3] = 0;			mat[7] = 0;			mat[11] = 0;			mat[15] = 0;
}

void FineFittingTracking::computeScalingGrad(double s, double *& mat)
{
	mat[0] = 1;			mat[4] = 0;			mat[8] = 0;				mat[12] = 0;
	mat[1] = 0;			mat[5] = 1;			mat[9] = 0;				mat[13] = 0;
	mat[2] = 0;			mat[6] = 0;			mat[10] = 1;			mat[14] = 0;
	mat[3] = 0;			mat[7] = 0;			mat[11] = 0;			mat[15] = 0;
}

void FineFittingTracking::computeHierarchy(int numLevels)
{
	//compute smallestMesh and hierarchy here
	int i, j;

	g_Part * currentMesh = (g_Part *) templateMesh;
	
	for(i = 0; i < numLevels; i++)
	{
		MeshSimplification * simplify = new MeshSimplification(currentMesh);
		currentMesh = simplify->getSimplifiedMesh((int)floor(templateMesh->getNumberOfNodes() / pow(2.0, i+1)), true);
		if(currentMesh->nodes().numberOfItems() == simplify->getOriginalMesh()->nodes().numberOfItems())
		{
			//In this case, the mesh cannot be reduced further. Every edge collapse is invalid.
			break;
		}
		else {
		  hierarchy.push_back(simplify);
		}
	}

	smallestMesh = new g_Part(*currentMesh);

	//Compute and store the normal vectors:
	meshNormals.resize(hierarchy.size()+1);
	for(i = 0; i < hierarchy.size()+1; i++)
	{
		if(i == 0)
		{
			templateMesh->calculateVertexNormals();
			for(j = 0; j < templateMesh->getNumberOfNodes(); j++)
				meshNormals[i].push_back(templateMesh->getNormalAtNode(j+1));
		}
		else
		{
			vector<int> indices = hierarchy[i-1]->getAdditionalIndices();
			for(j = 0; j < indices.size(); j++)
				meshNormals[i].push_back(meshNormals[i-1][indices[j]]);
		}
	}

	//Compute the cluster indices if they have been set:
	int k, l;
	if(clusterIndicesInHierarchy.size() > 0)
	{
		//Transfer to all hierarchy levels:
		for(i = 0; i < hierarchy.size(); i++)
		{
			vector<int> indices = hierarchy[i]->getAdditionalIndices();
			vector< vector<int> > clusterIndicesCurrentLevel(indices.size());

			for(j = 0; j < indices.size(); j++)
			{
				for(k = 0; k < clusterIndicesInHierarchy[i][indices[j]].size(); k++)
				{
					for(l = 0; l < indices.size(); l++) 
					{
						if(indices[l] == clusterIndicesInHierarchy[i][indices[j]][k])
							clusterIndicesCurrentLevel[j].push_back(l);
					}
				}
			}
			clusterIndicesInHierarchy.push_back(clusterIndicesCurrentLevel);
		}
	}

	//Compute the feature levels if they have been set:
	if(isFeatureInHierarchy.size() > 0)
	{
		//Transfer to all hierarchy levels
		for(i = 0; i < hierarchy.size(); i++)
		{
			vector<int> indices = hierarchy[i]->getAllIndices();
			vector<int> isFeature(indices.size());
			for(j = 0; j < indices.size(); j++)
				isFeature[j] = isFeatureInHierarchy[i][indices[j]];

			isFeatureInHierarchy.push_back(isFeature);
		}
	}
}

void FineFittingTracking::fitToData(g_Part * mesh, vector<g_Vector> normalVectors, vector< vector<int> > clusterIds, vector<int> * isFeature, bool onlyUseMarkers)
{
	int i;
	double energy, prevEnergy;

	//Check that everything is initialized
	if(!kdInitialized) return;
	for(i = 0; i < (int)allEdges.numberOfItems(); i++) delete allEdges[i];
	allEdges.clear();
	allEdges = mesh->uniqueEdges();

	//Initialize weights and transformations:
	nnWeight = 1.0;
	rigidWeight = 1000.0;
	accessWeight = 1000.0;
	mdWeight = 1.0;
	//Don't use regularization (we have access energy instead)
	regWeight = 0.0;
	//Don't use clustering (too slow)
	clusterWeight = 0.0;
	//Don't use isometry either (since that requirement is not crucial):
	isoWeight = 0.0;

	energy = -1;
	while(rigidWeight > 10.0)
	{
		prevEnergy = energy;
		energy = solveOptimization(mesh, normalVectors, clusterIds, isFeature, onlyUseMarkers);
		accessWeight = accessWeight / 2.0;
		rigidWeight = rigidWeight / 2.0;

		if(prevEnergy != -1)
		{
			if(fabs(prevEnergy - energy)/prevEnergy < ENERGY_TOLERANCE) 
				break;
		}
	}
}

double FineFittingTracking::solveOptimization(g_Part * mesh, vector<g_Vector> normalVectors, vector< vector<int> > clusterIds, vector<int> * isFeature, bool onlyUseMarkers)
{
	int i, dimension;
	double energy;

	dimension = (int)mesh->nodes().numberOfItems() * 12;

	if(nnWeight > 0) computeNearestNeighbors(mesh, normalVectors, isFeature, onlyUseMarkers);
	if(accessWeight > 0) computeAccessNeighbors(mesh);
	if(mdWeight > 0) computeMdNeighbors(mesh);

	//initialize the parameters for the nonlinear optimization
	vnl_vector<double> x( dimension, 0.0 );
	for(i = 0; i < dimension; i++) x[i] = transformationParameters[i];
	SolveOptimizationCostFunction costFct(this, mesh, normalVectors, clusterIds, isFeature, dimension);
	
	vnl_lbfgsb * minimizer = new vnl_lbfgsb(costFct);
	minimizer->set_cost_function_convergence_factor ( 1000000000000 );
	minimizer->set_projected_gradient_tolerance ( 0.001 ); 
	minimizer->set_max_function_evals(MAX_NUM_ITER);
	minimizer->minimize( x );

	energy = minimizer->get_end_error();
	
	//copy the result
	for(i = 0; i < dimension; i++) transformationParameters[i] = x[i];

	//delete stuff:
	delete minimizer;

	return energy;
}

void FineFittingTracking::SolveOptimizationCostFunction::compute(vnl_vector< double > const &x, double *f, vnl_vector< double > *g)
{
	int i, dimension;
	dimension = mesh->nodes().numberOfItems() * 12;
	
	double * helpGrad = new double[dimension];

	// Copy the result:
	for(i = 0; i < dimension; i++) tracker->transformationParameters[i] = x[i];

	// Compute Energy:
	*f = 0;
	if(tracker->isoWeight > 0) *f += tracker->isoWeight * tracker->isometricEnergy(mesh);
	if(tracker->nnWeight > 0) *f += tracker->nnWeight * tracker->nearestNeighborEnergy(mesh);
	if(tracker->regWeight > 0) *f += tracker->regWeight * tracker->regularizationEnergy(mesh);
	if(tracker->rigidWeight > 0) *f += tracker->rigidWeight * tracker->rigidEnergy(mesh);
	if((tracker->clusterWeight > 0) && (clusterIds.size() > 0)) *f += tracker->clusterWeight * tracker->clusterEnergy(mesh, clusterIds);
	if(tracker->accessWeight > 0) *f += tracker->accessWeight * tracker->accessEnergy(mesh);
	if(tracker->mdWeight > 0) *f += tracker->mdWeight * tracker->mdEnergy(mesh);

	// Compute Gradient:
	for(i = 0; i < dimension; i++) 
		(*g)[i] = 0;
	if(tracker->isoWeight > 0) 
	{
		tracker->isometricGradient(mesh, helpGrad);
		for(i = 0; i < dimension; i++) 
			(*g)[i] += tracker->isoWeight * helpGrad[i];
	}
	if(tracker->nnWeight > 0)
	{
		tracker->nearestNeighborGradient(mesh, helpGrad);
		for(i = 0; i < dimension; i++) 
			(*g)[i] += tracker->nnWeight * helpGrad[i];
	}
	if(tracker->regWeight > 0)
	{
		tracker->regularizationGradient(mesh, helpGrad);
		for(i = 0; i < dimension; i++) 
			(*g)[i] += tracker->regWeight * helpGrad[i];
	}
	if(tracker->rigidWeight > 0)
	{
		tracker->rigidGradient(mesh, helpGrad);
		for(i = 0; i < dimension; i++)
			(*g)[i] += tracker->rigidWeight * helpGrad[i];
	}
	if((tracker->clusterWeight > 0) && (clusterIds.size() > 0))
	{
		tracker->clusterGradient(mesh, helpGrad, clusterIds);
		for(i = 0; i < dimension; i++)
			(*g)[i] += tracker->clusterWeight * helpGrad[i];
	}
	if(tracker->accessWeight > 0)
	{
		tracker->accessGradient(mesh, helpGrad);
		for(i = 0; i < dimension; i++) 
			(*g)[i] += tracker->accessWeight * helpGrad[i];
	}
	if(tracker->mdWeight > 0)
	{
		tracker->mdGradient(mesh, helpGrad);
		for(i = 0; i < dimension; i++) 
			(*g)[i] += tracker->mdWeight * helpGrad[i];
	}

	delete [] helpGrad;
}

double FineFittingTracking::isometricEnergy(g_Part * mesh)
{
	int i, id1, id2;
	double dist1, dist2;
	double energy = 0;

	for(i = 0; i < (int)allEdges.numberOfItems(); i++)
	{
		id1 = allEdges[i]->firstNode().id()-1;
		id2 = allEdges[i]->lastNode().id()-1;

		dist1 = (mesh->nodes()[id1]->coordinate() - mesh->nodes()[id2]->coordinate()).SquaredLength();

		dist2 = (computeTransformedPoint(mesh, id1) - computeTransformedPoint(mesh, id2)).SquaredLength();

		energy += pow(dist2-dist1, 2.0);
	}

	return energy;
}

void FineFittingTracking::isometricGradient(g_Part * mesh, double *& g)
{
	int i, id1, id2;
	double dist1, dist2;
	
	for(i = 0; i < 12*(int)mesh->nodes().numberOfItems(); i++) g[i] = 0;

	for(i = 0; i < (int)allEdges.numberOfItems(); i++)
	{
		id1 = allEdges[i]->firstNode().id()-1;
		id2 = allEdges[i]->lastNode().id()-1;

		dist1 = (mesh->nodes()[id1]->coordinate() - mesh->nodes()[id2]->coordinate()).SquaredLength();

		g_Vector diff = computeTransformedPoint(mesh, id1) - computeTransformedPoint(mesh, id2);
		dist2 = (diff).SquaredLength();		

		g[12*id1]	 += 4.0 * (dist2-dist1) * diff.x() * mesh->nodes()[id1]->coordinate().x();
		g[12*id1+1]	 += 4.0 * (dist2-dist1) * diff.y() * mesh->nodes()[id1]->coordinate().x();
		g[12*id1+2]	 += 4.0 * (dist2-dist1) * diff.z() * mesh->nodes()[id1]->coordinate().x();
		g[12*id1+3]	 += 4.0 * (dist2-dist1) * diff.x() * mesh->nodes()[id1]->coordinate().y();
		g[12*id1+4]  += 4.0 * (dist2-dist1) * diff.y() * mesh->nodes()[id1]->coordinate().y();
		g[12*id1+5]  += 4.0 * (dist2-dist1) * diff.z() * mesh->nodes()[id1]->coordinate().y();
		g[12*id1+6]  += 4.0 * (dist2-dist1) * diff.x() * mesh->nodes()[id1]->coordinate().z();
		g[12*id1+7]	 += 4.0 * (dist2-dist1) * diff.y() * mesh->nodes()[id1]->coordinate().z();
		g[12*id1+8]	 += 4.0 * (dist2-dist1) * diff.z() * mesh->nodes()[id1]->coordinate().z();
		g[12*id1+9]	 += 4.0 * (dist2-dist1) * diff.x();
		g[12*id1+10] += 4.0 * (dist2-dist1) * diff.y();
		g[12*id1+11] += 4.0 * (dist2-dist1) * diff.z();

		g[12*id2]	 -= 4.0 * (dist2-dist1) * diff.x() * mesh->nodes()[id2]->coordinate().x();
		g[12*id2+1]	 -= 4.0 * (dist2-dist1) * diff.y() * mesh->nodes()[id2]->coordinate().x();
		g[12*id2+2]	 -= 4.0 * (dist2-dist1) * diff.z() * mesh->nodes()[id2]->coordinate().x();
		g[12*id2+3]	 -= 4.0 * (dist2-dist1) * diff.x() * mesh->nodes()[id2]->coordinate().y();
		g[12*id2+4]  -= 4.0 * (dist2-dist1) * diff.y() * mesh->nodes()[id2]->coordinate().y();
		g[12*id2+5]  -= 4.0 * (dist2-dist1) * diff.z() * mesh->nodes()[id2]->coordinate().y();
		g[12*id2+6]  -= 4.0 * (dist2-dist1) * diff.x() * mesh->nodes()[id2]->coordinate().z();
		g[12*id2+7]	 -= 4.0 * (dist2-dist1) * diff.y() * mesh->nodes()[id2]->coordinate().z();
		g[12*id2+8]	 -= 4.0 * (dist2-dist1) * diff.z() * mesh->nodes()[id2]->coordinate().z();
		g[12*id2+9]	 -= 4.0 * (dist2-dist1) * diff.x();
		g[12*id2+10] -= 4.0 * (dist2-dist1) * diff.y();
		g[12*id2+11] -= 4.0 * (dist2-dist1) * diff.z();
	}
}

double FineFittingTracking::clusterEnergy(g_Part * mesh, vector< vector<int> > clusterIds)
{
	int i, j, id1, id2;
	double dist1, dist2;
	double energy = 0;

	for(i = 0; i < (int)mesh->nodes().numberOfItems(); i++)
	{
		for(j = 0; j < (int)clusterIds[i].size(); j++)
		{
			id1 = i;
			id2 = clusterIds[i][j];

			dist1 = (mesh->nodes()[id1]->coordinate() - mesh->nodes()[id2]->coordinate()).SquaredLength();

			dist2 = (computeTransformedPoint(mesh, id1) - computeTransformedPoint(mesh, id2)).SquaredLength();

			energy += pow(dist2-dist1, 2.0);
		}
	}
	return energy;
}

void FineFittingTracking::clusterGradient(g_Part * mesh, double *& g, vector< vector<int> > clusterIds)
{
	int i, j, id1, id2;
	double dist1, dist2;
	
	for(i = 0; i < 12*(int)mesh->nodes().numberOfItems(); i++) g[i] = 0;

	for(i = 0; i < (int)mesh->nodes().numberOfItems(); i++)
	{
		for(j = 0; j < (int)clusterIds[i].size(); j++)
		{
			id1 = i;
			id2 = clusterIds[i][j];

			dist1 = (mesh->nodes()[id1]->coordinate() - mesh->nodes()[id2]->coordinate()).SquaredLength();

			g_Vector diff = computeTransformedPoint(mesh, id1) - computeTransformedPoint(mesh, id2);
			dist2 = (diff).SquaredLength();		

			g[12*id1]	 += 4.0 * (dist2-dist1) * diff.x() * mesh->nodes()[id1]->coordinate().x();
			g[12*id1+1]	 += 4.0 * (dist2-dist1) * diff.y() * mesh->nodes()[id1]->coordinate().x();
			g[12*id1+2]	 += 4.0 * (dist2-dist1) * diff.z() * mesh->nodes()[id1]->coordinate().x();
			g[12*id1+3]	 += 4.0 * (dist2-dist1) * diff.x() * mesh->nodes()[id1]->coordinate().y();
			g[12*id1+4]  += 4.0 * (dist2-dist1) * diff.y() * mesh->nodes()[id1]->coordinate().y();
			g[12*id1+5]  += 4.0 * (dist2-dist1) * diff.z() * mesh->nodes()[id1]->coordinate().y();
			g[12*id1+6]  += 4.0 * (dist2-dist1) * diff.x() * mesh->nodes()[id1]->coordinate().z();
			g[12*id1+7]	 += 4.0 * (dist2-dist1) * diff.y() * mesh->nodes()[id1]->coordinate().z();
			g[12*id1+8]	 += 4.0 * (dist2-dist1) * diff.z() * mesh->nodes()[id1]->coordinate().z();
			g[12*id1+9]	 += 4.0 * (dist2-dist1) * diff.x();
			g[12*id1+10] += 4.0 * (dist2-dist1) * diff.y();
			g[12*id1+11] += 4.0 * (dist2-dist1) * diff.z();

			g[12*id2]	 -= 4.0 * (dist2-dist1) * diff.x() * mesh->nodes()[id2]->coordinate().x();
			g[12*id2+1]	 -= 4.0 * (dist2-dist1) * diff.y() * mesh->nodes()[id2]->coordinate().x();
			g[12*id2+2]	 -= 4.0 * (dist2-dist1) * diff.z() * mesh->nodes()[id2]->coordinate().x();
			g[12*id2+3]	 -= 4.0 * (dist2-dist1) * diff.x() * mesh->nodes()[id2]->coordinate().y();
			g[12*id2+4]  -= 4.0 * (dist2-dist1) * diff.y() * mesh->nodes()[id2]->coordinate().y();
			g[12*id2+5]  -= 4.0 * (dist2-dist1) * diff.z() * mesh->nodes()[id2]->coordinate().y();
			g[12*id2+6]  -= 4.0 * (dist2-dist1) * diff.x() * mesh->nodes()[id2]->coordinate().z();
			g[12*id2+7]	 -= 4.0 * (dist2-dist1) * diff.y() * mesh->nodes()[id2]->coordinate().z();
			g[12*id2+8]	 -= 4.0 * (dist2-dist1) * diff.z() * mesh->nodes()[id2]->coordinate().z();
			g[12*id2+9]	 -= 4.0 * (dist2-dist1) * diff.x();
			g[12*id2+10] -= 4.0 * (dist2-dist1) * diff.y();
			g[12*id2+11] -= 4.0 * (dist2-dist1) * diff.z();
		}
	}
}

void FineFittingTracking::computeNearestNeighbors(g_Part * mesh, vector<g_Vector> normalVectors, vector<int> * isFeature, bool onlyUseMarkers,
												  bool symmetric)
{
	int i, j;
	double x, y, z;

	//Compute the nearest neighbors if not only markers are used
	if(onlyUseMarkers)
	{
		nearestNeighbors.resize((int)mesh->nodes().numberOfItems());
		for(i = 0; i < (int)mesh->nodes().numberOfItems(); i++) nearestNeighbors[i] = -1;
	}
	else
	{
		long int n = (long int) 3;
		long int * ipiv = new long int [n];
		long int info;
		double * work = new double[n];

		// Initialize kd-tree use
		ANNpoint			queryPt;
		ANNidxArray			nnIdx;	
		ANNdistArray		dists;
		queryPt = annAllocPt(3);			
		nnIdx = new ANNidx[1];	
		dists = new ANNdist[1];

		//Transform the normals:
		vector<g_Vector> transformedNormals((int)mesh->nodes().numberOfItems());
		double * currentTrafo = new double[9];
		for(i = 0; i < (int)mesh->nodes().numberOfItems(); i++)
		{
			currentTrafo[0] = transformationParameters[12*i];	currentTrafo[3] = transformationParameters[12*i+3]; currentTrafo[6] = transformationParameters[12*i+6];
			currentTrafo[1] = transformationParameters[12*i+1]; currentTrafo[4] = transformationParameters[12*i+4]; currentTrafo[7] = transformationParameters[12*i+7]; 
			currentTrafo[2] = transformationParameters[12*i+2]; currentTrafo[5] = transformationParameters[12*i+5]; currentTrafo[8] = transformationParameters[12*i+8];
			
			//Compute inverse:
			clapack::dgetrf_(&n, &n, currentTrafo, &n, ipiv, &info);
			if(info != 0) cout<<"Problem when transforming normals in computeNearestNeighbors(...)"<<endl;
			clapack::dgetri_(&n, currentTrafo, &n, ipiv, work, &n, &info);
			if(info != 0) cout<<"Problem when transforming normals in computeNearestNeighbors(...)"<<endl;
			
			//Compute transformed normal (use transposed currentTrafo) and store in transformedNormals:
			x = currentTrafo[0] * normalVectors[i].x() + currentTrafo[1] * normalVectors[i].y() + currentTrafo[2] * normalVectors[i].z();
			y = currentTrafo[3] * normalVectors[i].x() + currentTrafo[4] * normalVectors[i].y() + currentTrafo[5] * normalVectors[i].z();
			z = currentTrafo[6] * normalVectors[i].x() + currentTrafo[7] * normalVectors[i].y() + currentTrafo[8] * normalVectors[i].z();
			transformedNormals[i].Set(x, y, z);
		}

		//Copy the mesh and compute its normals
		TriangleMesh tempMesh;
		for(i = 0; i < (int)mesh->nodes().numberOfItems(); i++)
		{
			g_Node * node = new g_Node(computeTransformedPoint(mesh, i));
			tempMesh.node(node);
		}
		for(i = 0; i < (int)mesh->elements().numberOfItems(); i++)
		{
			g_Element * elem = new g_Element;
			elem->node(tempMesh.nodes()[mesh->elements()[i]->nodes()[0]->id()-1]);
			elem->node(tempMesh.nodes()[mesh->elements()[i]->nodes()[1]->id()-1]);
			elem->node(tempMesh.nodes()[mesh->elements()[i]->nodes()[2]->id()-1]);
			tempMesh.element(elem);
		}

		// Find the nearest neighbor for each transformed point
		nearestNeighbors.resize((int)mesh->nodes().numberOfItems());
		vector<double> distancesForSymmetry((int)mesh->nodes().numberOfItems());
		for(i = 0; i < (int)mesh->nodes().numberOfItems(); i++)
		{
			queryPt[0] = tempMesh.nodes()[i]->coordinate().x();
			queryPt[1] = tempMesh.nodes()[i]->coordinate().y();
			queryPt[2] = tempMesh.nodes()[i]->coordinate().z();

			kdTree->annkPriSearch(queryPt, 1, nnIdx, dists);
			distancesForSymmetry[i] = dists[0];
			if((sqrt(dists[0]) < largestDist) && (transformedNormals[i].AngleBetween(frameMesh->getNormalAtNode(nnIdx[0]+1)) < MAX_ANGLE/180.0*3.14159))
			{
				nearestNeighbors[i] = nnIdx[0];

				// If symmetric is true, only keep the smallest distance to each index
				if(symmetric)
				{
					for(j = 0; j < i; j++)
					{
						if(nearestNeighbors[j] == nearestNeighbors[i])
						{
							if(distancesForSymmetry[j] < distancesForSymmetry[i])
								nearestNeighbors[i] = -1;
							else nearestNeighbors[j] = -1;

							break;
						}
					}
				}
			}
			else nearestNeighbors[i] = -1;
		}

		// Free space
		annDeallocPt(queryPt);
		delete [] nnIdx;
		delete [] dists;

		delete [] ipiv;
		delete [] work;
	}

	//Set to the features in case it is not NULL
	if(isFeature != NULL)
	{
		//Check for dimension
		if(isFeature->size() != nearestNeighbors.size())
		{
			cout<<"Incorrect number of features in computeNearestNeighbors "<<isFeature->size()<<" "<<nearestNeighbors.size()<<endl;
			return;
		}

		//Set features
		for(i = 0; i < (int)(*isFeature).size(); i++)
		{
			if((*isFeature)[i] != -1) nearestNeighbors[i] = (*isFeature)[i];
		}
	}
}

double FineFittingTracking::nearestNeighborEnergy(g_Part * mesh)
{
	int i;
	double energy = 0;

	// Find the nearest neighbor for each transformed point and add the squared distance to the energy
	for(i = 0; i < (int)mesh->nodes().numberOfItems(); i++)
	{
		if(nearestNeighbors[i] != -1)
			energy += (computeTransformedPoint(mesh, i) - frameMesh->nodes()[nearestNeighbors[i]]->coordinate()).SquaredLength();
	}

	return energy;
}

void FineFittingTracking::nearestNeighborGradient(g_Part * mesh, double *& g)
{
	int i;
	g_Vector queryPt;

	for(i = 0; i < 12*(int)mesh->nodes().numberOfItems(); i++) g[i] = 0;

	for(i = 0; i < (int)mesh->nodes().numberOfItems(); i++)
	{
		if(nearestNeighbors[i] != -1)
		{
			queryPt = computeTransformedPoint(mesh, i);

			g[12*i]		 += 2.0 * (queryPt.x() - frameMesh->nodes()[nearestNeighbors[i]]->coordinate().x()) * mesh->nodes()[i]->coordinate().x();
			g[12*i+1]	 += 2.0 * (queryPt.y() - frameMesh->nodes()[nearestNeighbors[i]]->coordinate().y()) * mesh->nodes()[i]->coordinate().x();
			g[12*i+2]	 += 2.0 * (queryPt.z() - frameMesh->nodes()[nearestNeighbors[i]]->coordinate().z()) * mesh->nodes()[i]->coordinate().x();
			g[12*i+3]	 += 2.0 * (queryPt.x() - frameMesh->nodes()[nearestNeighbors[i]]->coordinate().x()) * mesh->nodes()[i]->coordinate().y();
			g[12*i+4]    += 2.0 * (queryPt.y() - frameMesh->nodes()[nearestNeighbors[i]]->coordinate().y()) * mesh->nodes()[i]->coordinate().y();
			g[12*i+5]    += 2.0 * (queryPt.z() - frameMesh->nodes()[nearestNeighbors[i]]->coordinate().z()) * mesh->nodes()[i]->coordinate().y();
			g[12*i+6]    += 2.0 * (queryPt.x() - frameMesh->nodes()[nearestNeighbors[i]]->coordinate().x()) * mesh->nodes()[i]->coordinate().z();
			g[12*i+7]	 += 2.0 * (queryPt.y() - frameMesh->nodes()[nearestNeighbors[i]]->coordinate().y()) * mesh->nodes()[i]->coordinate().z();
			g[12*i+8]	 += 2.0 * (queryPt.z() - frameMesh->nodes()[nearestNeighbors[i]]->coordinate().z()) * mesh->nodes()[i]->coordinate().z();
			g[12*i+9]	 += 2.0 * (queryPt.x() - frameMesh->nodes()[nearestNeighbors[i]]->coordinate().x());
			g[12*i+10]   += 2.0 * (queryPt.y() - frameMesh->nodes()[nearestNeighbors[i]]->coordinate().y());
			g[12*i+11]   += 2.0 * (queryPt.z() - frameMesh->nodes()[nearestNeighbors[i]]->coordinate().z());
		}
	}
}

double FineFittingTracking::distanceResolutionEnergy(g_Part * mesh, g_Part * reconstructedMesh, vector<int> indices)
{
	int i;

	double energy = 0.0;

	for(i = 0; i < indices.size(); i++)		
		energy += (computeTransformedPoint(mesh, indices[i]) - reconstructedMesh->nodes()[indices[i]]->coordinate()).SquaredLength();
	
	return energy;
}

void FineFittingTracking::distanceResolutionGradient(g_Part * mesh, g_Part * reconstructedMesh, vector<int> indices, double *& g)
{
	int i;
	g_Vector transformedPoint;

	for(i = 0; i < 12*indices.size(); i++) g[i] = 0;
	
	for(i = 0; i < indices.size(); i++)
	{
		transformedPoint = computeTransformedPoint(mesh, indices[i]);

		g[12*i] += 2 * (transformedPoint.x() - reconstructedMesh->nodes()[indices[i]]->coordinate().x()) * mesh->nodes()[indices[i]]->coordinate().x();
		g[12*i+1] += 2 * (transformedPoint.y() - reconstructedMesh->nodes()[indices[i]]->coordinate().y()) * mesh->nodes()[indices[i]]->coordinate().x();
		g[12*i+2] += 2 * (transformedPoint.z() - reconstructedMesh->nodes()[indices[i]]->coordinate().z()) * mesh->nodes()[indices[i]]->coordinate().x();
		g[12*i+3] += 2 * (transformedPoint.x() - reconstructedMesh->nodes()[indices[i]]->coordinate().x()) * mesh->nodes()[indices[i]]->coordinate().y();
		g[12*i+4] += 2 * (transformedPoint.y() - reconstructedMesh->nodes()[indices[i]]->coordinate().y()) * mesh->nodes()[indices[i]]->coordinate().y();
		g[12*i+5] += 2 * (transformedPoint.z() - reconstructedMesh->nodes()[indices[i]]->coordinate().z()) * mesh->nodes()[indices[i]]->coordinate().y();
		g[12*i+6] += 2 * (transformedPoint.x() - reconstructedMesh->nodes()[indices[i]]->coordinate().x()) * mesh->nodes()[indices[i]]->coordinate().z();
		g[12*i+7] += 2 * (transformedPoint.y() - reconstructedMesh->nodes()[indices[i]]->coordinate().y()) * mesh->nodes()[indices[i]]->coordinate().z();
		g[12*i+8] += 2 * (transformedPoint.z() - reconstructedMesh->nodes()[indices[i]]->coordinate().z()) * mesh->nodes()[indices[i]]->coordinate().z();
		g[12*i+9] += 2 * (transformedPoint.x() - reconstructedMesh->nodes()[indices[i]]->coordinate().x());
		g[12*i+10] += 2 * (transformedPoint.y() - reconstructedMesh->nodes()[indices[i]]->coordinate().y());
		g[12*i+11] += 2 * (transformedPoint.z() - reconstructedMesh->nodes()[indices[i]]->coordinate().z());
	}
}

double FineFittingTracking::regularizationEnergy(g_Part * mesh, vector<int> * indices)
{
	int i, j, id1, id2;
	double energy = 0;

	if(indices == NULL)
	{
		for(i = 0; i < (int)allEdges.numberOfItems(); i++)
		{
			id1 = allEdges[i]->firstNode().id()-1;
			id2 = allEdges[i]->lastNode().id()-1;
			for(j = 0; j < 12; j++)
				energy += pow(transformationParameters[id1*12+j] - transformationParameters[id2*12+j], 2.0);
		}
	}

	else
	{
		bool found;
		for(i = 0; i < (int)allEdges.numberOfItems(); i++)
		{
			id1 = allEdges[i]->firstNode().id()-1;
			id2 = allEdges[i]->lastNode().id()-1;
			found = false;
			for(j = 0; j < indices->size(); j++)
			{
				if(((*indices)[j] == id1) || ((*indices)[j] == id2)) 
				{
					found = true;
					break;
				}
			}

			if(found)
			{
				for(j = 0; j < 12; j++)
					energy += pow(transformationParameters[id1*12+j] - transformationParameters[id2*12+j], 2.0);
			}
		}
	}

	return energy;
}

void FineFittingTracking::regularizationGradient(g_Part * mesh, double *& g)
{
	int i, j, id1, id2;
	for(i = 0; i < 12*(int)mesh->nodes().numberOfItems(); i++) g[i] = 0;

	for(i = 0; i < (int)allEdges.numberOfItems(); i++)
	{
		id1 = allEdges[i]->firstNode().id()-1;
		id2 = allEdges[i]->lastNode().id()-1;
		for(j = 0; j < 12; j++)
		{
			g[id1*12+j] += 2.0 * (transformationParameters[id1*12+j] - transformationParameters[id2*12+j]);
			g[id2*12+j] += 2.0 * (transformationParameters[id2*12+j] - transformationParameters[id1*12+j]);
		}
	}
}

double FineFittingTracking::rigidEnergy(g_Part * mesh, vector<int> * indices)
{
	int i;
	double energy = 0;

	if(indices == NULL)
	{
		for(i = 0; i < (int)mesh->nodes().numberOfItems(); i++)
		{
			energy += pow(1.0 - pow(transformationParameters[12*i], 2) - pow(transformationParameters[12*i+1], 2) - pow(transformationParameters[12*i+2], 2), 2);
			energy += pow(1.0 - pow(transformationParameters[12*i+3], 2) - pow(transformationParameters[12*i+4], 2) - pow(transformationParameters[12*i+5], 2), 2);
			energy += pow(1.0 - pow(transformationParameters[12*i+6], 2) - pow(transformationParameters[12*i+7], 2) - pow(transformationParameters[12*i+8], 2), 2);
			energy += pow((transformationParameters[12*i]*transformationParameters[12*i+3] + transformationParameters[12*i+1]*transformationParameters[12*i+4] + 
				transformationParameters[12*i+2]*transformationParameters[12*i+5]), 2);
			energy += pow((transformationParameters[12*i+6]*transformationParameters[12*i+3] + transformationParameters[12*i+7]*transformationParameters[12*i+4] + 
				transformationParameters[12*i+8]*transformationParameters[12*i+5]), 2);
			energy += pow((transformationParameters[12*i]*transformationParameters[12*i+6] + transformationParameters[12*i+1]*transformationParameters[12*i+7] + 
				transformationParameters[12*i+2]*transformationParameters[12*i+8]), 2);
		}
	}

	else
	{
		for(i = 0; i < indices->size(); i++)
		{
			int id = (*indices)[i];
			energy += pow(1.0 - pow(transformationParameters[12*id], 2) - pow(transformationParameters[12*id+1], 2) - pow(transformationParameters[12*id+2], 2), 2);
			energy += pow(1.0 - pow(transformationParameters[12*id+3], 2) - pow(transformationParameters[12*id+4], 2) - pow(transformationParameters[12*id+5], 2), 2);
			energy += pow(1.0 - pow(transformationParameters[12*id+6], 2) - pow(transformationParameters[12*id+7], 2) - pow(transformationParameters[12*id+8], 2), 2);
			energy += pow((transformationParameters[12*id]*transformationParameters[12*id+3] + transformationParameters[12*id+1]*transformationParameters[12*id+4] + 
				transformationParameters[12*id+2]*transformationParameters[12*id+5]), 2);
			energy += pow((transformationParameters[12*id+6]*transformationParameters[12*id+3] + transformationParameters[12*id+7]*transformationParameters[12*id+4] + 
				transformationParameters[12*id+8]*transformationParameters[12*id+5]), 2);
			energy += pow((transformationParameters[12*id]*transformationParameters[12*id+6] + transformationParameters[12*id+1]*transformationParameters[12*id+7] + 
				transformationParameters[12*id+2]*transformationParameters[12*id+8]), 2);
		}
	}

	return energy;
}

void FineFittingTracking::rigidGradient(g_Part * mesh, double *& g)
{
	int i;
	double grad1, grad2, grad3, grad4, grad5, grad6;

	for(i = 0; i < 12*(int)mesh->nodes().numberOfItems(); i++) g[i] = 0;

	for(i = 0; i < (int)mesh->nodes().numberOfItems(); i++)
	{
		grad1 = 2 * (1.0 - pow(transformationParameters[12*i], 2) - pow(transformationParameters[12*i+1], 2) - pow(transformationParameters[12*i+2], 2));
		grad2 = 2 * (1.0 - pow(transformationParameters[12*i+3], 2) - pow(transformationParameters[12*i+4], 2) - pow(transformationParameters[12*i+5], 2)); 
		grad3 = 2 * (1.0 - pow(transformationParameters[12*i+6], 2) - pow(transformationParameters[12*i+7], 2) - pow(transformationParameters[12*i+8], 2));
		grad4 = 2 * (transformationParameters[12*i]*transformationParameters[12*i+3] + transformationParameters[12*i+1]*transformationParameters[12*i+4] + 
			transformationParameters[12*i+2]*transformationParameters[12*i+5]);
		grad5 = 2 * (transformationParameters[12*i+6]*transformationParameters[12*i+3] + transformationParameters[12*i+7]*transformationParameters[12*i+4] + 
			transformationParameters[12*i+8]*transformationParameters[12*i+5]);
		grad6 = 2 * (transformationParameters[12*i]*transformationParameters[12*i+6] + transformationParameters[12*i+1]*transformationParameters[12*i+7] + 
			transformationParameters[12*i+2]*transformationParameters[12*i+8]);

		g[12*i] += grad1 * (-2*transformationParameters[12*i]) + grad4 * transformationParameters[12*i+3] + grad6 * transformationParameters[12*i+6];
		g[12*i+1] += grad1 * (-2*transformationParameters[12*i+1]) + grad4 * transformationParameters[12*i+4] + grad6 * transformationParameters[12*i+7];
		g[12*i+2] += grad1 * (-2*transformationParameters[12*i+2]) + grad4 * transformationParameters[12*i+5] + grad6 * transformationParameters[12*i+8];
		g[12*i+3] += grad2 * (-2*transformationParameters[12*i+3]) + grad4 * transformationParameters[12*i] + grad5 * transformationParameters[12*i+6];
		g[12*i+4] += grad2 * (-2*transformationParameters[12*i+4]) + grad4 * transformationParameters[12*i+1] + grad5 * transformationParameters[12*i+7];
		g[12*i+5] += grad2 * (-2*transformationParameters[12*i+5]) + grad4 * transformationParameters[12*i+2] + grad5 * transformationParameters[12*i+8];
		g[12*i+6] += grad3 * (-2*transformationParameters[12*i+6]) + grad5 * transformationParameters[12*i+3] + grad6 * transformationParameters[12*i];
		g[12*i+7] += grad3 * (-2*transformationParameters[12*i+7]) + grad5 * transformationParameters[12*i+4] + grad6 * transformationParameters[12*i+1];
		g[12*i+8] += grad3 * (-2*transformationParameters[12*i+8]) + grad5 * transformationParameters[12*i+5] + grad6 * transformationParameters[12*i+2];
	}
}

void FineFittingTracking::computeAccessNeighbors(g_Part * mesh, vector<int> * indices)
{
	bool found;
	TriangleMesh tempMesh (*mesh);
	
	double radius = SIZE_NBHD * tempMesh.getMeshResolution();

	indicesForIntersection.clear();
	indicesForIntersection.resize(tempMesh.getNumberOfNodes());

	tempMesh.initGeodesicDistanceCalculation();

	//create kd-tree and find nearest neighbors
	int i, j, numNN;

	ANNpointArray		dataPtsForNN;
	ANNpoint			queryPtForNN;
	ANNidxArray			nnIdxForNN;	
	ANNdistArray		distsForNN;
	ANNkd_tree*			kdTreeForNN;	
	queryPtForNN = annAllocPt(3);				
	dataPtsForNN = annAllocPts(tempMesh.getNumberOfNodes(), 3);		
	nnIdxForNN = new ANNidx[tempMesh.getNumberOfNodes()];
	distsForNN = new ANNdist[tempMesh.getNumberOfNodes()];

	g_Vector transformedPoint;
	for(i = 0; i < tempMesh.getNumberOfNodes(); i++)
	{
		transformedPoint = computeTransformedPoint(mesh, i);
		dataPtsForNN[i][0] = transformedPoint.x();
		dataPtsForNN[i][1] = transformedPoint.y();
		dataPtsForNN[i][2] = transformedPoint.z();
	}

	kdTreeForNN = new ANNkd_tree(dataPtsForNN, tempMesh.getNumberOfNodes(), 3);

	//do spherical range search and store the list of neighbors:
	svector<double> geodesicDists(tempMesh.getNumberOfNodes());
	for(i = 0; i < tempMesh.getNumberOfNodes(); i++)
	{
		if(indices == NULL) found = true;
		else 
		{
			found = false;
			for(j = 0; j < indices->size(); j++)
			{
				if((*indices)[j] == i)
				{
					found = true;
					break;
				}
			}
		}
		if(found)
		{
			queryPtForNN[0] = dataPtsForNN[i][0];
			queryPtForNN[1] = dataPtsForNN[i][1];
			queryPtForNN[2] = dataPtsForNN[i][2];
			numNN = kdTreeForNN->annkFRSearch(queryPtForNN, radius*radius, tempMesh.getNumberOfNodes(), nnIdxForNN, distsForNN);
			indicesForIntersection[i].clear();
			
			tempMesh.getGeodesicDistance(i+1, radius, geodesicDists);

			for(j = 0; j < numNN; j++)
			{
				//intersect the Euclidean sphere with a geodesic circle
				if((nnIdxForNN[j] != i) && (geodesicDists[nnIdxForNN[j]] != -1) && (geodesicDists[nnIdxForNN[j]] < radius))
					indicesForIntersection[i].push_back(nnIdxForNN[j]);
			}
		}
		else indicesForIntersection[i].clear();
	}

	delete [] nnIdxForNN;
	delete [] distsForNN;
	annDeallocPt(queryPtForNN);
	annDeallocPts(dataPtsForNN);
	delete kdTreeForNN;
}

double FineFittingTracking::accessEnergy(g_Part * mesh)
{
	int i, j, k;
	double energy, trafoDifference, weight, dist;
	g_Vector queryPt1, queryPt2;

	TriangleMesh tempMesh (*mesh);
	
	double radius = SIZE_NBHD * tempMesh.getMeshResolution();

	energy = 0;

	//Go through the list of neighbors:
	for(i = 0; i < tempMesh.getNumberOfNodes(); i++)
	{
		queryPt1 = computeTransformedPoint(mesh, i);

		for(j = 0; j < indicesForIntersection[i].size(); j++)
		{
			queryPt2 = computeTransformedPoint(mesh, indicesForIntersection[i][j]);

			trafoDifference = 0;
			for(k = 0; k < 12; k++) trafoDifference += pow((transformationParameters[12*i+k]-transformationParameters[12*indicesForIntersection[i][j]+k]), 2.0);

			dist = (queryPt1-queryPt2).SquaredLength();

			weight = max(0.0, 1.0 - (dist / ((double)max(radius*radius, MD_DEFAULT_EPSILON))));
			energy += weight * trafoDifference;
		}
	}

	return energy;
}

void FineFittingTracking::accessGradient(g_Part * mesh, double *& g)
{
	int i, j, k, currentId;
	double factor, weight, trafoDifference;
	g_Vector queryPt1, queryPt2, diffVec;
	TriangleMesh tempMesh(*mesh);

	double radius = SIZE_NBHD * tempMesh.getMeshResolution();

	for(i = 0; i < 12*tempMesh.getNumberOfNodes(); i++) g[i] = 0.0;

	//compute the gradient using TriangleMesh indicesForIntersection
	for(i = 0; i < tempMesh.getNumberOfNodes(); i++)
	{
		queryPt1 = computeTransformedPoint(mesh, i);

		for(j = 0; j < indicesForIntersection[i].size(); j++)
		{	
			currentId = indicesForIntersection[i][j];
			queryPt2 = computeTransformedPoint(mesh, currentId);

			//Use the product rule: derivative of weight * trafoDifference + derivative of trafoDifference * weight
			diffVec = queryPt1 - queryPt2;
			weight = max(0.0, 1.0 - (diffVec.SquaredLength() / ((double)max(radius*radius, MD_DEFAULT_EPSILON))));
			trafoDifference = 0;
			for(k = 0; k < 12; k++) trafoDifference += pow((transformationParameters[12*i+k] - transformationParameters[12*currentId+k]), 2.0);

			factor = (-2.0/((double)max(radius*radius, MD_DEFAULT_EPSILON))) * trafoDifference; 

			g[12*i]				+= factor * diffVec.x() * mesh->nodes()[i]->coordinate().x();
			g[12*i+1]			+= factor * diffVec.y() * mesh->nodes()[i]->coordinate().x();
			g[12*i+2]			+= factor * diffVec.z() * mesh->nodes()[i]->coordinate().x();
			g[12*i+3]			+= factor * diffVec.x() * mesh->nodes()[i]->coordinate().y();
			g[12*i+4]			+= factor * diffVec.y() * mesh->nodes()[i]->coordinate().y();
			g[12*i+5]			+= factor * diffVec.z() * mesh->nodes()[i]->coordinate().y();
			g[12*i+6]			+= factor * diffVec.x() * mesh->nodes()[i]->coordinate().z();
			g[12*i+7]			+= factor * diffVec.y() * mesh->nodes()[i]->coordinate().z();
			g[12*i+8]			+= factor * diffVec.z() * mesh->nodes()[i]->coordinate().z();
			g[12*i+9]			+= factor * diffVec.x();
			g[12*i+10]			+= factor * diffVec.y();
			g[12*i+11]			+= factor * diffVec.z();	

			g[12*currentId]		-= factor * diffVec.x() * mesh->nodes()[currentId]->coordinate().x();
			g[12*currentId+1]	-= factor * diffVec.y() * mesh->nodes()[currentId]->coordinate().x();
			g[12*currentId+2]	-= factor * diffVec.z() * mesh->nodes()[currentId]->coordinate().x();
			g[12*currentId+3]	-= factor * diffVec.x() * mesh->nodes()[currentId]->coordinate().y();
			g[12*currentId+4]	-= factor * diffVec.y() * mesh->nodes()[currentId]->coordinate().y();
			g[12*currentId+5]	-= factor * diffVec.z() * mesh->nodes()[currentId]->coordinate().y();
			g[12*currentId+6]	-= factor * diffVec.x() * mesh->nodes()[currentId]->coordinate().z();
			g[12*currentId+7]	-= factor * diffVec.y() * mesh->nodes()[currentId]->coordinate().z();
			g[12*currentId+8]	-= factor * diffVec.z() * mesh->nodes()[currentId]->coordinate().z();
			g[12*currentId+9]	-= factor * diffVec.x();
			g[12*currentId+10]	-= factor * diffVec.y();
			g[12*currentId+11]	-= factor * diffVec.z();	

			factor = 2.0 * weight;
			for(k = 0; k < 12; k++) g[12*i+k] += factor * (transformationParameters[12*i+k] - transformationParameters[12*currentId+k]);	
			for(k = 0; k < 12; k++) g[12*currentId+k] -= factor * (transformationParameters[12*i+k] - transformationParameters[12*currentId+k]);	
		}
	}
}

void FineFittingTracking::computeMdNeighbors(g_Part * mesh, vector<int> * indices)
{
	int i, j, k, numNN, id;
	double radius, distance, threshold;
	double * barycenters = new double[6];
	bool found;
	set< vector<int> > triangleCandidates;
	g_Vector closePt1, closePt2;
	g_Element t1, t2;
	g_Node * elementNodes = new g_Node[6];
	t1.node(&(elementNodes[0])); t1.node(&(elementNodes[1])); t1.node(&(elementNodes[2]));
	t2.node(&(elementNodes[3])); t2.node(&(elementNodes[4])); t2.node(&(elementNodes[5]));
	TriangleMesh tempMesh (*mesh);	

	//Search within the maximum edge length (guarantee that no intersections are missed):
	radius = 0;
	g_PEdgeContainer allEdges = tempMesh.uniqueEdges();
	for(i = 0; i < allEdges.numberOfItems(); i++)
	{
		radius = max(radius, allEdges[i]->firstNode().coordinate().DistanceTo(allEdges[i]->lastNode().coordinate()));
		delete allEdges[i];
	}

	threshold = tempMesh.getMeshResolution() / 50.0;

	triangleIdsForMd.clear();
	barycentersForMd.clear();

	//create kd-tree and find nearest neighbors
	ANNpointArray		dataPtsForNN;
	ANNpoint			queryPtForNN;
	ANNidxArray			nnIdxForNN;	
	ANNdistArray		distsForNN;
	ANNkd_tree*			kdTreeForNN;	
	queryPtForNN = annAllocPt(3);				
	dataPtsForNN = annAllocPts(tempMesh.getNumberOfNodes(), 3);		
	nnIdxForNN = new ANNidx[tempMesh.getNumberOfNodes()];
	distsForNN = new ANNdist[tempMesh.getNumberOfNodes()];

	g_Vector transformedPoint;
	for(i = 0; i < tempMesh.getNumberOfNodes(); i++)
	{
		transformedPoint = computeTransformedPoint(mesh, i);
		dataPtsForNN[i][0] = transformedPoint.x();
		dataPtsForNN[i][1] = transformedPoint.y();
		dataPtsForNN[i][2] = transformedPoint.z();
	}

	kdTreeForNN = new ANNkd_tree(dataPtsForNN, tempMesh.getNumberOfNodes(), 3);

	//do spherical range search and look at all the triangle lists of close neighbors
	for(i = 0; i < tempMesh.getNumberOfNodes(); i++)
	{
		if(indices == NULL) found = true;
		else 
		{
			found = false;
			for(j = 0; j < indices->size(); j++)
			{
				if((*indices)[j] == i)
				{
					found = true;
					break;
				}
			}
		}
		if(found)
		{
			queryPtForNN[0] = dataPtsForNN[i][0];
			queryPtForNN[1] = dataPtsForNN[i][1];
			queryPtForNN[2] = dataPtsForNN[i][2];
			numNN = kdTreeForNN->annkFRSearch(queryPtForNN, radius*radius, tempMesh.getNumberOfNodes(), nnIdxForNN, distsForNN);

			//Find all the triangles to consider:
			g_ElementContainer neighborTriangles;
			g_PtrAdapter<g_Element *> adapter(neighborTriangles);
			for(j = 0; j < numNN; j++) 
			{
				if(nnIdxForNN[j] != i) neighborTriangles.append(tempMesh.nodes()[nnIdxForNN[j]]->elements());
			}
			adapter.removeDuplicates();

			//Store the triangle candidates (this avoids duplicates):
			for(j = 0; j < tempMesh.nodes()[i]->elements().numberOfItems(); j++)
			{
				for(k = 0; k < neighborTriangles.numberOfItems(); k++)
				{
					vector<int> currentCandidate;
					if(tempMesh.nodes()[i]->elements()[j]->id()-1 < neighborTriangles[k]->id()-1)
					{
						currentCandidate.push_back(tempMesh.nodes()[i]->elements()[j]->id()-1);
						currentCandidate.push_back(neighborTriangles[k]->id()-1);
					}
					else
					{
						currentCandidate.push_back(neighborTriangles[k]->id()-1);
						currentCandidate.push_back(tempMesh.nodes()[i]->elements()[j]->id()-1);				
					}
					triangleCandidates.insert(currentCandidate);
				}
			}
		}
	}

	//Go through all the triangle candidates:
	set< vector<int> >::iterator itOuter;
	for(itOuter = triangleCandidates.begin(); itOuter != triangleCandidates.end(); itOuter++)
	{
		id = tempMesh.elements()[(*itOuter)[0]]->nodes()[0]->id()-1;
		t1.nodes()[0]->coordinate(g_Vector(dataPtsForNN[id][0], dataPtsForNN[id][1], dataPtsForNN[id][2]));
		id = tempMesh.elements()[(*itOuter)[0]]->nodes()[1]->id()-1;
		t1.nodes()[1]->coordinate(g_Vector(dataPtsForNN[id][0], dataPtsForNN[id][1], dataPtsForNN[id][2]));
		id = tempMesh.elements()[(*itOuter)[0]]->nodes()[2]->id()-1;
		t1.nodes()[2]->coordinate(g_Vector(dataPtsForNN[id][0], dataPtsForNN[id][1], dataPtsForNN[id][2]));

		id = tempMesh.elements()[(*itOuter)[1]]->nodes()[0]->id()-1;
		t2.nodes()[0]->coordinate(g_Vector(dataPtsForNN[id][0], dataPtsForNN[id][1], dataPtsForNN[id][2]));
		id = tempMesh.elements()[(*itOuter)[1]]->nodes()[1]->id()-1;
		t2.nodes()[1]->coordinate(g_Vector(dataPtsForNN[id][0], dataPtsForNN[id][1], dataPtsForNN[id][2]));
		id = tempMesh.elements()[(*itOuter)[1]]->nodes()[2]->id()-1;
		t2.nodes()[2]->coordinate(g_Vector(dataPtsForNN[id][0], dataPtsForNN[id][1], dataPtsForNN[id][2]));

		if((t1.nodes()[0]->id() == t2.nodes()[0]->id()) || (t1.nodes()[0]->id() == t2.nodes()[1]->id()) || (t1.nodes()[0]->id() == t2.nodes()[2]->id()) ||
			(t1.nodes()[1]->id() == t2.nodes()[0]->id()) || (t1.nodes()[1]->id() == t2.nodes()[1]->id()) || (t1.nodes()[1]->id() == t2.nodes()[2]->id()) ||
			(t1.nodes()[2]->id() == t2.nodes()[0]->id()) || (t1.nodes()[2]->id() == t2.nodes()[1]->id()) || (t1.nodes()[2]->id() == t2.nodes()[2]->id())) 
			id = -1;
		else distance = MinimumDistance::minimumDistance(t1, t2, &id, &barycenters);

		//Only consider triangles that do not share a vertex and that are close:
		if((id != -1) && (distance < threshold))
		{
			triangleIdsForMd.push_back((*itOuter)[0]);
			triangleIdsForMd.push_back((*itOuter)[1]);
			barycentersForMd.push_back(barycenters[0]); barycentersForMd.push_back(barycenters[1]); barycentersForMd.push_back(barycenters[2]); 
			barycentersForMd.push_back(barycenters[3]); barycentersForMd.push_back(barycenters[4]); barycentersForMd.push_back(barycenters[5]);
		}
	}

	delete [] nnIdxForNN;
	delete [] distsForNN;
	annDeallocPt(queryPtForNN);
	annDeallocPts(dataPtsForNN);
	delete kdTreeForNN;
	delete [] elementNodes;
	delete [] barycenters;
}

double FineFittingTracking::mdEnergy(g_Part * mesh)
{
	int i;
	int id[6];
	double energy, threshold, weight;
	g_Vector p1, p2, p3, p4, p5, p6, closePt1, closePt2, diffVec;
	TriangleMesh tempMesh(*mesh);

	energy = 0;
	threshold = tempMesh.getMeshResolution() / 50.0;

	//compute the gradient using TriangleMesh indicesForIntersection
	for(i = 0; i < (int)(triangleIdsForMd.size()/2.0); i++)
	{
		id[0] = mesh->elements()[triangleIdsForMd[2*i]]->nodes()[0]->id()-1;
		p1 = computeTransformedPoint(mesh, id[0]);
		id[1] = mesh->elements()[triangleIdsForMd[2*i]]->nodes()[1]->id()-1;
		p2 = computeTransformedPoint(mesh, id[1]);
		id[2] = mesh->elements()[triangleIdsForMd[2*i]]->nodes()[2]->id()-1;
		p3 = computeTransformedPoint(mesh, id[2]);
		closePt1 = barycentersForMd[6*i] * p1 + barycentersForMd[6*i+1] * p2 + barycentersForMd[6*i+2] * p3;

		id[3] = mesh->elements()[triangleIdsForMd[2*i+1]]->nodes()[0]->id()-1;
		p4 = computeTransformedPoint(mesh, id[3]);
		id[4] = mesh->elements()[triangleIdsForMd[2*i+1]]->nodes()[1]->id()-1;
		p5 = computeTransformedPoint(mesh, id[4]);
		id[5] = mesh->elements()[triangleIdsForMd[2*i+1]]->nodes()[2]->id()-1;
		p6 = computeTransformedPoint(mesh, id[5]);
		closePt2 = barycentersForMd[6*i+3] * p4 + barycentersForMd[6*i+4] * p5 + barycentersForMd[6*i+5] * p6;

		//Use the product rule: derivative of weight * 1/distance + derivative of 1/distance * weight
		diffVec = closePt1 - closePt2;
		weight = max(0.0, 1.0 - (diffVec.SquaredLength() / ((double)max(threshold*threshold, MD_DEFAULT_EPSILON))));
		energy += weight * (1.0 / max(diffVec.SquaredLength(), MD_DEFAULT_EPSILON)); 
	}

	return energy;
}

void FineFittingTracking::mdGradient(g_Part * mesh, double *& g)
{
	int i, j;
	int id[6];
	double factor, weight, threshold;
	g_Vector p1, p2, p3, p4, p5, p6, closePt1, closePt2, diffVec;
	TriangleMesh tempMesh(*mesh);

	threshold = tempMesh.getMeshResolution() / 50.0;

	for(i = 0; i < 12*tempMesh.getNumberOfNodes(); i++) g[i] = 0.0;

	//compute the gradient using TriangleMesh indicesForIntersection
	for(i = 0; i < (int)(triangleIdsForMd.size()/2.0); i++)
	{
		id[0] = mesh->elements()[triangleIdsForMd[2*i]]->nodes()[0]->id()-1;
		p1 = computeTransformedPoint(mesh, id[0]);
		id[1] = mesh->elements()[triangleIdsForMd[2*i]]->nodes()[1]->id()-1;
		p2 = computeTransformedPoint(mesh, id[1]);
		id[2] = mesh->elements()[triangleIdsForMd[2*i]]->nodes()[2]->id()-1;
		p3 = computeTransformedPoint(mesh, id[2]);
		closePt1 = barycentersForMd[6*i] * p1 + barycentersForMd[6*i+1] * p2 + barycentersForMd[6*i+2] * p3;

		id[3] = mesh->elements()[triangleIdsForMd[2*i+1]]->nodes()[0]->id()-1;
		p4 = computeTransformedPoint(mesh, id[3]);
		id[4] = mesh->elements()[triangleIdsForMd[2*i+1]]->nodes()[1]->id()-1;
		p5 = computeTransformedPoint(mesh, id[4]);
		id[5] = mesh->elements()[triangleIdsForMd[2*i+1]]->nodes()[2]->id()-1;
		p6 = computeTransformedPoint(mesh, id[5]);
		closePt2 = barycentersForMd[6*i+3] * p4 + barycentersForMd[6*i+4] * p5 + barycentersForMd[6*i+5] * p6;

		//Use the product rule: derivative of weight * 1/distance + derivative of 1/distance * weight
		diffVec = closePt1 - closePt2;
		weight = max(0.0, 1.0 - (diffVec.SquaredLength() / ((double)max(threshold*threshold, MD_DEFAULT_EPSILON))));
		factor = (-2.0/((double)max(threshold*threshold, MD_DEFAULT_EPSILON))) * (1.0 / max(diffVec.SquaredLength(), MD_DEFAULT_EPSILON)); 

		//Gradient affects up to 6 points
		//First part of derivative
		for(j = 0; j < 3; j++)
		{
			if(barycentersForMd[6*i+j] > 0)
			{
				g[12*id[j]]				+= barycentersForMd[6*i+j] * factor * diffVec.x() * mesh->nodes()[id[j]]->coordinate().x();
				g[12*id[j]+1]			+= barycentersForMd[6*i+j] * factor * diffVec.y() * mesh->nodes()[id[j]]->coordinate().x();
				g[12*id[j]+2]			+= barycentersForMd[6*i+j] * factor * diffVec.z() * mesh->nodes()[id[j]]->coordinate().x();
				g[12*id[j]+3]			+= barycentersForMd[6*i+j] * factor * diffVec.x() * mesh->nodes()[id[j]]->coordinate().y();
				g[12*id[j]+4]			+= barycentersForMd[6*i+j] * factor * diffVec.y() * mesh->nodes()[id[j]]->coordinate().y();
				g[12*id[j]+5]			+= barycentersForMd[6*i+j] * factor * diffVec.z() * mesh->nodes()[id[j]]->coordinate().y();
				g[12*id[j]+6]			+= barycentersForMd[6*i+j] * factor * diffVec.x() * mesh->nodes()[id[j]]->coordinate().z();
				g[12*id[j]+7]			+= barycentersForMd[6*i+j] * factor * diffVec.y() * mesh->nodes()[id[j]]->coordinate().z();
				g[12*id[j]+8]			+= barycentersForMd[6*i+j] * factor * diffVec.z() * mesh->nodes()[id[j]]->coordinate().z();
				g[12*id[j]+9]			+= barycentersForMd[6*i+j] * factor * diffVec.x();
				g[12*id[j]+10]			+= barycentersForMd[6*i+j] * factor * diffVec.y();
				g[12*id[j]+11]			+= barycentersForMd[6*i+j] * factor * diffVec.z();	
			}
		}

		for(j = 3; j < 6; j++)
		{
			if(barycentersForMd[6*i+j] > 0)
			{
				g[12*id[j]]				-= barycentersForMd[6*i+j] * factor * diffVec.x() * mesh->nodes()[id[j]]->coordinate().x();
				g[12*id[j]+1]			-= barycentersForMd[6*i+j] * factor * diffVec.y() * mesh->nodes()[id[j]]->coordinate().x();
				g[12*id[j]+2]			-= barycentersForMd[6*i+j] * factor * diffVec.z() * mesh->nodes()[id[j]]->coordinate().x();
				g[12*id[j]+3]			-= barycentersForMd[6*i+j] * factor * diffVec.x() * mesh->nodes()[id[j]]->coordinate().y();
				g[12*id[j]+4]			-= barycentersForMd[6*i+j] * factor * diffVec.y() * mesh->nodes()[id[j]]->coordinate().y();
				g[12*id[j]+5]			-= barycentersForMd[6*i+j] * factor * diffVec.z() * mesh->nodes()[id[j]]->coordinate().y();
				g[12*id[j]+6]			-= barycentersForMd[6*i+j] * factor * diffVec.x() * mesh->nodes()[id[j]]->coordinate().z();
				g[12*id[j]+7]			-= barycentersForMd[6*i+j] * factor * diffVec.y() * mesh->nodes()[id[j]]->coordinate().z();
				g[12*id[j]+8]			-= barycentersForMd[6*i+j] * factor * diffVec.z() * mesh->nodes()[id[j]]->coordinate().z();
				g[12*id[j]+9]			-= barycentersForMd[6*i+j] * factor * diffVec.x();
				g[12*id[j]+10]			-= barycentersForMd[6*i+j] * factor * diffVec.y();
				g[12*id[j]+11]			-= barycentersForMd[6*i+j] * factor * diffVec.z();		
			}
		}

		//Second part of derivative
		factor = -2.0 * weight;

		for(j = 0; j < 3; j++)
		{
			if(barycentersForMd[6*i+j] > 0)
			{	
				g[12*id[j]]				+= barycentersForMd[6*i+j] * factor * (1.0 / max(diffVec.x(), MD_DEFAULT_EPSILON)) * mesh->nodes()[id[j]]->coordinate().x();
				g[12*id[j]+1]			+= barycentersForMd[6*i+j] * factor * (1.0 / max(diffVec.y(), MD_DEFAULT_EPSILON)) * mesh->nodes()[id[j]]->coordinate().x();
				g[12*id[j]+2]			+= barycentersForMd[6*i+j] * factor * (1.0 / max(diffVec.z(), MD_DEFAULT_EPSILON)) * mesh->nodes()[id[j]]->coordinate().x();
				g[12*id[j]+3]			+= barycentersForMd[6*i+j] * factor * (1.0 / max(diffVec.x(), MD_DEFAULT_EPSILON)) * mesh->nodes()[id[j]]->coordinate().y();
				g[12*id[j]+4]			+= barycentersForMd[6*i+j] * factor * (1.0 / max(diffVec.y(), MD_DEFAULT_EPSILON)) * mesh->nodes()[id[j]]->coordinate().y();
				g[12*id[j]+5]			+= barycentersForMd[6*i+j] * factor * (1.0 / max(diffVec.z(), MD_DEFAULT_EPSILON)) * mesh->nodes()[id[j]]->coordinate().y();
				g[12*id[j]+6]			+= barycentersForMd[6*i+j] * factor * (1.0 / max(diffVec.x(), MD_DEFAULT_EPSILON)) * mesh->nodes()[id[j]]->coordinate().z();
				g[12*id[j]+7]			+= barycentersForMd[6*i+j] * factor * (1.0 / max(diffVec.y(), MD_DEFAULT_EPSILON)) * mesh->nodes()[id[j]]->coordinate().z();
				g[12*id[j]+8]			+= barycentersForMd[6*i+j] * factor * (1.0 / max(diffVec.z(), MD_DEFAULT_EPSILON)) * mesh->nodes()[id[j]]->coordinate().z();
				g[12*id[j]+9]			+= barycentersForMd[6*i+j] * factor * (1.0 / max(diffVec.x(), MD_DEFAULT_EPSILON));
				g[12*id[j]+10]			+= barycentersForMd[6*i+j] * factor * (1.0 / max(diffVec.y(), MD_DEFAULT_EPSILON));
				g[12*id[j]+11]			+= barycentersForMd[6*i+j] * factor * (1.0 / max(diffVec.z(), MD_DEFAULT_EPSILON));	
			}
		}

		for(j = 3; j < 6; j++)
		{
			if(barycentersForMd[6*i+j] > 0)
			{	
				g[12*id[j]]				-= barycentersForMd[6*i+j] * factor * (1.0 / max(diffVec.x(), MD_DEFAULT_EPSILON)) * mesh->nodes()[id[j]]->coordinate().x();
				g[12*id[j]+1]			-= barycentersForMd[6*i+j] * factor * (1.0 / max(diffVec.y(), MD_DEFAULT_EPSILON)) * mesh->nodes()[id[j]]->coordinate().x();
				g[12*id[j]+2]			-= barycentersForMd[6*i+j] * factor * (1.0 / max(diffVec.z(), MD_DEFAULT_EPSILON)) * mesh->nodes()[id[j]]->coordinate().x();
				g[12*id[j]+3]			-= barycentersForMd[6*i+j] * factor * (1.0 / max(diffVec.x(), MD_DEFAULT_EPSILON)) * mesh->nodes()[id[j]]->coordinate().y();
				g[12*id[j]+4]			-= barycentersForMd[6*i+j] * factor * (1.0 / max(diffVec.y(), MD_DEFAULT_EPSILON)) * mesh->nodes()[id[j]]->coordinate().y();
				g[12*id[j]+5]			-= barycentersForMd[6*i+j] * factor * (1.0 / max(diffVec.z(), MD_DEFAULT_EPSILON)) * mesh->nodes()[id[j]]->coordinate().y();
				g[12*id[j]+6]			-= barycentersForMd[6*i+j] * factor * (1.0 / max(diffVec.x(), MD_DEFAULT_EPSILON)) * mesh->nodes()[id[j]]->coordinate().z();
				g[12*id[j]+7]			-= barycentersForMd[6*i+j] * factor * (1.0 / max(diffVec.y(), MD_DEFAULT_EPSILON)) * mesh->nodes()[id[j]]->coordinate().z();
				g[12*id[j]+8]			-= barycentersForMd[6*i+j] * factor * (1.0 / max(diffVec.z(), MD_DEFAULT_EPSILON)) * mesh->nodes()[id[j]]->coordinate().z();
				g[12*id[j]+9]			-= barycentersForMd[6*i+j] * factor * (1.0 / max(diffVec.x(), MD_DEFAULT_EPSILON));
				g[12*id[j]+10]			-= barycentersForMd[6*i+j] * factor * (1.0 / max(diffVec.y(), MD_DEFAULT_EPSILON));
				g[12*id[j]+11]			-= barycentersForMd[6*i+j] * factor * (1.0 / max(diffVec.z(), MD_DEFAULT_EPSILON));	
			}
		}
	}
}

g_Vector FineFittingTracking::computeTransformedPoint(g_Part * mesh, int index)
{
	g_Vector result;
	if(transformationParameters == NULL) return result;
	else
	{
		double x, y, z;
		x = transformationParameters[12*index] * mesh->nodes()[index]->coordinate().x() + transformationParameters[12*index+3] * mesh->nodes()[index]->coordinate().y() +
			transformationParameters[12*index+6] * mesh->nodes()[index]->coordinate().z() + transformationParameters[12*index+9];
		y = transformationParameters[12*index+1] * mesh->nodes()[index]->coordinate().x() + transformationParameters[12*index+4] * mesh->nodes()[index]->coordinate().y() +
			transformationParameters[12*index+7] * mesh->nodes()[index]->coordinate().z() + transformationParameters[12*index+10];
		z = transformationParameters[12*index+2] * mesh->nodes()[index]->coordinate().x() + transformationParameters[12*index+5] * mesh->nodes()[index]->coordinate().y() +
			transformationParameters[12*index+8] * mesh->nodes()[index]->coordinate().z() + transformationParameters[12*index+11];
		result.Set(x, y, z);
		return result;
	}
}
	
