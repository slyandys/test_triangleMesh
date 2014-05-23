#include "FineFittingTrackingLocalFrame.h"
#include "icosahedron.h"

#include <sstream>
#include <iostream>
#include <string>

// #define _DEBUG_FINE_FITTING_TRACKING_LOCAL_FRAME_

#include "distPoint3Line3.h"


FineFittingTrackingLocalFrame::FineFittingTrackingLocalFrame() : useProbeConstraint(false)
{
	smallestMeshNum = 1000;
	templateMesh = NULL;
	frameMesh = NULL;
	smallestMesh = NULL;
	hierarchy.clear();
	transformationParameters = NULL;
	kdInitialized = false;
	allEdges.clear();
	resultMesh = NULL;
	isFeatureInHierarchy.clear();
}

FineFittingTrackingLocalFrame::~FineFittingTrackingLocalFrame()
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
	if(resultMesh != NULL)
	{
		delete resultMesh;
		resultMesh = NULL;
	}
}

vector<int> FineFittingTrackingLocalFrame::getNearestNeighbors(bool excludeBoundary)
{
	computeNearestNeighbors(templateMesh, meshNormals[0], NULL, false, excludeBoundary); 

	return nearestNeighbors;
}

void FineFittingTrackingLocalFrame::exportTransformationParameters(char * filename)
{
	if(templateMesh == NULL) return;

	FILE * fp = fopen(filename, "w");
	for (int i = 0; i < (int)templateMesh->getNumberOfNodes() * 7; i++)
		fprintf(fp, "%f ", transformationParameters[i]);
	fclose(fp);
}

void FineFittingTrackingLocalFrame::setTemplate(TriangleMesh * mesh)
{
	templateMesh = mesh;
	largestDist = MAX_DIST * templateMesh->getMeshResolution(); 

	// JL: Moved calculation of vertex normals up here
	templateMesh->calculateVertexNormals(true);

/*/FOR START AFTER A FEW FRAMES
int numNodes = templateMesh->getNumberOfNodes();
int numLevels = (int) ceil(log((double)numNodes / (double)smallestMeshNum))+1;
if(hierarchy.size() == 0) computeHierarchy(numLevels);
numLevels = (int)hierarchy.size();
//END START AFTER A FEW FRAMES*/
}

void FineFittingTrackingLocalFrame::resetTemplateCoo(g_NodeContainer nodes, char * filename)
{
	int i, j;

	g_Part * currentMesh, * previousMesh;

	//Set new coordinates to the template
	if(templateMesh == NULL) return;
	else if(templateMesh->getNumberOfNodes() > (int)nodes.numberOfItems()) return;
	else
	{
	  const g_NodeContainer& mNodes = templateMesh->nodes();
	  i=0;
	  for ( g_NodeContainer::const_iterator nIt = mNodes.begin();
		nIt != mNodes.end(); ++nIt ) {
	    (*nIt)->coordinate(nodes[i]->coordinate());
	    ++i;
	  }
	}

	//Set new coordinates to the hierarchy
	if(hierarchy.size() > 0)
	{
		for(i = 0; i < (int)hierarchy.size()+1; i++)
		{
			if(i < hierarchy.size()) currentMesh = hierarchy[i]->getOriginalMesh();
			else currentMesh = smallestMesh;
			const g_NodeContainer& currentNodes = currentMesh->nodes();


			if(i == 0)
			{
				for(j = 0; j < (int)currentNodes.numberOfItems(); j++)
					currentNodes[j]->coordinate(nodes[j]->coordinate());
			}

			else
			{
				vector<int> indices = hierarchy[i-1]->getAdditionalIndices();
				for(j = 0; j < (int)currentNodes.numberOfItems(); j++)
					currentNodes[j]->coordinate(previousMesh->nodes()[indices[j]]->coordinate());
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

		//Compute and store the normal vectors:
		meshNormals.clear();
		meshNormals.resize(hierarchy.size()+1);
		for(i = 0; i < hierarchy.size()+1; i++)
		{
			if(i == 0)
			{
				templateMesh->calculateVertexNormals(true);
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
	}
}

void FineFittingTrackingLocalFrame::setFrame(TriangleMesh * mesh, g_Vector probePos, 
					     g_Vector probeOrient, double depth, double radius )
{
  // radius *= 4;
	frameMesh = mesh;
	//Initialize the kd-tree:
	if(kdInitialized)
	{
		delete kdTree;
		annDeallocPts(dataPts);
		annClose();	
		kdInitialized = false;
	}
	if (radius > 0 ) {
	  // Delete original nodes within probe
	  int nCnt = 0;
	  TriangleMesh* revMesh = new TriangleMesh();
	  g_NodeContainer frameNodes = frameMesh->nodes();
	  g_Vector notUsed, probeOrientUp = -1.0 * probeOrient;
	  for ( int id =0; id < frameMesh->getNumberOfNodes(); ++id ) {
	    g_Vector nodePos = frameNodes[id]->coordinate();
	    // Check if frame data away from probe
	    // if ( probePos.DistanceTo(nodePos) > 2.0*radius ) { // Go small 1.2
	    if (// ( probePos.DistanceTo(nodePos) > 1.2*radius ) &&
		( Wm5::distPoint3Segment3Squared( nodePos, probePos, probeOrientUp,
						  20.0, notUsed ) > 1.5*1.5*radius*radius ))
	    { 
	      revMesh->node(new g_Node(nodePos));
	      // get normal at node with id
	      g_Vector normal = frameMesh->getNormalAtNode(id+1);
	      revMesh->SetNormals(nCnt, 0, normal.x());
	      revMesh->SetNormals(nCnt, 1, normal.y());
	      revMesh->SetNormals(nCnt, 2, normal.z());
	      ++nCnt;
	    } // else drop node
	  }
	  // flip mesh
	  frameMesh = revMesh;
	  probeNodeBegin = nCnt;
	  // Extra nodes for probe -- 
	  // Icosahedron icosa(1,probePos,radius); // use default resolution
	  Icosahedron icosa(4,probePos-0.2*radius*probeOrient
			    ,1.25*radius); // high resolution
	  icosa.setRotation(probeOrient);
	  std::vector<g_Vector> probeNodes, probeNormals;
	  icosa.generateTesselation(probeNodes, probeNormals, (std::min)(depth,radius/2));
#define  _DEBUG_FINE_FITTING_TRACKING_LOCAL_FRAME_SPHERE
#ifdef  _DEBUG_FINE_FITTING_TRACKING_LOCAL_FRAME_SPHERE
	  // For debugging -- write out the nodes on the sphere separately
	  TriangleMesh probeSphereMesh;
	  int psmCnt = 0;
	  static int nameCnt = 0;
#endif
	  std::vector<g_Vector>::const_iterator nIter = probeNormals.begin();
	  for (std::vector<g_Vector>::const_iterator pIter=probeNodes.begin(); pIter != probeNodes.end() && nIter != probeNormals.end();
			  ++pIter,++nIter) {
		  frameMesh->node(new g_Node(*pIter));
#ifdef  _DEBUG_FINE_FITTING_TRACKING_LOCAL_FRAME
		  if ( probePos.DistanceTo(*pIter) > 1.1*radius ) {
			  cerr << "Probe Node at unexpected location: " << probePos << " <-> " << *pIter << endl;
		  }
#endif
		  // set normal at node with id
		  g_Vector normal = *nIter;
		  frameMesh->SetNormals(nCnt, 0, normal.x());
		  frameMesh->SetNormals(nCnt, 1, normal.y());
		  frameMesh->SetNormals(nCnt, 2, normal.z());
#ifdef  _DEBUG_FINE_FITTING_TRACKING_LOCAL_FRAME_SPHERE
		  probeSphereMesh.node(new g_Node(*pIter));
		  probeSphereMesh.SetNormals(psmCnt, 0, normal.x());
		  probeSphereMesh.SetNormals(psmCnt, 1, normal.y());
		  probeSphereMesh.SetNormals(psmCnt, 2, normal.z());
		  ++psmCnt;
#endif
		  ++nCnt;
	  }
#ifdef  _DEBUG_FINE_FITTING_TRACKING_LOCAL_FRAME_SPHERE
	  std::stringstream  strS;
	  strS << "sphere" << nameCnt << ".wrl" << std::endl;
	  std::string fName;
	  strS >> fName;
	  exportMeshWrapper( fName.c_str(), &probeSphereMesh );
	  strS.clear();
	  strS << "frame_modified" << nameCnt << ".wrl" << std::endl;
	  strS >> fName;
	  exportMeshWrapper( fName.c_str(), frameMesh );
	  ++nameCnt;
#endif 
	} else {
	  probeNodeBegin = frameMesh->getNumberOfNodes();
	}
	dataPts = annAllocPts(frameMesh->getNumberOfNodes(), 3);
	const g_NodeContainer& mNodes = frameMesh->nodes();
	int i=0;
	for ( g_NodeContainer::const_iterator nIt = mNodes.begin();
	      nIt != mNodes.end(); ++nIt ) {
	  g_Vector coord = (*nIt)->coordinate();
	  dataPts[i][0] = coord.x();
	  dataPts[i][1] = coord.y();
	  dataPts[i][2] = coord.z();
	  ++i;
	}
	kdTree = new ANNkd_tree(dataPts, frameMesh->getNumberOfNodes(), 3);
	kdInitialized = true;
}

bool FineFittingTrackingLocalFrame::setProbeConstraint( const g_Vector& probePos, const g_Vector& probeOrient ) {
  constraintNodesHierarchy.clear();
  if((templateMesh == NULL) || (frameMesh == NULL)) {
    useProbeConstraint = false;
    return useProbeConstraint;
  }
  std::vector<bool> constraintNodes( templateMesh->getNumberOfNodes(), false );
  //Find the closest point to the probe on the current template mesh 
  int closestIdx = -1; 
  double meshRes = templateMesh->getMeshResolution();
  double minDist = 10.0*meshRes;
  g_Vector probeUp(probeOrient);
  probeUp *= -1.0;
  const g_NodeContainer& tNodes = templateMesh->nodes();
  for(int i=0; i<templateMesh->getNumberOfNodes(); ++i) {
    double dist = probePos.DistanceTo(tNodes[i]->coordinate());
    if ((dist < minDist) && (probeUp.AngleBetween(templateMesh->getNormalAtNode(i+1)) < MAX_ANGLE/180.0*3.14159)) {
      minDist = dist;
      closestIdx = i;
    }
  }
  if (closestIdx >= 0) {
    cerr << "Closest index: " << closestIdx << endl;
    constraintNodes[closestIdx] = true;
    // Go over 1-ring neighbors of the closest vertex 
    std::set<int> indicesNeighbors;
    int numNeighboringElements = tNodes[closestIdx]->elements().numberOfItems();
    for(int j=0; j<numNeighboringElements; ++j) {
      g_Element* currentIndicentTriangle = tNodes[closestIdx]->elements()[j];
      // Add all neighbors to the set
#ifdef _DEBUG_FINE_FITTING_TRACKING_LOCAL_FRAME_
      cerr << "Adding: ";
#endif
      for ( int nCnt=0; nCnt<3; ++nCnt ) {
	if(currentIndicentTriangle->nodes()[nCnt]->id() != closestIdx+1) {
	  indicesNeighbors.insert(currentIndicentTriangle->nodes()[nCnt]->id()-1);
#ifdef _DEBUG_FINE_FITTING_TRACKING_LOCAL_FRAME_
	  cerr << currentIndicentTriangle->nodes()[nCnt]->id() << " ";
#endif
	}
      }
#ifdef _DEBUG_FINE_FITTING_TRACKING_LOCAL_FRAME_
      cerr << endl;
#endif
    }
    // Now loop over all added vertices and only keep those within the mesh resolution of the probe
    for (std::set<int>::iterator it=indicesNeighbors.begin(); it != indicesNeighbors.end(); ++it) {
    	// Use 3*mean edge length
      if ( probePos.DistanceTo(tNodes[*it]->coordinate()) <= 3.0*meshRes ) {
#ifdef _DEBUG_FINE_FITTING_TRACKING_LOCAL_FRAME_
    	  cerr << *it << endl;
#endif
    	  constraintNodes[*it] = true;
      }
    }
    useProbeConstraint = true;    
  }
// #ifdef _DEBUG_FINE_FITTING_TRACKING_LOCAL_FRAME_
  if ( useProbeConstraint ) {
	  cerr << "Probe: " << probePos << " Constraint nodes (" << constraintNodes.size() << "): " << endl;
	  int idx = 0;
	  for ( std::vector<bool>::const_iterator it=constraintNodes.begin();
		it!=constraintNodes.end(); ++it ) {
	    const g_NodeContainer& templateNodes = templateMesh->nodes();
	    if (*it) {
	      cerr << idx << ": " << templateNodes[idx]->coordinate() << endl;
	    }
	    ++idx;
	  }
  } else {
	  cerr << "No constraint node found!" << endl;
  }
// #endif
  // Do we need it for the whole hierarchy?
  constraintNodesHierarchy.push_back(constraintNodes);

  //Transfer to all hierarchy levels
  for(int i=0; i<hierarchy.size(); ++i) {
    vector<int> indices = hierarchy[i]->getAllIndices();
    vector<bool> isConstraint(indices.size());
    for(int j=0; j<indices.size(); ++j)
      isConstraint[j] = isFeatureInHierarchy[i][indices[j]];
    constraintNodesHierarchy.push_back(isConstraint);
  }
  return useProbeConstraint;
}

void FineFittingTrackingLocalFrame::setMarkers(char * markersToFirstFrame, int numMarkers, double scaleFactor)
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

	vector<int> isFeatureInTemplate (numNodes);
	for(i = 0; i < numNodes; i++) isFeatureInTemplate[i] = -1;
	
	//read the file		
	for(i = 0; i < numMarkers; i++)
	{
		fscanf(fp, "%f %f %f %f %f %f ", &(floatRead[0]), &(floatRead[1]), &(floatRead[2]), &(floatRead[3]), &(floatRead[4]), &(floatRead[5]));

		//Find the point on the current template corresponding to the feature
		currentPoint.Set(scaleFactor * floatRead[0], scaleFactor * floatRead[1], scaleFactor * floatRead[2]);
		const g_NodeContainer& templateNodes = templateMesh->nodes();
		for(j = 0; j < numNodes; j++)
		{
			dist = templateNodes[j]->coordinate().DistanceTo(currentPoint);
			if((j == 0) || (dist < minDist))
			{
				minDist = dist;
				minId1 = j;
			}
		}

		//Find the point on the data corresponding to the feature
		currentPoint.Set(scaleFactor * floatRead[3], scaleFactor * floatRead[4], scaleFactor * floatRead[5]);
		const g_NodeContainer& frameNodes = frameMesh->nodes();
		for(j = 0; j < numSamples; j++)
		{
			dist = frameNodes[j]->coordinate().DistanceTo(currentPoint);
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

void FineFittingTrackingLocalFrame::unsetMarkers()
{
	isFeatureInHierarchy.clear();
}

TriangleMesh * FineFittingTrackingLocalFrame::computeDeformedTemplate(char * filename, bool rigidAligment, bool * selfIntersectReturn, bool onlyUseMarkers)
{
	int i, j, k, numNodes, numNodesInCurrentLevel, numLevels, counter, resolutionIndex;
	g_Part * currentMesh;

	//First, align rigidly:
	if(rigidAligment) 
	{
		alignRigidlyToFrame(onlyUseMarkers);

		if(filename != NULL)
		{
			char buffer[200];
			sprintf(buffer, "%s_after_rigid_alignment.wrl", filename);
			exportMeshWrapper( buffer, templateMesh );
		}
	}

	//Second, compute multi-resolution non-rigid alignment:
	numNodes = templateMesh->getNumberOfNodes();
	numLevels = (int) ceil(log((double)numNodes / (double)smallestMeshNum))+1;

	if(hierarchy.size() == 0) computeHierarchy(numLevels);
	numLevels = (int)hierarchy.size();

	currentMesh = smallestMesh;
	TriangleMesh tempSmallest(*currentMesh);

	largestDist = MAX_DIST * tempSmallest.getMeshResolution();
	tempSmallest.calculateVertexNormals();

	numNodesInCurrentLevel = tempSmallest.getNumberOfNodes();

	if(transformationParameters != NULL) delete [] transformationParameters;
	transformationParameters = new double[numNodesInCurrentLevel * 7];

	for(i = 0; i < numNodesInCurrentLevel; i++)
	{
		transformationParameters[7*i] = 0.0;
		transformationParameters[7*i+1] = 0.0;
		transformationParameters[7*i+2] = 0.0;
		transformationParameters[7*i+3] = tempSmallest.getNormalAtNode(i+1).x();
		transformationParameters[7*i+4] = tempSmallest.getNormalAtNode(i+1).y();
		transformationParameters[7*i+5] = tempSmallest.getNormalAtNode(i+1).z();
		transformationParameters[7*i+6] = 0.0;
	}

	for(i = numLevels-1; i >= 0; i--)
	{
		vector<int> * isFeature;
		if(isFeatureInHierarchy.size() > 0) 
			isFeature = &(isFeatureInHierarchy[i+1]);
		else isFeature = NULL;
		
		vector< vector<int> > clusterIds;
		fitToData(currentMesh, meshNormals[i+1], clusterIds, isFeature, onlyUseMarkers);

		if(filename != NULL)
		{
			char buffer[200];
			TriangleMesh test(*currentMesh);	
			const g_NodeContainer& tNodes = test.nodes();
			int j=0;
			for ( g_NodeContainer::const_iterator nIt = tNodes.begin();
			      nIt != tNodes.end(); ++nIt ) {
			  (*nIt)->coordinate(computeTransformedPoint(&(transformationParameters[7*j]),
								     (*nIt)->coordinate()));
			  ++j;
			}
			sprintf(buffer, "%s_after_deformation_level_%d.wrl", filename, i);
			exportMeshWrapper( buffer, &test );
		}

		//Update both the current mesh and the transformation parameters:
		currentMesh = hierarchy[i]->getOriginalMesh();
		TriangleMesh tempMesh(*currentMesh);

		tempMesh.calculateVertexNormals();

		numNodesInCurrentLevel = tempMesh.getNumberOfNodes();

		double * newTrafoParams = new double[numNodesInCurrentLevel * 7];

		vector<int> indices = hierarchy[i]->getAdditionalIndices();
		
		for(j = 0; j < numNodesInCurrentLevel; j++)
		{
			newTrafoParams[7*j] = 0.0;
			newTrafoParams[7*j+1] = 0.0;
			newTrafoParams[7*j+2] = 0.0;
			newTrafoParams[7*j+3] = tempMesh.getNormalAtNode(j+1).x();
			newTrafoParams[7*j+4] = tempMesh.getNormalAtNode(j+1).y();
			newTrafoParams[7*j+5] = tempMesh.getNormalAtNode(j+1).z();
			newTrafoParams[7*j+6] = 0.0;
		}

		//Transformations from previous level:
		for(j = 0; j < indices.size(); j++)
		{
			for(k = 0; k < 7; k++) newTrafoParams[7*indices[j]+k] = transformationParameters[7*j+k];
		}

		delete [] transformationParameters;
		transformationParameters = newTrafoParams;

		smoothTransformations(i, filename);

		largestDist = MAX_DIST * tempMesh.getMeshResolution();

		if(filename != NULL)
		{
			char buffer[200];
			TriangleMesh test(*currentMesh);	
			const g_NodeContainer& tNodes = test.nodes();
			int j=0;
			for ( g_NodeContainer::const_iterator nIt = tNodes.begin();
			      nIt != tNodes.end(); ++nIt ) {
			  (*nIt)->coordinate(computeTransformedPoint(&(transformationParameters[7*j]),
								     (*nIt)->coordinate()));
			  ++j;
			}
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
	
	fitToData(currentMesh, meshNormals[0], clusterIds, isFeature, onlyUseMarkers);

	if(filename != NULL)
	{
		char buffer[200];
		TriangleMesh test(*currentMesh);	
		const g_NodeContainer& tNodes = test.nodes();
		int j=0;
		for ( g_NodeContainer::const_iterator nIt = tNodes.begin();
		      nIt != tNodes.end(); ++nIt ) {
		  (*nIt)->coordinate(computeTransformedPoint(&(transformationParameters[7*j]),
							     (*nIt)->coordinate()));
		  ++j;
		}
		sprintf(buffer, "%s_before_post_processing.wrl", filename);
		exportMeshWrapper( buffer, &test );
	}

	//NOW: Return the mesh to use with FEM
	if(selfIntersectReturn != NULL) 
	{
int t = (int)time(0);
		//Test for self-intersections:
		TriangleMesh checkIntersectionMesh(*currentMesh);	
		const g_NodeContainer& ciNodes = checkIntersectionMesh.nodes();
		int j=0;
		for ( g_NodeContainer::const_iterator nIt = ciNodes.begin();
		      nIt != ciNodes.end(); ++nIt ) {
		  (*nIt)->coordinate(computeTransformedPoint(&(transformationParameters[7*j]),
							     (*nIt)->coordinate()));
		  ++j;
		}
		g_ElementContainer tElements = checkIntersectionMesh.elements();
		g_NodeContainer tNodes = checkIntersectionMesh.nodes();
		MinimumDistance dist(&tElements, &tNodes);
		// MinimumDistance dist(&(checkIntersectionMesh.elements()), &(checkIntersectionMesh.nodes()));
		bool selfIntersect = dist.selfIntersection(0);

		*selfIntersectReturn = selfIntersect;
cout<<"Time for self-intersection test "<<(int)time(0)-t<<" seconds"<<endl;
	}

	//Post-processing after tracking to avoid "spikes":
	postProcess(currentMesh);

	// Return the result:
	if(resultMesh != NULL) delete resultMesh;
	resultMesh = new TriangleMesh(*currentMesh);

	//Set the result (modify the mesh):
	const g_NodeContainer& rNodes = resultMesh->nodes();
	j=0;
	for ( g_NodeContainer::const_iterator nIt = rNodes.begin();
	      nIt != rNodes.end(); ++nIt ) {
	  (*nIt)->coordinate(computeTransformedPoint(&(transformationParameters[7*j]),
						     (*nIt)->coordinate()));
	  ++j;
	}
	return resultMesh;
}

void FineFittingTrackingLocalFrame::postProcess(g_Part * currentMesh)
{
int t = (int)time(0);
	int i, j, k, l, currentId, neighborId, numNodes, numElementsForCurrent;
	double closestRelativeNeighbor;
	double averageTransformations[7];
	g_Vector currentNode, currentTransformedNode, neighborCentroid, currentNeighbor, transformedNeighbor, transformedNeighborCentroid, rotAxis;
	
	//If a vertex is far from the centroid of the neighborhood, move the vertex by averaging over the transformations in the neighborhood:
	

	const g_NodeContainer& currentNodes = currentMesh->nodes();
	numNodes = (int) currentNodes.numberOfItems();
	for(i = 0; i < numNodes; i++)
	{
		currentId = currentNodes[i]->id();

		currentNode = currentNodes[currentId-1]->coordinate();
		currentTransformedNode = computeTransformedPoint(&(transformationParameters[7*(currentId-1)]), currentNode);
		
		neighborCentroid.Set(0, 0, 0);
		transformedNeighborCentroid.Set(0, 0, 0);

		for(l = 0; l < 7; l++) averageTransformations[l] = 0;

		closestRelativeNeighbor = 100;

		numElementsForCurrent = (int)currentNodes[currentId-1]->elements().numberOfItems();

		for(j = 0; j < numElementsForCurrent; j++)
		{
			for(k = 0; k < 3; k++)
			{
				neighborId = currentNodes[currentId-1]->elements()[j]->nodes()[k]->id();
				if(currentId != neighborId)
				{
					currentNeighbor = currentNodes[neighborId-1]->coordinate();
					neighborCentroid = neighborCentroid + currentNeighbor;

					transformedNeighbor = computeTransformedPoint(&(transformationParameters[7*(neighborId-1)]), currentNeighbor);
					transformedNeighborCentroid = transformedNeighborCentroid + transformedNeighbor;
					
					for(l = 0; l < 7; l++) 
					{
						averageTransformations[l] += transformationParameters[7*(neighborId-1)+l];
					}

					closestRelativeNeighbor = min(closestRelativeNeighbor, transformedNeighbor.DistanceTo(currentTransformedNode) / (double)currentNeighbor.DistanceTo(currentNode));
				}
			}
		}

		transformedNeighborCentroid = 1.0 / (2.0 * (double)numElementsForCurrent) * transformedNeighborCentroid;
		neighborCentroid = 1.0 / (2.0 * (double)numElementsForCurrent) * neighborCentroid;
		
		for(l = 0; l < 7; l++) averageTransformations[l] = 1.0 / (2.0 * (double)numElementsForCurrent) * averageTransformations[l];
		rotAxis = g_Vector(averageTransformations[3], averageTransformations[4], averageTransformations[5]);
		rotAxis.Normalize();
		averageTransformations[3] = rotAxis.x();
		averageTransformations[4] = rotAxis.y();
		averageTransformations[5] = rotAxis.z();

		//Use a threshold here that controls how much deformation is allowed:
		if(closestRelativeNeighbor > 2.0) 
		{
			for(l = 0; l < 7; l++) 
				transformationParameters[7*(currentId-1)+l] = averageTransformations[l];
		}
	}
cout<<"Time for postProcess "<<(int)time(0)-t<<" seconds"<<endl;
}

void FineFittingTrackingLocalFrame::adjustTransformations(g_Part * targetMesh, vector<int> indicesModifiedWithFEM, bool * selfIntersectReturn)
{
	int i, j;

	g_Part * mesh = templateMesh;

	vector<int> filteredIndices;

	const g_NodeContainer& meshNodes = mesh->nodes();
	const g_NodeContainer& targetNodes = targetMesh->nodes();
	for(i = 0; i < (int)indicesModifiedWithFEM.size(); i++)
	{
		if(meshNodes[indicesModifiedWithFEM[i]]->coordinate().DistanceTo(targetNodes[indicesModifiedWithFEM[i]]->coordinate()) < largestDist)
			filteredIndices.push_back(indicesModifiedWithFEM[i]);
	}

	int dimension = (int)indicesModifiedWithFEM.size() * 7;

	vnl_vector< double > x (dimension, 0.0);
	for(i = 0; i < (int)indicesModifiedWithFEM.size(); i++)
	{
		for(j = 0; j < 7; j++) x[7*i+j] = transformationParameters[7*indicesModifiedWithFEM[i]+j];
	}

	// Set the weights used here:
	nnWeight = 1.0;
	accessWeight = 20.0;
	mdWeight = 0.0; 
		
	//Compute the access indices:
	if(accessWeight > 0) computeAccessNeighbors(mesh, &indicesModifiedWithFEM);
	if(mdWeight > 0) computeMdNeighbors(mesh, &indicesModifiedWithFEM);

	//perform the nonlinear optimization
	AdjustTransformationCostFunction costFct(this, mesh, targetMesh, indicesModifiedWithFEM, filteredIndices, dimension);
	vnl_lbfgsb * minimizer = new vnl_lbfgsb(costFct);
//	minimizer->set_cost_function_convergence_factor ( 10000000 ); 
//	minimizer->set_projected_gradient_tolerance ( 0.00001 ); 
	minimizer->set_max_function_evals(MAX_NUM_ITER);
	minimizer->minimize( x );

	//Copy the transformations:
	for(i = 0; i < (int)indicesModifiedWithFEM.size(); i++)
	{
		for(j = 0; j < 7; j++) transformationParameters[7*indicesModifiedWithFEM[i]+j] = x[7*i+j];
	}

	delete minimizer;

	postProcess(mesh);
	
	//Change the coordinates of the input
	for(i = 0; i < (int)meshNodes.numberOfItems(); i++)
		targetNodes[i]->coordinate(computeTransformedPoint(&(transformationParameters[7*i]), meshNodes[i]->coordinate()));

	//Check for self-intersections
	if(selfIntersectReturn != NULL) 
	{
		//Test for self-intersections:
	  g_ElementContainer tElements = targetMesh->elements();
	  g_NodeContainer tNodes = targetMesh->nodes();
	  MinimumDistance dist(&tElements, &tNodes);
	  // MinimumDistance dist(&(targetMesh->elements()), &(targetMesh->nodes())); // JL
		bool selfIntersect = dist.selfIntersection(0);

		*selfIntersectReturn = selfIntersect;
	}
}

void FineFittingTrackingLocalFrame::AdjustTransformationCostFunction::compute(vnl_vector< double > const &x, double *f, vnl_vector< double > *g)
{
	int i, j, dimension;

	dimension = (int)indices.size() * 7;

	double * helpGrad1 = new double[mesh->nodes().numberOfItems() * 7];
	double * helpGrad2 = new double[mesh->nodes().numberOfItems() * 7];

	// Update the transformation parameters:
	for(i = 0; i < (int)indices.size(); i++)
	{
		for(j = 0; j < 7; j++) tracker->transformationParameters[7*indices[i]+j] = x[7*i+j];
	} 

	// Compute Energy:
	*f = 0;
	if(tracker->accessWeight > 0) *f += tracker->accessWeight * tracker->accessEnergy(mesh);
	if(tracker->mdWeight > 0) *f += tracker->mdWeight * tracker->mdEnergy(mesh);
	if(tracker->nnWeight > 0) *f += tracker->nnWeight * tracker->distanceResolutionEnergy(mesh, targetMesh, indicesLowRes);

	// Compute Gradient:
	for(i = 0; i < dimension; i++) (*g)[i] = 0;
	int meshNumNodes = mesh->nodes().numberOfItems();
	for(i = 0; i < meshNumNodes * 7; i++) helpGrad2[i] = 0;
	if(tracker->accessWeight > 0)
	{
		tracker->accessGradient(mesh, helpGrad1);
		for(i = 0; i < meshNumNodes * 7; i++) 
			helpGrad2[i] += tracker->accessWeight * helpGrad1[i];
	}
	if(tracker->mdWeight > 0)
	{
		tracker->mdGradient(mesh, helpGrad1);
		for(i = 0; i < meshNumNodes * 7; i++) 
			helpGrad2[i] += tracker->mdWeight * helpGrad1[i];
	}
	if(tracker->nnWeight > 0) 
	{
		tracker->distanceResolutionGradient(mesh, targetMesh, indicesLowRes, helpGrad1);
		for(i = 0; i < (int)indicesLowRes.size(); i++)
		{
			for(j = 0; j < 7; j++) helpGrad2[7*indicesLowRes[i]+j] += tracker->nnWeight * helpGrad1[7*i+j];
		}
	}

	for(i = 0; i < (int)indices.size(); i++)
	{
		for(j = 0; j < 7; j++) (*g)[7*i+j] = helpGrad2[7*indices[i]+j];
	}

	delete [] helpGrad1;
	delete [] helpGrad2;
}

void FineFittingTrackingLocalFrame::smoothTransformations(int index, char * filename)
{
	int i, j, counterId; 

	g_Part * mesh = hierarchy[index]->getOriginalMesh();
	vector<int> indices = hierarchy[index]->getAdditionalIndices();
	vector<int> newIndices;

	int meshNumNodes = mesh->nodes().numberOfItems();
	int dimension = (int)(meshNumNodes - indices.size()) * 7;

	counterId = 0;
	for(i = 0; i < meshNumNodes; i++)
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
		for(j = 0; j < 7; j++) x[7*i+j] = transformationParameters[7*newIndices[i]+j];
	}

	// Set the weights used here:
	accessWeight = 1.0;
	mdWeight = 0.0;
	nnWeight = 0.0;

	// Compute the access indices
	if(accessWeight > 0) computeAccessNeighbors(mesh, &newIndices);
	if(mdWeight > 0) computeMdNeighbors(mesh, &newIndices);

	/*/ Compute the mean-value coordinates of the new points:
	g_NodeContainer nodeCont;
	for(i = 0; i < indices.size(); i++)
	{
		g_Node * node = new g_Node(computeTransformedPoint(&(transformationParameters[7*indices[i]]), mesh->nodes()[indices[i]]->coordinate()));
		nodeCont.insert(node);
	}
	g_Part * reconstructedMesh = hierarchy[index]->reconstructMesh(nodeCont, indices, true);

	if(filename != NULL)
	{
		char buffer[200];
		sprintf(buffer, "%s_reconstructed_mesh_%d.wrl", filename, index);
		TriangleMesh exportMesh(*reconstructedMesh);
		exportMeshWrapper( buffer, exportMesh );
	}

	for(i = 0; i < indices.size(); i++) delete nodeCont[i];
	*/
	//Much faster:
	g_Part * reconstructedMesh = NULL;

	//perform the nonlinear optimization
	SmoothTransformationCostFunction costFct(this, mesh, reconstructedMesh, newIndices, dimension);
	vnl_lbfgsb * minimizer = new vnl_lbfgsb(costFct);
//	minimizer->set_cost_function_convergence_factor ( 10000000 );
//	minimizer->set_projected_gradient_tolerance ( 0.00001 ); 
	minimizer->set_max_function_evals(MAX_NUM_ITER);
	minimizer->minimize( x );

	//Copy the transformations:
	for(i = 0; i < newIndices.size(); i++)
	{
		for(j = 0; j < 7; j++) transformationParameters[7*newIndices[i]+j] = x[7*i+j];
	}

	//delete stuff:
	delete minimizer;
}

void FineFittingTrackingLocalFrame::SmoothTransformationCostFunction::compute(vnl_vector< double > const &x, double *f, vnl_vector< double > *g)
{
	int i, j, dimension;

	dimension = newIndices.size() * 7;

	int meshNumNodes = mesh->nodes().numberOfItems();

	double * helpGrad1 = new double[meshNumNodes * 7];
	double * helpGrad2 = new double[meshNumNodes * 7];

	// Update the transformation parameters:
	for(i = 0; i < newIndices.size(); i++)
	{
		for(j = 0; j < 7; j++) tracker->transformationParameters[7*newIndices[i]+j] = x[7*i+j];
	}

	// Compute Energy:
	*f = 0;
	if(tracker->accessWeight > 0) *f += tracker->accessWeight * tracker->accessEnergy(mesh);
	if(tracker->mdWeight > 0) *f += tracker->mdWeight * tracker->mdEnergy(mesh);
	if((tracker->nnWeight > 0) && (reconstructedMesh != NULL)) *f += tracker->nnWeight * tracker->distanceResolutionEnergy(mesh, reconstructedMesh, newIndices);

	// Compute Gradient:
	for(i = 0; i < dimension; i++) (*g)[i] = 0;
	for(i = 0; i < meshNumNodes * 7; i++) helpGrad2[i] = 0;

	if(tracker->accessWeight > 0) 
	{
		tracker->accessGradient(mesh, helpGrad1);
		for(i = 0; i < meshNumNodes * 7; i++) 
			helpGrad2[i] += tracker->accessWeight * helpGrad1[i];
	}
	if(tracker->mdWeight > 0) 
	{
		tracker->mdGradient(mesh, helpGrad1);
		for(i = 0; i < meshNumNodes * 7; i++) 
			helpGrad2[i] += tracker->mdWeight * helpGrad1[i];
	}
	if((tracker->nnWeight > 0) && (reconstructedMesh != NULL)) 
	{
		tracker->distanceResolutionGradient(mesh, reconstructedMesh, newIndices, helpGrad1);
		for(i = 0; i < newIndices.size(); i++)
		{
			for(j = 0; j < 7; j++) helpGrad2[7*newIndices[i]+j] += tracker->nnWeight * helpGrad1[7*i+j];
		}
	}

	for(i = 0; i < newIndices.size(); i++)
	{
		for(j = 0; j < 7; j++) (*g)[7*i+j] = helpGrad2[7*newIndices[i]+j];
	}

	delete [] helpGrad1;
	delete [] helpGrad2;
}

void FineFittingTrackingLocalFrame::alignRigidlyToFrame(bool onlyUseMarkers)
{
	//Check that everything is initialized
	if(!kdInitialized) return;

	vnl_vector<double> x( 8, 0.0 );
	x[0] = x[4] = x[5] = x[6] = 1.0;
	x[1] = x[2] = x[3] = x[7] = 0.0;

	//initialize the parameters for the nonlinear optimization
	RigidAlignmentCostFunction costFct(this);
	
	vnl_lbfgsb * minimizer = new vnl_lbfgsb(costFct);
//	minimizer->set_cost_function_convergence_factor ( 10000000 );
//	minimizer->set_projected_gradient_tolerance ( 0.00001 ); 
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

	const g_NodeContainer& templateNodes = templateMesh->nodes();
	int templateMeshNumNodes = templateMesh->nodes().numberOfItems();

	for(int i = 0; i < templateMeshNumNodes; i++)
	{
		queryPt[0] = trafo[0]*templateNodes[i]->coordinate().x()+trafo[4]*templateNodes[i]->coordinate().y()+trafo[8]*templateNodes[i]->coordinate().z()+trafo[12];
		queryPt[1] = trafo[1]*templateNodes[i]->coordinate().x()+trafo[5]*templateNodes[i]->coordinate().y()+trafo[9]*templateNodes[i]->coordinate().z()+trafo[13];
		queryPt[2] = trafo[2]*templateNodes[i]->coordinate().x()+trafo[6]*templateNodes[i]->coordinate().y()+trafo[10]*templateNodes[i]->coordinate().z()+trafo[14];
		templateNodes[i]->coordinate(g_Vector(queryPt[0], queryPt[1], queryPt[2]));
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

void FineFittingTrackingLocalFrame::RigidAlignmentCostFunction::compute(vnl_vector< double > const &x, double *f, vnl_vector< double > *g)
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

	// tracker will only have templateMesh with normals
	// tracker->templateMesh->calculateVertexNormals();
	// Should not be called because the framemesh has no face information and calculate vertex normals will reset all normals to 0 0 0
	// tracker->frameMesh->calculateVertexNormals();

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
	const g_NodeContainer& templateNodes = tracker->templateMesh->nodes();
	int templateMeshNumNodes = tracker->templateMesh->nodes().numberOfItems();
	const g_NodeContainer& frameNodes = tracker->frameMesh->nodes();
	for(i = 0; i < templateMeshNumNodes; i++)
	{
		//deform point and compute nearest neighbor:
		queryPt[0] = trafo[0]*templateNodes[i]->coordinate().x()+trafo[4]*templateNodes[i]->coordinate().y()+
			trafo[8]*templateNodes[i]->coordinate().z()+trafo[12];
		queryPt[1] = trafo[1]*templateNodes[i]->coordinate().x()+trafo[5]*templateNodes[i]->coordinate().y()+
			trafo[9]*templateNodes[i]->coordinate().z()+trafo[13];
		queryPt[2] = trafo[2]*templateNodes[i]->coordinate().x()+trafo[6]*templateNodes[i]->coordinate().y()+
			trafo[10]*templateNodes[i]->coordinate().z()+trafo[14];

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
					if((templateNodes[i]->coordinate().DistanceTo(frameNodes[nnIdx[0]]->coordinate()) < tracker->largestDist) && 
						(tracker->templateMesh->getNormalAtNode(i+1).AngleBetween(tracker->frameMesh->getNormalAtNode(nnIdx[0]+1)) < MAX_ANGLE/180.0*3.14159))
						fixedNN.push_back(nnIdx[0]);
					else fixedNN.push_back(-1);
				}
			}
			else
			{
				tracker->kdTree->annkPriSearch(queryPt, 1, nnIdx, dists);
				if((templateNodes[i]->coordinate().DistanceTo(frameNodes[nnIdx[0]]->coordinate()) < tracker->largestDist) && 
					(tracker->templateMesh->getNormalAtNode(i+1).AngleBetween(tracker->frameMesh->getNormalAtNode(nnIdx[0]+1)) < MAX_ANGLE/180.0*3.14159))
					fixedNN.push_back(nnIdx[0]);
				else fixedNN.push_back(-1);
			}
		}

		if(fixedNN[i] != -1)
		{
			g_Vector diffVec(g_Vector(queryPt[0], queryPt[1], queryPt[2])-frameNodes[fixedNN[i]]->coordinate());
			*f += diffVec.SquaredLength();

			//compute gradient:
			for(j = 0; j < 8; j++)
			{
				queryPt[0] = trafoGrad[16*j+0]*templateNodes[i]->coordinate().x()+trafoGrad[16*j+4]*templateNodes[i]->coordinate().y()+
					trafoGrad[16*j+8]*templateNodes[i]->coordinate().z()+trafoGrad[16*j+12];
				queryPt[1] = trafoGrad[16*j+1]*templateNodes[i]->coordinate().x()+trafoGrad[16*j+5]*templateNodes[i]->coordinate().y()+
					trafoGrad[16*j+9]*templateNodes[i]->coordinate().z()+trafoGrad[16*j+13];
				queryPt[2] = trafoGrad[16*j+2]*templateNodes[i]->coordinate().x()+trafoGrad[16*j+6]*templateNodes[i]->coordinate().y()+
					trafoGrad[16*j+10]*templateNodes[i]->coordinate().z()+trafoGrad[16*j+14];
	
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

void FineFittingTrackingLocalFrame::computeTranslationMat(double tx, double ty, double tz, double *& mat)
{
	mat[0] = 1;		mat[4] = 0;		mat[8] = 0;			mat[12] = tx;
	mat[1] = 0;		mat[5] = 1;		mat[9] = 0;			mat[13] = ty;
	mat[2] = 0;		mat[6] = 0;		mat[10] = 1;		mat[14] = tz;
	mat[3] = 0;		mat[7] = 0;		mat[11] = 0;		mat[15] = 1;
}

void FineFittingTrackingLocalFrame::computeRotationMat(double tx, double ty, double tz, double angle, double *& mat)
{
	double quaternion[4];
	g_Vector crossProd(tx, ty, tz);
	crossProd.Normalize();
	crossProd = crossProd * sin(angle/2.0);
	quaternion[0] = cos(angle/2.0);		quaternion[1] = crossProd.x();		quaternion[2] = crossProd.y();		quaternion[3] = crossProd.z();
		
	mat[0] = pow(quaternion[0],2) + pow(quaternion[1],2) - pow(quaternion[2],2) - pow(quaternion[3],2);		
	mat[4] = 2*quaternion[1]*quaternion[2] - 2*quaternion[0]*quaternion[3];			
	mat[8] = 2*quaternion[1]*quaternion[3] + 2*quaternion[0]*quaternion[2];				
	mat[12] = 0;
	mat[1] = 2*quaternion[1]*quaternion[2] + 2*quaternion[0]*quaternion[3];	
	mat[5] = pow(quaternion[0],2) - pow(quaternion[1],2) + pow(quaternion[2],2) - pow(quaternion[3],2);	
	mat[9] = 2*quaternion[2]*quaternion[3] - 2*quaternion[0]*quaternion[1];				
	mat[13] = 0;
	mat[2] = 2*quaternion[1]*quaternion[3] - 2*quaternion[0]*quaternion[2];	
	mat[6] = 2*quaternion[2]*quaternion[3] + 2*quaternion[0]*quaternion[1];	
	mat[10] = pow(quaternion[0],2) - pow(quaternion[1],2) - pow(quaternion[2],2) + pow(quaternion[3],2);	
	mat[14] = 0;
	mat[3] = 0;			mat[7] = 0;			mat[11] = 0;			mat[15] = 1;	
}

void FineFittingTrackingLocalFrame::computeScalingMat(double scale, double *& mat)
{
	mat[0] = scale;		mat[4] = 0;			mat[8] = 0;				mat[12] = 0;
	mat[1] = 0;			mat[5] = scale;		mat[9] = 0;				mat[13] = 0;
	mat[2] = 0;			mat[6] = 0;			mat[10] = scale;		mat[14] = 0;
	mat[3] = 0;			mat[7] = 0;			mat[11] = 0;			mat[15] = 1;
}

void FineFittingTrackingLocalFrame::computeTranslationGrad(double tx, double ty, double tz, double *& mat, int index)
{
	for(int i = 0; i < 16; i++) mat[i] = 0;
	if(index == 0) mat[12] = 1;
	else if(index == 1) mat[13] = 1;
	else mat[14] = 1;
}

void FineFittingTrackingLocalFrame::computeRotationGrad(double tx, double ty, double tz, double alpha, double *& mat, int index)
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

void FineFittingTrackingLocalFrame::computeScalingGrad(double s, double *& mat)
{
	mat[0] = 1;			mat[4] = 0;			mat[8] = 0;				mat[12] = 0;
	mat[1] = 0;			mat[5] = 1;			mat[9] = 0;				mat[13] = 0;
	mat[2] = 0;			mat[6] = 0;			mat[10] = 1;			mat[14] = 0;
	mat[3] = 0;			mat[7] = 0;			mat[11] = 0;			mat[15] = 0;
}

void FineFittingTrackingLocalFrame::computeHierarchy(int numLevels)
{
	//compute smallestMesh and hierarchy here
	int i, j;

	cerr << "Template #nodes: " << templateMesh->getNumberOfNodes() << endl;
	g_Part * currentMesh = (g_Part *) templateMesh;
	cerr << "Template #nodes: " << currentMesh->getNumberOfNodes() << endl;
	
	for(i = 0; i < numLevels; i++)
	{
		MeshSimplification * simplify = new MeshSimplification(currentMesh);
		currentMesh = simplify->getSimplifiedMesh((int)floor(templateMesh->getNumberOfNodes() / pow(2.0, i+1)));
		cerr << "Level " << i << " #nodes: " << currentMesh->getNumberOfNodes() << endl;
		if(currentMesh->nodes().numberOfItems() == simplify->getOriginalMesh()->nodes().numberOfItems())
		{
			//In this case, the mesh cannot be reduced further. Every edge collapse is invalid.
			break;
		}
		else {
		  std::ostringstream os;
		  os << "hierarchy_" << i << ".wrl";
		  TriangleMesh tmpMesh(*currentMesh);
		  exportMeshWrapper(os.str(), &tmpMesh);
		  hierarchy.push_back(simplify);
		}
	}

	smallestMesh = new g_Part(*currentMesh);
	cerr << "Smallest mesh #nodes: " << smallestMesh->getNumberOfNodes() << endl;

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

void FineFittingTrackingLocalFrame::fitToData(g_Part * mesh, vector<g_Vector> normalVectors, vector< vector<int> > clusterIds, vector<int> * isFeature, bool onlyUseMarkers)
{
int t = (int)time(0);
	int i;
	double energy, prevEnergy;

	//Check that everything is initialized
	if(!kdInitialized) return;
	for(i = 0; i < (int)allEdges.numberOfItems(); i++) delete allEdges[i];
	allEdges.clear();
	allEdges = mesh->uniqueEdges();

	cerr << "Computing neighbors ... ";
	computeAccessNeighbors(mesh);
	cerr << "done." << endl;

	//Initialize weights and transformations:
	// JL 
	// accessWeight = 1.0;	
	nnWeight = 1.0;
	// Apr. 5: Starting with nnWeight = accessWeight = 100 seems good.
	// nnWeight = 1.0;	
	accessWeight = 10.0;
	mdWeight = 0.0; 

	energy = -1;
	while(accessWeight > 1.0) // 20
	  // while(accessWeight > 0.99)
	{
		prevEnergy = energy;
		energy = solveOptimization(mesh, normalVectors, clusterIds, isFeature, onlyUseMarkers);
		cerr << "Access Weight: " << accessWeight << " Tracking energy: " << energy << " <- " << prevEnergy << endl;
		accessWeight = accessWeight / 2.0;

		if(prevEnergy != -1)
		{
			if(fabs(prevEnergy - energy)/prevEnergy < ENERGY_TOLERANCE) 
				break;
		}
	}
cout<<"Time for fitToData "<<(int)time(0)-t<<" seconds"<<endl;
}

double FineFittingTrackingLocalFrame::solveOptimization(g_Part * mesh, vector<g_Vector> normalVectors, vector< vector<int> > clusterIds, vector<int> * isFeature, bool onlyUseMarkers)
{
	int i, dimension;
	double energy;

	dimension = (int)mesh->nodes().numberOfItems() * 7;

	if(nnWeight > 0) computeNearestNeighbors(mesh, normalVectors, isFeature, onlyUseMarkers, true);
	if(mdWeight > 0) computeMdNeighbors(mesh);

	//initialize the parameters for the nonlinear optimization
	vnl_vector<double> x( dimension, 0.0 );
	for(i = 0; i < dimension; i++) x[i] = transformationParameters[i];

	SolveOptimizationCostFunction costFct(this, mesh, normalVectors, clusterIds, isFeature, dimension);
	
	vnl_lbfgsb * minimizer = new vnl_lbfgsb(costFct);
//	minimizer->set_cost_function_convergence_factor ( 100000000000 );
//	minimizer->set_projected_gradient_tolerance ( 0.0001 ); 
	minimizer->set_max_function_evals(MAX_NUM_ITER);
	minimizer->minimize( x );

	energy = minimizer->get_end_error();
	
	//copy the result
	for(i = 0; i < dimension; i++) transformationParameters[i] = x[i];

	//delete stuff:
	delete minimizer;

	return energy;
}

void FineFittingTrackingLocalFrame::SolveOptimizationCostFunction::compute(vnl_vector< double > const &x, double *f, vnl_vector< double > *g)
{
	int i, dimension;
	dimension = mesh->nodes().numberOfItems() * 7;

	bool isCurrentFeature;
	if(isFeature != NULL) isCurrentFeature = true;
	else isCurrentFeature = false;
	
	double * helpGrad = new double[dimension];
	// Copy the result:
	for(i = 0; i < dimension; i++) tracker->transformationParameters[i] = x[i];
	// Compute Energy:
	*f = 0;
	if(tracker->nnWeight > 0) *f += tracker->nnWeight * tracker->nearestNeighborEnergy(mesh, isCurrentFeature);
	if(tracker->accessWeight > 0) *f += tracker->accessWeight * tracker->accessEnergy(mesh);
	if(tracker->mdWeight > 0) *f += tracker->mdWeight * tracker->mdEnergy(mesh);

	// Compute Gradient:
	for(i = 0; i < dimension; i++) 
		(*g)[i] = 0;
	if(tracker->nnWeight > 0)
	{
		tracker->nearestNeighborGradient(mesh, helpGrad);
		for(i = 0; i < dimension; i++) 
			(*g)[i] += tracker->nnWeight * helpGrad[i];
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

void FineFittingTrackingLocalFrame::computeNearestNeighbors(g_Part * mesh, vector<g_Vector> normalVectors, vector<int> * isFeature, bool onlyUseMarkers,
												  bool excludeBoundary)
{
	int i, j, numNodes, numElements;
	double x, y, z;

	//Compute the nearest neighbors if not only markers are used
	if(onlyUseMarkers)
	{
		numNodes = (int)mesh->nodes().numberOfItems();
		nearestNeighbors.resize(numNodes);
		approxNearestPointPlaneNeighbors.resize(numNodes);
		for(i = 0; i < numNodes; i++) nearestNeighbors[i] = -1;
	}
	else
	{
	        const long int _n = 3;
		long int n = _n;
		long int ipiv[_n+1];
		long int info;
                double work[_n];

		// Initialize kd-tree use
		ANNpoint			queryPt;
		ANNidxArray			nnIdx;	
		ANNdistArray		dists;
		queryPt = annAllocPt(3);			
		nnIdx = new ANNidx[1];	
		dists = new ANNdist[1];

		//Transform the normals:
		numNodes = (int)mesh->nodes().numberOfItems();
		vector<g_Vector> transformedNormals(numNodes);
		double* currentTrafo = new double[16];
		for(i = 0; i < numNodes; i++)
		{
			computeRotationMat(transformationParameters[7*i+3], transformationParameters[7*i+4], transformationParameters[7*i+5], transformationParameters[7*i+6], currentTrafo); 
			// copy the 4x4 result to the 3x3 subpart
			// currentTrafo[0] = currentTrafo[0];
			// currentTrafo[1] = currentTrafo[1];
			// currentTrafo[2] = currentTrafo[2];
			currentTrafo[3] = currentTrafo[4];
			currentTrafo[4] = currentTrafo[5]; 
			currentTrafo[5] = currentTrafo[6]; 
			currentTrafo[6] = currentTrafo[8];
			currentTrafo[7] = currentTrafo[9]; 
			currentTrafo[8] = currentTrafo[10];

			double debugTrafo[9];
			for (int k=0; k<n*n; ++k) 
			  debugTrafo[k] = currentTrafo[k];
			for (int k=0; k<n; ++k) ipiv[k] = -1;
			info = -1;
			//Compute inverse:
			clapack::dgetrf_(&n, &n, currentTrafo, &n, ipiv, &info);
			if(static_cast<int>(info) != 0) {
			  cerr<<"dgetrf: Problem when transforming normals in computeNearestNeighbors(...) : info "<<info<<endl;
			  cerr << "n: " << n << endl; 
			  for (int k=0;k<n*n;++k) {
			    if (k%n==0 && k!=0) cerr <<endl;
			    cerr << debugTrafo[k] << " ";
			  }
			  cerr << endl << "=======================" << endl;
			  for (int k=0;k<n*n;++k) {
			    if (k%n==0 && k!=0) cerr <<endl;
			    cerr << currentTrafo[k] << " ";
			  }
			  cerr << endl;
			  for (int k=0;k<=n;++k) cerr << ipiv[k] << " ";
			  cerr << endl;
			}
			clapack::dgetri_(&n, currentTrafo, &n, ipiv, work, &n, &info);
			if(static_cast<int>(info) != 0) {
			  cerr<<"dgetri: Problem when transforming normals in computeNearestNeighbors(...) : info "<<info<<endl;
			  cerr << "n: " << n << endl; 
			  for (int k=0;k<n*n;++k) {
			    if (k%n==0 && k!=0) cerr <<endl;
			    cerr << debugTrafo[k] << " ";
			  }
			  cerr << endl << "=======================" << endl;
			  for (int k=0;k<n*n;++k) {
			    if (k%n==0 && k!=0) cerr <<endl;
			    cerr << currentTrafo[k] << " ";
			  }
			  cerr << endl;
			}

			
			//Compute transformed normal (use transposed currentTrafo) and store in transformedNormals:
			x = currentTrafo[0] * normalVectors[i].x() + currentTrafo[1] * normalVectors[i].y() + currentTrafo[2] * normalVectors[i].z();
			y = currentTrafo[3] * normalVectors[i].x() + currentTrafo[4] * normalVectors[i].y() + currentTrafo[5] * normalVectors[i].z();
			z = currentTrafo[6] * normalVectors[i].x() + currentTrafo[7] * normalVectors[i].y() + currentTrafo[8] * normalVectors[i].z();
			transformedNormals[i].Set(x, y, z);
		}
		delete [] currentTrafo;

		//Copy the mesh and compute its normals
		TriangleMesh tempMesh;
		for(i = 0; i < numNodes; i++)
		{
			g_Node * node = new g_Node(computeTransformedPoint(&(transformationParameters[7*i]), mesh->nodes()[i]->coordinate()));
			tempMesh.node(node);
		}

		const g_NodeContainer& tempNodes = tempMesh.nodes();
		const g_ElementContainer& meshElements = mesh->elements();
		numElements = (int)meshElements.numberOfItems();
		for(i = 0; i < numElements; i++)
		{
			g_Element * elem = new g_Element;
			elem->node(tempNodes[meshElements[i]->nodes()[0]->id()-1]);
			elem->node(tempNodes[meshElements[i]->nodes()[1]->id()-1]);
			elem->node(tempNodes[meshElements[i]->nodes()[2]->id()-1]);
			tempMesh.element(elem);
		}

		// Find the nearest neighbor for each transformed point
		nearestNeighbors.resize(numNodes);
		approxNearestPointPlaneNeighbors.resize(numNodes);
		for(i = 0; i < numNodes; i++)
		{
			queryPt[0] = tempNodes[i]->coordinate().x();
			queryPt[1] = tempNodes[i]->coordinate().y();
			queryPt[2] = tempNodes[i]->coordinate().z();

			kdTree->annkPriSearch(queryPt, 1, nnIdx, dists);
#ifdef _DEBUG_FINE_FITTING_TRACKING_LOCAL_FRAME_
			cerr << i << ": " << transformedNormals[i] << " vs " << nnIdx[0] << ": " << frameMesh->getNormalAtNode(nnIdx[0]+1) << endl;
#endif
			if ( nnIdx[0] >= probeNodeBegin ) {
				cerr << "Considering probeNode: " << i << "<->" << nnIdx[0] << " with " << sqrt(dists[0]) << endl;
				cerr << transformedNormals[i] << " <-> " << frameMesh->getNormalAtNode(nnIdx[0]+1) << " alpha=" <<
						transformedNormals[i].AngleBetween(frameMesh->getNormalAtNode(nnIdx[0]+1))/3.14159*180.0 << endl;
			}
			if((sqrt(dists[0]) < largestDist) && (transformedNormals[i].AngleBetween(frameMesh->getNormalAtNode(nnIdx[0]+1)) < MAX_ANGLE/180.0*3.14159)) {
				nearestNeighbors[i] = nnIdx[0];
				if ( nnIdx[0] >= probeNodeBegin ) {
					cerr << "Match to probeNode: " << i << "<->" << nnIdx[0] << endl;
				}
			}
			else nearestNeighbors[i] = -1;
		}

		// If excludeBoundary is true, only keep matches that do not match to the boundary of the input point cloud (heuristic)
		if(excludeBoundary)
		{
			vector<int> counter ((int)frameMesh->getNumberOfNodes());
			for(i = 0; i < (int)frameMesh->getNumberOfNodes(); i++) counter[i] = 0;
			for(i = 0; i < numNodes; i++) 
			{
				if(nearestNeighbors[i] != -1)
					counter[nearestNeighbors[i]]++;
			}

			double averageCount = 0;
			int numInterest = 0;

			for(i = 0; i < (int)frameMesh->getNumberOfNodes(); i++)
			{
				if(counter[i] > 0)
				{
					averageCount += counter[i];
					numInterest++;
				}
			}

			averageCount = averageCount / (double)numInterest;

			for(i = 0; i < numNodes; i++)
			{
				if(nearestNeighbors[i] != -1)
				{
					if(counter[nearestNeighbors[i]] > 2.0 * averageCount)
						nearestNeighbors[i] = -1;
				}
			}
		}

		// Free space
		annDeallocPt(queryPt);
		delete [] nnIdx;
		delete [] dists;
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

double FineFittingTrackingLocalFrame::nearestNeighborEnergy(g_Part * mesh, bool useFeatures)
{
	int i;
	double signedPointPlaneDist;
	g_Vector currentPointCloudNormal, transformedPoint, targetLocation;
	double energy = 0;

	const g_NodeContainer& meshNodes = mesh->nodes();
	const g_NodeContainer& frameNodes = frameMesh->nodes();
	int meshNumNodes = mesh->nodes().numberOfItems();
	// Find the nearest neighbor for each transformed point and add the squared distance to the energy
	for(i = 0; i < meshNumNodes; i++)
	{
		if(nearestNeighbors[i] != -1)
		{
			if(useFeatures)
			{
				//in this case, use the point-to-point distance as this is for rough alignment, and as the points can be far away and far in normal angles:
				transformedPoint = computeTransformedPoint(&(transformationParameters[7*i]), meshNodes[i]->coordinate());

				targetLocation = frameNodes[nearestNeighbors[i]]->coordinate();
				approxNearestPointPlaneNeighbors[i].Set(targetLocation.x(), targetLocation.y(), targetLocation.z());

				energy += (transformedPoint-targetLocation).SquaredLength();
			}
			else
			{
				//compute the target point based on the point-plane distance:
				currentPointCloudNormal = frameMesh->getNormalAtNode(nearestNeighbors[i]+1);
				currentPointCloudNormal.Normalize();

				transformedPoint = computeTransformedPoint(&(transformationParameters[7*i]), meshNodes[i]->coordinate());

				signedPointPlaneDist = currentPointCloudNormal.Dot(frameNodes[nearestNeighbors[i]]->coordinate() - transformedPoint);
				
				targetLocation = transformedPoint + signedPointPlaneDist * currentPointCloudNormal;
				approxNearestPointPlaneNeighbors[i].Set(targetLocation.x(), targetLocation.y(), targetLocation.z());

				if (useProbeConstraint && constraintNodesHierarchy[0][i]) {
#ifdef _DEBUG_FINE_FITTING_TRACKING_LOCAL_FRAME_
				  double pEnergy = mesh->getNumberOfNodes() * pow(signedPointPlaneDist, 2);
				  cerr << "Probe constraint energy " << i << " :" << pEnergy << " sppd: " << signedPointPlaneDist << endl; 
				  energy += pEnergy;
#else
				  energy += mesh->getNumberOfNodes() * pow(signedPointPlaneDist, 2);
#endif
				} else {
				  energy += pow(signedPointPlaneDist, 2);
				}
				// add energy for hard constraint
			}
		}
	}

	return energy;
}

void FineFittingTrackingLocalFrame::nearestNeighborGradient(g_Part * mesh, double *& g)
{
	int i, numNodes;

	g_Vector differenceVec, partialGrad;

	numNodes = (int)mesh->nodes().numberOfItems();
	const g_NodeContainer& meshNodes = mesh->nodes();

	for(i = 0; i < 7*numNodes; i++) g[i] = 0;

	for(i = 0; i < numNodes; i++)
	{
		if(nearestNeighbors[i] != -1)
		{
			//Compute gradient of nearest neighbor energy with respect to the transformation matrix 
			differenceVec = 2.0 * (computeTransformedPoint(&(transformationParameters[7*i]), meshNodes[i]->coordinate()) -  approxNearestPointPlaneNeighbors[i]);
			if ( useProbeConstraint && constraintNodesHierarchy[0][i]) {
				differenceVec *= mesh->getNumberOfNodes();
			}
			//Compute the gradient w.r.t. the translation
			partialGrad = computeTransformedPoint(&(transformationParameters[7*i]), meshNodes[i]->coordinate(), 0);
			g[7*i] += differenceVec.Dot(partialGrad);

			partialGrad = computeTransformedPoint(&(transformationParameters[7*i]), meshNodes[i]->coordinate(), 1);
			g[7*i+1] += differenceVec.Dot(partialGrad);

			partialGrad = computeTransformedPoint(&(transformationParameters[7*i]), meshNodes[i]->coordinate(), 2);
			g[7*i+2] += differenceVec.Dot(partialGrad);
			
			//Compute the gradient w.r.t. the rotation
			partialGrad = computeTransformedPoint(&(transformationParameters[7*i]), meshNodes[i]->coordinate(), 3);
			g[7*i+3] += differenceVec.Dot(partialGrad);

			partialGrad = computeTransformedPoint(&(transformationParameters[7*i]), meshNodes[i]->coordinate(), 4);
			g[7*i+4] += differenceVec.Dot(partialGrad);

			partialGrad = computeTransformedPoint(&(transformationParameters[7*i]), meshNodes[i]->coordinate(), 5);
			g[7*i+5] += differenceVec.Dot(partialGrad);

			partialGrad = computeTransformedPoint(&(transformationParameters[7*i]), meshNodes[i]->coordinate(), 6);
			g[7*i+6] += differenceVec.Dot(partialGrad);
		}
	}
}

double FineFittingTrackingLocalFrame::distanceResolutionEnergy(g_Part * mesh, g_Part * reconstructedMesh, vector<int> indices)
{
	int i;

	double energy = 0.0;

	const g_NodeContainer& meshNodes = mesh->nodes();
	const g_NodeContainer& reconstructedNodes = reconstructedMesh->nodes();
	for(i = 0; i < indices.size(); i++)	
	  energy += (computeTransformedPoint(&(transformationParameters[7*indices[i]]), meshNodes[indices[i]]->coordinate()) - reconstructedNodes[indices[i]]->coordinate()).SquaredLength();
	
	return energy;
}

void FineFittingTrackingLocalFrame::distanceResolutionGradient(g_Part * mesh, g_Part * reconstructedMesh, vector<int> indices, double *& g)
{
	int i, numNodes;

	g_Vector differenceVec, partialGrad;

	numNodes = (int)mesh->nodes().numberOfItems();

	for(i = 0; i < 7*numNodes; i++) g[i] = 0;
	const g_NodeContainer& meshNodes = mesh->nodes();
	const g_NodeContainer& reconstructedNodes = reconstructedMesh->nodes();
	for(i = 0; i < indices.size(); i++)
	{
		//Compute gradient of nearest neighbor energy with respect to the transformation matrix 
		differenceVec = 2.0 * (computeTransformedPoint(&(transformationParameters[7*indices[i]]), meshNodes[indices[i]]->coordinate()) -  reconstructedNodes[indices[i]]->coordinate());

		//Compute the gradient w.r.t. the translation
		partialGrad = computeTransformedPoint(&(transformationParameters[7*indices[i]]), meshNodes[indices[i]]->coordinate(), 0);
		g[7*i] += differenceVec.Dot(partialGrad);

		partialGrad = computeTransformedPoint(&(transformationParameters[7*indices[i]]), meshNodes[indices[i]]->coordinate(), 1);
		g[7*i+1] += differenceVec.Dot(partialGrad);

		partialGrad = computeTransformedPoint(&(transformationParameters[7*indices[i]]), meshNodes[indices[i]]->coordinate(), 2);
		g[7*i+2] += differenceVec.Dot(partialGrad);
		
		//Compute the gradient w.r.t. the rotation
		partialGrad = computeTransformedPoint(&(transformationParameters[7*indices[i]]), meshNodes[indices[i]]->coordinate(), 3);
		g[7*i+3] += differenceVec.Dot(partialGrad);

		partialGrad = computeTransformedPoint(&(transformationParameters[7*indices[i]]), meshNodes[indices[i]]->coordinate(), 4);
		g[7*i+4] += differenceVec.Dot(partialGrad);

		partialGrad = computeTransformedPoint(&(transformationParameters[7*indices[i]]), meshNodes[indices[i]]->coordinate(), 5);
		g[7*i+5] += differenceVec.Dot(partialGrad);

		partialGrad = computeTransformedPoint(&(transformationParameters[7*indices[i]]), meshNodes[indices[i]]->coordinate(), 6);
		g[7*i+6] += differenceVec.Dot(partialGrad);
	}
}

void FineFittingTrackingLocalFrame::computeAccessNeighbors(g_Part * mesh, vector<int> * indices)
{
	int i, j, numNN, numNodes, numNeighboringElements;
	bool found;
	double radius, resolution;
	g_Element * currentIndicentTriangle;
	TriangleMesh tempMesh (*mesh);
	
	resolution = tempMesh.getMeshResolution();
	radius = SIZE_NBHD * resolution;

	numNodes = tempMesh.getNumberOfNodes();

	indicesForIntersection.clear();
	indicesForIntersection.resize(numNodes);

	indicesForSimpleRegularization.clear();
	indicesForSimpleRegularization.resize(numNodes);

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

	const g_NodeContainer& tempNodes = tempMesh.nodes();
	for(i = 0; i < tempMesh.getNumberOfNodes(); i++)
	{
		dataPtsForNN[i][0] = tempNodes[i]->coordinate().x();
		dataPtsForNN[i][1] = tempNodes[i]->coordinate().y();
		dataPtsForNN[i][2] = tempNodes[i]->coordinate().z();
	}

	kdTreeForNN = new ANNkd_tree(dataPtsForNN, tempMesh.getNumberOfNodes(), 3);

	//do spherical range search and store the list of neighbors:
	cerr << "Spherical range search: " << endl;
	bool tmpUseProbe = useProbeConstraint;
	useProbeConstraint = false;
	for(i = 0; i < tempMesh.getNumberOfNodes(); i++)
	{
		indicesForIntersection[i].clear();
		indicesForSimpleRegularization[i].clear();
		// assume hierarchy level 0 ???
		if ( useProbeConstraint && constraintNodesHierarchy[0][i] ) continue;
	  // cerr << i << " ";
	  cerr.flush();
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

		if(found) {
			queryPtForNN[0] = dataPtsForNN[i][0];
			queryPtForNN[1] = dataPtsForNN[i][1];
			queryPtForNN[2] = dataPtsForNN[i][2];
			numNN = kdTreeForNN->annkFRSearch(queryPtForNN, radius*radius, tempMesh.getNumberOfNodes(), nnIdxForNN, distsForNN);

			//Compute the indices for the smoothness energy (Euclidean distance > radius):
			for(j = 0; j < numNN; j++)
			{
				if(nnIdxForNN[j] != i)
					indicesForSimpleRegularization[i].push_back(nnIdxForNN[j]);
			}

			//Also add all 1-ring neighbors of a vertex to the indices for the simple regularization:
			std::set<int> indicesNeighbors;
			numNeighboringElements = tempNodes[i]->elements().numberOfItems();

			for(j = 0; j < numNeighboringElements; j++)
			{
				currentIndicentTriangle = tempNodes[i]->elements()[j];
				assert(currentIndicentTriangle->nodes().numberOfItems() == 3);
				for ( int nCnt=0; nCnt<3; ++nCnt) {
				  if(currentIndicentTriangle->nodes()[nCnt]->id() != i+1)
					indicesNeighbors.insert(currentIndicentTriangle->nodes()[nCnt]->id()-1);
				}
			}

			std::set<int>::iterator it;

			for(it = indicesNeighbors.begin(); it != indicesNeighbors.end(); it++)
			{
				//Recall that the ones closer to p_i than radius have already been added
				if(tempNodes[i]->coordinate().DistanceTo(tempNodes[*it]->coordinate()) > radius)
					indicesForSimpleRegularization[i].push_back(*it);
			}
		}
	}
	cerr << endl;
	useProbeConstraint = tmpUseProbe;

	delete [] nnIdxForNN;
	delete [] distsForNN;
	annDeallocPt(queryPtForNN);
	annDeallocPts(dataPtsForNN);
	delete kdTreeForNN;
}

double FineFittingTrackingLocalFrame::accessEnergy(g_Part * mesh)
{
	int i, j, k, id1, id2, numNodes, numAccessNeighbors, numSimpleNeighbors, numNeighbors;
	double energy, trafoDifference, weight, dotProdNormals, dotProdAxes, radius;
	g_Vector queryPt1, queryPt2, rotAxis1, rotAxis2, normal1, normal2;

	char trans = 'N';
	long int m = 4;
	double alpha = 1.0;
	double beta = 0.0;

	double * partX1 = new double[7];
	double * partX2 = new double[7];

	TriangleMesh tempMesh (*mesh);
	tempMesh.calculateVertexNormals();
	radius = SIZE_NBHD * tempMesh.getMeshResolution();
	numNodes = tempMesh.getNumberOfNodes();
	const g_NodeContainer& tempNodes = tempMesh.nodes();


	energy = 0;	

	//First, transform angles to the interval from 0 to 2*PI (CAREFUL: make copies):
	double * rotationAngles = new double[numNodes];
	for(i = 0; i < numNodes; i++)
	{
		rotationAngles[i] = transformationParameters[7*i+6];
		while(rotationAngles[i] < 0)
			rotationAngles[i] += 2.0*PI;
		while(rotationAngles[i] >= 2.0*PI)
			rotationAngles[i] -= 2.0*PI;
	}

	//Go through the list of neighbors:
	for(i = 0; i < numNodes; i++)
	{
		id1 = i;
		queryPt1 = tempNodes[id1]->coordinate();

		numAccessNeighbors = indicesForIntersection[i].size();
		numSimpleNeighbors = indicesForSimpleRegularization[i].size();
		numNeighbors =  numAccessNeighbors+ numSimpleNeighbors;
#ifdef _DEBUG_FINE_FITTING_TRACKING_LOCAL_FRAME_
		bool _print_debug_info_ = false;
		if (useProbeConstraint && constraintNodesHierarchy[0][i]) {
		  cerr << "Probe constraint access " << i << " :" << numAccessNeighbors << " + " << numSimpleNeighbors << endl; 
		  _print_debug_info_ = true;
		}
#endif
		for(j = 0; j < numNeighbors; j++)
		{
			//In the first case, weigh by the distance from the point, in the second one, simply give weight one
			if(j < numAccessNeighbors)
			{
				id2 = indicesForIntersection[i][j];
				queryPt2 = tempNodes[id2]->coordinate();
				weight = max(0.0, 1.0 - ((queryPt1-queryPt2).SquaredLength() / 
							 ((double)max(radius*radius, MD_DEFAULT_EPSILON))));
			}
			else 
			{
				id2 = indicesForSimpleRegularization[i][j-numAccessNeighbors];
				queryPt2 = tempNodes[id2]->coordinate();
				weight = 1.0;
			}
#ifdef _DEBUG_FINE_FITTING_TRACKING_LOCAL_FRAME_
			if ( _print_debug_info_ ) {
			  cerr << "Weight: " << weight << endl;
			}
#endif
			//Make sure the distance function is truncated at zero:
			if(weight > 0)
			{
				for(k = 0; k < 6; k++) partX1[k] = transformationParameters[7*id1+k];
				partX1[6] = rotationAngles[id1];
				for(k = 0; k < 6; k++) partX2[k] = transformationParameters[7*id2+k];
				partX2[6] = rotationAngles[id2];

				//Compute the difference in translation
				trafoDifference = 0.0;
				for(k = 0; k < 3; k++) trafoDifference += pow(partX1[k] - partX2[k], 2);

				/*/Compute the difference in rotation axis
				normal1 = tempMesh.getNormalAtNode(id1+1);
				normal2 = tempMesh.getNormalAtNode(id2+1);
				
				//Normalize the rotation axes
				rotAxis1.Set(partX1[3], partX1[4], partX1[5]);
				rotAxis1.Normalize();

				rotAxis2.Set(partX2[3], partX2[4], partX2[5]);
				rotAxis2.Normalize();

				//Add information related to the angle between the rotation axes compared to the angle between the outer normals
				dotProdNormals = normal1.Dot(normal2);
				dotProdAxes = rotAxis1.Dot(rotAxis2);
				trafoDifference += max(0.0, pow(dotProdAxes, 2) - pow(dotProdNormals, 2));
				*/
				//Compute the difference in rotation angle
				trafoDifference += min(pow(partX1[6] - partX2[6], 2), pow(2.0*PI - fabs(partX1[6] - partX2[6]), 2));

				//Normalize per point
				trafoDifference = trafoDifference / (double)numNeighbors;
#ifdef _DEBUG_FINE_FITTING_TRACKING_LOCAL_FRAME_
				if ( _print_debug_info_ ) {
				  cerr << "Energy added: " << weight * trafoDifference << endl;
				}
#endif
				energy += weight * trafoDifference;
			}
		}
	}

	delete [] partX1;
	delete [] partX2;

	delete [] rotationAngles;

	return energy;
}

void FineFittingTrackingLocalFrame::accessGradient(g_Part * mesh, double *& g)
{	
	int i, j, k, id1, id2, numNodes, numAccessNeighbors, numSimpleNeighbors, numNeighbors;
	double radius, weight, gradMultiplier, sqLength1, sqLength2, dotProd, diff;
	g_Vector queryPt1, queryPt2, rotAxis1, rotAxis2, normal1, normal2, partialDerivative;

	char trans = 'N';
	long int m = 4;
	double alpha = 1.0;
	double beta = 0.0;

	double * partX1 = new double[7];
	double * partX2 = new double[7];

	TriangleMesh tempMesh (*mesh);
	tempMesh.calculateVertexNormals();
	radius = SIZE_NBHD * tempMesh.getMeshResolution();
	numNodes = (int)tempMesh.getNumberOfNodes();
	const g_NodeContainer& tempNodes = tempMesh.nodes();

	for(i = 0; i < 7*numNodes; i++) g[i] = 0;

	//First, transform angles to the interval from 0 to 2*PI (CAREFUL: make copies):
	double * rotationAngles = new double[numNodes];
	for(i = 0; i < numNodes; i++)
	{
		rotationAngles[i] = transformationParameters[7*i+6];
		while(rotationAngles[i] < 0)
			rotationAngles[i] += 2.0*PI;
		while(rotationAngles[i] >= 2.0*PI)
			rotationAngles[i] -= 2.0*PI;
	}

	for(i = 0; i < numNodes; i++)
	{
		id1 = i;
		queryPt1 = tempNodes[id1]->coordinate();

		numAccessNeighbors = indicesForIntersection[i].size();
		numSimpleNeighbors = indicesForSimpleRegularization[i].size();
		numNeighbors = numAccessNeighbors + numSimpleNeighbors;

		for(j = 0; j < numNeighbors; j++)
		{
			//In the first case, weigh by the distance from the point, in the second one, simply give weight one
			if(j < numAccessNeighbors)
			{
				id2 = indicesForIntersection[i][j];
				queryPt2 = tempNodes[id2]->coordinate();
				weight = max(0.0, 1.0 - ((queryPt1-queryPt2).SquaredLength() / 
							 ((double)max(radius*radius, MD_DEFAULT_EPSILON))));
			}
			else
			{
				id2 = indicesForSimpleRegularization[i][j-numAccessNeighbors];
				queryPt2 = tempNodes[id2]->coordinate();
				weight = 1.0;
			}

			//Make sure the distance function is truncated at zero:
			if(weight > 0)
			{	
				gradMultiplier = 1.0 / (double)numNeighbors * weight;

				for(k = 0; k < 6; k++) partX1[k] = transformationParameters[7*id1+k];
				partX1[6] = rotationAngles[id1];
				for(k = 0; k < 6; k++) partX2[k] = transformationParameters[7*id2+k];
				partX2[6] = rotationAngles[id2];

				//Derviative w.r.t. translation:
				g[7*id1]   += gradMultiplier * 2.0 * (partX1[0] - partX2[0]);
				g[7*id1+1] += gradMultiplier * 2.0 * (partX1[1] - partX2[1]);
				g[7*id1+2] += gradMultiplier * 2.0 * (partX1[2] - partX2[2]);

				g[7*id2]   -= gradMultiplier * 2.0 * (partX1[0] - partX2[0]);
				g[7*id2+1] -= gradMultiplier * 2.0 * (partX1[1] - partX2[1]);
				g[7*id2+2] -= gradMultiplier * 2.0 * (partX1[2] - partX2[2]);

				//Derivative w.r.t. rotation parameters:
				/*/Compute derivatives w.r.t. rotation axes angle's deviation from the normal angle
				normal1 = tempMesh.getNormalAtNode(id1+1);
				normal2 = tempMesh.getNormalAtNode(id2+1);

				rotAxis1.Set(partX1[3], partX1[4], partX1[5]);
				rotAxis2.Set(partX2[3], partX2[4], partX2[5]);

				if(pow(rotAxis1.Dot(rotAxis2)/(rotAxis1.Length()*rotAxis2.Length()), 2) - pow(normal1.Dot(normal2), 2) > 0)
				{
					sqLength1 = rotAxis1.SquaredLength();
					sqLength2 = rotAxis2.SquaredLength();
					dotProd = rotAxis1.Dot(rotAxis2);

					partialDerivative = (2.0 / (pow(sqLength1 * sqLength2, 2))) * 
						((dotProd * sqLength1 * sqLength2) * rotAxis2 - (pow(dotProd, 2) * sqLength2) * rotAxis1);
					g[7*id1+3] += gradMultiplier * partialDerivative.x();
					g[7*id1+4] += gradMultiplier * partialDerivative.y();
					g[7*id1+5] += gradMultiplier * partialDerivative.z();

					partialDerivative = (2.0 / (pow(sqLength1 * sqLength2, 2))) * 
						((dotProd * sqLength1 * sqLength2) * rotAxis1 - (pow(dotProd, 2) * sqLength1) * rotAxis2);
					g[7*id2+3] += gradMultiplier * partialDerivative.x();
					g[7*id2+4] += gradMultiplier * partialDerivative.y();
					g[7*id2+5] += gradMultiplier * partialDerivative.z();
				}
				*/
				//Derivative w.r.t. angle difference
				if(partX1[6] >= partX2[6])
				{
					diff = partX1[6] - partX2[6];
					if(diff <= (2.0*PI - diff))
					{
						g[7*id1+6] += gradMultiplier * 2.0 * diff;
						g[7*id2+6] -= gradMultiplier * 2.0 * diff;
					}
					else
					{
						g[7*id1+6] -= gradMultiplier * 2.0 * (2.0*PI - diff);
						g[7*id2+6] += gradMultiplier * 2.0 * (2.0*PI - diff);
					}
				}
				else
				{
					diff = partX2[6] - partX1[6];
					if(diff <= (2.0*PI - diff))
					{
						g[7*id1+6] -= gradMultiplier * 2.0 * diff;
						g[7*id2+6] += gradMultiplier * 2.0 * diff;
					}
					else
					{
						g[7*id1+6] += gradMultiplier * 2.0 * (2.0*PI - diff);
						g[7*id2+6] -= gradMultiplier * 2.0 * (2.0*PI - diff);
					}
				}
			}
		}
	}
	
	delete [] partX1;
	delete [] partX2;

	delete [] rotationAngles;
}

void FineFittingTrackingLocalFrame::computeMdNeighbors(g_Part * mesh, vector<int> * indices)
{
	int i, j, numNN;
	bool found, isNeighbor;
	g_Vector transformedPoint;
	TriangleMesh tempMesh (*mesh);
	
	double radius = tempMesh.getMeshResolution();

	tempMesh.initGeodesicDistanceCalculation();

	indicesForMd.clear();
	indicesForMd.resize(tempMesh.getNumberOfNodes());

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

	const g_NodeContainer& tempNodes = tempMesh.nodes();
	for(i = 0; i < tempMesh.getNumberOfNodes(); i++)
	{
		transformedPoint = computeTransformedPoint(&(transformationParameters[7*i]), tempNodes[i]->coordinate());
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
			indicesForMd[i].clear();

			tempMesh.getGeodesicDistance(i+1, radius, geodesicDists);

			for(j = 0; j < numNN; j++)
			{
				//intersect the Euclidean sphere with a geodesic circle
				isNeighbor = false;
				for(int k = 0; k < (int)indicesForIntersection[i].size(); k++)
				{
					if((geodesicDists[nnIdxForNN[j]] != -1) && (geodesicDists[nnIdxForNN[j]] < radius))
						isNeighbor = true;
				}
				if(nnIdxForNN[j] == i) isNeighbor = true;

				if(!isNeighbor)
					indicesForMd[i].push_back(nnIdxForNN[j]);
			}
		}
		else indicesForMd[i].clear();
	}

	delete [] nnIdxForNN;
	delete [] distsForNN;
	annDeallocPt(queryPtForNN);
	annDeallocPts(dataPtsForNN);
	delete kdTreeForNN;
}

double FineFittingTrackingLocalFrame::mdEnergy(g_Part * mesh)
{
	int i, j, numPickedNeighbors;
	int id;
	double length, squaredLength, energy, perPtEnergy, weight, threshold, thresholdMin;
	g_Vector closePt1, closePt2, diffVec;
	TriangleMesh tempMesh(*mesh);

	energy = 0;

	//Consider neighbors within a disk of radius threshold:
	threshold = tempMesh.getMeshResolution() / MD_THRESHOLD_DIVIDE;
	//Let the maximum scale of one energy term be around the same size as one energy term in nearestNeighborEnergy:
	thresholdMin = 0.0; // 1.0 / ((double) MAX_DIST * tempMesh.getMeshResolution());

	numUsedNeighbors.resize((int)tempMesh.getNumberOfNodes());
	const g_NodeContainer& tempNodes = tempMesh.nodes();

	//compute the gradient using TriangleMesh indicesForIntersection
	for(i = 0; i < (int)tempMesh.getNumberOfNodes(); i++)
	{
		numPickedNeighbors = 0;
		perPtEnergy = 0;

		for(j = 0; j < (int)indicesForMd[i].size(); j++)
		{
			id = indicesForMd[i][j];
			closePt1 = computeTransformedPoint(&(transformationParameters[7*i]), tempNodes[i]->coordinate());
			closePt2 = computeTransformedPoint(&(transformationParameters[7*id]), tempNodes[id]->coordinate());
		
			diffVec = closePt1 - closePt2;

			length = diffVec.Length();

			if(length < thresholdMin) length = thresholdMin;

			if(length < threshold)
			{
				numPickedNeighbors++;

				squaredLength = pow(length, 2);
				weight = 1.0 - squaredLength / (double)max(threshold*threshold, MD_DEFAULT_EPSILON);

				perPtEnergy += weight * (1.0 / (double)max(squaredLength, MD_DEFAULT_EPSILON));
			}
		}

		if(numPickedNeighbors > 0) energy += perPtEnergy / (double)numPickedNeighbors;

		numUsedNeighbors[i] = numPickedNeighbors;
	}

	return energy;
}

void FineFittingTrackingLocalFrame::mdGradient(g_Part * mesh, double *& g)
{
	int i, j, k;
	int id;
	double length, squaredLength, weight, invDist, threshold, thresholdMin, gradWeight, gradInvDist, grad;
	g_Vector closePt1, closePt2, diffVec, partialGrad;
	TriangleMesh tempMesh(*mesh);

	//Truncate at threshold:
	threshold = tempMesh.getMeshResolution() / MD_THRESHOLD_DIVIDE;
	//Let the maximum scale of one energy term be around the same size as one energy term in nearestNeighborEnergy:
	thresholdMin = 0.0; // 1.0 / ((double) MAX_DIST * tempMesh.getMeshResolution());

	for(i = 0; i < 7*(int)tempMesh.getNumberOfNodes(); i++) g[i] = 0.0;

	const g_NodeContainer& tempNodes = tempMesh.nodes();
	//compute the gradient using TriangleMesh indicesForIntersection
	for(i = 0; i < (int)tempMesh.getNumberOfNodes(); i++)
	{
		for(j = 0; j < (int)indicesForMd[i].size(); j++)
		{
			id = indicesForMd[i][j];
			closePt1 = computeTransformedPoint(&(transformationParameters[7*i]), tempNodes[i]->coordinate());
			closePt2 = computeTransformedPoint(&(transformationParameters[7*id]), tempNodes[id]->coordinate());
			
			diffVec = closePt1 - closePt2;

			length = diffVec.Length();

			if(length < thresholdMin)
			{
				length = thresholdMin;
				diffVec.Normalize();
				diffVec = length * diffVec;
			}

			if(length < threshold)
			{
				squaredLength = pow(length, 2);
				weight = 1.0 - squaredLength / (double)max(threshold*threshold, MD_DEFAULT_EPSILON);
				invDist = 1.0 / (double)max(squaredLength, MD_DEFAULT_EPSILON);

				//Gradient affects two points
				for(k = 0; k < 7; k++)
				{
					//Derivative w.r.t. transformation parameters of first point
					partialGrad = computeTransformedPoint(&(transformationParameters[7*i]), tempNodes[i]->coordinate(), k);
					
					gradWeight = -1.0 / (double)max(threshold*threshold, MD_DEFAULT_EPSILON) * 2.0 * diffVec.Dot(partialGrad);
					gradInvDist = -1.0 / (double)max(squaredLength*squaredLength, MD_DEFAULT_EPSILON) * 2.0 * diffVec.Dot(partialGrad);

					grad = gradWeight * invDist + weight * gradInvDist;

					g[7*i + k] += grad / (double)numUsedNeighbors[i];
					
					//Derivative w.r.t. transformation parameters of second point
					partialGrad = computeTransformedPoint(&(transformationParameters[7*id]), tempNodes[id]->coordinate(), k);
					
					gradWeight = -1.0 / (double)max(threshold*threshold, MD_DEFAULT_EPSILON) * 2.0 * diffVec.Dot(partialGrad);
					gradInvDist = -1.0 / (double)max(squaredLength*squaredLength, MD_DEFAULT_EPSILON) * 2.0 * diffVec.Dot(partialGrad);

					grad = gradWeight * invDist + weight * gradInvDist;

					g[7*id + k] -= grad / (double)numUsedNeighbors[i];
				}
			}
		}
	}
}

void FineFittingTrackingLocalFrame::computeMdNeighborsFixedBarycenters(g_Part * mesh, vector<int> * indices)
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

	/*/Search within the maximum edge length (guarantee that no intersections are missed):
	radius = 0;
	g_PEdgeContainer allEdges = tempMesh.uniqueEdges();
	for(i = 0; i < allEdges.numberOfItems(); i++)
	{
		radius = max(radius, allEdges[i]->firstNode().coordinate().DistanceTo(allEdges[i]->lastNode().coordinate()));
		delete allEdges[i];
	}
	//Much faster:*/
	radius = templateMesh->getMeshResolution() / 2.0;
	
	threshold = templateMesh->getMeshResolution() / MD_THRESHOLD_DIVIDE; // 2.0;

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
	const g_NodeContainer& tempNodes = tempMesh.nodes();
	for(i = 0; i < tempMesh.getNumberOfNodes(); i++)
	{
		transformedPoint = computeTransformedPoint(&(transformationParameters[7*i]), tempNodes[i]->coordinate());
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
				if(nnIdxForNN[j] != i) neighborTriangles.append(tempNodes[nnIdxForNN[j]]->elements());
			}
			adapter.removeDuplicates();

			//Store the triangle candidates (this avoids duplicates):
			for(j = 0; j < tempNodes[i]->elements().numberOfItems(); j++)
			{
				for(k = 0; k < neighborTriangles.numberOfItems(); k++)
				{
					vector<int> currentCandidate;
					if(tempNodes[i]->elements()[j]->id()-1 < neighborTriangles[k]->id()-1)
					{
						currentCandidate.push_back(tempNodes[i]->elements()[j]->id()-1);
						currentCandidate.push_back(neighborTriangles[k]->id()-1);
					}
					else
					{
						currentCandidate.push_back(neighborTriangles[k]->id()-1);
						currentCandidate.push_back(tempNodes[i]->elements()[j]->id()-1);				
					}
					triangleCandidates.insert(currentCandidate);
				}
			}
		}
	}

	//Go through all the triangle candidates:
	set< vector<int> >::iterator itOuter;
	const g_ElementContainer& tempElements = tempMesh.elements();
	for(itOuter = triangleCandidates.begin(); itOuter != triangleCandidates.end(); itOuter++)
	{
		id = tempElements[(*itOuter)[0]]->nodes()[0]->id()-1;
		t1.nodes()[0]->coordinate(g_Vector(dataPtsForNN[id][0], dataPtsForNN[id][1], dataPtsForNN[id][2]));
		t1.nodes()[0]->id(id);
		id = tempElements[(*itOuter)[0]]->nodes()[1]->id()-1;
		t1.nodes()[1]->coordinate(g_Vector(dataPtsForNN[id][0], dataPtsForNN[id][1], dataPtsForNN[id][2]));
		t1.nodes()[1]->id(id);
		id = tempElements[(*itOuter)[0]]->nodes()[2]->id()-1;
		t1.nodes()[2]->coordinate(g_Vector(dataPtsForNN[id][0], dataPtsForNN[id][1], dataPtsForNN[id][2]));
		t1.nodes()[2]->id(id);

		id = tempElements[(*itOuter)[1]]->nodes()[0]->id()-1;
		t2.nodes()[0]->coordinate(g_Vector(dataPtsForNN[id][0], dataPtsForNN[id][1], dataPtsForNN[id][2]));
		t2.nodes()[0]->id(id);
		id = tempElements[(*itOuter)[1]]->nodes()[1]->id()-1;
		t2.nodes()[1]->coordinate(g_Vector(dataPtsForNN[id][0], dataPtsForNN[id][1], dataPtsForNN[id][2]));
		t2.nodes()[1]->id(id);
		id = tempElements[(*itOuter)[1]]->nodes()[2]->id()-1;
		t2.nodes()[2]->coordinate(g_Vector(dataPtsForNN[id][0], dataPtsForNN[id][1], dataPtsForNN[id][2]));
		t2.nodes()[2]->id(id);

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

double FineFittingTrackingLocalFrame::mdEnergyFixedBarycenters(g_Part * mesh)
{
	int i;
	int id[6];
	double length, squaredLength, energy, weight, threshold;
	g_Vector p1, p2, p3, p4, p5, p6, closePt1, closePt2, diffVec;
	TriangleMesh tempMesh(*mesh);

	energy = 0;

	//Truncate at threshold:
	threshold = templateMesh->getMeshResolution() / MD_THRESHOLD_DIVIDE;

	const g_ElementContainer& tempElements = tempMesh.elements();
	const g_NodeContainer& tempNodes = tempMesh.nodes();
	//compute the gradient using TriangleMesh indicesForIntersection
	for(i = 0; i < (int)(triangleIdsForMd.size()/2.0); i++)
	{
		id[0] = tempElements[triangleIdsForMd[2*i]]->nodes()[0]->id()-1;
		p1 = computeTransformedPoint(&(transformationParameters[7*id[0]]), tempNodes[id[0]]->coordinate());
		id[1] = tempElements[triangleIdsForMd[2*i]]->nodes()[1]->id()-1;
		p2 = computeTransformedPoint(&(transformationParameters[7*id[1]]), tempNodes[id[1]]->coordinate());
		id[2] = tempElements[triangleIdsForMd[2*i]]->nodes()[2]->id()-1;
		p3 = computeTransformedPoint(&(transformationParameters[7*id[2]]), tempNodes[id[2]]->coordinate());
		closePt1 = barycentersForMd[6*i] * p1 + barycentersForMd[6*i+1] * p2 + barycentersForMd[6*i+2] * p3;

		id[3] = tempElements[triangleIdsForMd[2*i+1]]->nodes()[0]->id()-1;
		p4 = computeTransformedPoint(&(transformationParameters[7*id[3]]), tempNodes[id[3]]->coordinate());
		id[4] = tempElements[triangleIdsForMd[2*i+1]]->nodes()[1]->id()-1;
		p5 = computeTransformedPoint(&(transformationParameters[7*id[4]]), tempNodes[id[4]]->coordinate());
		id[5] = tempElements[triangleIdsForMd[2*i+1]]->nodes()[2]->id()-1;
		p6 = computeTransformedPoint(&(transformationParameters[7*id[5]]), tempNodes[id[5]]->coordinate());
		closePt2 = barycentersForMd[6*i+3] * p4 + barycentersForMd[6*i+4] * p5 + barycentersForMd[6*i+5] * p6;

		//Use the product rule: derivative of weight * 1/distance + derivative of 1/distance * weight
		diffVec = closePt1 - closePt2;

		length = diffVec.Length();

		if(length < threshold)
		{
			squaredLength = pow(length, 2);
			weight = 1.0 - squaredLength / (double)max(threshold*threshold, MD_DEFAULT_EPSILON);
			energy += weight * (1.0 / (double)max(squaredLength, MD_DEFAULT_EPSILON)); 
		}
	}

	return energy;
}

void FineFittingTrackingLocalFrame::mdGradientFixedBarycenters(g_Part * mesh, double *& g)
{
	int i, j, k;
	int id[6];
	double length, squaredLength, weight, invDist, threshold, gradWeight, gradInvDist, grad;
	g_Vector p1, p2, p3, p4, p5, p6, closePt1, closePt2, diffVec, partialGrad;
	TriangleMesh tempMesh(*mesh);

	//Truncate at threshold:
	threshold = templateMesh->getMeshResolution() / MD_THRESHOLD_DIVIDE;

	for(i = 0; i < 7*tempMesh.getNumberOfNodes(); i++) g[i] = 0.0;

	const g_ElementContainer& tempElements = tempMesh.elements();
	const g_NodeContainer& tempNodes = tempMesh.nodes();
	//compute the gradient using TriangleMesh indicesForIntersection
	for(i = 0; i < (int)(triangleIdsForMd.size()/2.0); i++)
	{
		id[0] = tempElements[triangleIdsForMd[2*i]]->nodes()[0]->id()-1;
		p1 = computeTransformedPoint(&(transformationParameters[7*id[0]]), tempNodes[id[0]]->coordinate());
		id[1] = tempElements[triangleIdsForMd[2*i]]->nodes()[1]->id()-1;
		p2 = computeTransformedPoint(&(transformationParameters[7*id[1]]), tempNodes[id[1]]->coordinate());
		id[2] = tempElements[triangleIdsForMd[2*i]]->nodes()[2]->id()-1;
		p3 = computeTransformedPoint(&(transformationParameters[7*id[2]]), tempNodes[id[2]]->coordinate());
		closePt1 = barycentersForMd[6*i] * p1 + barycentersForMd[6*i+1] * p2 + barycentersForMd[6*i+2] * p3;

		id[3] = tempElements[triangleIdsForMd[2*i+1]]->nodes()[0]->id()-1;
		p4 = computeTransformedPoint(&(transformationParameters[7*id[3]]), tempNodes[id[3]]->coordinate());
		id[4] = tempElements[triangleIdsForMd[2*i+1]]->nodes()[1]->id()-1;
		p5 = computeTransformedPoint(&(transformationParameters[7*id[4]]), tempNodes[id[4]]->coordinate());
		id[5] = tempElements[triangleIdsForMd[2*i+1]]->nodes()[2]->id()-1;
		p6 = computeTransformedPoint(&(transformationParameters[7*id[5]]), tempNodes[id[5]]->coordinate());
		closePt2 = barycentersForMd[6*i+3] * p4 + barycentersForMd[6*i+4] * p5 + barycentersForMd[6*i+5] * p6;

		diffVec = closePt1 - closePt2;

		length = diffVec.Length();

		if(length < threshold)
		{
			squaredLength = pow(length, 2);
			weight = 1.0 - squaredLength / (double)max(threshold*threshold, MD_DEFAULT_EPSILON);
			invDist = 1.0 / (double)max(squaredLength, MD_DEFAULT_EPSILON);

			//Gradient affects up to 6 points
			for(j = 0; j < 6; j++)
			{
				if(barycentersForMd[6*i+j] > 0)
				{	
					for(k = 0; k < 7; k++)
					{
						//Derivative w.r.t. transformation parameters
						partialGrad = computeTransformedPoint(&(transformationParameters[7*id[j]]), tempNodes[id[j]]->coordinate(), k);
						
						gradWeight = -1.0 / (double)max(threshold*threshold, MD_DEFAULT_EPSILON) * 2.0 * diffVec.Dot(partialGrad);
						gradInvDist = -1.0 / (double)max(squaredLength*squaredLength, MD_DEFAULT_EPSILON) * 2.0 * diffVec.Dot(partialGrad);

						grad = barycentersForMd[6*i+j] * (gradWeight * invDist + weight * gradInvDist);

						if(j < 3) g[7*id[j] + k] += grad;
						else g[7*id[j] + k] -= grad;
					}
				}
			}
		}
	}
}

g_Vector FineFittingTrackingLocalFrame::computeTransformedPoint(double * X, g_Vector point, int index)
{
	int i;

	double * transformation1 = new double[16];
	double * transformation2 = new double[16];
	double * transformation3 = new double[16];
	double * transformation4 = new double[16];
	double * partialResult1 = new double[16];
	double * partialResult2 = new double[16];

	// 1. Translate the point to the origin
	computeTranslationMat(-point.x(), -point.y(), -point.z(), transformation1);

	// 2. Translate the point by X[0], X[1], X[2]
	if((index  < 0) || (index > 2)) computeTranslationMat(X[0], X[1], X[2], transformation2);
	else computeTranslationGrad(X[0], X[1], X[2], transformation2, index);
	
	// 3. Rotate the point using quaternions
	if((index < 3) || (index > 6)) computeRotationMat(X[3], X[4], X[5], X[6], transformation3);
	else computeRotationGrad(X[3], X[4], X[5], X[6], transformation3, index-2);

	// 4. Translate the point back
	computeTranslationMat(point.x(), point.y(), point.z(), transformation4);

	// 5. Compute composite transformation
	char trans = 'N';
	char transY = 'T';
	long int m = 4;
	double alpha = 1.0;
	double beta = 0.0;

	clapack::dgemm_(&trans, &trans, &m, &m, &m, &alpha, transformation4, &m, transformation3, &m, &beta, partialResult1, &m);
	clapack::dgemm_(&trans, &trans, &m, &m, &m, &alpha, transformation2, &m, transformation1, &m, &beta, partialResult2, &m);
	clapack::dgemm_(&trans, &trans, &m, &m, &m, &alpha, partialResult1, &m, partialResult2, &m, &beta, transformation1, &m);	

	// 6. Transform the point
	g_Vector resultPoint;
	double resultX, resultY, resultZ;

	resultX = transformation1[0]*point.x() + transformation1[4]*point.y() + transformation1[8]*point.z() + transformation1[12];
	resultY = transformation1[1]*point.x() + transformation1[5]*point.y() + transformation1[9]*point.z() + transformation1[13];
	resultZ = transformation1[2]*point.x() + transformation1[6]*point.y() + transformation1[10]*point.z() + transformation1[14];

	resultPoint.Set(resultX, resultY, resultZ);

	delete [] transformation1;
	delete [] transformation2;
	delete [] transformation3;
	delete [] transformation4;
	delete [] partialResult1;
	delete [] partialResult2;

	return resultPoint;
}
