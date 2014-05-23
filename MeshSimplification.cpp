#include "MeshSimplification.h"

#ifdef USE_MOELLER
#include "tritri.h"
#endif




bool operator < (const priorityEdge c1, const priorityEdge c2)
{
	return c1.weight < c2.weight;
}
bool operator > (const priorityEdge c1, const priorityEdge c2)
{
	return c1.weight > c2.weight;
}

MeshSimplification::MeshSimplification()
{
	originalMesh = NULL;
	workingMesh = NULL;
	finalMesh = NULL;
	geomCost = true;
}

MeshSimplification::MeshSimplification(g_Part * mesh)
{
	originalMesh = mesh;
	resolution = (TriangleMesh(*originalMesh)).getMeshResolution();
	workingMesh = NULL;
	finalMesh = NULL;
	geomCost = true;

	initOriginal();
}

MeshSimplification::~MeshSimplification()
{
	if(workingMesh != NULL)
	{
		delete workingMesh;
		workingMesh = NULL;
	}
	if(finalMesh != NULL)
	{
		delete finalMesh;
		finalMesh = NULL;
	}
}

void MeshSimplification::initNewMesh(g_Part * mesh)
{
	if(workingMesh != NULL)
	{
		delete workingMesh;
		workingMesh = NULL;
	}
	if(finalMesh != NULL)
	{
		delete finalMesh;
		finalMesh = NULL;
	}
	originalMesh = mesh;
	resolution = (TriangleMesh(*originalMesh)).getMeshResolution();
	initOriginal();
}

g_Part * MeshSimplification::getSimplifiedMesh(int numNodes, bool checkIntersection)
{
	int numberCurrentNodes = (int)originalMesh->nodes().numberOfItems();
	cerr << numberCurrentNodes << " simplify to " << numNodes << " requested" << endl;

	int i, numEdges, id1, previouslyDeleted, previouslyCollapsed;
	double edgeWeight;
	vector<int> newNeighbors;
	priorityEdge prioEdge;
	g_PEdge * edge;
	g_Node firstNode, lastNode;
	g_PEdgeContainer allEdges;

	//build the tree structure:
	treeStructure.resize((int)originalMesh->nodes().numberOfItems());
	for(i = 0; i < (int)originalMesh->nodes().numberOfItems(); i++)
	{
		treeNode treeNd;
		treeNd.id = i+1;
		treeNd.parent = i+1;
		treeNd.child2.clear();
		treeStructure[i] = treeNd;
	}

	//collapse edges one by one
	previouslyDeleted = -1;
	priority_queue<priorityEdge, vector<priorityEdge>, greater<priorityEdge> > edgesToCollapse;
	while(numberCurrentNodes > numNodes)
	{
		//find edgesToCollapse and store them in a priority queue:
		if(previouslyDeleted == -1)
		{
			allEdges = workingMesh->uniqueEdges();
			numEdges = (int)allEdges.numberOfItems();
			cerr << "Adding " << numEdges << " to priority queue" << endl;
			for(i = 0; i < numEdges; i++)
			{
				firstNode = allEdges[i]->firstNode();
				lastNode = allEdges[i]->lastNode();
				edge = new g_PEdge(lastNode, firstNode);
				allEdges.insert(edge);
			}
			// each half-edge is now inserted
		}
		else
		{
#ifdef DEBUG_MESH_SIMPLIFICATION
		        cerr << "allEdges clear" << endl;
#endif
			allEdges.clear();
			for(i = 0; i < (int)newNeighbors.size(); i++)
			{
#ifdef DEBUG_MESH_SIMPLIFICATION
			  cerr << "Looping over neighbors: " << i << endl;
#endif
				edge = new g_PEdge(*(workingMesh->nodes()[previouslyCollapsed]), *(workingMesh->nodes()[newNeighbors[i]]));
				allEdges.insert(edge);
				edge = new g_PEdge(*(workingMesh->nodes()[newNeighbors[i]]),*(workingMesh->nodes()[previouslyCollapsed]));
				allEdges.insert(edge);
			}
		}
#ifdef DEBUG_MESH_SIMPLIFICATION
		cerr << "# edges: " << allEdges.size() << endl;
#endif
		numEdges = (int)allEdges.numberOfItems();
		for(i = 0; i < numEdges; i++)
		{
			edge = allEdges[i];
			id1 = edge->firstNode().id()-1;
			if(geomCost) edgeWeight = computeHalfEdgeCostGeom(edge, checkIntersection);
			else edgeWeight = computeHalfEdgeCostVolume(edge, checkIntersection);
#ifdef DEBUG_MESH_SIMPLIFICATION
			cerr << "Edge " << i << " : " << edgeWeight << " (" << id1 << ")" <<endl;
#endif
			if(edgeWeight == -1) delete edge;
			else
			{
				prioEdge.firstNode = workingMesh->nodes()[edge->firstNode().id()-1];
				prioEdge.lastNode = workingMesh->nodes()[edge->lastNode().id()-1];
				prioEdge.weight = edgeWeight;
				edgesToCollapse.push(prioEdge);
				delete edge;
			}
		}
		//collapse the edges:
		if(edgesToCollapse.empty()) break;
		else
		{
			edge = findBestHalfEdge(edgesToCollapse, checkIntersection);
			if(edge == NULL) break;
			else
			{
#ifdef DEBUG_MESH_SIMPLIFICATION
			  cerr << "Half edge collapse: " <<
			    edge->firstNode().id()-1 << " " << edge->lastNode().id()-1 << "...";
#endif
				halfEdgeCollapse(edge->firstNode().id()-1, edge->lastNode().id()-1, newNeighbors);
				previouslyDeleted = edge->firstNode().id()-1;
				previouslyCollapsed = edge->lastNode().id()-1;
				numberCurrentNodes--;
#ifdef DEBUG_MESH_SIMPLIFICATION
				cerr << "done." << endl;
#endif
				delete edge;
				edge = NULL;
			}
		}
	}
	cerr << "Final mesh creation: " << endl;

	vector<int> indicesToKeep;
	createFinalMesh(indicesToKeep, false, true);
	for ( vector<int>::iterator it = indicesToKeep.begin();
	      it != indicesToKeep.end(); ++it ) {
	  cerr << *it << " ";
	}
	cerr << endl;
	return finalMesh;
}

g_Part * MeshSimplification::getSimplifiedMesh(int numNodes, vector<bool> importantToCollapse, bool checkIntersection)
{
	int numberCurrentNodes = (int)originalMesh->nodes().numberOfItems();

	int i, numEdges, id1, previouslyDeleted, previouslyCollapsed;
	double edgeWeight, maxPriority;
	vector<int> newNeighbors;
	priorityEdge prioEdge;
	g_PEdge * edge;
	g_Node firstNode, lastNode;
	g_PEdgeContainer allEdges;

	//build the tree structure:
	treeStructure.resize((int)originalMesh->nodes().numberOfItems());
	for(i = 0; i < (int)originalMesh->nodes().numberOfItems(); i++)
	{
		treeNode treeNd;
		treeNd.id = i+1;
		treeNd.parent = i+1;
		treeNd.child2.clear();
		treeStructure[i] = treeNd;
	}

	//compute the maximum priority
	maxPriority = 0;
	allEdges = workingMesh->uniqueEdges();
	numEdges = (int)allEdges.numberOfItems();
	for(i = 0; i < numEdges; i++)
	{
		firstNode = allEdges[i]->firstNode();
		lastNode = allEdges[i]->lastNode();
		edge = new g_PEdge(lastNode, firstNode);
		allEdges.insert(edge);
	}
	numEdges = (int)allEdges.numberOfItems();
	for(i = 0; i < numEdges; i++)
	{
		edge = allEdges[i];
		if(geomCost) edgeWeight = computeHalfEdgeCostGeom(edge, checkIntersection);
		else edgeWeight = computeHalfEdgeCostVolume(edge, checkIntersection);
		if(edgeWeight > maxPriority) maxPriority = edgeWeight;
	}

	//collapse edges one by one
	previouslyDeleted = -1;
	priority_queue<priorityEdge, vector<priorityEdge>, greater<priorityEdge> > edgesToCollapse;
	while(numberCurrentNodes > numNodes)
	{
		//find edgesToCollapse and store them in a priority queue:
		if(previouslyDeleted != -1)
		{
			allEdges.clear();
			for(i = 0; i < (int)newNeighbors.size(); i++)
			{
				edge = new g_PEdge(*(workingMesh->nodes()[previouslyCollapsed]), *(workingMesh->nodes()[newNeighbors[i]]));
				allEdges.insert(edge);
				edge = new g_PEdge(*(workingMesh->nodes()[newNeighbors[i]]),*(workingMesh->nodes()[previouslyCollapsed]));
				allEdges.insert(edge);
			}
		}
		numEdges = (int)allEdges.numberOfItems();
		for(i = 0; i < numEdges; i++)
		{
			edge = allEdges[i];
			id1 = edge->firstNode().id()-1;
			if(geomCost) edgeWeight = computeHalfEdgeCostGeom(edge, checkIntersection);
			else edgeWeight = computeHalfEdgeCostVolume(edge, checkIntersection);
			if(edgeWeight == -1) delete edge;
			else
			{
				if(importantToCollapse[id1]) edgeWeight += 3*maxPriority;
				prioEdge.firstNode = workingMesh->nodes()[edge->firstNode().id()-1];
				prioEdge.lastNode = workingMesh->nodes()[edge->lastNode().id()-1];
				prioEdge.weight = edgeWeight;
				edgesToCollapse.push(prioEdge);
				delete edge;
			}
		}
		//collapse the edges:
		if(edgesToCollapse.empty()) break;
		else
		{
			edge = findBestHalfEdge(edgesToCollapse, checkIntersection);
			if(edge == NULL) break;
			else
			{
				halfEdgeCollapse(edge->firstNode().id()-1, edge->lastNode().id()-1, newNeighbors);
				previouslyDeleted = edge->firstNode().id()-1;
				previouslyCollapsed = edge->lastNode().id()-1;
				numberCurrentNodes--;
				delete edge;
				edge = NULL;
			}
		}
	}

	vector<int> indicesToKeep;
	createFinalMesh(indicesToKeep, false, true);
	return finalMesh;
}

g_Part * MeshSimplification::getSimplifiedMesh(int numNodes, vector<int> & indicesToKeep, bool checkIntersection)
{
	int numberCurrentNodes = (int)originalMesh->nodes().numberOfItems();

	int i, j, numEdges, id1, previouslyDeleted, previouslyCollapsed;
	double edgeWeight;
	bool id1Present;
	vector<int> newNeighbors;
	priorityEdge prioEdge;
	g_PEdge * edge;
	g_Node firstNode, lastNode;
	g_PEdgeContainer allEdges;

	//build the tree structure:
	treeStructure.resize((int)originalMesh->nodes().numberOfItems());
	for(i = 0; i < (int)originalMesh->nodes().numberOfItems(); i++)
	{
		treeNode treeNd;
		treeNd.id = i+1;
		treeNd.parent = i+1;
		treeNd.child2.clear();
		treeStructure[i] = treeNd;
	}

	//collapse edges one by one
	previouslyDeleted = -1;
	priority_queue<priorityEdge, vector<priorityEdge>, greater<priorityEdge> > edgesToCollapse;
	while(numberCurrentNodes > numNodes)
	{
		//find edgesToCollapse and store them in a priority queue:
		if(previouslyDeleted == -1)
		{
			allEdges = workingMesh->uniqueEdges();
			numEdges = (int)allEdges.numberOfItems();
			for(i = 0; i < numEdges; i++)
			{
				firstNode = allEdges[i]->firstNode();
				lastNode = allEdges[i]->lastNode();
				edge = new g_PEdge(lastNode, firstNode);
				allEdges.insert(edge);
			}
		}
		else
		{
			allEdges.clear();
			for(i = 0; i < (int)newNeighbors.size(); i++)
			{
				edge = new g_PEdge(*(workingMesh->nodes()[previouslyCollapsed]), *(workingMesh->nodes()[newNeighbors[i]]));
				allEdges.insert(edge);
				edge = new g_PEdge(*(workingMesh->nodes()[newNeighbors[i]]),*(workingMesh->nodes()[previouslyCollapsed]));
				allEdges.insert(edge);
			}
		}
		numEdges = (int)allEdges.numberOfItems();
		for(i = 0; i < numEdges; i++)
		{
			edge = allEdges[i];
			id1 = edge->firstNode().id()-1;
			id1Present = false; 
			for(j = 0; j < (int)indicesToKeep.size(); j++)
			{
				if(id1 == indicesToKeep[j]) 
				{
					id1Present = true;
					break;
				}
			}
			if(id1Present) delete edge;
			else 
			{
				if(geomCost) edgeWeight = computeHalfEdgeCostGeom(edge, checkIntersection);
				else edgeWeight = computeHalfEdgeCostVolume(edge, checkIntersection);
				if(edgeWeight == -1) delete edge;
				else
				{
					prioEdge.firstNode = workingMesh->nodes()[edge->firstNode().id()-1];
					prioEdge.lastNode = workingMesh->nodes()[edge->lastNode().id()-1];
					prioEdge.weight = edgeWeight;
					edgesToCollapse.push(prioEdge);
					delete edge;
				}
			}
		}
		//collapse the edges:
		if(edgesToCollapse.empty()) break;
		else
		{
			edge = findBestHalfEdge(edgesToCollapse, checkIntersection);
			if(edge == NULL) break;
			else
			{
				halfEdgeCollapse(edge->firstNode().id()-1, edge->lastNode().id()-1, newNeighbors);
				previouslyDeleted = edge->firstNode().id()-1;
				previouslyCollapsed = edge->lastNode().id()-1;
				numberCurrentNodes--;
				delete edge;
				edge = NULL;
			}
		}
	}

	createFinalMesh(indicesToKeep, false, true);
	return finalMesh;
}

g_Part * MeshSimplification::getSimplifiedMesh(vector<int> & indicesToKeep, bool checkIntersection, 
											   bool keepNonCollapsables)
{
	int i, j, numEdges, id1, previouslyDeleted, previouslyCollapsed;
	bool id1Present, keepGoing, allDeleted;
	double edgeWeight;
	vector<int> newNeighbors;
	priorityEdge prioEdge;
	g_PEdge * edge;
	g_Node firstNode, lastNode;
	g_PEdgeContainer allEdges;
	allDeleted = true;

	//build the tree structure:
	treeStructure.resize((int)originalMesh->nodes().numberOfItems());
	for(i = 0; i < (int)originalMesh->nodes().numberOfItems(); i++)
	{
		treeNode treeNd;
		treeNd.id = i+1;
		treeNd.parent = i+1;
		treeNd.child2.clear();
		treeStructure[i] = treeNd;
	}

	//collapse edges one by one
	keepGoing = true;
	previouslyDeleted = -1;
	priority_queue<priorityEdge, vector<priorityEdge>, greater<priorityEdge> > edgesToCollapse;
	do
	{
		//find edgesToCollapse and store them in a priority queue:
		if(previouslyDeleted == -1)
		{
			allEdges = workingMesh->uniqueEdges();
			numEdges = (int)allEdges.numberOfItems();
			for(i = 0; i < numEdges; i++)
			{
				firstNode = allEdges[i]->firstNode();
				lastNode = allEdges[i]->lastNode();
				edge = new g_PEdge(lastNode, firstNode);
				allEdges.insert(edge);
			}
		}
		else
		{
			allEdges.clear();
			for(i = 0; i < (int)newNeighbors.size(); i++)
			{
				edge = new g_PEdge(*(workingMesh->nodes()[previouslyCollapsed]), *(workingMesh->nodes()[newNeighbors[i]]));
				allEdges.insert(edge);
				edge = new g_PEdge(*(workingMesh->nodes()[newNeighbors[i]]),*(workingMesh->nodes()[previouslyCollapsed]));
				allEdges.insert(edge);
			}
		}
		numEdges = (int)allEdges.numberOfItems();
		for(i = 0; i < numEdges; i++)
		{
			edge = allEdges[i];
			id1 = edge->firstNode().id()-1;
			id1Present = false; 
			for(j = 0; j < (int)indicesToKeep.size(); j++)
			{
				if(id1 == indicesToKeep[j]) 
				{
					id1Present = true;
					break;
				}
			}
			if(id1Present) delete edge;
			else 
			{
				if(geomCost) edgeWeight = computeHalfEdgeCostGeom(edge, checkIntersection);
				else edgeWeight = computeHalfEdgeCostVolume(edge, checkIntersection);
				if(edgeWeight == -1) delete edge;
				else
				{
					prioEdge.firstNode = workingMesh->nodes()[edge->firstNode().id()-1];
					prioEdge.lastNode = workingMesh->nodes()[edge->lastNode().id()-1];
					prioEdge.weight = edgeWeight;
					edgesToCollapse.push(prioEdge);
					delete edge;
				}
			}
		}
		//collapse the edges:
		if(edgesToCollapse.empty()) keepGoing = false;
		else
		{
			edge = findBestHalfEdge(edgesToCollapse, checkIntersection);
			if(edge == NULL) 
			{
				keepGoing = false;
				allDeleted = false;
			}
			else
			{
				halfEdgeCollapse(edge->firstNode().id()-1, edge->lastNode().id()-1, newNeighbors);
				previouslyDeleted = edge->firstNode().id()-1;
				previouslyCollapsed = edge->lastNode().id()-1;
				delete edge;
				edge = NULL;
				allDeleted = true;
			}
		}
	}while (keepGoing);

	createFinalMesh(indicesToKeep, allDeleted, keepNonCollapsables);
	return finalMesh;
}

void MeshSimplification::initOriginal()
{
	int i;
	g_Node * node;
	g_Element * element;

	//copy originalMesh to workingMesh
	if(workingMesh != NULL)
		delete workingMesh;
	workingMesh = new g_Part;
	for(i = 0; i < (int)originalMesh->nodes().numberOfItems(); i++)
	{
		node = new g_Node(originalMesh->nodes()[i]->coordinate());
		workingMesh->node(node);
	}
	for(i = 0; i < (int)originalMesh->elements().numberOfItems(); i++)
	{
		element = new g_Element();
		element->node(workingMesh->nodes()[originalMesh->elements()[i]->nodes()[0]->id()-1]);
		element->node(workingMesh->nodes()[originalMesh->elements()[i]->nodes()[1]->id()-1]);
		element->node(workingMesh->nodes()[originalMesh->elements()[i]->nodes()[2]->id()-1]);
		workingMesh->element(element);
	}
	storeOriginalMeshWeights();
}

void MeshSimplification::storeOriginalMeshWeights()
{
	int i;

	//build the tree structure:
	originalMeshWeights.resize((int)originalMesh->nodes().numberOfItems());
	for(i = 0; i < (int)originalMesh->nodes().numberOfItems(); i++)
	{
		treeNode treeNd;
		treeNd.id = i+1;
		treeNd.parent = i+1;
		treeNd.child2.clear();
		originalMeshWeights[i] = treeNd;
	}
	for(i = 0; i < (int)originalMesh->nodes().numberOfItems(); i++) storeMeanValueWeights(i, originalMeshWeights);
}

bool MeshSimplification::isMeshValid(g_Part * mesh)
{
	int i, j, numberElements, id10, id11, id12, id20, id21, id22;
	//make sure that the same element never occurs more than once:
	numberElements = (int)mesh->elements().numberOfItems();
	for(i = 0; i < numberElements; i++)
	{
		if((int)mesh->elements()[i]->nodes().numberOfItems() > 0)
		{
			id10 = mesh->elements()[i]->nodes()[0]->id();
			id11 = mesh->elements()[i]->nodes()[1]->id();
			id12 = mesh->elements()[i]->nodes()[2]->id();
			for(j = 0; j < i; j++)
			{
				if((int)mesh->elements()[j]->nodes().numberOfItems() > 0)
				{
					id20 = mesh->elements()[j]->nodes()[0]->id();
					id21 = mesh->elements()[j]->nodes()[1]->id();
					id22 = mesh->elements()[j]->nodes()[2]->id();
					if(((id10 == id20) || (id10 == id21) || (id10 == id22)) && 
						((id11 == id20) || (id11 == id21) || (id11 == id22)) &&
						((id12 == id20) || (id12 == id21) || (id12 == id22)))
					{
						cout<<"Mesh invalid "<<i<<" ("<<id10<<","<<id11<<","<<id12<<") "<<
							j<<" ("<<id20<<","<<id21<<","<<id22<<")"<<endl;
						return false;
					}
				}
			}
		}
	}
	return true;
}

//orients half edge such that the first vertex should be collapsed into the second vertex
g_PEdge * MeshSimplification::findBestHalfEdge(priority_queue<priorityEdge, vector<priorityEdge>, greater<priorityEdge> > &
											   edgesToCollapse, bool checkIntersection)
{
	g_PEdge * edge = NULL;
	priorityEdge prioEdge;
	bool validEdge = false;
	while((!edgesToCollapse.empty()) && (!validEdge))
	{
		prioEdge = edgesToCollapse.top();
		edgesToCollapse.pop();
		edge = new g_PEdge(*(prioEdge.firstNode), *(prioEdge.lastNode));
		validEdge = isLegalMove(edge, checkIntersection);
		if(!validEdge) 
		{
			delete edge;
			edge = NULL;
		}
	}

	return edge;
}

// Cost is set to -1 if the collapse is illegal (see paper by Hoppe for conditions)
// Otherwise, the cost gives the change in volume:
double MeshSimplification::computeHalfEdgeCostVolume(g_PEdge * halfEdge, bool checkIntersection)
{
	if(!isLegalMove(halfEdge, false)) return -1;

	int i, id1, id2, localId, numElements;
	double cost, angleCost, minAngleAfter, angle;
	g_Element * element;
	g_Vector node1, node2, triNode1, triNode2, helpVec1, helpVec2;
	
	id1 = halfEdge->firstNode().id()-1;
	id2 = halfEdge->lastNode().id()-1;
	numElements = workingMesh->nodes()[id1]->elements().numberOfItems();
	node1 = workingMesh->nodes()[id1]->coordinate();
	node2 = workingMesh->nodes()[id2]->coordinate();
	cost = 0;
	angleCost = 0;

	for(i = 0; i < numElements; i++)
	{
		element = workingMesh->nodes()[id1]->elements()[i];
		if((element->nodes()[0]->id()-1 == id2) || (element->nodes()[1]->id()-1 == id2) || 
			(element->nodes()[2]->id()-1 == id2)) continue;
		if(element->nodes()[0]->id()-1 == id1) localId = 0;
		else if(element->nodes()[1]->id()-1 == id1) localId = 1;
		else if(element->nodes()[2]->id()-1 == id1) localId = 2;
		triNode1 = element->nodes()[(localId+1)%3]->coordinate();
		triNode2 = element->nodes()[(localId+2)%3]->coordinate();
		//compute and add the change in volume to the cost
		//Formula given by Lindstrom and Turk
		cost += fabs(-node2.Dot(triNode1.Cross(triNode2)) + node1.Dot(triNode1.Cross(triNode2)) - node1.Dot(node2.Cross(triNode2)) + node1.Dot(node2.Cross(triNode1)));

		//Compute and add the inverse minimum angle:
		angle = (node2-triNode1).AngleBetween(triNode2-node2);
		if(i == 0) minAngleAfter = angle;
		else minAngleAfter = min(angle, minAngleAfter);
		angle = (triNode1-triNode2).AngleBetween(node2-triNode1);
		minAngleAfter = min(angle, minAngleAfter);
		angle = (triNode2-node2).AngleBetween(triNode1-triNode2);
		minAngleAfter = min(angle, minAngleAfter);
	}

	cost = cost/(max(pow(resolution, 3.0), MDS_THRESHOLD_DIVIDE0));
	angleCost = 1.0/(max(minAngleAfter, MDS_THRESHOLD_DIVIDE0));

//	cost = cost + 0.001 * angleCost;

	return cost;
}

//Cost is set to -1 is the collapse is illegal. Otherwise, the cost is according to Garland and Heckbert's geometry criterion (related to curvature)
double MeshSimplification::computeHalfEdgeCostGeom(g_PEdge * halfEdge, bool checkIntersection)
{
	if(!isLegalMove(halfEdge, false)) return -1;

	int i, id1, id2, numElements;
	double a, b, c, d, cost;
	double Q1 [16];
	g_Vector normal, node1, node2;
	g_Element * element;
	
	id1 = halfEdge->firstNode().id()-1;
	id2 = halfEdge->lastNode().id()-1;
	node1 = workingMesh->nodes()[id1]->coordinate();
	node2 = workingMesh->nodes()[id2]->coordinate();

	//Compute the matrix for the first vertex:
	numElements = workingMesh->nodes()[id1]->elements().numberOfItems();
	for(i = 0; i < 16; i++) Q1[i] = 0.0;
	for(i = 0; i < numElements; i++)
	{
		element = workingMesh->nodes()[id1]->elements()[i];
		normal = (element->nodes()[1]->coordinate()-element->nodes()[0]->coordinate()).Cross((element->nodes()[2]->coordinate()-element->nodes()[0]->coordinate()));
		normal.Normalize();
		a = normal.x();
		b = normal.y();
		c = normal.z();
		d = -(a * element->nodes()[0]->coordinate().x() + b * element->nodes()[0]->coordinate().y() + c * element->nodes()[0]->coordinate().z());

		Q1[0] += a*a; Q1[1] += a*b; Q1[2] += a*c; Q1[3] += a*d;
		Q1[4] += b*a; Q1[5] += b*b; Q1[6] += b*c; Q1[7] += b*d;
		Q1[8] += c*a; Q1[9] += c*b; Q1[10] += c*c; Q1[11] += c*d;
		Q1[12] += d*a; Q1[13] += d*b; Q1[14] += d*c; Q1[15] += d*d;
	}

	//Compute the cost of keeping node2 instead of node1:
	Q1[0] = node2.x() * Q1[0] + node2.y() * Q1[4] + node2.z() * Q1[8] + Q1[12];
	Q1[1] = node2.x() * Q1[1] + node2.y() * Q1[5] + node2.z() * Q1[9] + Q1[13];
	Q1[2] = node2.x() * Q1[2] + node2.y() * Q1[6] + node2.z() * Q1[10] + Q1[14];
	Q1[3] = node2.x() * Q1[3] + node2.y() * Q1[7] + node2.z() * Q1[11] + Q1[15];
	cost = node2.x() * Q1[0] + node2.y() * Q1[1] + node2.z() * Q1[2] + Q1[3];

	return cost;
}

//compute whether or not the edge collapse is a legal move
bool MeshSimplification::isLegalMove(g_PEdge * halfEdge, bool checkIntersection)
{
	int i, j, numElements, id00, id01, id02;
	bool returnVal, isBoundVert1, isBoundVert2, triangleFound;
	g_Node * node1, * node2;
	g_NodeContainer neighbors1, neighbors2, intersection, triangle;

	numElements = workingMesh->nodes().numberOfItems();
	returnVal = true;

	//test if the edge is still in the current mesh:
	id00 = halfEdge->firstNode().id();
	if(treeStructure[id00-1].parent != id00) {
#ifdef DEBUG_MESH_SIMPLIFICATION
	  cerr << "No parent 0" << endl;
#endif
	  return false;
	}
	id00 = halfEdge->lastNode().id();
	if(treeStructure[id00-1].parent != id00) {
#ifdef DEBUG_MESH_SIMPLIFICATION
	  cerr << "No parent 1" << endl;
#endif
	  return false;
	}
	node1 = workingMesh->nodes()[halfEdge->firstNode().id()-1];
	node2 = workingMesh->nodes()[halfEdge->lastNode().id()-1];
	g_PtrAdapter<g_Node*> adapter1(neighbors1);
	for(i = 0; i < (int)node1->elements().numberOfItems(); i++)
	{
#ifdef DEBUG_MESH_SIMPLIFICATION
	         if ( node1->id() == 1156 ) {
		   cerr << "1156 : " << i << endl;
		   cerr << node1->elements()[i]->nodes()[0]->id() << " ";
		   cerr << node1->elements()[i]->nodes()[1]->id() << " ";
		   cerr << node1->elements()[i]->nodes()[2]->id() << endl;
		 }
#endif
		if(node1->elements()[i]->nodes()[0]->id() != node1->id())neighbors1.insert(node1->elements()[i]->nodes()[0]);
		if(node1->elements()[i]->nodes()[1]->id() != node1->id())neighbors1.insert(node1->elements()[i]->nodes()[1]);
		if(node1->elements()[i]->nodes()[2]->id() != node1->id())neighbors1.insert(node1->elements()[i]->nodes()[2]);
	}
	adapter1.removeDuplicates();
	g_PtrAdapter<g_Node*> adapter2(neighbors2);
	for(i = 0; i < (int)node2->elements().numberOfItems(); i++)
	{
		if(node2->elements()[i]->nodes()[0]->id() != node2->id())neighbors2.insert(node2->elements()[i]->nodes()[0]);
		if(node2->elements()[i]->nodes()[1]->id() != node2->id())neighbors2.insert(node2->elements()[i]->nodes()[1]);
		if(node2->elements()[i]->nodes()[2]->id() != node2->id())neighbors2.insert(node2->elements()[i]->nodes()[2]);
	}
	adapter2.removeDuplicates();

	//enough triangles in mesh
	if(numElements < 5) 
	{
		returnVal = false;
	}
	else
	{
		//all shared nodes are part of a triangle
		for(i = 0; i < (int)neighbors1.numberOfItems(); i++)
		{
			for(j = 0; j < (int)neighbors2.numberOfItems(); j++)
			{
				if(neighbors1[i]->id() == neighbors2[j]->id()) intersection.insert(neighbors1[i]);
			}
		}
		
		for(i = 0; i < (int)intersection.numberOfItems(); i++)
		{
			triangleFound = false;
			for(j = 0; j < (int)intersection[i]->elements().numberOfItems(); j++)
			{
				id00 = intersection[i]->elements()[j]->nodes()[0]->id();
				id01 = intersection[i]->elements()[j]->nodes()[1]->id();
				id02 = intersection[i]->elements()[j]->nodes()[2]->id();
				if(((id00 == node1->id()) || (id01 == node1->id()) || (id02 == node1->id())) &&
					((id00 == node2->id()) || (id01 == node2->id()) || (id02 == node2->id())) &&
					((id00 == intersection[i]->id()) || (id01 == intersection[i]->id()) || (id02 == intersection[i]->id())))
				{
					triangleFound = true;
					break;
				}
			}
			if(!triangleFound) {
#ifdef DEBUG_MESH_SIMPLIFICATION
			  cerr << "Triangle not found" << endl;
#endif
			  returnVal = false;
			}
		}
	}

	//if two nodes are boundary nodes, then edge is a boundary edge
	if(returnVal)
	{
		if(halfEdge->elements().numberOfItems() > 1)
		{
			isBoundVert1 = isBoundaryNode(node1->id()-1);
			isBoundVert2 = isBoundaryNode(node2->id()-1);
			if(isBoundVert1 && isBoundVert2) 
			{
			  cerr << "Bondary nodes" << endl;
			  returnVal = false;
			}
		}
	}	

	if(leadsProjectionProblems(halfEdge)) {
#ifdef DEBUG_MESH_SIMPLIFICATION
	  cerr << "projection" << endl;
#endif
	  returnVal = false;
	}
	//test if the collapse would lead to self-intersecting triangles:
	if(checkIntersection)
	{ 
	  if(leadsSelfIntersection(halfEdge)) {
#ifdef DEBUG_MESH_SIMPLIFICATION
	    cerr << "Self intersection" << endl;
#endif
	    returnVal = false;
	  }
	}

	return returnVal;
}

bool MeshSimplification::leadsProjectionProblems(g_PEdge * halfEdge)
{
	int i, j;
	bool id0, id2;
	g_Node node = halfEdge->firstNode();
	g_Vector projectedNode;
	g_NodeContainer elementNeighbors;
	g_NodeContainer projectedNeighbors;
	vector<bool> reversed ((int)node.elements().numberOfItems());
#ifdef DEBUG_MESH_SIMPLIFICATION
	cerr << "First node #neighbors: " << node.elements().numberOfItems() << endl;
#endif
	//Check if the neighborhood projections of the first node's neighborhood yield a self intersecting boundary:
	for(i = 0; i < (int)node.elements().numberOfItems(); i++)
	{
		id0 = id2 = false;
		if(node.elements()[i]->nodes()[0]->id() != node.id())
		{
			elementNeighbors.insert(node.elements()[i]->nodes()[0]);
			id0 = true;
		}
		if(node.elements()[i]->nodes()[1]->id() != node.id())
			elementNeighbors.insert(node.elements()[i]->nodes()[1]);
		if(node.elements()[i]->nodes()[2]->id() != node.id())
		{
			elementNeighbors.insert(node.elements()[i]->nodes()[2]);
			id2 = true;
		}
		if(id0 && id2) reversed[i] = true;
		else reversed[i] = false;
	}

	//find tangent plane and project everything to it:
	g_Vector l;
	l.Set(0, 0, 0);

	g_NodeContainer uniqueNeighbors;
	g_PtrAdapter <g_Node *> adapter(uniqueNeighbors);
	
#ifdef DEBUG_MESH_SIMPLIFICATION
	cerr << "# Element Neighbors: " << elementNeighbors.numberOfItems() << endl;
#endif
	if((int)elementNeighbors.numberOfItems() > 0)
	{
		for(i = 0; i < (int)elementNeighbors.numberOfItems(); i++) uniqueNeighbors.insert(elementNeighbors[i]);
		adapter.removeDuplicates();
		for(i = 0; i < (int)uniqueNeighbors.numberOfItems(); i++)
			l = l + uniqueNeighbors[i]->coordinate();
		l = l * (1.0/(double)uniqueNeighbors.numberOfItems());
	}
#ifdef DEBUG_MESH_SIMPLIFICATION
	cerr << "# Unique Neighbors: " << uniqueNeighbors.numberOfItems() << endl;
#endif
	g_Vector normal;
	normal.Set(0, 0, 0);
	for(i = 0; i < (int)node.elements().numberOfItems(); i++)
		normal = normal + (elementNeighbors[2*i]->coordinate()-l).Cross(elementNeighbors[2*i+1]->coordinate()-l);		
#ifdef DEBUG_MESH_SIMPLIFICATION
	cerr << "Normal " << normal << endl;
#endif
	if(normal.Length() < MDS_THRESHOLD_DIVIDE0) return true;

	normal.Normalize();


	double d = 0.0;
	if((int)elementNeighbors.numberOfItems() > 0)
	{
		for(i = 0; i < (int)elementNeighbors.numberOfItems(); i++)
			d += normal.Dot(elementNeighbors[i]->coordinate());
		d = -d/(double)elementNeighbors.numberOfItems();
		for(i = 0; i < (int)elementNeighbors.numberOfItems(); i++)
		{
			g_Node * projNode = new g_Node;
			projNode->coordinate(elementNeighbors[i]->coordinate()-(d+elementNeighbors[i]->coordinate().Dot(normal))*normal);
			projectedNeighbors.insert(projNode);
		}
	}
	
	bool intersection = false;
	for(i = 0; i < (int)node.elements().numberOfItems(); i++)
	{
		for(j = i+1; j < (int)node.elements().numberOfItems(); j++)
		{
			if((elementNeighbors[2*i]->id() == elementNeighbors[2*j]->id()) || (elementNeighbors[2*i+1]->id() == elementNeighbors[2*j]->id()) || 
				(elementNeighbors[2*i]->id() == elementNeighbors[2*j+1]->id()) || (elementNeighbors[2*i+1]->id() == elementNeighbors[2*j+1]->id()))
			{
				g_Vector vec1, vec2;
				if((elementNeighbors[2*i]->id() == elementNeighbors[2*j]->id()))
				{
					vec1 = elementNeighbors[2*i]->coordinate() - elementNeighbors[2*i+1]->coordinate();
					vec2 = elementNeighbors[2*j+1]->coordinate() - elementNeighbors[2*j]->coordinate();
				}
				else if((elementNeighbors[2*i+1]->id() == elementNeighbors[2*j]->id()))
				{
					vec1 = elementNeighbors[2*i+1]->coordinate() - elementNeighbors[2*i]->coordinate();
					vec2 = elementNeighbors[2*j+1]->coordinate() - elementNeighbors[2*j]->coordinate();
				}
				else if((elementNeighbors[2*i]->id() == elementNeighbors[2*j+1]->id()))
				{
					vec1 = elementNeighbors[2*i]->coordinate() - elementNeighbors[2*i+1]->coordinate();
					vec2 = elementNeighbors[2*j]->coordinate() - elementNeighbors[2*j+1]->coordinate();
				}
				else if((elementNeighbors[2*i+1]->id() == elementNeighbors[2*j+1]->id()))
				{
					vec1 = elementNeighbors[2*i+1]->coordinate() - elementNeighbors[2*i]->coordinate();
					vec2 = elementNeighbors[2*j]->coordinate() - elementNeighbors[2*j+1]->coordinate();
				}
				vec1.Normalize();
				vec2.Normalize();

				if(fabs(vec1.Dot(vec2) + 1.0) < MDS_THRESHOLD_DIVIDE0)
				{
					intersection = true;
					break;
				}
			}
			else
			{
				double dist = MinimumDistance::minimumDistance(projectedNeighbors[2*i]->coordinate(), projectedNeighbors[2*i+1]->coordinate() ,
					projectedNeighbors[2*j]->coordinate(), projectedNeighbors[2*j+1]->coordinate(), false);
				if(dist < MDS_THRESHOLD_DIVIDE0)
				{
					intersection = true;
					break;
				}
			}
		}
	}
#ifdef DEBUG_MESH_SIMPLIFICATION
	if(intersection)
	  cerr << "Intersection" << endl;
	else
	  cerr << "No interesection" << endl;
#endif
	if(intersection)
	{
		for(i = 0; i < (int)elementNeighbors.numberOfItems(); i++) delete projectedNeighbors[i];
		projectedNeighbors.clear();
		return true;
	}

	//Check if the projection of the first node is contained in the CH of the projected boundary:
	projectedNode = node.coordinate()-(d+node.coordinate().Dot(normal))*normal;
	bool inside = false;
	for(i = 1; i < (int)node.elements().numberOfItems(); i++)
	{
		if((elementNeighbors[0]->id() == elementNeighbors[2*i]->id()) || (elementNeighbors[0]->id() == elementNeighbors[2*i+1]->id())) continue;
		g_Vector turn1 = (projectedNeighbors[2*i]->coordinate() - projectedNeighbors[0]->coordinate()).Cross(projectedNode - projectedNeighbors[2*i]->coordinate());
		g_Vector turn2 = (projectedNeighbors[2*i+1]->coordinate() - projectedNeighbors[2*i]->coordinate()).Cross(projectedNode - projectedNeighbors[2*i+1]->coordinate());
		g_Vector turn3 = (projectedNeighbors[0]->coordinate() - projectedNeighbors[2*i+1]->coordinate()).Cross(projectedNode - projectedNeighbors[0]->coordinate());

		if((turn1.Dot(turn2) > 0) && (turn1.Dot(turn3) > 0)) 
		{
			inside = true;
			break;
		}
	}
#ifdef DEBUG_MESH_SIMPLIFICATION
	if(inside)
	  cerr << "Inside" << endl;
	else
	  cerr << "Not inside" << endl;
#endif
	if(!inside) 
	{
		for(i = 0; i < (int)elementNeighbors.numberOfItems(); i++) delete projectedNeighbors[i];
			projectedNeighbors.clear();
		return true;
	}

	//Check if projection is in star-shaped kernel:
	bool inKernel = true;
	g_Vector baseTurn, currentTurn;
	for(i = 0; i < (int)node.elements().numberOfItems(); i++)
	{
		if(i == 0) 
		{
			if(reversed[i]) baseTurn = (projectedNeighbors[2*i]->coordinate() - projectedNode).Cross(projectedNeighbors[2*i+1]->coordinate() - projectedNeighbors[2*i]->coordinate());
			else baseTurn = (projectedNeighbors[2*i+1]->coordinate() - projectedNode).Cross(projectedNeighbors[2*i]->coordinate() - projectedNeighbors[2*i+1]->coordinate());
		}
		else
		{
			if(reversed[i]) currentTurn = (projectedNeighbors[2*i]->coordinate() - projectedNode).Cross(projectedNeighbors[2*i+1]->coordinate() - projectedNeighbors[2*i]->coordinate());
			else currentTurn = (projectedNeighbors[2*i+1]->coordinate() - projectedNode).Cross(projectedNeighbors[2*i]->coordinate() - projectedNeighbors[2*i+1]->coordinate());
			if(baseTurn.Dot(currentTurn) < 0) 
			{
				inKernel = false;
				break;
			}
		}
	}

	for(i = 0; i < (int)elementNeighbors.numberOfItems(); i++) delete projectedNeighbors[i];
		projectedNeighbors.clear();
#ifdef DEBUG_MESH_SIMPLIFICATION
	if(inKernel)
	  cerr << "inKernel" << endl;
	else
	  cerr << "Not inKernel" << endl;
#endif

	if(!inKernel) return true;
	else return false;
}


bool MeshSimplification::leadsSelfIntersection(g_PEdge * halfEdge)
{
#if USE_ANN	
	//Test the nearest neighbors in the kd-tree (neighbors found based on vertex distance)
	//Build a kd-tree for the nodes
	int i, j, k, numNodes, numLocalElements, kValueElement;
	bool node1Present, node2Present, returnVal;
	g_Node * node1, * node2;
	g_NodeContainer elem1Cont;
	g_ElementContainer neighborCont;
#ifdef USE_MOELLER
	float * a, * b, * c, * d, * e, * f;
	a = new float[3];
	b = new float[3];
	c = new float[3];
	d = new float[3];
	e = new float[3];
	f = new float[3];
#else
	tetgenmesh * tetmesh = new tetgenmesh();

	double * a, * b, * c, * d, * e, * f;
	a = new double[3];
	b = new double[3];
	c = new double[3];
	d = new double[3];
	e = new double[3];
	f = new double[3];
#endif

	returnVal = false;
	numNodes = workingMesh->nodes().numberOfItems();
	ANNpointArray		dataPts;
	ANNpoint			queryPt;
	ANNidxArray			nnIdx;	
	ANNdistArray		dists;
	ANNkd_tree*			kdTree;	
	queryPt = annAllocPt(3);				
	dataPts = annAllocPts(numNodes, 3);		
	nnIdx = new ANNidx[numNodes];	
	dists = new ANNdist[numNodes];
	for(i = 0; i < numNodes; i++)
	{
		dataPts[i][0] = workingMesh->nodes()[i]->coordinate().x();
		dataPts[i][1] = workingMesh->nodes()[i]->coordinate().y();
		dataPts[i][2] = workingMesh->nodes()[i]->coordinate().z();
	}
	kdTree = new ANNkd_tree(dataPts, numNodes, 3);

	//build the changed elements:
	node1 = workingMesh->nodes()[halfEdge->firstNode().id()-1];
	node2 = workingMesh->nodes()[halfEdge->lastNode().id()-1];
	neighborCont = node1->elements();
	neighborCont.append(node2->elements());
	numLocalElements = neighborCont.numberOfItems();

	//Only do intersection test for the changed elements:
	for(i = 0; i < numLocalElements; i++)
	{
		if(neighborCont[i]->nodes().numberOfItems() != 3) continue;
		node1Present = node2Present = false;
		elem1Cont.clear();
		for(j = 0; j < 3; j++)
		{
			if(neighborCont[i]->nodes()[j] == node1) 
			{
				node1Present = true;
				elem1Cont.insert(node2);
			}
			else if(neighborCont[i]->nodes()[j] == node2)
			{
				node2Present = true;
				elem1Cont.insert(node2);
			}
			else elem1Cont.insert(neighborCont[i]->nodes()[j]);
		}
		if((!node1Present) || (node1Present && node2Present)) continue;

		g_Element helpElem;
		for(j = 0; j < 3; j++) helpElem.node(elem1Cont[j]);

		// Check angles of helpElem: if the angle is too small, then reject the edge collapse
		double angle = (helpElem.nodes()[1]->coordinate()-helpElem.nodes()[0]->coordinate()).AngleBetween(
			helpElem.nodes()[2]->coordinate()-helpElem.nodes()[0]->coordinate());
		angle = min(angle, (helpElem.nodes()[0]->coordinate()-helpElem.nodes()[1]->coordinate()).AngleBetween(
			helpElem.nodes()[2]->coordinate()-helpElem.nodes()[1]->coordinate()));
		angle = min(angle, (helpElem.nodes()[0]->coordinate()-helpElem.nodes()[2]->coordinate()).AngleBetween(
			helpElem.nodes()[1]->coordinate()-helpElem.nodes()[2]->coordinate()));
		if(angle < 0.1745)
		{
			for(j = 0; j < 3; j++) elem1Cont[j]->removeElement(helpElem);
			returnVal = true;
			goto INTERSECT_EXIT;
		}

		//find the k nearest neighbors
		kValueElement = elem1Cont[0]->elements().numberOfItems() + elem1Cont[1]->elements().numberOfItems() + 
									elem1Cont[2]->elements().numberOfItems() - 3 + 1;
		g_ElementContainer contToCheck;
		g_PtrAdapter<g_Element *> adapt(contToCheck);

		for(j = 0; j < 3; j++)
		{
			queryPt[0] = elem1Cont[j]->coordinate().x();
			queryPt[1] = elem1Cont[j]->coordinate().y();
			queryPt[2] = elem1Cont[j]->coordinate().z();
			kdTree->annkPriSearch(queryPt, kValueElement, nnIdx, dists);	
			for(k = 0; k < kValueElement; k++)
				contToCheck.append(workingMesh->nodes()[nnIdx[k]]->elements());
		}
		adapt.removeDuplicates();
		int contNum = contToCheck.numberOfItems();

		//Find nearest neighbor
		for(j = 0; j < contNum; j++)
		{
			if(contToCheck[j]->nodes().numberOfItems() != 3) continue;

			//Test using tetgen:
			if((contToCheck[j]->nodes()[0] == node2) || (contToCheck[j]->nodes()[1] == node2) ||
			   (contToCheck[j]->nodes()[2] == node2))
				continue;

			g_Element testElement;
			g_Node * node = new g_Node(elem1Cont[0]->coordinate());
			testElement.node(node);
			node = new g_Node(elem1Cont[1]->coordinate());
			testElement.node(node);
			node = new g_Node(elem1Cont[2]->coordinate());
			testElement.node(node);

			a[0] = testElement.nodes()[0]->coordinate().x(); 
			a[1] = testElement.nodes()[0]->coordinate().y(); 
			a[2] = testElement.nodes()[0]->coordinate().z();

			b[0] = testElement.nodes()[1]->coordinate().x(); 
			b[1] = testElement.nodes()[1]->coordinate().y(); 
			b[2] = testElement.nodes()[1]->coordinate().z();

			c[0] = testElement.nodes()[2]->coordinate().x(); 
			c[1] = testElement.nodes()[2]->coordinate().y(); 
			c[2] = testElement.nodes()[2]->coordinate().z();

			d[0] = contToCheck[j]->nodes()[0]->coordinate().x();
			d[1] = contToCheck[j]->nodes()[0]->coordinate().y();
			d[2] = contToCheck[j]->nodes()[0]->coordinate().z();

			e[0] = contToCheck[j]->nodes()[1]->coordinate().x();
			e[1] = contToCheck[j]->nodes()[1]->coordinate().y();
			e[2] = contToCheck[j]->nodes()[1]->coordinate().z();

			f[0] = contToCheck[j]->nodes()[2]->coordinate().x();
			f[1] = contToCheck[j]->nodes()[2]->coordinate().y();
			f[2] = contToCheck[j]->nodes()[2]->coordinate().z();

#ifdef USE_MOELLER
			if ( !NoDivTriTriIsect(a, b, c, d, e, f))
#else
			  tetgenmesh::interresult inter = static_cast<tetgenmesh::interresult>(tetmesh->tri_tri_inter(a, b, c, d, e, f));
			  // if(inter == tetgenmesh::interresult::INTERSECT)
			  // JL: Does VC++ need the enum in-between?
			if(inter == tetgenmesh::INTERSECT)
#endif
			{
				for(k = 0; k < 3; k++) elem1Cont[k]->removeElement(helpElem);
				for(k = 0; k < 3; k++) delete testElement.nodes()[k];
				returnVal = true;
				goto INTERSECT_EXIT;
			}
			for(k = 0; k < 3; k++) delete testElement.nodes()[k];
		}
		for(j = 0; j < 3; j++) elem1Cont[j]->removeElement(helpElem);
	}

	INTERSECT_EXIT:

	//Delete kd-tree
	delete [] nnIdx;
	delete [] dists;
	annDeallocPt(queryPt);
	annDeallocPts(dataPts);
	delete kdTree;
	annClose();	

#ifndef USE_MOELLER
	delete tetmesh;
#endif
	//Delete other variables
	delete [] a;
	delete [] b;
	delete [] c;
	delete [] d;
	delete [] e;
	delete [] f;

	return returnVal;

#else
	int i, j, numElements, numLocalElements;
	bool node1Present, node2Present;
	g_Node * node1, * node2;
	g_NodeContainer elem1Cont;
	g_ElementContainer neighborCont;

	node1 = workingMesh->nodes()[halfEdge->firstNode().id()-1];
	node2 = workingMesh->nodes()[halfEdge->lastNode().id()-1];
	neighborCont = node1->elements();
	neighborCont.append(node2->elements());

	numElements = workingMesh->elements().numberOfItems();
	numLocalElements = neighborCont.numberOfItems();
	
	for(i = 0; i < numLocalElements; i++)
	{
		if(neighborCont[i]->nodes().numberOfItems() != 3) continue;
		node1Present = node2Present = false;
		elem1Cont.clear();
		for(j = 0; j < 3; j++)
		{
			if(neighborCont[i]->nodes()[j] == node1) 
			{
				node1Present = true;
				elem1Cont.insert(node2);
			}
			else if(neighborCont[i]->nodes()[j] == node2)
			{
				node2Present = true;
				elem1Cont.insert(node2);
			}
			else elem1Cont.insert(neighborCont[i]->nodes()[j]);
		}
		if(node1Present && node2Present) continue;
		g_Element helpElem;
		for(j = 0; j < 3; j++) helpElem.node(elem1Cont[j]);

		for(j = 0; j < numElements; j++)
		{
			if(workingMesh->elements()[j]->nodes().numberOfItems() != 3) continue;
			if((workingMesh->elements()[j]->nodes()[0] == node1) || (workingMesh->elements()[j]->nodes()[0] == node2) || 
				(workingMesh->elements()[j]->nodes()[1] == node1) || (workingMesh->elements()[j]->nodes()[1] == node2) || 
				(workingMesh->elements()[j]->nodes()[2] == node1) || (workingMesh->elements()[j]->nodes()[2] == node2))
				continue;
			if((workingMesh->elements()[j]->nodes()[0] == elem1Cont[0]) || 
				(workingMesh->elements()[j]->nodes()[0] == elem1Cont[1]) ||
				(workingMesh->elements()[j]->nodes()[0] == elem1Cont[2]) ||
				(workingMesh->elements()[j]->nodes()[1] == elem1Cont[0]) ||
				(workingMesh->elements()[j]->nodes()[1] == elem1Cont[1]) ||
				(workingMesh->elements()[j]->nodes()[1] == elem1Cont[2]) ||
				(workingMesh->elements()[j]->nodes()[2] == elem1Cont[0]) ||
				(workingMesh->elements()[j]->nodes()[2] == elem1Cont[1]) ||
				(workingMesh->elements()[j]->nodes()[2] == elem1Cont[2]))
				continue;
			if(MinimumDistance::intersect(helpElem, *(workingMesh->elements()[j]))) 
			{
				for(k = 0; k < 3; k++) elem1Cont[k]->removeElement(helpElem);
				return true;
			}
		}

		for(j = 0; j < 3; j++) elem1Cont[j]->removeElement(helpElem);
	}

	return false;
#endif
}

bool MeshSimplification::isBoundaryNode(int nodeIdx)
{
	int i, j;
	bool isBoundVert = false;
	g_Node * node = workingMesh->nodes()[nodeIdx];
	for(i = 0; i < (int)node->elements().numberOfItems(); i++)
	{
		if(isBoundVert) break;
		g_PEdgeContainer tempCont = node->elements()[i]->pEdges();
		for(j = 0; j < (int)tempCont.numberOfItems(); j++)
		{
			if((tempCont[j]->elements().numberOfItems() == 1) && 
				((tempCont[j]->firstNode().id() == node->id()) || (tempCont[j]->lastNode().id() == node->id())))
				isBoundVert = true;
			delete tempCont[j];
			tempCont[j] = NULL;
		}
		tempCont.clear();
	}
	return isBoundVert;
}

//collapses the first vertex into the second vertex
void MeshSimplification::halfEdgeCollapse(int id1, int id2, vector<int> & newNeighbors)
{
	int i, j, numElements, localId;
	bool nodeInThere;
	vector<int> commonNodes;
	vector<int> replacedNodes;
	g_Node * nodeToDelete, * collapseNode;
	g_ElementContainer elCont;
	g_Element * element;
	nodeToDelete = workingMesh->nodes()[id1];
	collapseNode = workingMesh->nodes()[id2];
	elCont = nodeToDelete->elements();
	numElements = elCont.numberOfItems();
	
#ifdef DEBUG_MESH_SIMPLIFICATION
	nodeToDelete->print(cerr);
	collapseNode->print(cerr);
#endif

	//update the tree structure:
	treeStructure[nodeToDelete->id()-1].parent = collapseNode->id();
	treeStructure[collapseNode->id()-1].child2.push_back(nodeToDelete->id());
	storeMeanValueWeights(nodeToDelete->id()-1, treeStructure);

	for(i = 0; i < numElements; i++)
	{
		element = elCont[i];
		if(element->nodes()[0]->id()-1 == id1) localId = 0;
		else if(element->nodes()[1]->id()-1 == id1) localId = 1;
		else if(element->nodes()[2]->id()-1 == id1) localId = 2;
		//delete triangles along the edge by removing them from their vertices and by removing their vertex lists 
		if((element->nodes()[0]->id()-1 == id2) || (element->nodes()[1]->id()-1 == id2) || 
			(element->nodes()[2]->id()-1 == id2))
		{
#ifdef DEBUG_MESH_SIMPLIFICATION
		        cerr << "Removing: " << element->id() << endl;
#endif
			if(element->nodes()[(localId+1)%3]->id()-1 == id2) commonNodes.push_back((localId+2)%3);
			else if(element->nodes()[(localId+2)%3]->id()-1 == id2) commonNodes.push_back((localId+1)%3);
#ifdef DEBUG_MESH_SIMPLIFICATION
			cerr << "from: " << element->nodes()[0]->id() << endl;
#endif
			element->nodes()[0]->removeElement(element);
#ifdef DEBUG_MESH_SIMPLIFICATION
			cerr << "from: " << element->nodes()[1]->id() << endl;
#endif
			element->nodes()[1]->removeElement(element);
#ifdef DEBUG_MESH_SIMPLIFICATION
			cerr << "from: " << element->nodes()[2]->id() << endl;
#endif
			element->nodes()[2]->removeElement(element);
			element->emptyNodeList();
		}
		//Change the other triangles to use collapseNode instead of nodeToDelete
		else
		{
			element->replaceNodeAt(localId, collapseNode);			
			nodeToDelete->removeElement(element);
			//check if the element is already in the mesh:
			nodeInThere = false;
			int id00, id01, id02, id10, id11, id12;
			id00 = element->nodes()[0]->id();
			id01 = element->nodes()[1]->id();
			id02 = element->nodes()[2]->id();
			for(j = 0; j < (int)collapseNode->elements().numberOfItems(); j++)
			{
				id10 = collapseNode->elements()[j]->nodes()[0]->id();
				id11 = collapseNode->elements()[j]->nodes()[1]->id();
				id12 = collapseNode->elements()[j]->nodes()[2]->id();
				if(((id00 == id10) || (id00 == id11) || (id00 == id12)) &&
					((id01 == id10) || (id01 == id11) || (id01 == id12)) &&
					((id02 == id10) || (id02 == id11) || (id02 == id12)))
				{
					nodeInThere = true;
					break;
				}
			}
			//If so, simply delete the old element
			if(nodeInThere)
			{
				element->nodes()[(localId+1)%3]->removeElement(element);
				element->nodes()[(localId+2)%3]->removeElement(element);
#ifdef DEBUG_MESH_SIMPLIFICATION
				cerr << "Empty node list: " << element->id() << endl;
#endif
				element->emptyNodeList();
			}
			else
			{
				replacedNodes.push_back(element->nodes()[(localId+1)%3]->id()-1);
				replacedNodes.push_back(element->nodes()[(localId+2)%3]->id()-1);
#ifdef DEBUG_MESH_SIMPLIFICATION
				cerr << collapseNode->id() << " adding " << element->id() << endl;
#endif
				collapseNode->element(*element);
			}
		}
	}

	newNeighbors.clear();
	//clean up duplicates from newNeighbors:
	for(i = 0; i < (int)replacedNodes.size(); i++)
	{
		nodeInThere = false;
		for(j = 0; j < (int)newNeighbors.size(); j++)
		{
			if(replacedNodes[i] == newNeighbors[j])
			{
				nodeInThere = true;
				break;
			}
		}
		if(nodeInThere) continue;
		for(j = 0; j < (int)commonNodes.size(); j++)
		{
			if(replacedNodes[i] == commonNodes[j])
			{
				nodeInThere = true;
				break;
			}
		}
		if(!nodeInThere) newNeighbors.push_back(replacedNodes[i]);
	}
#ifdef DEBUG_MESH_SIMPLIFICATION
	collapseNode->print(cerr);
#endif
}

void MeshSimplification::storeMeanValueWeights(int index, vector<treeNode> & structure)
{
	//find the neighbors of the node index:
	int i, j, localId;
	double sumOmegaDash, c, d, length, angle, omegaDash;
	g_Vector normal, l, projNode, proj1, proj2;
	g_Node * node;
	g_PEdge edge;
	g_NodeContainer neighbors;
	g_ElementContainer neighborsOfEdge;
	g_PtrAdapter<g_Node*> adapter(neighbors);
	vector<double> omegaDashes;

	//clear the information in the treeNode:
	structure[index].neighborIndices.clear();
	structure[index].elementIndices.clear();
	structure[index].omegas.clear();
	structure[index].bs.clear();

	node = workingMesh->nodes()[index];	
	for(i = 0; i < (int)node->elements().numberOfItems(); i++)
	{
		if(node->elements()[i]->nodes()[0]->id() != node->id())
		{
			neighbors.insert(node->elements()[i]->nodes()[0]);
			structure[index].elementIndices.push_back(node->elements()[i]->nodes()[0]->id());
		}
		if(node->elements()[i]->nodes()[1]->id() != node->id())
		{
			neighbors.insert(node->elements()[i]->nodes()[1]);
			structure[index].elementIndices.push_back(node->elements()[i]->nodes()[1]->id());
		}
		if(node->elements()[i]->nodes()[2]->id() != node->id())
		{
			neighbors.insert(node->elements()[i]->nodes()[2]);
			structure[index].elementIndices.push_back(node->elements()[i]->nodes()[2]->id());
		}
	}
	adapter.removeDuplicates();

	for(i = 0; i < (int)neighbors.numberOfItems(); i++)
		structure[index].neighborIndices.push_back(neighbors[i]->id());

	//find tangent plane and project everything to it:
	l.Set(0, 0, 0);
	if((int)neighbors.numberOfItems() > 0)
	{
		for(i = 0; i < (int)neighbors.numberOfItems(); i++)
			l = l + neighbors[i]->coordinate();
		l = l * (1.0/(double)neighbors.numberOfItems());
	}
	normal.Set(0, 0, 0);
	for(i = 0; i < (int)node->elements().numberOfItems(); i++)
	{
		normal = normal + (workingMesh->nodes()[structure[index].elementIndices[2*i]-1]->coordinate()-l).Cross(
			workingMesh->nodes()[structure[index].elementIndices[2*i+1]-1]->coordinate()-l);
	}
	normal.Normalize();

	d = 0.0;
	if((int)neighbors.numberOfItems() > 0)
	{
		for(i = 0; i < (int)neighbors.numberOfItems(); i++)
			d += normal.Dot(neighbors[i]->coordinate());
		d = -d/(double)neighbors.numberOfItems();
	}
	projNode = node->coordinate()-(d+node->coordinate().Dot(normal))*normal;
	//compute weights and store them:
	omegaDashes.clear();

	for(i = 0; i < (int)neighbors.numberOfItems(); i++)
	{
		edge.nodes(*node, *(neighbors[i]));
		neighborsOfEdge = edge.elements();
		proj1 = neighbors[i]->coordinate()-(d+neighbors[i]->coordinate().Dot(normal))*normal;
		length = projNode.DistanceTo(proj1);
		omegaDash = 0;
		for(j = 0; j < (int)neighborsOfEdge.numberOfItems(); j++)
		{
			if((neighborsOfEdge[j]->nodes()[0]->id() != node->id()) && (neighborsOfEdge[j]->nodes()[0]->id() != neighbors[i]->id())) 
				localId = 0;
			if((neighborsOfEdge[j]->nodes()[1]->id() != node->id()) && (neighborsOfEdge[j]->nodes()[1]->id() != neighbors[i]->id()))
				localId = 1;
			if((neighborsOfEdge[j]->nodes()[2]->id() != node->id()) && (neighborsOfEdge[j]->nodes()[2]->id() != neighbors[i]->id()))
				localId = 2;

			proj2 = neighborsOfEdge[j]->nodes()[localId]->coordinate()-(d+neighborsOfEdge[j]->nodes()[localId]->coordinate().Dot(normal))*normal;
			angle = (proj1-projNode).SquaredLength() + (proj2-projNode).SquaredLength() - (proj1-proj2).SquaredLength();
			angle = angle / max((2.0 * (proj1-projNode).Length() * (proj2-projNode).Length()), MDS_THRESHOLD_DIVIDE0);
			if(angle < -1.0) angle = -1.0;
			if(angle > 1.0) angle = 1.0;
			angle = acos(angle);
			angle = angle/2.0;
			omegaDash += tan(angle);
		}
		omegaDash = omegaDash/max(length, MDS_THRESHOLD_DIVIDE0);		
		omegaDashes.push_back(omegaDash);
	}
	sumOmegaDash = 0;
	for(i = 0; i < (int)neighbors.numberOfItems(); i++)
		sumOmegaDash += omegaDashes[i];
	structure[index].omegas.clear();
	for(i = 0; i < (int)neighbors.numberOfItems(); i++)
		structure[index].omegas.push_back(omegaDashes[i]/max(sumOmegaDash, MDS_THRESHOLD_DIVIDE0));

	for(i = 0; i < (int)neighbors.numberOfItems(); i++)
	{
		c = (node->coordinate()-neighbors[i]->coordinate()).Dot(normal) / max((node->coordinate().DistanceTo(neighbors[i]->coordinate())), MDS_THRESHOLD_DIVIDE0);
		structure[index].bs.push_back(c/max(sqrt(1-pow(c, 2)), MDS_THRESHOLD_DIVIDE0));
	}
}

g_Vector MeshSimplification::getMeanValueCoordinate(g_Part * deformedMesh, int index)
{
	g_Vector result;

	g_Part * meshToStore = workingMesh;
	workingMesh = deformedMesh;

	result = computeMeanValueVertexPos(index, originalMeshWeights);

	workingMesh = meshToStore;

	return result;
}

g_Vector MeshSimplification::computeMeanValueVertexPos(int index, vector<treeNode> & structure)
{
	int i, j;
	g_Vector l, normal, vectorSum, result;
	double * matrixN = new double[9];

	l.Set(0, 0, 0);
	for(i = 0; i < (int)structure[index].neighborIndices.size(); i++)		
		l = l + workingMesh->nodes()[structure[index].neighborIndices[i]-1]->coordinate();
	l = l / (int)structure[index].neighborIndices.size();
	normal.Set(0, 0, 0);
	for(i = 0; i < (int)structure[index].elementIndices.size()/2; i++)
		normal = normal + (workingMesh->nodes()[structure[index].elementIndices[2*i]-1]->coordinate()-l).Cross(
			workingMesh->nodes()[structure[index].elementIndices[2*i+1]-1]->coordinate()-l);		
	normal.Normalize();

	matrixN[0] = 1-normal.x()*normal.x(); matrixN[1] = -normal.y()*normal.x(); matrixN[2] = -normal.z()*normal.x();
	matrixN[3] = -normal.y()*normal.x(); matrixN[4] = 1-normal.y()*normal.y(); matrixN[5] = -normal.z()*normal.y();
	matrixN[6] = -normal.z()*normal.x(); matrixN[7] = -normal.y()*normal.z(); matrixN[8] = 1-normal.z()*normal.z();

	result.Set(0, 0, 0);
	for(i = 0; i < (int)structure[index].neighborIndices.size(); i++)
	{
		vectorSum.Set(0, 0, 0);
		for(j = 0; j < (int)structure[index].neighborIndices.size(); j++)
			vectorSum = vectorSum + structure[index].omegas[j]*
				(workingMesh->nodes()[structure[index].neighborIndices[j]-1]->coordinate()-
				workingMesh->nodes()[structure[index].neighborIndices[i]-1]->coordinate());
		result = result + 
			structure[index].omegas[i]* (workingMesh->nodes()[structure[index].neighborIndices[i]-1]->coordinate() +
			sqrt(pow(matrixN[0]*vectorSum.x() + matrixN[1]*vectorSum.y() + matrixN[2]*vectorSum.z(), 2) + 
			pow(matrixN[3]*vectorSum.x() + matrixN[4]*vectorSum.y() + matrixN[5]*vectorSum.z(), 2) + 
			pow(matrixN[6]*vectorSum.x() + matrixN[7]*vectorSum.y() + matrixN[8]*vectorSum.z(), 2)) * 
			structure[index].bs[i]*normal);
	}

	delete [] matrixN;

	return result;
}

bool MeshSimplification::minimizeEnergyKilian(vector<int> & indicesToKeep)
{
	int i, numberNodes, dimension, count1, count2;
	
	numberNodes = (int)workingMesh->nodes().numberOfItems() - (int)indicesToKeep.size();
	dimension = numberNodes * 3;
	
	int * indexToKeep = (int *) malloc((int)workingMesh->nodes().numberOfItems() * sizeof(int));
	vnl_vector<double> x(dimension, 0.0);

	//intialize the variables
	count1 = 0;
	count2 = 0;
	for(i = 0; i < (int)workingMesh->nodes().numberOfItems(); i++)	
	{
		if(count1 < (int)indicesToKeep.size())
		{
			if(indicesToKeep[count1] == i) 
			{
				count1++;
				indexToKeep[i] = -1;
			}
			else indexToKeep[i] = count2;
		}
		else indexToKeep[i] = count2;
		if(indexToKeep[i] != -1)
		{
			x[3*count2] = workingMesh->nodes()[i]->coordinate().x();
			x[3*count2+1] = workingMesh->nodes()[i]->coordinate().y();
			x[3*count2+2] = workingMesh->nodes()[i]->coordinate().z();	
			count2++;
		}
	}

	//perform the nonlinear optimization
	KilianCostFunction costFct(this, indexToKeep, dimension);
	vnl_lbfgsb * minimizer = new vnl_lbfgsb(costFct);
	minimizer->set_cost_function_convergence_factor ( 100000000 );
	minimizer->set_projected_gradient_tolerance ( 0.0001 ); 
	minimizer->set_max_function_evals(1000);
	minimizer->minimize( x );

	
	//compute and copy the result:
	for(i = 0; i < (int)workingMesh->nodes().numberOfItems(); i++)	
	{
		if(indexToKeep[i] != -1) 
			workingMesh->nodes()[i]->coordinate(g_Vector(x[3*indexToKeep[i]], x[3*indexToKeep[i]+1], x[3*indexToKeep[i]+2]));
	}

	//delete stuff:
	delete minimizer;
	free(indexToKeep);
	
	return true;
}

void MeshSimplification::KilianCostFunction::compute(vnl_vector< double > const &x, double *f, vnl_vector< double > *g)
{
	int i, dimension, id1, id2;
	g_Vector coo1, coo2;
	double wantedLength, actualLength, commonFactor;
	dimension = x.size();

	g_PEdgeContainer edges = simplify->workingMesh->uniqueEdges();

	// evaluate stress function at x and gradient vector of stress function at x:
	*f = 0;
	for(i = 0; i < dimension; i++) (*g)[i] = 0.0;
	for(i = 0; i < (int)edges.numberOfItems(); i++) 
	{
		id1 = edges[i]->firstNode().id()-1;
		id2 = edges[i]->lastNode().id()-1;
		wantedLength = simplify->originalMesh->nodes()[id1]->coordinate().DistanceTo(simplify->originalMesh->nodes()[id2]->coordinate());
		if(indexToKeep[id1] == -1) coo1 = simplify->workingMesh->nodes()[id1]->coordinate();
		else coo1 = g_Vector(x[3*indexToKeep[id1]], x[3*indexToKeep[id1]+1], x[3*indexToKeep[id1]+2]);
		if(indexToKeep[id2] == -1) coo2 = simplify->workingMesh->nodes()[id2]->coordinate();
		else coo2 = g_Vector(x[3*indexToKeep[id2]], x[3*indexToKeep[id2]+1], x[3*indexToKeep[id2]+2]);
		actualLength = coo1.DistanceTo(coo2);
		*f += pow(actualLength-wantedLength, 2);
		commonFactor = 2.0*(actualLength-wantedLength)/(max(GRAD_TOL, actualLength));
		if(indexToKeep[id1] != -1) 
		{
			(*g)[3*indexToKeep[id1]] += commonFactor*(coo1-coo2).x();
			(*g)[3*indexToKeep[id1]+1] += commonFactor*(coo1-coo2).y();
			(*g)[3*indexToKeep[id1]+2] += commonFactor*(coo1-coo2).z();
		}
		if(indexToKeep[id2] != -1) 
		{
			(*g)[3*indexToKeep[id2]] += commonFactor*(coo2-coo1).x();
			(*g)[3*indexToKeep[id2]+1] += commonFactor*(coo2-coo1).y();
			(*g)[3*indexToKeep[id2]+2] += commonFactor*(coo2-coo1).z();
		}
	}

	for(i = 0; i < (int)edges.numberOfItems(); i++)
		delete edges[i];
	edges.clear();
}

bool MeshSimplification::minimizeEnergy(vector<int> & indicesToKeep)
{
	//initialize the paramteters for the nonlinear optimization
	int i, numberNodes, dimension, count1, count2;
	numberNodes = (int)workingMesh->nodes().numberOfItems() - (int)indicesToKeep.size();
	dimension = numberNodes * 3;

	vnl_vector<double> x( dimension, 0.0 );

	perturbConst = max(0.5 * resolution, MDS_THRESHOLD_DIVIDE0);

	int * indexToKeep = (int *) malloc((int)workingMesh->nodes().numberOfItems() * sizeof(int));
	
	//intialize the variables
	count1 = 0;
	count2 = 0;
	for(i = 0; i < (int)workingMesh->nodes().numberOfItems(); i++)	
	{
		if(count1 < (int)indicesToKeep.size())
		{
			if(indicesToKeep[count1] == i) 
			{
				count1++;
				indexToKeep[i] = -1;
			}
			else indexToKeep[i] = count2;
		}
		else indexToKeep[i] = count2;
		if(indexToKeep[i] != -1)
		{
			x[3*count2] = workingMesh->nodes()[i]->coordinate().x();
			x[3*count2+1] = workingMesh->nodes()[i]->coordinate().y();
			x[3*count2+2] = workingMesh->nodes()[i]->coordinate().z();	
			count2++;
		}
	}

	//perform the nonlinear optimization
	MeanValueCostFunction costFct(this, indexToKeep, dimension);
	vnl_lbfgsb * minimizer = new vnl_lbfgsb(costFct);
	minimizer->set_cost_function_convergence_factor ( 100000000 );
	minimizer->set_projected_gradient_tolerance ( 0.0001 ); 
	minimizer->set_max_function_evals(1000);
	minimizer->minimize( x );

	//copy the result:
	for(i = 0; i < (int)workingMesh->nodes().numberOfItems(); i++)	
	{
		if(indexToKeep[i] != -1) 
			workingMesh->nodes()[i]->coordinate(g_Vector(x[3*indexToKeep[i]], x[3*indexToKeep[i]+1], x[3*indexToKeep[i]+2]));
	}

	//delete stuff:
	delete minimizer;
	free(indexToKeep);

	return true;
}

void MeshSimplification::MeanValueCostFunction::compute(vnl_vector< double > const &x, double *f, vnl_vector< double > *g)
{
	int i, dimension, id1, id2;
	g_Vector coo;
	dimension = x.size();

	g_PEdgeContainer edges = simplify->workingMesh->uniqueEdges();

	// evaluate stress function at x and gradient vector of stress function at x:
	*f = 0;
	for(i = 0; i < dimension; i++) (*g)[i] = 0.0;
	//copy the result:
	for(i = 0; i < (int)simplify->workingMesh->nodes().numberOfItems(); i++)	
	{
		if(indexToKeep[i] != -1) 
			simplify->workingMesh->nodes()[i]->coordinate(g_Vector(x[3*indexToKeep[i]], x[3*indexToKeep[i]+1], x[3*indexToKeep[i]+2]));
	}

	for(i = 0; i < (int)simplify->workingMesh->nodes().numberOfItems(); i++) 
	{
		if(indexToKeep[i] != -1) 
		{
			coo = simplify->computeMeanValueVertexPos(i, simplify->originalMeshWeights);
			*f += 0.5*(simplify->workingMesh->nodes()[i]->coordinate()-coo).SquaredLength();
			(*g)[3*indexToKeep[i]] = (simplify->workingMesh->nodes()[i]->coordinate()-coo).x();
			(*g)[3*indexToKeep[i]+1] = (simplify->workingMesh->nodes()[i]->coordinate()-coo).y();
			(*g)[3*indexToKeep[i]+2] = (simplify->workingMesh->nodes()[i]->coordinate()-coo).z();
		}
	}
	//NOTE that the numerical derivative is a MATRIX
	for(i = 0; i < (int)edges.numberOfItems(); i++) 
	{
		id1 = edges[i]->firstNode().id()-1;
		id2 = edges[i]->lastNode().id()-1;
		if((indexToKeep[id1] != -1) && (indexToKeep[id2] != -1)) 
		{
			(*g)[3*indexToKeep[id1]] += simplify->getPartialDerivative(id1, id2, 0);
			(*g)[3*indexToKeep[id1]+1] += simplify->getPartialDerivative(id1, id2, 1);
			(*g)[3*indexToKeep[id1]+2] += simplify->getPartialDerivative(id1, id2, 2);

			(*g)[3*indexToKeep[id2]] += simplify->getPartialDerivative(id2, id1, 0);
			(*g)[3*indexToKeep[id2]+1] += simplify->getPartialDerivative(id2, id1, 1);
			(*g)[3*indexToKeep[id2]+2] += simplify->getPartialDerivative(id2, id1, 2);
		}
	}

	for(i = 0; i < (int)edges.numberOfItems(); i++)
		delete edges[i];
	edges.clear();
}

double MeshSimplification::getPartialDerivative(int id1, int id2, int index)
{
	g_Vector coo1, coo2, coo3, perturbation;

	if(index == 0) perturbation.Set(perturbConst, 0, 0);
	else if(index == 1) perturbation.Set(0, perturbConst, 0);
	else perturbation.Set(0, 0, perturbConst);

	coo3 = workingMesh->nodes()[id2]->coordinate()-computeMeanValueVertexPos(id2, originalMeshWeights);			
	//perturb the vertex at id1 slightly:
	workingMesh->nodes()[id1]->coordinate(workingMesh->nodes()[id1]->coordinate() + perturbation);
	coo1 = computeMeanValueVertexPos(id2, originalMeshWeights);
	workingMesh->nodes()[id1]->coordinate(workingMesh->nodes()[id1]->coordinate() - 2*perturbation);
	coo2 = computeMeanValueVertexPos(id2, originalMeshWeights);
	coo1 = (1.0/(2.0*perturbConst))*(coo1-coo2);
	workingMesh->nodes()[id1]->coordinate(workingMesh->nodes()[id1]->coordinate() + perturbation);

	return -(coo3.Dot(coo1));
}

bool MeshSimplification::minimizeEnergy(g_Part *& newMesh, vector<int> & indicesToKeep)
{
	int i;

	if(workingMesh->nodes().numberOfItems() != newMesh->nodes().numberOfItems()) return false;
	//After this check, we assume the mesh topology is OK:
	for(i = 0; i < (int)workingMesh->nodes().numberOfItems(); i++) workingMesh->nodes()[i]->coordinate(newMesh->nodes()[i]->coordinate());

	minimizeEnergy(indicesToKeep);
	for(i = 0; i < (int)workingMesh->nodes().numberOfItems(); i++) newMesh->nodes()[i]->coordinate(workingMesh->nodes()[i]->coordinate());

	return true;
}

g_Part * MeshSimplification::reconstructMesh(g_NodeContainer cont, vector<int> & indicesToKeep,
											 bool minimizeKilianEnergy)
{
	int i, j, counter, childIdx;
	g_Node * node;
	g_Element * element;
	bool decoding, ready;

	//initialize the mesh:
	if(workingMesh != NULL)
		delete workingMesh;
	workingMesh = new g_Part;
	counter = 0;
	for(i = 0; i < (int)originalMesh->nodes().numberOfItems(); i++)
	{
		if(counter < (int)indicesToKeep.size())
		{
			if(indicesToKeep[counter] == i)
			{
				node = new g_Node(cont[counter]->coordinate());
				counter++;
			}
			else node = new g_Node(originalMesh->nodes()[i]->coordinate());
		}
		else node = new g_Node();
		workingMesh->node(node);
	}
	for(i = 0; i < (int)originalMesh->elements().numberOfItems(); i++)
	{
		element = new g_Element();
		element->node(workingMesh->nodes()[originalMesh->elements()[i]->nodes()[0]->id()-1]);
		element->node(workingMesh->nodes()[originalMesh->elements()[i]->nodes()[1]->id()-1]);
		element->node(workingMesh->nodes()[originalMesh->elements()[i]->nodes()[2]->id()-1]);
		workingMesh->element(element);
	}

	//recompute the coordinates one by one:
	do
	{
		decoding = false;
		for(i = 0; i < (int)treeStructure.size(); i++)
		{
			if((treeStructure[i].parent == treeStructure[i].id) && (treeStructure[i].child2.size() > 0))
			{
				//test if the child is ready to be decoded: all of its neighbors need to be roots
				childIdx = treeStructure[i].child2.back();
				ready = true;
				for(j = 0; j < (int)treeStructure[childIdx-1].neighborIndices.size(); j++)
				{
					if(treeStructure[treeStructure[childIdx-1].neighborIndices[j]-1].parent !=
						treeStructure[treeStructure[childIdx-1].neighborIndices[j]-1].id)				
					{
						ready = false;
						break;
					}
				}
				if(ready)
				{
					//set the node coordinate to the result:
					workingMesh->nodes()[childIdx-1]->coordinate(computeMeanValueVertexPos(childIdx-1, treeStructure));
					treeStructure[childIdx-1].parent = treeStructure[childIdx-1].id;
					treeStructure[i].child2.pop_back();
					decoding = true;
				}
			}
		}
	}while(decoding);
	if(minimizeKilianEnergy)
	{
		if(!minimizeEnergy(indicesToKeep))cout<<"Problem when minimizing energy during recontruction"<<endl;
		if(!minimizeEnergyKilian(indicesToKeep))cout<<"Problem when minimizing energy during recontruction"<<endl;
	}

	//return the resulting mesh:
	return workingMesh;
}

g_Part * MeshSimplification::reconstructMesh(g_NodeContainer cont, vector<int> & indicesToKeep, vector<int> & indicesForMinimization, bool minKilian)
{
	reconstructMesh(cont, indicesToKeep, false);
	if(!minimizeEnergy(indicesToKeep))cout<<"Problem when minimizing energy during recontruction"<<endl;
	if(minKilian) 
	{
		if(!minimizeEnergyKilian(indicesForMinimization))cout<<"Problem when minimizing energy during recontruction"<<endl;
	}

	//return the resulting mesh:
	return workingMesh;
}

void MeshSimplification::createFinalMesh(vector<int> & indicesToKeep, bool allDeleted, bool keepNonCollapsables)
{
	int i, j, number, numberDelVertices;
	bool found;
	g_Node * node;
	g_Element * element;
	int * helpVec;

	additionalIndicesKept.clear();
	allIndicesKept.clear();
	finalMesh = new g_Part;
	number = (int)workingMesh->nodes().numberOfItems();
	numberDelVertices = 0;
	//delete nodes that are not kept and that were not deleted before due to only using collapse:
	if(!allDeleted)
	{
		for(i = 0; i < number; i++)
		{
			found = false;
			if((int)workingMesh->nodes()[i]->elements().numberOfItems() > 0)
			{
				for(j = 0; j < (int)indicesToKeep.size(); j++) 
				{
					if(indicesToKeep[j] == workingMesh->nodes()[i]->id()-1)
					{
						found = true; 
						allIndicesKept.push_back(workingMesh->nodes()[i]->id()-1);
						break;
					}
				}
				if(!found)
				{
					if(!keepNonCollapsables)
					{
						numberDelVertices++;
						while((int)workingMesh->nodes()[i]->elements().numberOfItems()>0)
						{
							element = workingMesh->nodes()[i]->elements()[0];
							element->nodes()[0]->removeElement(element);
							element->nodes()[1]->removeElement(element);
							element->nodes()[2]->removeElement(element);
							element->emptyNodeList();
						}
					}
					else 
					{
						additionalIndicesKept.push_back(workingMesh->nodes()[i]->id()-1);
						allIndicesKept.push_back(workingMesh->nodes()[i]->id()-1);
					}
				}
			}
		}
		if(!keepNonCollapsables) cout<<"WE NEED TO DELETE VERTICES "<<numberDelVertices<<endl;
	}
	//Add valid nodes and elements to the mesh:
	helpVec = new int[number];
#pragma omp parallel for private(i)
	for(i = 0; i < number; i++) helpVec[i] = -1;
	for(i = 0; i < number; i++)
	{
		found = false;
		for(j = 0; j < (int)indicesToKeep.size(); j++) 
		{
			if(indicesToKeep[j] == workingMesh->nodes()[i]->id()-1)
			{
				found = true; 
				break;
			}
		}
		if(keepNonCollapsables)
		{
			for(j = 0; j < (int)additionalIndicesKept.size(); j++)
			{
				if(additionalIndicesKept[j] == workingMesh->nodes()[i]->id()-1)
				{
					found = true;
					break;
				}
			}
		}
		if(found)
		{
			node = new g_Node(workingMesh->nodes()[i]->coordinate());
			finalMesh->node(node);
			helpVec[i] = node->id()-1;
		}
	}
	number = workingMesh->elements().numberOfItems();
	for(i = 0; i < number; i++)
	{
		if(workingMesh->elements()[i]->nodes().numberOfItems() == 3)
		{
			if((helpVec[workingMesh->elements()[i]->nodes()[0]->id()-1] != -1) &&
				(helpVec[workingMesh->elements()[i]->nodes()[1]->id()-1] != -1) &&
				(helpVec[workingMesh->elements()[i]->nodes()[2]->id()-1] != -1))
			{
				element = new g_Element();
				element->node(finalMesh->nodes()[helpVec[workingMesh->elements()[i]->nodes()[0]->id()-1]]);
				element->node(finalMesh->nodes()[helpVec[workingMesh->elements()[i]->nodes()[1]->id()-1]]);
				element->node(finalMesh->nodes()[helpVec[workingMesh->elements()[i]->nodes()[2]->id()-1]]);
				finalMesh->element(element);
			}
		}
	}
	delete [] helpVec;
	//Delete the workingMesh:
	delete workingMesh;
	workingMesh = NULL;
}
