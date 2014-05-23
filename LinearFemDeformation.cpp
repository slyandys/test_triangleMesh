#include "LinearFemDeformation.h"

LinearFemDeformation::LinearFemDeformation()
{
	initialize();
}

LinearFemDeformation::LinearFemDeformation(TriangleMesh * inMesh, TriangleMesh * tempMesh, char * exportTet, int mode, g_Node * contactPoint)
{
	initialize();

	//Copy the surface mesh:
	mesh = new TriangleMesh(*inMesh);

	//Simplify the mesh: 	
	if(tempMesh == NULL)
	{
		lowResMesh = tempMesh;
		if(intersectionTetgen(inMesh, lowResMesh, mode))
		{
			cout<<"Problem simplifying mesh"<<endl;
			return;
		}
	}
	else lowResMesh = new TriangleMesh(*tempMesh);

	if(exportTet != NULL)
	{
		char buffer[1024];
		sprintf(buffer, "%s.wrl", exportTet);
		exportMeshWrapper(buffer, lowResMesh );
	}

	//Compute tetrahedralization of the mesh
	int i;

	//Compute tetrahedralization of the mesh
	tetgenio * in, * out;
	in = NULL;
	createTetgenMesh(in, lowResMesh);
	out = new tetgenio();

	// Tetrahedralize the PLC. Switches are chosen to read a PLC (p)
	tetgenbehavior b;
	if(exportTet != NULL)
	{
		b.parse_commandline(ARG_O); 
		for(i = 0; i < 1024; i++) 
		{
			b.outfilename[i] = exportTet[i];
			if(exportTet[i] == '\0') break;
		}
	}
	else b.parse_commandline(ARG); 
	tetrahedralize(&b, in, out);

	//Copy the information to the class variables:
	numNodes = out->numberofpoints;
	numTetrahedra = out->numberoftetrahedra;

	//Store the tetrahedral mesh: nodes
	allNodes.clear();
	for(i = 0; i < numNodes; i++)
	{
		g_Node * newNode = new g_Node(g_Vector(out->pointlist[3*i], out->pointlist[3*i+1], out->pointlist[3*i+2])); 
		allNodes.insert(newNode);
	}

	//Tetrahedra:
	allTetrahedra.clear();
	for(i = 0; i < numTetrahedra; i++)
	{
		vector<int> tetrahedron;
		tetrahedron.push_back(out->tetrahedronlist[4*i]-1);
		tetrahedron.push_back(out->tetrahedronlist[4*i+1]-1);
		tetrahedron.push_back(out->tetrahedronlist[4*i+2]-1);
		tetrahedron.push_back(out->tetrahedronlist[4*i+3]-1);
		allTetrahedra.push_back(tetrahedron);
	}

	//Compute and store the matching information:
	CorrespondTemplates * corres = new CorrespondTemplates();
	corres->setAlignedSurface(inMesh);
	corres->setTetModel(allNodes, allTetrahedra);
	corres->computeMappings(NULL,NULL, contactPoint);

	mappingIndices = corres->getMappingIndicesTet();
	mappingWeights = corres->getMappingWeightsTet();
	mappingOffsets = corres->getMappingOffsetsTet();
	invMappingIndices = corres->getMappingIndicesSurf();
	invMappingWeights = corres->getMappingWeightsSurf();
	invMappingOffsets = corres->getMappingOffsetsSurf();
	contactNodes = corres->getContactPoints();
	delete corres;

	delete in;
	delete out;
}

void LinearFemDeformation::initialize()
{
	matrixK = NULL;
	offsets = NULL;
	forces = NULL;
	mesh = NULL;
	lowResMesh = NULL;
	contactPointAndForce = NULL;
	numContactPoints = 0;
	allNodes.clear();
	allTetrahedra.clear();
	mappingIndices.clear();
	updatedPositions.clear();
	validNodes.clear();

	youngsModulus = 1000000;
	poissonRatio = 0.45;

	lengthRBFOffsets.clear();
}

LinearFemDeformation::~LinearFemDeformation()
{
	if(matrixK != NULL)
	{
		delete [] matrixK;
		matrixK = NULL;
	}

	if(offsets != NULL)
	{
		delete [] offsets;
		offsets = NULL;
	}

	if(forces != NULL)
	{
		delete [] forces;
		forces = NULL;
	}

	if(contactPointAndForce != NULL)
	{
		delete [] contactPointAndForce;
		contactPointAndForce = NULL;
	}

	if(mesh != NULL) 
	{
		delete mesh;
		mesh = NULL;
	}

	if(lowResMesh != NULL) 
	{
		delete lowResMesh;
		lowResMesh = NULL;
	}

	int i;
	for(i = 0; i < numNodes; i++)
		delete allNodes[i];
	allNodes.clear();

	for(i = 0; i < (int)updatedPositions.numberOfItems(); i++)
		delete updatedPositions[i];
	updatedPositions.clear();

}

void LinearFemDeformation::setParameters(double youngsModulus, double poissonRatio)
{
	this->youngsModulus = youngsModulus;
	this->poissonRatio = poissonRatio;

	computeMatrixK();
}

void LinearFemDeformation::setContactPtAndForce(double * contactPointAndForceIn, int numContactPoints)
{
	int i;

	this->numContactPoints = numContactPoints;
	if(contactPointAndForce != NULL) delete [] contactPointAndForce;
	contactPointAndForce = new double[numContactPoints * 4];
	for(i = 0; i < numContactPoints; i++)
	{
		contactPointAndForce[4*i + 0] = (double)findClosestIndex(g_Vector(contactPointAndForceIn[6*i + 0], contactPointAndForceIn[6*i + 1], contactPointAndForceIn[6*i + 2]), allNodes)[0];
		contactPointAndForce[4*i + 1] = contactPointAndForceIn[6*i + 3];
		contactPointAndForce[4*i + 2] = contactPointAndForceIn[6*i + 4];
		contactPointAndForce[4*i + 3] = contactPointAndForceIn[6*i + 5];
		
		cout<<"Original Force: "<<g_Vector(contactPointAndForce[4*i + 1],contactPointAndForce[4*i + 2],contactPointAndForce[4*i + 3])<<endl;

#ifdef _FORCE_DIRECTION_NORMAL_TO_TRI
		cerr << ">>>>> Correcting force " << contactPointAndForce[4*i + 0] << endl;

		g_Vector frameF = g_Vector(contactPointAndForceIn[6*i + 3], contactPointAndForceIn[6*i + 4], contactPointAndForceIn[6*i + 5]);

		g_Vector a = g_Vector(allNodes[contactNodes[1]]->coordinate() - allNodes[contactNodes[0]]->coordinate());
		g_Vector b = g_Vector(allNodes[contactNodes[2]]->coordinate() - allNodes[contactNodes[0]]->coordinate());
		g_Vector dir  =  a.Cross(b);
		if(dir.Dot(frameF) < 0) //if the direction have more than 90 deg angle with original force vector
			dir *= -1;
			  
		double mag = frameF.length();
		dir.normalize();
		dir = mag * dir;

		contactPointAndForce[4*i + 1] = dir.x();
		contactPointAndForce[4*i + 2] = dir.y();
		contactPointAndForce[4*i + 3] = dir.z();
		
		cout<<"Corrected Force: "<<g_Vector(contactPointAndForce[4*i + 1],contactPointAndForce[4*i + 2],contactPointAndForce[4*i + 3])<<endl;
#endif
	}
}

void LinearFemDeformation::setNewPositions(g_NodeContainer newPositions, vector<int> validNodes)
{
	int i;

	//Copy the new positions and the valid node indices:
	for(i = 0; i < (int)updatedPositions.numberOfItems(); i++)
		delete updatedPositions[i];
	updatedPositions.clear();

	for(i = 0; i < (int)newPositions.numberOfItems(); i++)
	{
		g_Node * node = new g_Node(*(newPositions[i]));
		updatedPositions.insert(node);
	}

	this->validNodes = validNodes;
}

bool LinearFemDeformation::intersectionTetgen(TriangleMesh * inMesh, TriangleMesh *& tempMesh, int mode)
{
	bool returnVal;

	//Simplify the mesh: 

	if(mode == 0)
	{
		//Case 1: Use MeshSimplification class (explicitly disallow self-intersections)
		if(tempMesh != NULL) delete tempMesh;

		MeshSimplification simplify(inMesh);
//		simplify.useVolumeCost();
		tempMesh = new TriangleMesh(*(simplify.getSimplifiedMesh(SIMPLIFIED_NUM, true)));
		exportMeshWrapper( "lowRes.wrl", tempMesh ); 
	}

	else if(mode == 1)
	{
		//Case 2 (faster): Use MeshLab
		if(tempMesh != NULL) delete tempMesh;

		exportMeshWrapper("highRes.wrl", inMesh );

		system("\"C:\\Program Files (x86)\\VCG\\MeshLab\\meshlabserver\" -i highRes.wrl -o lowRes.wrl -s meshlab_compress_script.mlx");

		tempMesh = new TriangleMesh();
		importMeshWrapper("lowRes.wrl", tempMesh);
	}

	//Case 3 (even faster) assumes the mesh is given.

	//Compute tetrahedralization of the mesh
	tetgenio * in, * out;
	in = NULL;
	createTetgenMesh(in, tempMesh);
	out = new tetgenio();

	// Check for intersections
	tetgenbehavior b;
	b.parse_commandline("dQ");
	tetrahedralize(&b, in, out);
	
	//In -d mode, points are only reported if there are self-intersections
	if(out->numberofpoints > 0) returnVal = true;
	else
	{
		//Reset the meshes (just in case tetgen changes anything here)
		createTetgenMesh(in, tempMesh);
		delete out;
		out = new tetgenio();

		tetgenbehavior bOther;
		bOther.parse_commandline(ARG_Q);
		tetrahedralize(&bOther, in, out);

		if(out->numberofpoints > MAX_NUM_VERT_IN_TET) returnVal = true;
		else returnVal = false;
	}

	delete in;
	delete out;

	return returnVal;
}

void LinearFemDeformation::createTetgenMesh(tetgenio *& in, TriangleMesh * tempMesh)
{
	int i, numVertices, numTriangles;
	tetgenio::facet *f;
	tetgenio::polygon *p;

	if(in != NULL) delete in;
	in = new tetgenio();

	// All indices start from 1.
	in->firstnumber = 1;

	numVertices = tempMesh->getNumberOfNodes();
	numTriangles = tempMesh->getNumberOfTriangles();

	// Vertices
	in->numberofpoints = numVertices;
	in->pointlist = new REAL[numVertices * 3];
	for(i = 0; i < numVertices; i++)
	{
		in->pointlist[3*i] = tempMesh->nodes()[i]->coordinate().x(); 
		in->pointlist[3*i+1] = tempMesh->nodes()[i]->coordinate().y(); 
		in->pointlist[3*i+2] = tempMesh->nodes()[i]->coordinate().z(); 
	}

	// Triangles
	in->numberoffacets = numTriangles;
	in->facetlist = new tetgenio::facet[numTriangles];
	in->facetmarkerlist = new int[numTriangles];
	for(i = 0; i < numTriangles; i++)
	{
		f = &in->facetlist[i];
		f->numberofpolygons = 1;
		f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
		f->numberofholes = 0;
		f->holelist = NULL;
		p = &f->polygonlist[0];
		p->numberofvertices = 3;
		p->vertexlist = new int[p->numberofvertices];
		p->vertexlist[0] = tempMesh->elements()[i]->nodes()[0]->id();
		p->vertexlist[1] = tempMesh->elements()[i]->nodes()[1]->id();
		p->vertexlist[2] = tempMesh->elements()[i]->nodes()[2]->id();
		in->facetmarkerlist[i] = 0;
	}
}

vector<int> LinearFemDeformation::getPredictedIndices()
{
	vector<int> returnIndices;

	for(int i = 0; i < (int)updatedPositions.numberOfItems(); i++)
	{
		if(validNodes[i] == 2) 
			returnIndices.push_back(i);
	}

	return returnIndices;
}

bool LinearFemDeformation::computeSimulation(bool updateSurfaceNodes, bool useNormalOffsets, double * forceVector)
{
	//Check that things are initialized properly:
	if(contactPointAndForce == NULL) return false;
	if(updatedPositions.numberOfItems() == 0) return false;
	if(validNodes.size() == 0) return false;
	if(youngsModulus < 0) return false;
	if((poissonRatio < 0) || (poissonRatio >= 0.5)) return false;

	//Compute the offsets using the mapping
	double * modifiedK = new double[3*numNodes*3*numNodes];
	double * modifiedF = new double[3*numNodes];

	if(!computeSimulationInfo(modifiedK, modifiedF, useNormalOffsets, forceVector)) return false;

	//Do FEM deformation by solving a linear system of equations
	long int dimension, nrhs, info;
	char trans;
	dimension = 3*numNodes;
	nrhs = 1;
	trans = 'N';

	long int lwork = 2 * dimension;
	double * work = new double[lwork];

	clapack::dgels_(&trans, &dimension, &dimension, &nrhs, modifiedK, &dimension, modifiedF, &dimension, work, &lwork, &info);

	if(static_cast<int>(info) != 0)
	{
		//Throw an exception:
		cout<<"Problem with dgels in FEM in computeSimulation, info = "<<info<<endl;
		throw "Problem with dgels in FEM in computeSimulation";
	}

	delete [] work;

	//Set the result to the offset vector of the invalid nodes
	updateOffsets(modifiedF);

	if(updateSurfaceNodes) 
	{
		if(!computeInverseOfMapping()) return false;
	}

	delete [] modifiedK;
	delete [] modifiedF;

	return true;
}

bool LinearFemDeformation::computeSimulationInfo(double *& modifiedK, double *& modifiedF, bool useNormalOffsets, double * forceVector)
{
	int i, j;

	if(mappingIndices.size() == 0) return false;

	if(computeMatrixK() == false) return false;
	
	for(i = 0; i < 3*numNodes * 3*numNodes; i++) 
		modifiedK[i] = matrixK[i];

	if(forceVector != NULL)
	{
		if(forces == NULL) 
			forces = new double[3*numNodes];
		for(i = 0; i < 3*numNodes; i++) 
			forces[i] = forceVector[i];
	}
	else 
	{
		if(forces == NULL) 
			setForces();
	}

	for(i = 0; i < 3*numNodes; i++) 
		modifiedF[i] = forces[i];

	vector<bool> validOffsets = setOffsets(useNormalOffsets);

	//Modify matrixK and modifiedF to account for the fixed nodes:
	for(i = 0; i < numNodes; i++)
	{
		if(validOffsets[i])
		{
			//Set modifiedF to offsets
			modifiedF[3*i] = offsets[3*i];
			modifiedF[3*i+1] = offsets[3*i+1];
			modifiedF[3*i+2] = offsets[3*i+2];

			//Modify the rest of the force vector
			for(j = 0; j < numNodes; j++)
			{
				//If invalid, modify entry of force vector
				if(!validOffsets[j])
				{
					modifiedF[3*j] -= modifiedF[3*i] * matrixK[3*i*(3*numNodes) + 3*j];
					modifiedF[3*j+1] -= modifiedF[3*i] * matrixK[3*i*(3*numNodes) + 3*j+1];
					modifiedF[3*j+2] -= modifiedF[3*i] * matrixK[3*i*(3*numNodes) + 3*j+2];

					modifiedF[3*j] -= modifiedF[3*i+1] * matrixK[(3*i+1)*(3*numNodes) + 3*j];
					modifiedF[3*j+1] -= modifiedF[3*i+1] * matrixK[(3*i+1)*(3*numNodes) + 3*j+1];
					modifiedF[3*j+2] -= modifiedF[3*i+1] * matrixK[(3*i+1)*(3*numNodes) + 3*j+2];

					modifiedF[3*j] -= modifiedF[3*i+2] * matrixK[(3*i+2)*(3*numNodes) + 3*j];
					modifiedF[3*j+1] -= modifiedF[3*i+2] * matrixK[(3*i+2)*(3*numNodes) + 3*j+1];
					modifiedF[3*j+2] -= modifiedF[3*i+2] * matrixK[(3*i+2)*(3*numNodes) + 3*j+2];
				}
			}

			//Set corresponding rows / colummns to zero
			for(j = 0; j < 3*numNodes; j++)
			{
				modifiedK[3*i*(3*numNodes) + j] = 0;
				modifiedK[(3*i+1)*(3*numNodes) + j] = 0;
				modifiedK[(3*i+2)*(3*numNodes) + j] = 0;
				modifiedK[j*(3*numNodes) + 3*i] = 0;
				modifiedK[j*(3*numNodes) + (3*i+1)] = 0;
				modifiedK[j*(3*numNodes) + (3*i+2)] = 0;
			}

			//Set corresponding diagonal elements to one
			modifiedK[3*i*(3*numNodes) + 3*i] = 1;
			modifiedK[(3*i+1)*(3*numNodes) + (3*i+1)] = 1;
			modifiedK[(3*i+2)*(3*numNodes) + (3*i+2)] = 1;
		}	
	}

	return true;
}

void LinearFemDeformation::updateOffsets(double *& modifiedF)
{
	int i;
	double scaleFactor, lengthOffset;

	//Compute the valid nodes
	vector<bool> valid(numNodes);
	for(i = 0; i < numNodes; i++)
	{
		valid[i] = false;
		//Check if the node has a valid mapping and if at least one of the mapped nodes had a valid observation during tracking
		if(mappingIndices[3*i] != -1)
		{
			if((validNodes[mappingIndices[3*i]] != -1) || (validNodes[mappingIndices[3*i+1]] != -1) || (validNodes[mappingIndices[3*i+2]] != -1))
				valid[i] = true;
		}
	}

	for(i = 0; i < numNodes; i++)
	{
		if(!valid[i])
		{
			//If the modifiedF are much larger than the offsets estimated using RBF, rescale them:
			if(lengthRBFOffsets.size() < numNodes) scaleFactor = 1.0;
			else
			{
				lengthOffset = sqrt(pow(modifiedF[3*i],2) + pow(modifiedF[3*i+1],2) + pow(modifiedF[3*i+2],2));
				if((lengthOffset < THRESHOLD_DIV_ZERO) || (lengthOffset <= (MAX_OFFSET * lengthRBFOffsets[i]))) scaleFactor = 1.0;
				else
				{
					scaleFactor = MAX_OFFSET * lengthRBFOffsets[i] / lengthOffset;
				}
			}

			offsets[3*i] = scaleFactor * modifiedF[3*i];
			offsets[3*i+1] = scaleFactor * modifiedF[3*i+1];
			offsets[3*i+2] = scaleFactor * modifiedF[3*i+2];
		}
	}
}

bool LinearFemDeformation::optimizeAndSimulate(bool useNormalOffsets)
{
	int i;

	//Check that things are initialized properly:
	if(contactPointAndForce == NULL) return false;
	if(updatedPositions.numberOfItems() == 0) return false;
	if(validNodes.size() == 0) return false;
	
	//FIRST: OPTIMIZE W.R.T. Young's modulus and Poisson ratio:
	vector<bool> validForces = setForces();
	computeOffsetsRBF(useNormalOffsets);

	vector<bool> forcesForEstimate(numNodes);
	for(i = 0; i < numNodes; i++)
		forcesForEstimate[i] = validForces[i];

	int numIter = 0;
	while(numIter < NUM_ITERATIONS)
	{
		numIter++;

		//initialize the parameters for the optimization
		FemCostFunction costFct(this, forcesForEstimate);

		vnl_lbfgsb * minimizer = new vnl_lbfgsb(costFct);

	#if ACCURACY
		//Settings for moderate accuracy:
		minimizer->set_cost_function_convergence_factor ( 10000000 );
		minimizer->set_projected_gradient_tolerance ( 0.00001 ); 
		minimizer->set_max_function_evals(1000);
	#else
		//Settings for low accuracy:
		minimizer->set_cost_function_convergence_factor ( 1000000000000 );
		minimizer->set_projected_gradient_tolerance ( 0.001 ); 
		minimizer->set_max_function_evals(10);
	#endif
		

		//Set up the vectors
		vnl_vector<double> x(2, 0.0);
		x[0] = poissonRatio;
		x[1] = youngsModulus;

		//Set upper and lower bounds:
		vnl_vector<long> nbd(2, 0);
		nbd[0] = 2;
		nbd[1] = 1;

		vnl_vector<double> u(2, 0.0);
		u[0] = 0.499;

		vnl_vector<double> l(2, 0.0);
		l[0] = 0.001;
		l[1] = 0.001;

		minimizer->set_bound_selection(nbd);
		minimizer->set_upper_bound(u);
		minimizer->set_lower_bound(l);

	#if OUTPUT_OPT_INFO
		//Output (only use during debugging)
		minimizer->set_trace(true);
	#endif

		minimizer->minimize( x );

		delete minimizer;

		//SECOND: UPDATE INFORMATION:
		cout<<"Parameters "<<x[1]<<" "<<x[0]<<endl;
		setParameters(x[1], x[0]);
		
		//Update the unknown forces based on the estimated parameters and offsets
		long int dimension, nrhs;
		double alpha, beta;
		char trans;
		dimension = 3*numNodes;
		nrhs = 1;
		trans = 'N';
		alpha = 1.0;
		beta = 0.0;

		double * computedForces = new double[dimension];
		clapack::dgemm_(&trans, &trans, &dimension, &nrhs, &dimension, &alpha, matrixK, &dimension, offsets, &dimension, &beta, computedForces, &dimension);

		for(i = 0; i < numNodes; i++)
		{
			if(!validForces[i])
			{
				forces[3*i] = computedForces[3*i];
				forces[3*i+1] = computedForces[3*i+1];
				forces[3*i+2] = computedForces[3*i+2];
				
				forcesForEstimate[i] = true;
			}
		}

		delete [] computedForces;
	}

	//THIRD: DO SIMULATION:
	try
	{
		computeSimulation(true, useNormalOffsets);
	}
	catch(char * exceptionString)
	{
		throw exceptionString;
	}

	return true;
}

void LinearFemDeformation::computeOffsetsRBF(bool useNormalOffsets)
{
	int i, countSurfNodes;

	setOffsets(useNormalOffsets, false);

	//Look at the valid offsets and use them to learn an RBF
	countSurfNodes = 0;
	for(i = 0; i < numNodes; i++)
	{
		if(mappingIndices[3*i] != -1) countSurfNodes++;
	}

	double * surfacePts = new double[3*countSurfNodes];
	double * updatedSurfacePts = new double[3*countSurfNodes];
	double * newPt = new double[3];
	double * updatedNewPt = new double[3];

	countSurfNodes = 0;

	for(i = 0; i < numNodes; i++)
	{
		if(mappingIndices[3*i] != -1)
		{
			//Set the surface nodes and the updated surface nodes
			surfacePts[3*countSurfNodes] = allNodes[i]->coordinate().x();
			surfacePts[3*countSurfNodes+1] = allNodes[i]->coordinate().y();
			surfacePts[3*countSurfNodes+2] = allNodes[i]->coordinate().z();

			updatedSurfacePts[3*countSurfNodes] = surfacePts[3*countSurfNodes] + offsets[3*i];
			updatedSurfacePts[3*countSurfNodes+1] = surfacePts[3*countSurfNodes+1] + offsets[3*i+1];
			updatedSurfacePts[3*countSurfNodes+2] = surfacePts[3*countSurfNodes+2] + offsets[3*i+2];

			countSurfNodes++;
		}
	}

	//Compute an RBF mapping
	RBFDeform rbf(countSurfNodes, 3, surfacePts, updatedSurfacePts);

	//Store the lengths of the RBF offsets (to avoid excessive deformation later)
	lengthRBFOffsets.clear();
	lengthRBFOffsets.resize(numNodes);
	
	//Set the offsets at the interior nodes:
	for(i = 0; i < numNodes; i++)
	{
		if(mappingIndices[3*i] == -1)
		{
			//Compute an offset using RBF
			newPt[0] = allNodes[i]->coordinate().x();
			newPt[1] = allNodes[i]->coordinate().y();
			newPt[2] = allNodes[i]->coordinate().z();

			rbf.evaluate(newPt, updatedNewPt);

			offsets[3*i] = updatedNewPt[0] - allNodes[i]->coordinate().x();
			offsets[3*i+1] = updatedNewPt[1] - allNodes[i]->coordinate().y();
			offsets[3*i+2] = updatedNewPt[2] - allNodes[i]->coordinate().z();
		}

		lengthRBFOffsets[i] = sqrt(pow(offsets[3*i],2) + pow(offsets[3*i+1],2) + pow(offsets[3*i+2],2));
	}

	delete [] surfacePts;
	delete [] updatedSurfacePts;
	delete [] newPt;
	delete [] updatedNewPt;
}

vector<bool> LinearFemDeformation::setForces()
{
	int i, j, isContactPoint;
	double weight;

	vector<bool> validForces(numNodes);

	if(forces == NULL) forces = new double[3*numNodes];
	for(i = 0; i < 3*numNodes; i++) forces[i] = 0;

	//Set the forces for the contact points:
	for(i = 0; i < numNodes; i++)
	{
		validForces[i] = false;

		if(mappingIndices[3*i] != -1)
		{
			isContactPoint = -1;
			for(j = 0; j < numContactPoints; j++)
			{
				if((int)contactPointAndForce[4*j + 0] == i) isContactPoint = j;
			}

			//Forces at observed nodes (except possibly contact point) are zero
			if((validNodes[mappingIndices[3*i]] != -1) || (validNodes[mappingIndices[3*i+1]] != -1) || (validNodes[mappingIndices[3*i+2]] != -1))
				validForces[i] = true;
		}
		else 
		{
			isContactPoint = -1;

			//Forces at interior nodes are zero
			validForces[i] = true;
		}

		if(isContactPoint >= 0)
		{
			validForces[i] = true;

			//Optionally weigh the forces at contact points to make the lengths of the two vectors more comparable
			weight = 1.0; // (double)numNodes/(double)numContactPoints;
			forces[3*i] = weight * contactPointAndForce[4*isContactPoint + 1];
			forces[3*i+1] = weight * contactPointAndForce[4*isContactPoint + 2];
			forces[3*i+2] = weight * contactPointAndForce[4*isContactPoint + 3];
		}
	}

	return validForces;
}

vector<bool> LinearFemDeformation::setOffsets(bool useNormalOffsets, bool useOnlyObservedNodes)
{
	int i;
	g_Vector mappedPosition, barycentric, normal;

	if(offsets == NULL) offsets = new double[3*numNodes];
	for(i = 0; i < 3*numNodes; i++) offsets[i] = 0;

	//Compute the valid nodes
	vector<bool> validOffsets(numNodes);
	for(i = 0; i < numNodes; i++)
	{
		validOffsets[i] = false;
		//Check if the node has a valid mapping and if at least one of the mapped nodes had a valid observation during tracking
		if(mappingIndices[3*i] != -1)
		{
			if(useOnlyObservedNodes)
			{
				if((validNodes[mappingIndices[3*i]] != -1) || (validNodes[mappingIndices[3*i+1]] != -1) || (validNodes[mappingIndices[3*i+2]] != -1))
					validOffsets[i] = true;
			}
			else validOffsets[i] = true;
		}
	}

	//Find offsets based on given positions and mapping. Also modify the matrixK and forces
	for(i = 0; i < numNodes; i++)
	{
		//Check if the node has a valid mapping
		if(validOffsets[i])
		{
			barycentric = mappingWeights[3*i] * updatedPositions[mappingIndices[3*i]]->coordinate() + 
				mappingWeights[3*i+1] * updatedPositions[mappingIndices[3*i+1]]->coordinate() + 
				mappingWeights[3*i+2] * updatedPositions[mappingIndices[3*i+2]]->coordinate();
			normal = (updatedPositions[mappingIndices[3*i+1]]->coordinate()-updatedPositions[mappingIndices[3*i]]->coordinate()).Cross(
				(updatedPositions[mappingIndices[3*i+2]]->coordinate()-updatedPositions[mappingIndices[3*i]]->coordinate()));
			normal.Normalize();
			
			if(useNormalOffsets) mappedPosition = barycentric + mappingOffsets[i] * normal;
			else mappedPosition = barycentric;

			//Set offsets
			offsets[3*i] = mappedPosition.x() - allNodes[i]->coordinate().x();
			offsets[3*i+1] = mappedPosition.y() - allNodes[i]->coordinate().y();
			offsets[3*i+2] = mappedPosition.z() - allNodes[i]->coordinate().z();
		}
	}

	return validOffsets;
}

vector<int> LinearFemDeformation::findClosestIndex(g_Vector vec, g_NodeContainer closeNodes, double contactRadius)
{
	int i, bestId;
	vector<int> ids;
	double dist, minDist;

	ids.clear();

	if(contactRadius <= 0)
	{
		for(i = 0; i < (int)closeNodes.numberOfItems(); i++)
		{
			dist = vec.DistanceTo(closeNodes[i]->coordinate());
			if((i == 0) || (dist < minDist))
			{
				minDist = dist;
				bestId = i;
			}
		}
		ids.push_back(bestId);
	}
	else
	{
		for(i = 0; i < (int)closeNodes.numberOfItems(); i++)
		{
			dist = vec.DistanceTo(closeNodes[i]->coordinate());
			if(dist < contactRadius)
				ids.push_back(i);
		}
	}

	return ids;
}

void LinearFemDeformation::FemCostFunction::compute(vnl_vector< double > const &x, double *f, vnl_vector< double > *g)
{
	double * gradient = new double[2];
	*f = getEnergyAndGradient(x[1], x[0], gradient);
	(*g)[0] = gradient[0];
	(*g)[1] = gradient[1];
	std::cerr << "Cost: " << *f << " at " << x[0] << " " << x[1] << " grad: " << gradient[0] << " " << gradient[1] << std::endl;
	delete [] gradient;
}

void LinearFemDeformation::updateCoordinates(g_NodeContainer * deformedPositions, char * exportName)
{
	if(offsets == NULL) return;

	int i;
	for(i = 0; i < numNodes; i++)
		allNodes[i]->coordinate(allNodes[i]->coordinate() + g_Vector(offsets[3*i], offsets[3*i+1], offsets[3*i+2]));

	for(i = 0; i < numNodes; i++)
	{
		if(mappingIndices[3*i] != -1)
		{
			g_Vector barycentric = mappingWeights[3*i] * (*deformedPositions)[mappingIndices[3*i]]->coordinate() + mappingWeights[3*i+1] * (*deformedPositions)[mappingIndices[3*i+1]]->coordinate() + 
						 mappingWeights[3*i+2] * (*deformedPositions)[mappingIndices[3*i+2]]->coordinate();

			allNodes[i]->coordinate(barycentric);
		}
	}

	// EXPORT
	if(exportName != NULL)
	{
		FILE * fpDeformed = fopen(exportName, "w");
		fprintf(fpDeformed, "OFF\n %d 0 0\n", numNodes);
		for(i = 0; i < numNodes; i++)
			fprintf(fpDeformed, "%f %f %f\n", allNodes[i]->coordinate().x(), allNodes[i]->coordinate().y(), allNodes[i]->coordinate().z());
		fclose(fpDeformed);
	}
	//END EXPORT
}

void LinearFemDeformation::setCoordinates(char * filename)
{
	int i, numberNodesCurrentFile, intRead1, intRead2;
	char * charRead = new char[5];
	float floatRead1, floatRead2, floatRead3;

	FILE * fp = fopen(filename, "r");
	fscanf(fp, "%s %d %d %d ", charRead, &numberNodesCurrentFile, &intRead1, &intRead2);

	for(i = 0; i < numberNodesCurrentFile; i++)
	{
		fscanf(fp, "%f %f %f ", &floatRead1, &floatRead2, &floatRead3);
		g_Vector coo(floatRead1, floatRead2, floatRead3);
		allNodes[i]->coordinate(coo);
	}

	fclose(fp);

	delete charRead;
}

double LinearFemDeformation::FemCostFunction::getEnergyAndGradient(double youngsModulusTest, double poissonRatioTest, double *& gradient)
{
	int i, j;

	//Evaluate the energy at youngsModulusTest, poissonRatioTest:
	double youngsModTemp = fem->youngsModulus;
	double poissonRatioTemp = fem->poissonRatio;

	fem->youngsModulus = youngsModulusTest;
	fem->poissonRatio = poissonRatioTest;

	fem->computeMatrixK();

	//Compute energy
	long int dimension, nrhs;
	double alpha, beta;
	char trans;
	dimension = 3*fem->numNodes;
	nrhs = 1;
	trans = 'N';
	alpha = 1.0;
	beta = 0.0;

	double * computedForces = new double[dimension];
	clapack::dgemm_(&trans, &trans, &dimension, &nrhs, &dimension, &alpha, fem->matrixK, &dimension, fem->offsets, &dimension, &beta, computedForces, &dimension);

	double energy = 0;
	for(i = 0; i < fem->numNodes; i++)
	{
		if(validForces[i])
		{
			computedForces[3*i] -= fem->forces[3*i];
			computedForces[3*i+1] -= fem->forces[3*i+1];
			computedForces[3*i+2] -= fem->forces[3*i+2];
			energy += pow(computedForces[3*i], 2);
			energy += pow(computedForces[3*i+1], 2);
			energy += pow(computedForces[3*i+2], 2);
		}
		else
		{
			//Set to zero if invalid in order to facilitate computation of gradient
			computedForces[3*i] = computedForces[3*i+1] = computedForces[3*i+2] = 0;
		}
	}

	//Gradient w.r.t. Poission ratio
	double delta = 0.001;
	double * modifiedK = new double[3*fem->numNodes*3*fem->numNodes];

	//Compute analytic gradient:
	fem->computeMatrixK(0);

	for(i = 0; i < 3*fem->numNodes*3*fem->numNodes; i++) modifiedK[i] = fem->matrixK[i];

	gradient[0] = 0;
	for(i = 0; i < 3*fem->numNodes; i++)
	{
		for(j = 0; j < 3*fem->numNodes; j++)
			gradient[0] += 2.0*computedForces[i]*fem->offsets[j]*modifiedK[j*3*fem->numNodes+i];
	}

	//Gradient w.r.t. Young's modulus
	//Compute analytic gradient:
	fem->computeMatrixK(1);

	for(i = 0; i < 3*fem->numNodes*3*fem->numNodes; i++) modifiedK[i] = fem->matrixK[i];


	gradient[1] = 0;
	for(i = 0; i < 3*fem->numNodes; i++)
	{
		for(j = 0; j < 3*fem->numNodes; j++)
			gradient[1] += 2.0*computedForces[i]*fem->offsets[j]*modifiedK[j*3*fem->numNodes+i];
	}

	delete [] computedForces;
	delete [] modifiedK;

	//Revert the temporaty change:
	fem->youngsModulus = youngsModTemp;
	fem->poissonRatio = poissonRatioTemp;

	return energy;
}

bool LinearFemDeformation::computeInverseOfMapping(bool useNormalOffsets)
{
	if(mappingIndices.size() == 0) return false;

	//use offsets to compute the new positions on the previously invalid nodes that have a matching vertex
	int i, numVerticesInSurfaceMesh, id1, id2, id3;
	g_Vector node1, node2, node3, corresPos, barycentric, normal;

	numVerticesInSurfaceMesh = (int)updatedPositions.numberOfItems();

	for(i = 0; i < numVerticesInSurfaceMesh; i++)
	{
		//If we were not sure about the position before and if there exists a valid match, set the new coordinate:
		if((validNodes[i] == -1) && (invMappingIndices[3*i] != -1))
		{
			id1 = invMappingIndices[3*i];
			id2 = invMappingIndices[3*i+1];
			id3 = invMappingIndices[3*i+2];

			node1 = allNodes[id1]->coordinate() + g_Vector(offsets[3*id1], offsets[3*id1+1], offsets[3*id1+2]);
			node2 = allNodes[id2]->coordinate() + g_Vector(offsets[3*id2], offsets[3*id2+1], offsets[3*id2+2]);
			node3 = allNodes[id3]->coordinate() + g_Vector(offsets[3*id3], offsets[3*id3+1], offsets[3*id3+2]);

			normal = (node2 - node1).Cross(node3 - node1);
			normal.Normalize();
			barycentric = invMappingWeights[3*i] * node1 + invMappingWeights[3*i+1] * node2 + invMappingWeights[3*i+2] * node3;

			if(useNormalOffsets) corresPos = barycentric + invMappingOffsets[i] * normal;
			else corresPos = barycentric;
			
			updatedPositions[i]->coordinate(corresPos);
			validNodes[i] = 2;
		}
	}

	return true;
}

bool LinearFemDeformation::computeMatrixK(int derivative)
{
	int i, j, k, mappedJ, mappedK;

	if((allNodes.numberOfItems() == 0) || (allTetrahedra.size() == 0)) return false;
	if(youngsModulus < 0) return false;
	if((poissonRatio < 0) || (poissonRatio >= 0.5)) return false;

	if(matrixK == NULL) matrixK = new double[3*numNodes * 3*numNodes];
	for(i = 0; i < 3*numNodes * 3*numNodes; i++) matrixK[i] = 0;
	double * localK = new double[12*12];
	for(i = 0; i < allTetrahedra.size(); i++)
	{
		computeLocalMatrixK(localK, i, derivative);
		for(j = 0; j < 4; j++)
		{
			mappedJ = allTetrahedra[i][j];
			for(k = 0; k < 4; k++)
			{
				mappedK = allTetrahedra[i][k];

				matrixK[(3*mappedK)*3*numNodes + (3*mappedJ)] += localK[(3*k)*12 + (3*j)];
				matrixK[(3*mappedK)*3*numNodes + (3*mappedJ+1)] += localK[(3*k)*12 + (3*j+1)];
				matrixK[(3*mappedK)*3*numNodes + (3*mappedJ+2)] += localK[(3*k)*12 + (3*j+2)];
				
				matrixK[(3*mappedK+1)*3*numNodes + (3*mappedJ)] += localK[(3*k+1)*12 + (3*j)];
				matrixK[(3*mappedK+1)*3*numNodes + (3*mappedJ+1)] += localK[(3*k+1)*12 + (3*j+1)];
				matrixK[(3*mappedK+1)*3*numNodes + (3*mappedJ+2)] += localK[(3*k+1)*12 + (3*j+2)];
				
				matrixK[(3*mappedK+2)*3*numNodes + (3*mappedJ)] += localK[(3*k+2)*12 + (3*j)];
				matrixK[(3*mappedK+2)*3*numNodes + (3*mappedJ+1)] += localK[(3*k+2)*12 + (3*j+1)];
				matrixK[(3*mappedK+2)*3*numNodes + (3*mappedJ+2)] += localK[(3*k+2)*12 + (3*j+2)];
			}
		}
	}
	delete [] localK;

	return true;
}

void LinearFemDeformation::computeLocalMatrixK(double *& localK, int tetrahedronId, int derivative)
{
	if((tetrahedronId < 0) || (tetrahedronId >= allTetrahedra.size())) return;
	int i, offset;
	double volume, mu, lambda, det, helper;
	double * matrixC, * matrixB, * helpMat;

	matrixC = new double[6*6];
	matrixB = new double[6*12];
	helpMat = new double[6*12];
	vector<double> x, y, z;

	for(i = 0; i < 12*12; i++) localK[i] = 0;

	//Compute the constant matrix C:
	
	if(derivative == 0)
	{
		//Derivative of Lame parameters w.r.t. poissonRatio
		lambda = (youngsModulus * ((1.0+poissonRatio)*(1.0-2.0*poissonRatio) + poissonRatio*(1.0+4.0*poissonRatio))) / (pow((1.0+poissonRatio)*(1.0-2.0*poissonRatio), 2));
		mu = -youngsModulus / (2.0*pow((1.0+poissonRatio), 2));
	}
	else if(derivative == 1)
	{
		//Derivative of Lame parameters w.r.t. youngsModulus
		lambda = poissonRatio / ((1.0+poissonRatio)*(1.0-2.0*poissonRatio));
		mu = 1.0 / (2.0*(1.0+poissonRatio));
	}
	else
	{
		//Convert to Lame parameters (formulas, see http://www.engr.uconn.edu/~lanbo/G229Lect08031RockMech2.pdf)
		lambda = (youngsModulus*poissonRatio) / ((1.0+poissonRatio)*(1.0-2.0*poissonRatio));
		mu = youngsModulus / (2.0*(1.0+poissonRatio)); 
	}

	for(i = 0; i < 6*6; i++) matrixC[i] = 0;
	matrixC[0] = matrixC[7] = matrixC[14] = lambda + 2*mu;
	matrixC[21] = matrixC[28] = matrixC[35] = mu;
	matrixC[1] = matrixC[2] = matrixC[6] = matrixC[8] = matrixC[12] = matrixC[13] = lambda;

	//Compute the matrix B:
	x.push_back(allNodes[allTetrahedra[tetrahedronId][0]]->coordinate().x());
	x.push_back(allNodes[allTetrahedra[tetrahedronId][1]]->coordinate().x());
	x.push_back(allNodes[allTetrahedra[tetrahedronId][2]]->coordinate().x());
	x.push_back(allNodes[allTetrahedra[tetrahedronId][3]]->coordinate().x());

	y.push_back(allNodes[allTetrahedra[tetrahedronId][0]]->coordinate().y());
	y.push_back(allNodes[allTetrahedra[tetrahedronId][1]]->coordinate().y());
	y.push_back(allNodes[allTetrahedra[tetrahedronId][2]]->coordinate().y());
	y.push_back(allNodes[allTetrahedra[tetrahedronId][3]]->coordinate().y());

	z.push_back(allNodes[allTetrahedra[tetrahedronId][0]]->coordinate().z());
	z.push_back(allNodes[allTetrahedra[tetrahedronId][1]]->coordinate().z());
	z.push_back(allNodes[allTetrahedra[tetrahedronId][2]]->coordinate().z());
	z.push_back(allNodes[allTetrahedra[tetrahedronId][3]]->coordinate().z());

	//1. Compute the volume
	for(i = 0; i < 6*12; i++) helpMat[i] = 0;
	helpMat[0] = helpMat[4] = helpMat[8] = helpMat[12] = 1.0;
	helpMat[1] = x[0]; helpMat[2] = y[0]; helpMat[3] = z[0];
	helpMat[5] = x[1]; helpMat[6] = y[1]; helpMat[7] = z[1];
	helpMat[9] = x[2]; helpMat[10] = y[2]; helpMat[11] = z[2];
	helpMat[13] = x[3]; helpMat[14] = y[3]; helpMat[15] = z[3];

	volume = computeDeterminant(helpMat, 4) / 6.0;

	//If the volume is zero, just give zero back (no need to compute everything since it will be multiplied by zero in the end anyways)
	if(volume == 0)
	{
		for(i = 0; i < 12*12; i++) localK[i] = 0;
		
		delete [] matrixC;
		delete [] matrixB;
		delete [] helpMat;
	}

	else
	{
		//2. Fill the matrix
		for(i = 0; i < 6*12; i++) matrixB[i] = 0;

		//b1
		for(i = 0; i < 6*12; i++) helpMat[i] = 0;
		helpMat[0] = helpMat[3] = helpMat[6] = 1.0;
		helpMat[1] = y[1]; helpMat[2] = z[1]; helpMat[4] = y[2]; helpMat[5] = z[2]; helpMat[7] = y[3]; helpMat[8] = z[3];
		det = -computeDeterminant(helpMat, 3);
		matrixB[0] = matrixB[9] = matrixB[17] = det;
		//c1
		for(i = 0; i < 6*12; i++) helpMat[i] = 0;
		helpMat[0] = helpMat[3] = helpMat[6] = 1.0;
		helpMat[1] = x[1]; helpMat[2] = z[1]; helpMat[4] = x[2]; helpMat[5] = z[2]; helpMat[7] = x[3]; helpMat[8] = z[3];
		det = computeDeterminant(helpMat, 3);
		matrixB[3] = matrixB[7] = matrixB[16] = det;
		//d1
		for(i = 0; i < 6*12; i++) helpMat[i] = 0;
		helpMat[0] = helpMat[3] = helpMat[6] = 1.0;
		helpMat[1] = x[1]; helpMat[2] = y[1]; helpMat[4] = x[2]; helpMat[5] = y[2]; helpMat[7] = x[3]; helpMat[8] = y[3];
		det = -computeDeterminant(helpMat, 3);
		matrixB[5] = matrixB[10] = matrixB[14] = det;
		
		//b2
		offset = 18;
		for(i = 0; i < 6*12; i++) helpMat[i] = 0;
		helpMat[0] = helpMat[3] = helpMat[6] = 1.0;
		helpMat[1] = y[2]; helpMat[2] = z[2]; helpMat[4] = y[3]; helpMat[5] = z[3]; helpMat[7] = y[0]; helpMat[8] = z[0];
		det = computeDeterminant(helpMat, 3);
		matrixB[offset+0] = matrixB[offset+9] = matrixB[offset+17] = det;
		//c2
		for(i = 0; i < 6*12; i++) helpMat[i] = 0;
		helpMat[0] = helpMat[3] = helpMat[6] = 1.0;
		helpMat[1] = x[2]; helpMat[2] = z[2]; helpMat[4] = x[3]; helpMat[5] = z[3]; helpMat[7] = x[0]; helpMat[8] = z[0];
		det = -computeDeterminant(helpMat, 3);
		matrixB[offset+3] = matrixB[offset+7] = matrixB[offset+16] = det;
		//d2
		for(i = 0; i < 6*12; i++) helpMat[i] = 0;
		helpMat[0] = helpMat[3] = helpMat[6] = 1.0;
		helpMat[1] = x[2]; helpMat[2] = y[2]; helpMat[4] = x[3]; helpMat[5] = y[3]; helpMat[7] = x[0]; helpMat[8] = y[0];
		det = computeDeterminant(helpMat, 3);
		matrixB[offset+5] = matrixB[offset+10] = matrixB[offset+14] = det;

		//b3
		offset = 36;
		for(i = 0; i < 6*12; i++) helpMat[i] = 0;
		helpMat[0] = helpMat[3] = helpMat[6] = 1.0;
		helpMat[1] = y[3]; helpMat[2] = z[3]; helpMat[4] = y[0]; helpMat[5] = z[0]; helpMat[7] = y[1]; helpMat[8] = z[1];
		det = -computeDeterminant(helpMat, 3);
		matrixB[offset+0] = matrixB[offset+9] = matrixB[offset+17] = det;
		//c3
		for(i = 0; i < 6*12; i++) helpMat[i] = 0;
		helpMat[0] = helpMat[3] = helpMat[6] = 1.0;
		helpMat[1] = x[3]; helpMat[2] = z[3]; helpMat[4] = x[0]; helpMat[5] = z[0]; helpMat[7] = x[1]; helpMat[8] = z[1];
		det = computeDeterminant(helpMat, 3);
		matrixB[offset+3] = matrixB[offset+7] = matrixB[offset+16] = det;
		//d3
		for(i = 0; i < 6*12; i++) helpMat[i] = 0;
		helpMat[0] = helpMat[3] = helpMat[6] = 1.0;
		helpMat[1] = x[3]; helpMat[2] = y[3]; helpMat[4] = x[0]; helpMat[5] = y[0]; helpMat[7] = x[1]; helpMat[8] = y[1];
		det = -computeDeterminant(helpMat, 3);
		matrixB[offset+5] = matrixB[offset+10] = matrixB[offset+14] = det;

		//b4
		offset = 54;
		for(i = 0; i < 6*12; i++) helpMat[i] = 0;
		helpMat[0] = helpMat[3] = helpMat[6] = 1.0;
		helpMat[1] = y[0]; helpMat[2] = z[0]; helpMat[4] = y[1]; helpMat[5] = z[1]; helpMat[7] = y[2]; helpMat[8] = z[2];
		det = computeDeterminant(helpMat, 3);
		matrixB[offset+0] = matrixB[offset+9] = matrixB[offset+17] = det;
		//c4
		for(i = 0; i < 6*12; i++) helpMat[i] = 0;
		helpMat[0] = helpMat[3] = helpMat[6] = 1.0;
		helpMat[1] = x[0]; helpMat[2] = z[0]; helpMat[4] = x[1]; helpMat[5] = z[1]; helpMat[7] = x[2]; helpMat[8] = z[2];
		det = -computeDeterminant(helpMat, 3);
		matrixB[offset+3] = matrixB[offset+7] = matrixB[offset+16] = det;
		//d4
		for(i = 0; i < 6*12; i++) helpMat[i] = 0;
		helpMat[0] = helpMat[3] = helpMat[6] = 1.0;
		helpMat[1] = x[0]; helpMat[2] = y[0]; helpMat[4] = x[1]; helpMat[5] = y[1]; helpMat[7] = x[2]; helpMat[8] = y[2];
		det = computeDeterminant(helpMat, 3);
		matrixB[offset+5] = matrixB[offset+10] = matrixB[offset+14] = det;

		//Multiply the matrices: V_e B_e^T C B_e
		long int dim1, dim2;
		char transY, transN;
		double alpha, beta;
		transY = 'T';
		transN = 'N';
		dim1 = 6; 
		dim2 = 12;
		alpha = 1.0;
		beta = 0.0;

		for(i = 0; i < 6*12; i++) helpMat[i] = 0;

		clapack::dgemm_(&transY, &transN, &dim2, &dim1, &dim1, &alpha, matrixB, &dim1, matrixC, &dim1, &beta, helpMat, &dim2);
		clapack::dgemm_(&transN, &transN, &dim2, &dim2, &dim1, &alpha, helpMat, &dim2, matrixB, &dim1, &beta, localK, &dim2);

		for(i = 0; i < 12*12; i++) 
		{
			helper = localK[i] / (36.0*volume); 
			if((fem_isnan(helper) == 0) && (fem_isinf(helper) == 0))
				localK[i] = helper;
			else
				localK[i] = localK[i] / THRESHOLD_DIV_ZERO;
		}

		delete [] matrixC;
		delete [] matrixB;
		delete [] helpMat;
	}
}

double LinearFemDeformation::computeDeterminant(double * matrix, int matrixDimension)
{
	//Check for validity of dimension
	if((matrixDimension < 1) || (matrixDimension > 4)) 
	{
		cout<<"Invalid dimensions in computeDeterminant"<<endl;
		return -1;
	}
	else if(matrixDimension == 1) return matrix[0];
	else if(matrixDimension == 2) return matrix[0]*matrix[3]-matrix[1]*matrix[2];
	else if(matrixDimension == 3) 
	{
		double det = matrix[0]*matrix[4]*matrix[8] + matrix[1]*matrix[5]*matrix[6] + matrix[2]*matrix[3]*matrix[7] - 
			   matrix[0]*matrix[5]*matrix[7] - matrix[1]*matrix[3]*matrix[8] - matrix[2]*matrix[4]*matrix[6];
		return det;
	}
	else 
	{
		double * helpMat = new double[9];
		double det = 0;

		//Use Laplace formula
		helpMat[0] = matrix[5]; helpMat[1] = matrix[6]; helpMat[2] = matrix[7];
		helpMat[3] = matrix[9]; helpMat[4] = matrix[10]; helpMat[5] = matrix[11];
		helpMat[6] = matrix[13]; helpMat[7] = matrix[14]; helpMat[8] = matrix[15];


		det += matrix[0] * computeDeterminant(helpMat, 3);

		helpMat[0] = matrix[1]; helpMat[1] = matrix[2]; helpMat[2] = matrix[3];
		helpMat[3] = matrix[9]; helpMat[4] = matrix[10]; helpMat[5] = matrix[11];
		helpMat[6] = matrix[13]; helpMat[7] = matrix[14]; helpMat[8] = matrix[15];

		det -= matrix[4] * computeDeterminant(helpMat, 3);

		helpMat[0] = matrix[1]; helpMat[1] = matrix[2]; helpMat[2] = matrix[3];
		helpMat[3] = matrix[5]; helpMat[4] = matrix[6]; helpMat[5] = matrix[7];
		helpMat[6] = matrix[13]; helpMat[7] = matrix[14]; helpMat[8] = matrix[15];

		det += matrix[8] * computeDeterminant(helpMat, 3);

		helpMat[0] = matrix[1]; helpMat[1] = matrix[2]; helpMat[2] = matrix[3];
		helpMat[3] = matrix[5]; helpMat[4] = matrix[6]; helpMat[5] = matrix[7];
		helpMat[6] = matrix[9]; helpMat[7] = matrix[10]; helpMat[8] = matrix[11];

		det -= matrix[12] * computeDeterminant(helpMat, 3);

		delete [] helpMat;

		return det;
	}
}
