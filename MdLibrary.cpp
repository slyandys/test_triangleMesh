#include "MdLibrary.h"

double MinimumDistance::minimumDistance(g_Vector &p1, g_Vector &p2, g_Vector &p3, g_Vector &p4, bool bSquared, double ** barycenters)
{
	g_Vector u, v;
	double d2121;

	u = p2 - p1;
	v = p4 - p3;
	d2121 = u.Dot(u);
	return minimumDistance(p1, u, p3, v, d2121, bSquared, barycenters);
}

double MinimumDistance::minimumDistance(g_Vector &p1, g_Vector &u, g_Vector &p3, g_Vector &v, double d2121, bool bSquared, double ** barycenters)
{
	double dMD = 0.0;
	g_Vector w;
	double d4343, d4321, d1343, d1321;
	double mu_a, mu_b, cl_a, cl_b;
	w = p1 - p3;

	d4343 = v.Dot(v);
	d4321 = v.Dot(u);
	d1343 = w.Dot(v);
	d1321 = w.Dot(u);

	mu_a = (d1343 * d4321 - d1321 * d4343) / (d2121 * d4343 - d4321 * d4321);
	mu_b = (d1343 + mu_a * d4321) / d4343;

	cl_a = max(0.0, mu_a);
	cl_a = min(1.0, cl_a);
	cl_b = max(0.0, mu_b);
	cl_b = min(1.0, cl_b);

	if(barycenters != NULL)
	{
		(*barycenters)[0] = 1.0 - cl_a;
		(*barycenters)[1] = cl_a;
		(*barycenters)[2] = 1.0 - cl_b;
		(*barycenters)[3] = cl_b;
	}

	if (bSquared)
		dMD = ((p1 + u * cl_a) - (p3 + v * cl_b)).SquaredLength();
	else
		dMD = ((p1 + u * cl_a) - (p3 + v * cl_b)).Length();
	return dMD;
}

double MinimumDistance::minimumDistance(g_Element & t1, g_Element & t2, int * returnSmallest, double ** barycenters)
{
	double dist, minDist;
	double * tempBarycenters = new double[6];

	g_NodeContainer nodes1, nodes2;
	nodes1 = t1.nodes();
	nodes2 = t2.nodes();
	g_Plane plane1(nodes1[0]->coordinate(), nodes1[1]->coordinate(), nodes1[2]->coordinate());
	g_Plane plane2(nodes2[0]->coordinate(), nodes2[1]->coordinate(), nodes2[2]->coordinate());

	dist = -1;
	if(barycenters != NULL) (*barycenters)[0] = (*barycenters)[1] = (*barycenters)[2] = (*barycenters)[3] = (*barycenters)[4] = (*barycenters)[5] = 0;
	if(returnSmallest != NULL) *returnSmallest = -1;
	if((nodes1[0] == nodes2[0]) || (nodes1[0] == nodes2[1]) || (nodes1[0] == nodes2[2]) ||
		(nodes1[1] == nodes2[0]) || (nodes1[1] == nodes2[1]) || (nodes1[1] == nodes2[2]) ||
		(nodes1[2] == nodes2[0]) || (nodes1[2] == nodes2[1]) || (nodes1[2] == nodes2[2])) 
		minDist = dist = 0.0;
	if(dist == -1)
	{
		// Segment-segment distances:
		minDist = minimumDistance(nodes1[0]->coordinate(), nodes1[1]->coordinate(), 
			nodes2[0]->coordinate(), nodes2[1]->coordinate(), false, &tempBarycenters);
		if(returnSmallest != NULL) *returnSmallest = 0;
		if(barycenters != NULL)
		{
			(*barycenters)[0] = tempBarycenters[0]; (*barycenters)[1] = tempBarycenters[1]; (*barycenters)[2] = 0;
			(*barycenters)[3] = tempBarycenters[2]; (*barycenters)[4] = tempBarycenters[3]; (*barycenters)[5] = 0; 
		}

		dist = minimumDistance(nodes1[0]->coordinate(), nodes1[1]->coordinate(), 
			nodes2[1]->coordinate(), nodes2[2]->coordinate(), false, &tempBarycenters);
		if(dist < minDist)
		{
			minDist = dist;
			if(returnSmallest != NULL) *returnSmallest = 1;
			if(barycenters != NULL)
			{
				(*barycenters)[0] = tempBarycenters[0]; (*barycenters)[1] = tempBarycenters[1]; (*barycenters)[2] = 0;
				(*barycenters)[3] = 0; (*barycenters)[4] = tempBarycenters[2]; (*barycenters)[5] = tempBarycenters[3]; 
			}
		}

		dist = minimumDistance(nodes1[0]->coordinate(), nodes1[1]->coordinate(), 
			nodes2[2]->coordinate(), nodes2[0]->coordinate(), false, &tempBarycenters);
		if(dist < minDist)
		{
			minDist = dist;
			if(returnSmallest != NULL) *returnSmallest = 2;
			if(barycenters != NULL)
			{
				(*barycenters)[0] = tempBarycenters[0]; (*barycenters)[1] = tempBarycenters[1]; (*barycenters)[2] = 0;
				(*barycenters)[3] = tempBarycenters[3]; (*barycenters)[4] = 0; (*barycenters)[5] = tempBarycenters[2]; 
			}
		}

		dist = minimumDistance(nodes1[1]->coordinate(), nodes1[2]->coordinate(), 
			nodes2[0]->coordinate(), nodes2[1]->coordinate(), false, &tempBarycenters);
		if(dist < minDist)
		{
			minDist = dist;
			if(returnSmallest != NULL) *returnSmallest = 3;
			if(barycenters != NULL)
			{
				(*barycenters)[0] = 0; (*barycenters)[1] = tempBarycenters[0]; (*barycenters)[2] = tempBarycenters[1];
				(*barycenters)[3] = tempBarycenters[2]; (*barycenters)[4] = tempBarycenters[3]; (*barycenters)[5] = 0; 
			}
		}

		dist = minimumDistance(nodes1[1]->coordinate(), nodes1[2]->coordinate(), 
			nodes2[1]->coordinate(), nodes2[2]->coordinate(), false, &tempBarycenters);
		if(dist < minDist)
		{
			minDist = dist;
			if(returnSmallest != NULL) *returnSmallest = 4;
			if(barycenters != NULL)
			{
				(*barycenters)[0] = 0; (*barycenters)[1] = tempBarycenters[0]; (*barycenters)[2] = tempBarycenters[1];
				(*barycenters)[3] = 0; (*barycenters)[4] = tempBarycenters[2]; (*barycenters)[5] = tempBarycenters[3]; 
			}
		}

		dist = minimumDistance(nodes1[1]->coordinate(), nodes1[2]->coordinate(), 
			nodes2[2]->coordinate(), nodes2[0]->coordinate(), false, &tempBarycenters);
		if(dist < minDist)
		{
			minDist = dist;
			if(returnSmallest != NULL) *returnSmallest = 5;
			if(barycenters != NULL)
			{
				(*barycenters)[0] = 0; (*barycenters)[1] = tempBarycenters[0]; (*barycenters)[2] = tempBarycenters[1];
				(*barycenters)[3] = tempBarycenters[3]; (*barycenters)[4] = 0; (*barycenters)[5] = tempBarycenters[2];
			}
		}

		dist = minimumDistance(nodes1[2]->coordinate(), nodes1[0]->coordinate(), 
			nodes2[0]->coordinate(), nodes2[1]->coordinate(), false, &tempBarycenters);
		if(dist < minDist)
		{
			minDist = dist;
			if(returnSmallest != NULL) *returnSmallest = 6;
			if(barycenters != NULL)
			{
				(*barycenters)[0] = tempBarycenters[1]; (*barycenters)[1] = 0; (*barycenters)[2] = tempBarycenters[0];
				(*barycenters)[3] = tempBarycenters[2]; (*barycenters)[4] = tempBarycenters[3]; (*barycenters)[5] = 0;
			}
		}

		dist = minimumDistance(nodes1[2]->coordinate(), nodes1[0]->coordinate(), 
			nodes2[1]->coordinate(), nodes2[2]->coordinate(), false, &tempBarycenters);
		if(dist < minDist)
		{
			minDist = dist;
			if(returnSmallest != NULL) *returnSmallest = 7;
			if(barycenters != NULL)
			{
				(*barycenters)[0] = tempBarycenters[1]; (*barycenters)[1] = 0; (*barycenters)[2] = tempBarycenters[0];
				(*barycenters)[3] = 0; (*barycenters)[4] = tempBarycenters[2]; (*barycenters)[5] = tempBarycenters[3];
			}
		}

		dist = minimumDistance(nodes1[2]->coordinate(), nodes1[0]->coordinate(), 
			nodes2[2]->coordinate(), nodes2[0]->coordinate(), false, &tempBarycenters);
		if(dist < minDist)
		{
			minDist = dist;
			if(returnSmallest != NULL) *returnSmallest = 8;
			if(barycenters != NULL)
			{
				(*barycenters)[0] = tempBarycenters[1]; (*barycenters)[1] = 0; (*barycenters)[2] = tempBarycenters[0];
				(*barycenters)[3] = tempBarycenters[3]; (*barycenters)[4] = 0; (*barycenters)[5] = tempBarycenters[2];
			}
		}

		// Point-plane distances:
		if(projectionInTriangle(t1, nodes2[0]->coordinate(), &tempBarycenters))
		{
			dist = plane1.distanceToPoint(nodes2[0]->coordinate());
			if(dist < minDist)
			{
				minDist = dist;
				if(returnSmallest != NULL) *returnSmallest = 9;
				if(barycenters != NULL)
				{
					(*barycenters)[0] = tempBarycenters[0]; (*barycenters)[1] = tempBarycenters[1]; (*barycenters)[2] = tempBarycenters[2];
					(*barycenters)[3] = 1; (*barycenters)[4] = 0; (*barycenters)[5] = 0;
				}
			}
		}
		if(projectionInTriangle(t1, nodes2[1]->coordinate(), &tempBarycenters))
		{
			dist = plane1.distanceToPoint(nodes2[1]->coordinate());
			if(dist < minDist)
			{
				minDist = dist;
				if(returnSmallest != NULL) *returnSmallest = 10;
				if(barycenters != NULL)
				{
					(*barycenters)[0] = tempBarycenters[0]; (*barycenters)[1] = tempBarycenters[1]; (*barycenters)[2] = tempBarycenters[2];
					(*barycenters)[3] = 0; (*barycenters)[4] = 1; (*barycenters)[5] = 0;
				}
			}
		}
		if(projectionInTriangle(t1, nodes2[2]->coordinate(), &tempBarycenters))
		{
			dist = plane1.distanceToPoint(nodes2[2]->coordinate());
			if(dist < minDist)
			{
				minDist = dist;
				if(returnSmallest != NULL) *returnSmallest = 11;
				if(barycenters != NULL)
				{
					(*barycenters)[0] = tempBarycenters[0]; (*barycenters)[1] = tempBarycenters[1]; (*barycenters)[2] = tempBarycenters[2];
					(*barycenters)[3] = 0; (*barycenters)[4] = 0; (*barycenters)[5] = 1;
				}
			}
		}
		if(projectionInTriangle(t2, nodes1[0]->coordinate(), &tempBarycenters))
		{
			dist = plane2.distanceToPoint(nodes1[0]->coordinate());
			if(dist < minDist)
			{
				minDist = dist;
				if(returnSmallest != NULL) *returnSmallest = 12;
				if(barycenters != NULL)
				{
					(*barycenters)[0] = 1; (*barycenters)[1] = 0; (*barycenters)[2] = 0;
					(*barycenters)[3] = tempBarycenters[0]; (*barycenters)[4] = tempBarycenters[1]; (*barycenters)[5] = tempBarycenters[2];
				}
			}
		}
		if(projectionInTriangle(t2, nodes1[1]->coordinate(), &tempBarycenters))
		{
			dist = plane2.distanceToPoint(nodes1[1]->coordinate());
			if(dist < minDist)
			{
				minDist = dist;
				if(returnSmallest != NULL) *returnSmallest = 13;
				if(barycenters != NULL)
				{
					(*barycenters)[0] = 0; (*barycenters)[1] = 1; (*barycenters)[2] = 0;
					(*barycenters)[3] = tempBarycenters[0]; (*barycenters)[4] = tempBarycenters[1]; (*barycenters)[5] = tempBarycenters[2];
				}
			}
		}
		if(projectionInTriangle(t2, nodes1[2]->coordinate(), &tempBarycenters))
		{
			dist = plane2.distanceToPoint(nodes1[2]->coordinate());
			if(dist < minDist)
			{
				minDist = dist;
				if(returnSmallest != NULL) *returnSmallest = 14;
				if(barycenters != NULL)
				{
					(*barycenters)[0] = 0; (*barycenters)[1] = 0; (*barycenters)[2] = 1;
					(*barycenters)[3] = tempBarycenters[0]; (*barycenters)[4] = tempBarycenters[1]; (*barycenters)[5] = tempBarycenters[2];
				}
			}
		}
	}

	delete [] tempBarycenters;

	return minDist;
}

bool MinimumDistance::projectionInTriangle(g_Element t, g_Vector p, double ** barycenters)
{
	int turn1, turn2, turn3;
	double dot;

	g_Plane plane(t.nodes()[0]->coordinate(), t.nodes()[1]->coordinate(), t.nodes()[2]->coordinate());
	g_Vector projection = plane.project(p);

	dot = (t.nodes()[1]->coordinate()-t.nodes()[0]->coordinate()).Cross(projection-t.nodes()[1]->coordinate()).Dot(plane.normal());
	if(dot < 0) turn1 = -1;
	else turn1 = 1;
	dot = (t.nodes()[2]->coordinate()-t.nodes()[1]->coordinate()).Cross(projection-t.nodes()[2]->coordinate()).Dot(plane.normal());
	if(dot < 0) turn2 = -1;
	else turn2 = 1;
	dot = (t.nodes()[0]->coordinate()-t.nodes()[2]->coordinate()).Cross(projection-t.nodes()[0]->coordinate()).Dot(plane.normal());
	if(dot < 0) turn3 = -1;
	else turn3 = 1;

	if(turn1 == turn2 && turn2 == turn3) 
	{
		if(barycenters != NULL) 
		{
			//compute barycentric coordinates of projection in t:
			double areaABC, areaPBC, areaPCA;	
			g_Vector normal = (t.nodes()[1]->coordinate()-t.nodes()[0]->coordinate()).Cross(t.nodes()[2]->coordinate()-t.nodes()[0]->coordinate());
			normal.Normalize();
			// Compute twice area of triangle ABC
			areaABC = normal.Dot((t.nodes()[1]->coordinate()-t.nodes()[0]->coordinate()).Cross(t.nodes()[2]->coordinate()-t.nodes()[0]->coordinate()));
			// Compute barycentrics
			areaPBC = normal.Dot((t.nodes()[1]->coordinate()-projection).Cross(t.nodes()[2]->coordinate()-projection));
			(*barycenters)[0] = areaPBC / areaABC;
			areaPCA = normal.Dot((t.nodes()[2]->coordinate()-projection).Cross(t.nodes()[0]->coordinate()-projection));
			(*barycenters)[1] = areaPCA / areaABC;
			(*barycenters)[2] = 1.0 - (*barycenters)[0] - (*barycenters)[1];
		}

		return true;
	}

	else return false;
}

bool MinimumDistance::intersect(g_Element & t1, g_Element & t2, double threshold)
{
	g_NodeContainer nodes1, nodes2;
	g_Vector intersection;
	double tVal;
	bool computationOK;
	nodes1 = t1.nodes();
	nodes2 = t2.nodes();
	g_Plane plane1(nodes1[0]->coordinate(), nodes1[1]->coordinate(), nodes1[2]->coordinate());
	g_Plane plane2(nodes2[0]->coordinate(), nodes2[1]->coordinate(), nodes2[2]->coordinate());

	if((nodes1[0] == nodes2[0]) || (nodes1[0] == nodes2[1]) || (nodes1[0] == nodes2[2]) ||
		(nodes1[1] == nodes2[0]) || (nodes1[1] == nodes2[1]) || (nodes1[1] == nodes2[2]) ||
		(nodes1[2] == nodes2[0]) || (nodes1[2] == nodes2[1]) || (nodes1[2] == nodes2[2])) 
		return true;
	computationOK = plane1.line_plane_intersection(nodes2[1]->coordinate()-nodes2[0]->coordinate(),nodes2[0]->coordinate(),intersection,tVal);
	if(computationOK)
	{
		if((0.0-threshold<=tVal) && (tVal<=1.0+threshold))
		{
			if(projectionInTriangle(t1, intersection)) 
				return true;
		}
	}
	else
	{
		if(intersectInPlane(t1, nodes2[0]->coordinate(), nodes2[1]->coordinate(), threshold))
			return true;
	}
	computationOK = plane1.line_plane_intersection(nodes2[2]->coordinate()-nodes2[1]->coordinate(),nodes2[1]->coordinate(),intersection,tVal);
	if(computationOK)
	{
		if((0.0-threshold<=tVal) && (tVal<=1.0+threshold))
		{
			if(projectionInTriangle(t1, intersection)) 
				return true;
		}
	}
	else
	{
		if(intersectInPlane(t1, nodes2[1]->coordinate(), nodes2[2]->coordinate(), threshold))
			return true;
	}
	computationOK = plane1.line_plane_intersection(nodes2[0]->coordinate()-nodes2[2]->coordinate(),nodes2[2]->coordinate(),intersection,tVal);
	if(computationOK)
	{
		if((0.0-threshold<=tVal) && (tVal<=1.0+threshold))
		{
			if(projectionInTriangle(t1, intersection)) 
				return true;
		}
	}
	else
	{
		if(intersectInPlane(t1, nodes2[2]->coordinate(), nodes2[0]->coordinate(), threshold))
			return true;
	}
	if(plane1.distanceToPoint(nodes2[0]->coordinate()) < threshold)
	{
		if(projectionInTriangle(t1, nodes2[0]->coordinate())) 
			return true;
	}
	if(plane1.distanceToPoint(nodes2[1]->coordinate()) < threshold)
	{
		if(projectionInTriangle(t1, nodes2[1]->coordinate())) 
			return true;
	}
	if(plane1.distanceToPoint(nodes2[2]->coordinate()) < threshold)
	{
		if(projectionInTriangle(t1, nodes2[2]->coordinate())) 
			return true;
	}
	computationOK = plane2.line_plane_intersection(nodes1[1]->coordinate()-nodes1[0]->coordinate(),nodes1[0]->coordinate(),intersection,tVal);
	if(computationOK)
	{
		if((0.0-threshold<=tVal) && (tVal<=1.0+threshold))
		{
			if(projectionInTriangle(t2, intersection)) 
				return true;
		}
	}
	else
	{
		if(intersectInPlane(t2, nodes1[0]->coordinate(), nodes1[1]->coordinate(), threshold))
			return true;
	}
	computationOK = plane2.line_plane_intersection(nodes1[2]->coordinate()-nodes1[1]->coordinate(),nodes1[1]->coordinate(),intersection,tVal);
	if(computationOK)
	{
		if((0.0-threshold<=tVal) && (tVal<=1.0+threshold))
		{
			if(projectionInTriangle(t2, intersection)) 
				return true;
		}
	}
	else
	{
		if(intersectInPlane(t2, nodes1[1]->coordinate(), nodes1[2]->coordinate(), threshold))
			return true;
	}
	computationOK = plane2.line_plane_intersection(nodes1[0]->coordinate()-nodes1[2]->coordinate(),nodes1[2]->coordinate(),intersection,tVal);
	if(computationOK)
	{
		if((0.0-threshold<=tVal) && (tVal<=1.0+threshold))
		{
			if(projectionInTriangle(t2, intersection)) 
				return true;
		}
	}
	else
	{
		if(intersectInPlane(t2, nodes1[2]->coordinate(), nodes1[0]->coordinate(), threshold))
			return true;
	}
	if(plane2.distanceToPoint(nodes1[0]->coordinate()) < threshold)
	{
		if(projectionInTriangle(t2, nodes1[0]->coordinate())) 
			return true;
	}
	if(plane2.distanceToPoint(nodes1[1]->coordinate()) < threshold)
	{
		if(projectionInTriangle(t2, nodes1[1]->coordinate())) 
			return true;
	}
	if(plane2.distanceToPoint(nodes1[2]->coordinate()) < threshold)
	{
		if(projectionInTriangle(t2, nodes1[2]->coordinate())) 
			return true;
	}
	return false;
}

bool MinimumDistance::intersectInPlane(g_Element & t, g_Vector & p1, g_Vector & p2, double threshold)
{
	if(minimumDistance(t.nodes()[0]->coordinate(), t.nodes()[1]->coordinate(), p1, p2, false) < threshold)
		return true;
	if(minimumDistance(t.nodes()[1]->coordinate(), t.nodes()[2]->coordinate(), p1, p2, false) < threshold)
		return true;
	if(minimumDistance(t.nodes()[2]->coordinate(), t.nodes()[0]->coordinate(), p1, p2, false) < threshold)
		return true;
	return false;
}

double MinimumDistance::boundingBoxDist(g_Vector& p1, g_Vector& p2, g_Vector& p3, g_Vector& p4)
{
	double dist, minDist;
	bool intersect = true;

	if(max(p1.x(), p2.x()) < min(p3.x(), p4.x())) intersect = false;
	else if(max(p3.x(), p4.x()) < min(p1.x(), p2.x())) intersect = false;
	else if(max(p1.y(), p2.y()) < min(p3.y(), p4.y())) intersect = false;
	else if(max(p3.y(), p4.y()) < min(p1.y(), p2.y())) intersect = false;
	else if(max(p1.z(), p2.z()) < min(p3.z(), p4.z())) intersect = false;
	else if(max(p3.z(), p4.z()) < min(p1.z(), p2.z())) intersect = false;

	if(intersect) return 0.0;
	
	minDist = -1;
	dist = fabs(min(p3.x(), p4.x()) - max(p1.x(), p2.x()));
	minDist = dist;
	dist = fabs(min(p1.x(), p2.x()) - max(p3.x(), p4.x()));
	if(dist < minDist) minDist = dist;

	dist = fabs(min(p3.y(), p4.y()) - max(p1.y(), p2.y()));
	if(dist < minDist) minDist = dist;
	dist = fabs(min(p1.y(), p2.y()) - max(p3.y(), p4.y()));
	if(dist < minDist) minDist = dist;
	
	dist = fabs(min(p3.z(), p4.z()) - max(p1.z(), p2.z()));
	if(dist < minDist) minDist = dist;
	dist = fabs(min(p1.z(), p2.z()) - max(p3.z(), p4.z()));
	if(dist < minDist) minDist = dist;

	return minDist;
}

void MinimumDistance::mdGrad_p1(g_Vector &p1, g_Vector &p2, g_Vector &p3, g_Vector &p4, g_Vector &grad)
{
	g_Vector u, v, w, z;
	double d2121;
	double d4343, d4321, d1343, d1321;
	double numer_a, denom_a;
	double rcpDenom, rcpD4343;
	double mu_a, mu_b;
	double cl_a, cl_b;
	g_Vector del_numer, del_denom;
	g_Vector del_mu_a, del_mu_b;
	g_Vector del_cl_a, del_cl_b;
	g_Vector ones(1.0, 1.0, 1.0);
	g_Vector displace;
	double dMD, rcpDMD;

	u = p2 - p1;
	v = p4 - p3;
	w = p1 - p3;
	z = p2 + p3 - (p1 * 2);

	d2121 = u.Dot(u);
	d4343 = v.Dot(v);
	d4321 = v.Dot(u);
	d1343 = w.Dot(v);
	d1321 = w.Dot(u);
	if(fabs(d4343) < MD_DEFAULT_EPSILON) d4343 = MD_DEFAULT_EPSILON;
	rcpD4343 = 1.0 / d4343;

	numer_a = d1343 * d4321 - d1321 * d4343;
	denom_a = d2121 * d4343 - d4321 * d4321;
	if(fabs(denom_a) < MD_DEFAULT_EPSILON) denom_a = MD_DEFAULT_EPSILON;
	rcpDenom = 1.0 / denom_a;
	mu_a = numer_a * rcpDenom;
	mu_b = (d1343 + mu_a * d4321) * rcpD4343;

	cl_a = max(0.0, mu_a);
	cl_a = min(1.0, cl_a);
	cl_b = max(0.0, mu_b);
	cl_b = min(1.0, cl_b);
	displace = (p1 + u * cl_a) - (p3 + v * cl_b);
	dMD = displace.Length();
	if (dMD < MD_DEFAULT_EPSILON)
		dMD = MD_DEFAULT_EPSILON;
	rcpDMD = 1.0 / dMD;

	del_numer = (v * v.Dot(z)) - z * d4343;
	del_denom = (v * d4321 - u * d4343) * 2;

	del_mu_a = (del_numer * denom_a - del_denom * numer_a) * (rcpDenom * rcpDenom);
	del_mu_b = (v * (1.0 - mu_a) + del_mu_a * d4321) * rcpD4343;

	if (mu_a <= 0.0 || mu_a >= 1.0)
		del_cl_a.Set(0.0, 0.0, 0.0);
	else
		del_cl_a = del_mu_a;
	if (mu_b <= 0.0 || mu_b >= 1.0)
		del_cl_b.Set(0.0, 0.0, 0.0);
	else
		del_cl_b = del_mu_b;

	double * matrix = new double[9];
	double * addMatrix = new double[9];
	initializeIdentity(matrix);
	matrixMultiply(del_cl_a, u, addMatrix);
	matrixAdd(matrix, addMatrix, matrix);
	initializeIdentity(addMatrix);
	scalarMultiply(-cl_a, addMatrix);
	matrixAdd(matrix, addMatrix, matrix);
	matrixMultiply(del_cl_b, v, addMatrix);
	scalarMultiply(-1.0, addMatrix);
	matrixAdd(matrix, addMatrix, matrix);
	scalarMultiply(rcpDMD, addMatrix);
	grad = matrixMultiply(matrix, displace);

	delete [] matrix;
	delete [] addMatrix;
}


void  MinimumDistance::initializeIdentity(double *& matrix)
{
	matrix[0] = matrix[4] = matrix[8] = 1.0;
	matrix[1] = matrix[2] = matrix[3] = matrix[5] = matrix[6] = matrix[7] = 0.0;
}

void MinimumDistance::matrixMultiply(g_Vector vec1, g_Vector vec2, double *& result)
{
	result[0] = vec1.x()*vec2.x();
	result[1] = vec1.y()*vec2.x();
	result[2] = vec1.z()*vec2.x();
	result[3] = vec1.x()*vec2.y();
	result[4] = vec1.y()*vec2.y();
	result[5] = vec1.z()*vec2.y();
	result[6] = vec1.x()*vec2.z();
	result[7] = vec1.y()*vec2.z();
	result[8] = vec1.z()*vec2.z();
}

g_Vector MinimumDistance::matrixMultiply(double *& matrix, g_Vector vec)
{
	double x, y, z;
	
	x = matrix[0]*vec.x() + matrix[3] * vec.y() + matrix[6] * vec.z();
	y = matrix[1]*vec.x() + matrix[4] * vec.y() + matrix[7] * vec.z();
	z = matrix[2]*vec.x() + matrix[5] * vec.y() + matrix[8] * vec.z();
	
	g_Vector result(x, y, z);
	return result;
}

void MinimumDistance::transposeMatrix(double *& matrix)
{
	double swap;
	swap = matrix[1];
	matrix[1] = matrix[3];
	matrix[3] = swap;
	swap = matrix[2];
	matrix[2] = matrix[6];
	matrix[6] = swap;
	swap = matrix[5];
	matrix[5] = matrix[7];
	matrix[7] = swap;
}

void MinimumDistance::scalarMultiply(double scalarMultiplier, double *& result)
{
	int i;
	for(i = 0; i < 9; i++) result[i] = scalarMultiplier * result[i];
}

void MinimumDistance::matrixAdd(double *& mat1, double *& mat2, double *& result)
{
	int i;
	for(i = 0; i < 9; i++) result[i] = mat1[i] + mat2[i];
}

double MinimumDistance::evaluateSegments()
{
	if(edges == NULL) return 0;
	else
	{
		double stress, dMD2;
		g_Node p1, p2, p3, p4;
		int i, j, numEdges;
		stress = 0.0;
		//compute MD energy for all pairs of edges
		numEdges = edges->numberOfItems();
		for(i = 0; i < numEdges; i++)
		{
			p1 = (*edges)[i]->firstNode();
			p2 = (*edges)[i]->lastNode();
			for(j = i+1; j < numEdges; j++)
			{
				p3 = (*edges)[j]->firstNode();
				p4 = (*edges)[j]->lastNode();
				if((p1 == p3) || (p1 == p4) || (p2 == p3) || (p2 == p4)) continue;
				dMD2 = minimumDistance(p1.coordinate(), p2.coordinate(), p3.coordinate(), p4.coordinate());
				if (dMD2 < MD_DEFAULT_EPSILON)
					dMD2 = MD_DEFAULT_EPSILON;
				stress += 1.0 / dMD2;
			}
		}
		return stress;
	}
}
	
void MinimumDistance::gradientSegments(double*& gradient)
{
	if(edges == NULL) return;
	else
	{
		int i, j, numEdges, sampleSize, targetDim;
		double dist3;
		g_Node p1, p2, p3, p4;
		g_Vector grad;

		sampleSize = nodes->numberOfItems();
		targetDim = 3;

		//compute gradient of MD energy for all pairs of edges
		for(i = 0; i < sampleSize * targetDim; i++) gradient[i] = 0.0;

		numEdges = edges->numberOfItems();
		for(i = 0; i < numEdges; i++)
		{
			p1 = (*edges)[i]->firstNode();
			p2 = (*edges)[i]->lastNode();
			for(j = i+1; j < numEdges; j++)
			{
				p3 = (*edges)[j]->firstNode();
				p4 = (*edges)[j]->lastNode();
				if((p1 == p3) || (p1 == p4) || (p2 == p3) || (p2 == p4)) continue;
				MinimumDistance::mdGrad_p1(p1.coordinate(), p2.coordinate(), p3.coordinate(), p4.coordinate(), grad);
				dist3 = (max(pow(MinimumDistance::minimumDistance(p1.coordinate(), p2.coordinate(), 
					p3.coordinate(), p4.coordinate(), false), 3.0),MD_DEFAULT_EPSILON));
				grad = -2.0 * grad/dist3;
				if((grad.x() != grad.x()) || (grad.y() != grad.y()) || (grad.z() != grad.z()))
				{
					cout<<"PROBLEM IN "<<i<<" "<<j<<endl;
				}
				gradient[(p1.id()-1)*targetDim] += grad.x();
				gradient[(p1.id()-1)*targetDim+1] += grad.y();
				gradient[(p1.id()-1)*targetDim+2] += grad.z();

				MinimumDistance::mdGrad_p1(p2.coordinate(), p1.coordinate(), p3.coordinate(), p4.coordinate(), grad);
				grad = -2.0 * grad/dist3;
				if((grad.x() != grad.x()) || (grad.y() != grad.y()) || (grad.z() != grad.z()))
				{
					cout<<"PROBLEM IN "<<i<<" "<<j<<endl;
				}
				gradient[(p2.id()-1)*targetDim] += grad.x();
				gradient[(p2.id()-1)*targetDim+1] += grad.y();
				gradient[(p2.id()-1)*targetDim+2] += grad.z();

				MinimumDistance::mdGrad_p1(p3.coordinate(), p4.coordinate(), p1.coordinate(), p2.coordinate(), grad);
				grad = -2.0 * grad/dist3;
				if((grad.x() != grad.x()) || (grad.y() != grad.y()) || (grad.z() != grad.z()))
				{
					cout<<"PROBLEM IN "<<i<<" "<<j<<endl;
				}
				gradient[(p3.id()-1)*targetDim] += grad.x();
				gradient[(p3.id()-1)*targetDim+1] += grad.y();
				gradient[(p3.id()-1)*targetDim+2] += grad.z();

				MinimumDistance::mdGrad_p1(p4.coordinate(), p3.coordinate(), p1.coordinate(), p2.coordinate(), grad);
				grad = -2.0 * grad/dist3;
				if((grad.x() != grad.x()) || (grad.y() != grad.y()) || (grad.z() != grad.z()))
				{
					cout<<"PROBLEM IN "<<i<<" "<<j<<endl;
				}
				gradient[(p4.id()-1)*targetDim] += grad.x();
				gradient[(p4.id()-1)*targetDim+1] += grad.y();
				gradient[(p4.id()-1)*targetDim+2] += grad.z();
			}
		}
	}
}

double MinimumDistance::evaluateSegments(int id)
{
	if(edges == NULL) return 0;
	else 
	{
		int i, j, numEdges;
		double stress, dMD2;
		g_PEdgeContainer localCont;
		g_ElementContainer elems;
		g_PtrAdapter<g_PEdge *> adapt(localCont);
		g_Node p1, p2, p3, p4;
		g_Node * current;
		stress = 0.0;

		current = (*nodes)[id];
		//find the edges adjacent to the node id:
		elems = current->elements();
		for(i = 0; i < (int)elems.numberOfItems(); i++)
		{
			g_PEdgeContainer edgeCont = elems[i]->pEdges();
			if((edgeCont[0]->firstNode().id() == current->id()) || (edgeCont[0]->lastNode().id() == current->id()))
				localCont.insert(edgeCont[0]);
			else delete edgeCont[0];
			if((edgeCont[1]->firstNode().id() == current->id()) || (edgeCont[1]->lastNode().id() == current->id()))
				localCont.insert(edgeCont[1]);
			else delete edgeCont[1];
			if((edgeCont[2]->firstNode().id() == current->id()) || (edgeCont[2]->lastNode().id() == current->id()))
				localCont.insert(edgeCont[2]);
			else delete edgeCont[2];
		}
		adapt.removeDuplicates();
		//compute MD energy for all pairs of edges
		numEdges = edges->numberOfItems();
		for(i = 0; i < (int)localCont.numberOfItems(); i++)
		{
			p1 = localCont[i]->firstNode();
			p2 = localCont[i]->lastNode();
			for(j = 0; j < numEdges; j++)
			{
				p3 = (*edges)[j]->firstNode();
				p4 = (*edges)[j]->lastNode();
				if((p1 == p3) || (p1 == p4) || (p2 == p3) || (p2 == p4)) continue;
				dMD2 = MinimumDistance::minimumDistance(p1.coordinate(), p2.coordinate(), p3.coordinate(), p4.coordinate());
				if (dMD2 < MD_DEFAULT_EPSILON)
					dMD2 = MD_DEFAULT_EPSILON;
				stress += 1.0 / dMD2;
			}
		}

		for(i = 0; i < (int)localCont.numberOfItems(); i++) delete localCont[i];
		localCont.clear();

		return stress;
	}
}

void MinimumDistance::gradientSegments(int id, double*& gradient)
{
	if(edges == NULL) return;
	else
	{
		int i, j, numEdges, targetDim;
		double dist3;
		g_PEdgeContainer localCont;
		g_ElementContainer elems;
		g_Node * current;
		g_PtrAdapter<g_PEdge *> adapt(localCont);
		g_Node p1, p2, p3, p4;
		g_Vector grad;

		targetDim =  3;

		//compute gradient of MD energy for all pairs of edges
		for(i = 0; i < targetDim; i++) gradient[i] = 0.0;
		
		current = (*nodes)[id];
		//find the edges adjacent to the node id:
		elems = current->elements();
		for(i = 0; i < (int)elems.numberOfItems(); i++)
		{
			g_PEdgeContainer edgeCont = elems[i]->pEdges();
			if((edgeCont[0]->firstNode().id() == current->id()) || (edgeCont[0]->lastNode().id() == current->id()))
				localCont.insert(edgeCont[0]);
			else delete edgeCont[0];
			if((edgeCont[1]->firstNode().id() == current->id()) || (edgeCont[1]->lastNode().id() == current->id()))
				localCont.insert(edgeCont[1]);
			else delete edgeCont[1];
			if((edgeCont[2]->firstNode().id() == current->id()) || (edgeCont[2]->lastNode().id() == current->id()))
				localCont.insert(edgeCont[2]);
			else delete edgeCont[2];
		}
		adapt.removeDuplicates();
		//compute MD energy for all pairs of edges
		numEdges = edges->numberOfItems();
		for(i = 0; i < (int)localCont.numberOfItems(); i++)
		{
			p1 = localCont[i]->firstNode();
			p2 = localCont[i]->lastNode();
			if(p2.id() == current->id()) 
			{
				p3 = p1;
				p1 = p2; 
				p2 = p3;
			}
			for(j = 0; j < numEdges; j++)
			{
				p3 = (*edges)[j]->firstNode();
				p4 = (*edges)[j]->lastNode();
				if((p1 == p3) || (p1 == p4) || (p2 == p3) || (p2 == p4)) continue;
				MinimumDistance::mdGrad_p1(p1.coordinate(), p2.coordinate(), p3.coordinate(), p4.coordinate(), grad);
				dist3 = (max(pow(MinimumDistance::minimumDistance(p1.coordinate(), p2.coordinate(), 
					p3.coordinate(), p4.coordinate(), false), 3.0),MD_DEFAULT_EPSILON));
				grad = -2.0 * grad/dist3;
				gradient[0] += grad.x();
				gradient[1] += grad.y();
				gradient[2] += grad.z();
			}
		}

		for(i = 0; i < (int)localCont.numberOfItems(); i++) delete localCont[i];
		localCont.clear();
	}
}

double MinimumDistance::getMaxStepSizeSegments()
{
	if(edges == NULL) return 0;
	else
	{
		double minimum, dMD2;
		g_Node p1, p2, p3, p4;
		int i, j, numEdges;
		minimum = -1.0;
		//compute MD energy for all pairs of edges
		numEdges = edges->numberOfItems();
		for(i = 0; i < numEdges; i++)
		{
			p1 = (*edges)[i]->firstNode();
			p2 = (*edges)[i]->lastNode();
#pragma omp parallel for private(j, p3, p4, dMD2)
			for(j = i+1; j < numEdges; j++)
			{
				p3 = (*edges)[j]->firstNode();
				p4 = (*edges)[j]->lastNode();
				if((p1 == p3) || (p1 == p4) || (p2 == p3) || (p2 == p4)) continue;
				if((minimum == -1) || (boundingBoxDist(p1.coordinate(), p2.coordinate(), p3.coordinate(), p4.coordinate()) < minimum))
				{
					dMD2 = minimumDistance(p1.coordinate(), p2.coordinate(), p3.coordinate(), p4.coordinate(), false);
#pragma omp critical
					{
						if((minimum == -1) || (dMD2 < minimum)) minimum = dMD2;
					}
				}
			}
		}
		return minimum/3.0;
	}
}

double MinimumDistance::getMaxStepSizeSegments(int id)
{
	if(edges == NULL) return 0;
	else
	{
		double minimum, dMD2;
		g_Node p1, p2, p3, p4;
		int i, j, numEdges, numElems;
		g_NodeContainer neighbors;
		g_PtrAdapter<g_Node *> adapter(neighbors);
		g_ElementContainer elements;
		minimum = -1.0;
		//compute MD energy for all pairs of edges
		numEdges = edges->numberOfItems();
		p1 = *((*nodes)[id]);
		elements = p1.elements();
		numElems = elements.numberOfItems();
		for(i = 0; i < numElems; i++)
		{
			if(*(elements[i]->nodes()[0]) != p1) neighbors.insert(elements[i]->nodes()[0]);
			if(*(elements[i]->nodes()[1]) != p1) neighbors.insert(elements[i]->nodes()[1]);
			if(*(elements[i]->nodes()[2]) != p1) neighbors.insert(elements[i]->nodes()[2]);
		}
		adapter.removeDuplicates();
		numElems = neighbors.numberOfItems();
		for(i = 0; i < numElems; i++)
		{
			p2 = *(neighbors[i]);
#pragma omp parallel for private(j, p3, p4, dMD2)
			for(j = 0; j < numEdges; j++)
			{
				p3 = (*edges)[j]->firstNode();
				p4 = (*edges)[j]->lastNode();
				if((p1 == p3) || (p1 == p4) || (p2 == p3) || (p2 == p4)) continue;
				if((minimum == -1) || (boundingBoxDist(p1.coordinate(), p2.coordinate(), p3.coordinate(), p4.coordinate()) < minimum))
				{
					dMD2 = minimumDistance(p1.coordinate(), p2.coordinate(), p3.coordinate(), p4.coordinate(), false);
#pragma omp critical
					{
						if((minimum == -1) || (dMD2 < minimum)) minimum = dMD2;
					}
				}
			}
		}
		return minimum/3.0;
	}
}


double MinimumDistance::evaluateTriangles()
{
	if(meshElements == NULL) return 0;
	else
	{
		int i, j, numElems;
		double stress;

		numElems = meshElements->numberOfItems();
		stress = 0.0;
		for(i = 0; i < numElems; i++)
		{
			for(j = i+1; j < numElems; j++)
			{
				if(((*meshElements)[i]->nodes()[0] == (*meshElements)[j]->nodes()[0]) || 
					((*meshElements)[i]->nodes()[0] == (*meshElements)[j]->nodes()[1]) || 
					((*meshElements)[i]->nodes()[0] == (*meshElements)[j]->nodes()[2]) || 
					((*meshElements)[i]->nodes()[1] == (*meshElements)[j]->nodes()[0]) || 
					((*meshElements)[i]->nodes()[1] == (*meshElements)[j]->nodes()[1]) || 
					((*meshElements)[i]->nodes()[1] == (*meshElements)[j]->nodes()[2]) || 
					((*meshElements)[i]->nodes()[2] == (*meshElements)[j]->nodes()[0]) || 
					((*meshElements)[i]->nodes()[2] == (*meshElements)[j]->nodes()[1]) || 
					((*meshElements)[i]->nodes()[2] == (*meshElements)[j]->nodes()[2]))
					continue;
				stress += 1.0/max(pow(minimumDistance((*(*meshElements)[i]), (*(*meshElements)[j])), 2), MD_DEFAULT_EPSILON);
			}
		}
		return stress;
	}
}

bool MinimumDistance::selfIntersection(double threshold)
{
	if(meshElements == NULL) return false;
	else
	{
		int i, j, numElems;
		bool returnVal = false;

		numElems = meshElements->numberOfItems();
		for(i = 0; i < numElems; i++)
		{
			for(j = i+1; j < numElems; j++)
			{
				if(((*meshElements)[i]->nodes()[0] == (*meshElements)[j]->nodes()[0]) || 
					((*meshElements)[i]->nodes()[0] == (*meshElements)[j]->nodes()[1]) || 
					((*meshElements)[i]->nodes()[0] == (*meshElements)[j]->nodes()[2]) || 
					((*meshElements)[i]->nodes()[1] == (*meshElements)[j]->nodes()[0]) || 
					((*meshElements)[i]->nodes()[1] == (*meshElements)[j]->nodes()[1]) || 
					((*meshElements)[i]->nodes()[1] == (*meshElements)[j]->nodes()[2]) || 
					((*meshElements)[i]->nodes()[2] == (*meshElements)[j]->nodes()[0]) || 
					((*meshElements)[i]->nodes()[2] == (*meshElements)[j]->nodes()[1]) || 
					((*meshElements)[i]->nodes()[2] == (*meshElements)[j]->nodes()[2]))
					continue;
				else if(intersect((*(*meshElements)[i]), (*(*meshElements)[j])), threshold) 
				{
					returnVal = true;
					return returnVal;
				}
			}
		}
		return returnVal;
	}
}

void MinimumDistance::gradientTriangles(double*& gradient)
{
	if(meshElements == NULL) return;
	else
	{
		int i, j, numElems, nNodes, index;
		nNodes = (int) nodes->numberOfItems();
		numElems = meshElements->numberOfItems();
		for(i = 0; i < nNodes*3; i++) gradient[i] = 0.0;
		for(i = 0; i < numElems; i++)
		{
			for(j = i+1; j < numElems; j++)
			{
				if(((*meshElements)[i]->nodes()[0] == (*meshElements)[j]->nodes()[0]) || 
					((*meshElements)[i]->nodes()[0] == (*meshElements)[j]->nodes()[1]) || 
					((*meshElements)[i]->nodes()[0] == (*meshElements)[j]->nodes()[2]) || 
					((*meshElements)[i]->nodes()[1] == (*meshElements)[j]->nodes()[0]) || 
					((*meshElements)[i]->nodes()[1] == (*meshElements)[j]->nodes()[1]) || 
					((*meshElements)[i]->nodes()[1] == (*meshElements)[j]->nodes()[2]) || 
					((*meshElements)[i]->nodes()[2] == (*meshElements)[j]->nodes()[0]) || 
					((*meshElements)[i]->nodes()[2] == (*meshElements)[j]->nodes()[1]) || 
					((*meshElements)[i]->nodes()[2] == (*meshElements)[j]->nodes()[2]))
					continue;
				minimumDistance((*(*meshElements)[i]), (*(*meshElements)[j]), &index);
				switch(index)
				{
				case 0:
					addSegmentGrad((*meshElements)[i]->nodes()[0], (*meshElements)[i]->nodes()[1], 
						(*meshElements)[j]->nodes()[0], (*meshElements)[j]->nodes()[1], gradient);
					break;
				case 1:
					addSegmentGrad((*meshElements)[i]->nodes()[0], (*meshElements)[i]->nodes()[1], 
						(*meshElements)[j]->nodes()[1], (*meshElements)[j]->nodes()[2], gradient);
					break;
				case 2:
					addSegmentGrad((*meshElements)[i]->nodes()[0], (*meshElements)[i]->nodes()[1], 
						(*meshElements)[j]->nodes()[2], (*meshElements)[j]->nodes()[0], gradient);
					break;
				case 3:
					addSegmentGrad((*meshElements)[i]->nodes()[1], (*meshElements)[i]->nodes()[2], 
						(*meshElements)[j]->nodes()[0], (*meshElements)[j]->nodes()[1], gradient);
					break;
				case 4:
					addSegmentGrad((*meshElements)[i]->nodes()[1], (*meshElements)[i]->nodes()[2], 
						(*meshElements)[j]->nodes()[1], (*meshElements)[j]->nodes()[2], gradient);
					break;
				case 5:
					addSegmentGrad((*meshElements)[i]->nodes()[1], (*meshElements)[i]->nodes()[2], 
						(*meshElements)[j]->nodes()[2], (*meshElements)[j]->nodes()[0], gradient);
					break;
				case 6:
					addSegmentGrad((*meshElements)[i]->nodes()[2], (*meshElements)[i]->nodes()[0], 
						(*meshElements)[j]->nodes()[0], (*meshElements)[j]->nodes()[1], gradient);
					break;
				case 7:
					addSegmentGrad((*meshElements)[i]->nodes()[2], (*meshElements)[i]->nodes()[0], 
						(*meshElements)[j]->nodes()[1], (*meshElements)[j]->nodes()[2], gradient);
					break;
				case 8:
					addSegmentGrad((*meshElements)[i]->nodes()[2], (*meshElements)[i]->nodes()[0], 
						(*meshElements)[j]->nodes()[2], (*meshElements)[j]->nodes()[0], gradient);
					break;
				case 9:
					addPointPlaneGrad((*meshElements)[i], (*meshElements)[j]->nodes()[0], gradient);
					break;
				case 10:
					addPointPlaneGrad((*meshElements)[i], (*meshElements)[j]->nodes()[1], gradient);
					break;
				case 11:
					addPointPlaneGrad((*meshElements)[i], (*meshElements)[j]->nodes()[2], gradient);
					break;
				case 12:
					addPointPlaneGrad((*meshElements)[j], (*meshElements)[i]->nodes()[0], gradient);
					break;
				case 13:
					addPointPlaneGrad((*meshElements)[j], (*meshElements)[i]->nodes()[1], gradient);
					break;
				case 14:
					addPointPlaneGrad((*meshElements)[j], (*meshElements)[i]->nodes()[2], gradient);
					break;
				}
			}
		}
	}
}

void MinimumDistance::addPointPlaneGrad(g_Element * t, g_Node * p, double *& gradient, int id)
{
	int index;
	double dist3;
	g_Vector grad, helpVec;
	double * matrix = new double[9];

	g_NodeContainer nodes = t->nodes();
	g_Vector normal = (nodes[1]->coordinate()-nodes[0]->coordinate()).Cross(
		nodes[2]->coordinate()-nodes[0]->coordinate());
	if(normal.Length() < MD_DEFAULT_EPSILON) 
		return;
	g_Plane plane(normal, nodes[0]->coordinate());
	dist3 = max(pow(plane.distanceToPoint(p->coordinate()), 3), MD_DEFAULT_EPSILON);

	if(id == -1 || id == 0)
	{
		grad = normal;
		grad = -2.0 * grad/dist3;
		if((grad.x() != grad.x()) || (grad.y() != grad.y()) || (grad.z() != grad.z()))
		{
			cout<<"PROBLEM"<<endl;
		}
		if(id == -1) index = p->id()-1;
		else index = 0;
		gradient[index*3] += grad.x();
		gradient[index*3+1] += grad.y();
		gradient[index*3+2] += grad.z();
	}

	if(id == -1 || id == 1)
	{
		helpVec = p->coordinate() - nodes[0]->coordinate();
		matrix[0] = 0.0; matrix[4] = 0.0; matrix[8] = 0.0;
		matrix[1] = nodes[1]->coordinate().z()-nodes[2]->coordinate().z(); matrix[3] = -matrix[1];
		matrix[2] = nodes[2]->coordinate().y()-nodes[1]->coordinate().y(); matrix[6] = -matrix[2];
		matrix[5] = nodes[1]->coordinate().x()-nodes[2]->coordinate().x(); matrix[7] = -matrix[5];
		grad = matrixMultiply(matrix, helpVec);
		grad = grad-normal;
		grad = -2.0 * grad/dist3;
		if((grad.x() != grad.x()) || (grad.y() != grad.y()) || (grad.z() != grad.z()))
		{
			cout<<"PROBLEM"<<endl;
		}
		if(id == -1) index = nodes[0]->id()-1;
		else index = 0;
		gradient[index*3] += grad.x();
		gradient[index*3+1] += grad.y();
		gradient[index*3+2] += grad.z();
	}

	if(id == -1 || id == 2)
	{
		helpVec = p->coordinate() - nodes[0]->coordinate();
		matrix[0] = 0.0; matrix[4] = 0.0; matrix[8] = 0.0;
		matrix[1] = nodes[2]->coordinate().z()-nodes[0]->coordinate().z(); matrix[3] = -matrix[1];
		matrix[2] = nodes[0]->coordinate().y()-nodes[2]->coordinate().y(); matrix[6] = -matrix[2];
		matrix[5] = nodes[2]->coordinate().x()-nodes[0]->coordinate().x(); matrix[7] = -matrix[5];
		grad = matrixMultiply(matrix, helpVec);
		grad = -2.0 * grad/dist3;
		if((grad.x() != grad.x()) || (grad.y() != grad.y()) || (grad.z() != grad.z()))
		{
			cout<<"PROBLEM"<<endl;
		}
		if(id == -1) index = nodes[1]->id()-1;
		else index = 0;
		gradient[index*3] += grad.x();
		gradient[index*3+1] += grad.y();
		gradient[index*3+2] += grad.z();
	}

	if(id == -1 || id == 3)
	{
		helpVec = p->coordinate() - nodes[0]->coordinate();
		matrix[0] = 0.0; matrix[4] = 0.0; matrix[8] = 0.0;
		matrix[1] = nodes[0]->coordinate().z()-nodes[1]->coordinate().z(); matrix[3] = -matrix[1];
		matrix[2] = nodes[1]->coordinate().y()-nodes[0]->coordinate().y(); matrix[6] = -matrix[2];
		matrix[5] = nodes[0]->coordinate().x()-nodes[1]->coordinate().x(); matrix[7] = -matrix[5];
		grad = matrixMultiply(matrix, helpVec);
		grad = -2.0 * grad/dist3;
		if((grad.x() != grad.x()) || (grad.y() != grad.y()) || (grad.z() != grad.z()))
		{
			cout<<"PROBLEM"<<endl;
		}
		if(id == -1) index = nodes[2]->id()-1;
		else index = 0;
		gradient[index*3] += grad.x();
		gradient[index*3+1] += grad.y();
		gradient[index*3+2] += grad.z();
	}

	delete [] matrix;
}

void MinimumDistance::addSegmentGrad(g_Node * p1, g_Node * p2, g_Node * p3, g_Node * p4, double *& gradient, int id)
{
	int index;
	double dist3;
	g_Vector grad;

	if((p1 == p3) || (p1 == p4) || (p2 == p3) || (p2 == p4)) return;
	dist3 = (max(pow(MinimumDistance::minimumDistance(p1->coordinate(), p2->coordinate(), 
		p3->coordinate(), p4->coordinate(), false), 3.0),MD_DEFAULT_EPSILON));

	if(id == -1 || id == 0)
	{
		MinimumDistance::mdGrad_p1(p1->coordinate(), p2->coordinate(), p3->coordinate(), p4->coordinate(), grad);	
		grad = -2.0 * grad/dist3;
		if((grad.x() != grad.x()) || (grad.y() != grad.y()) || (grad.z() != grad.z()))
		{
			cout<<"PROBLEM"<<endl;
		}
		if(id == -1) index = p1->id()-1;
		else index = 0;
		gradient[index*3] += grad.x();
		gradient[index*3+1] += grad.y();
		gradient[index*3+2] += grad.z();
	}

	if(id == -1 || id == 1)
	{
		MinimumDistance::mdGrad_p1(p2->coordinate(), p1->coordinate(), p3->coordinate(), p4->coordinate(), grad);
		grad = -2.0 * grad/dist3;
		if((grad.x() != grad.x()) || (grad.y() != grad.y()) || (grad.z() != grad.z()))
		{
			cout<<"PROBLEM"<<endl;
		}
		if(id == -1) index = p2->id()-1;
		else index = 0;
		gradient[index*3] += grad.x();
		gradient[index*3+1] += grad.y();
		gradient[index*3+2] += grad.z();
	}

	if(id == -1 || id == 2)
	{
		MinimumDistance::mdGrad_p1(p3->coordinate(), p4->coordinate(), p1->coordinate(), p2->coordinate(), grad);
		grad = -2.0 * grad/dist3;
		if((grad.x() != grad.x()) || (grad.y() != grad.y()) || (grad.z() != grad.z()))
		{
			cout<<"PROBLEM"<<endl;
		}
		if(id == -1) index = p3->id()-1;
		else index = 0;
		gradient[index*3] += grad.x();
		gradient[index*3+1] += grad.y();
		gradient[index*3+2] += grad.z();
	}

	if(id == -1 || id == 3)
	{
		MinimumDistance::mdGrad_p1(p4->coordinate(), p3->coordinate(), p1->coordinate(), p2->coordinate(), grad);
		grad = -2.0 * grad/dist3;
		if((grad.x() != grad.x()) || (grad.y() != grad.y()) || (grad.z() != grad.z()))
		{
			cout<<"PROBLEM"<<endl;
		}
		if(id == -1) index = p4->id()-1;
		else index = 0;
		gradient[index*3] += grad.x();
		gradient[index*3+1] += grad.y();
		gradient[index*3+2] += grad.z();
	}
}

double MinimumDistance::evaluateTriangles(int id)
{
	if(meshElements == NULL) return 0;
	else
	{
		int i, j, numElems;
		double stress;
		g_ElementContainer localElems = (*nodes)[id]->elements();
		numElems = meshElements->numberOfItems();
		stress = 0.0;
		for(i = 0; i < (int)localElems.numberOfItems(); i++)
		{
			for(j = 0; j < numElems; j++)
			{
				if(((*meshElements)[j]->nodes()[0] == localElems[i]->nodes()[0]) || 
					((*meshElements)[j]->nodes()[0] == localElems[i]->nodes()[1]) || 
					((*meshElements)[j]->nodes()[0] == localElems[i]->nodes()[2]) || 
					((*meshElements)[j]->nodes()[1] == localElems[i]->nodes()[0]) || 
					((*meshElements)[j]->nodes()[1] == localElems[i]->nodes()[1]) || 
					((*meshElements)[j]->nodes()[1] == localElems[i]->nodes()[2]) || 
					((*meshElements)[j]->nodes()[2] == localElems[i]->nodes()[0]) || 
					((*meshElements)[j]->nodes()[2] == localElems[i]->nodes()[1]) || 
					((*meshElements)[j]->nodes()[2] == localElems[i]->nodes()[2]))
					continue;
				stress += 1.0/max(pow(minimumDistance((*(*meshElements)[j]), (*localElems[i])), 2), MD_DEFAULT_EPSILON);
			}
		}
		return stress;
	}
}

void MinimumDistance::gradientTriangles(int id, double*& gradient)
{
	if(meshElements == NULL) return;
	else
	{
		int i, j, numElems, index, localIdx, helpIdx;
		gradient[0] = 0.0; gradient[1] = 0.0; gradient[2] = 0.0;
		g_ElementContainer localElems = (*nodes)[id]->elements();
		numElems = meshElements->numberOfItems();
		for(i = 0; i < (int)localElems.numberOfItems(); i++)
		{
			if(localElems[i]->nodes()[0] == (*nodes)[id]) localIdx = 0;
			else if(localElems[i]->nodes()[1] == (*nodes)[id]) localIdx = 1;
			else localIdx = 2;

			for(j = 0; j < numElems; j++)
			{
				if(((*meshElements)[j]->nodes()[0] == localElems[i]->nodes()[0]) || 
					((*meshElements)[j]->nodes()[0] == localElems[i]->nodes()[1]) || 
					((*meshElements)[j]->nodes()[0] == localElems[i]->nodes()[2]) || 
					((*meshElements)[j]->nodes()[1] == localElems[i]->nodes()[0]) || 
					((*meshElements)[j]->nodes()[1] == localElems[i]->nodes()[1]) || 
					((*meshElements)[j]->nodes()[1] == localElems[i]->nodes()[2]) || 
					((*meshElements)[j]->nodes()[2] == localElems[i]->nodes()[0]) || 
					((*meshElements)[j]->nodes()[2] == localElems[i]->nodes()[1]) || 
					((*meshElements)[j]->nodes()[2] == localElems[i]->nodes()[2]))
					continue;
				minimumDistance((*(*meshElements)[j]), (*localElems[i]), &index);
				switch(index)
				{
				case 0:
					if(localIdx != 2)
					{
						if(localIdx == 0) helpIdx = 2;
						if(localIdx == 1) helpIdx = 3;
						addSegmentGrad((*meshElements)[j]->nodes()[0], (*meshElements)[j]->nodes()[1], 
							localElems[i]->nodes()[0], localElems[i]->nodes()[1], gradient, helpIdx);
					}
					break;
				case 1:
					if(localIdx != 0)
					{
						if(localIdx == 1) helpIdx = 2;
						if(localIdx == 2) helpIdx = 3;
						addSegmentGrad((*meshElements)[j]->nodes()[0], (*meshElements)[j]->nodes()[1], 
							localElems[i]->nodes()[1], localElems[i]->nodes()[2], gradient, helpIdx);
					}
					break;
				case 2:
					if(localIdx != 1)
					{
						if(localIdx == 2) helpIdx = 2;
						if(localIdx == 0) helpIdx = 3;
						addSegmentGrad((*meshElements)[j]->nodes()[0], (*meshElements)[j]->nodes()[1], 
							localElems[i]->nodes()[2], localElems[i]->nodes()[0], gradient, helpIdx);
					}
					break;
				case 3:
					if(localIdx != 2)
					{
						if(localIdx == 0) helpIdx = 2;
						if(localIdx == 1) helpIdx = 3;
						addSegmentGrad((*meshElements)[j]->nodes()[1], (*meshElements)[j]->nodes()[2], 
							localElems[i]->nodes()[0], localElems[i]->nodes()[1], gradient, helpIdx);
					}
					break;
				case 4:
					if(localIdx != 0)
					{
						if(localIdx == 1) helpIdx = 2;
						if(localIdx == 2) helpIdx = 3;
						addSegmentGrad((*meshElements)[j]->nodes()[1], (*meshElements)[j]->nodes()[2], 
							localElems[i]->nodes()[1], localElems[i]->nodes()[2], gradient, helpIdx);
					}
					break;
				case 5:
					if(localIdx != 1)
					{
						if(localIdx == 2) helpIdx = 2;
						if(localIdx == 0) helpIdx = 3;
						addSegmentGrad((*meshElements)[j]->nodes()[1], (*meshElements)[j]->nodes()[2], 
							localElems[i]->nodes()[2], localElems[i]->nodes()[0], gradient, helpIdx);
					}
					break;
				case 6:
					if(localIdx != 2)
					{
						if(localIdx == 0) helpIdx = 2;
						if(localIdx == 1) helpIdx = 3;
						addSegmentGrad((*meshElements)[j]->nodes()[2], (*meshElements)[j]->nodes()[0], 
							localElems[i]->nodes()[0], localElems[i]->nodes()[1], gradient, helpIdx);
					}
					break;
				case 7:
					if(localIdx != 0)
					{
						if(localIdx == 1) helpIdx = 2;
						if(localIdx == 2) helpIdx = 3;
						addSegmentGrad((*meshElements)[j]->nodes()[2], (*meshElements)[j]->nodes()[0], 
							localElems[i]->nodes()[1], localElems[i]->nodes()[2], gradient, helpIdx);
					}
					break;
				case 8:
					if(localIdx != 1)
					{
						if(localIdx == 2) helpIdx = 2;
						if(localIdx == 0) helpIdx = 3;
						addSegmentGrad((*meshElements)[j]->nodes()[2], (*meshElements)[j]->nodes()[0], 
							localElems[i]->nodes()[2], localElems[i]->nodes()[0], gradient, helpIdx);
					}
					break;
				case 9:
					if(localIdx == 0)
						addPointPlaneGrad((*meshElements)[j], localElems[i]->nodes()[0], gradient, 0);
					break;
				case 10:
					if(localIdx == 1)
						addPointPlaneGrad((*meshElements)[j], localElems[i]->nodes()[1], gradient, 0);
					break;
				case 11:
					if(localIdx == 2)
						addPointPlaneGrad((*meshElements)[j], localElems[i]->nodes()[2], gradient, 0);
					break;
				case 12:
					addPointPlaneGrad(localElems[i], (*meshElements)[j]->nodes()[0], gradient, localIdx+1);
					break;
				case 13:
					addPointPlaneGrad(localElems[i], (*meshElements)[j]->nodes()[1], gradient, localIdx+1);
					break;
				case 14:
					addPointPlaneGrad(localElems[i], (*meshElements)[j]->nodes()[2], gradient, localIdx+1);
					break;
				}
			}
		}
	}
}

double MinimumDistance::getMaxStepSizeTriangles()
{
	if(meshElements == NULL) return 0;
	else
	{
		int i, j, numElems;
		double minimum, dist;

		numElems = meshElements->numberOfItems();
		minimum = -1.0;
		for(i = 0; i < numElems; i++)
		{
#pragma omp parallel for private (j, dist)
			for(j = i+1; j < numElems; j++)
			{
				if(((*meshElements)[j]->nodes()[0] == (*meshElements)[i]->nodes()[0]) || 
					((*meshElements)[j]->nodes()[0] == (*meshElements)[i]->nodes()[1]) || 
					((*meshElements)[j]->nodes()[0] == (*meshElements)[i]->nodes()[2]) || 
					((*meshElements)[j]->nodes()[1] == (*meshElements)[i]->nodes()[0]) || 
					((*meshElements)[j]->nodes()[1] == (*meshElements)[i]->nodes()[1]) || 
					((*meshElements)[j]->nodes()[1] == (*meshElements)[i]->nodes()[2]) || 
					((*meshElements)[j]->nodes()[2] == (*meshElements)[i]->nodes()[0]) || 
					((*meshElements)[j]->nodes()[2] == (*meshElements)[i]->nodes()[1]) || 
					((*meshElements)[j]->nodes()[2] == (*meshElements)[i]->nodes()[2]))
					continue;

				dist = minimumDistance((*(*meshElements)[i]), (*(*meshElements)[j]));
#pragma omp critical
				{
					if((minimum == -1.0) || (dist < minimum)) minimum = dist;
				}
			}
		}
		return minimum / 3.0;
	}
}

double MinimumDistance::getMaxStepSizeTriangles(int id)
{
	if(meshElements == NULL) return 0;
	else
	{
		int i, j, numElems;
		double minimum, dist;
		g_ElementContainer localElems = (*nodes)[id]->elements();

		numElems = meshElements->numberOfItems();
		minimum = -1.0;
		for(i = 0; i < (int)localElems.numberOfItems(); i++)
		{
#pragma omp parallel for private(j, dist)
			for(j = 0; j < numElems; j++)
			{
				if(((*meshElements)[j]->nodes()[0] == localElems[i]->nodes()[0]) || 
					((*meshElements)[j]->nodes()[0] == localElems[i]->nodes()[1]) || 
					((*meshElements)[j]->nodes()[0] == localElems[i]->nodes()[2]) || 
					((*meshElements)[j]->nodes()[1] == localElems[i]->nodes()[0]) || 
					((*meshElements)[j]->nodes()[1] == localElems[i]->nodes()[1]) || 
					((*meshElements)[j]->nodes()[1] == localElems[i]->nodes()[2]) || 
					((*meshElements)[j]->nodes()[2] == localElems[i]->nodes()[0]) || 
					((*meshElements)[j]->nodes()[2] == localElems[i]->nodes()[1]) || 
					((*meshElements)[j]->nodes()[2] == localElems[i]->nodes()[2]))
					continue;
				dist = minimumDistance((*(*meshElements)[j]), (*localElems[i]));
#pragma omp critical
				{
					if((minimum == -1) || (dist < minimum)) minimum = dist;
				}
			}
		}
		return minimum / 3.0;
	}
}

void MinimumDistance::updateNN()
{
	if(meshElements == NULL) return;
	else
	{
#if MEASURE_TIME
		int t = (int)time(0);
#endif
		int i, j, numElems;
		double distance; 

		nearestNeighbor.clear();
		distNN.clear();
		numElems = meshElements->numberOfItems();
		nearestNeighbor.resize(numElems);
		distNN.resize(numElems);
#pragma omp parallel for private(i)
		for(i = 0; i < numElems; i++) nearestNeighbor[i] = -1;
#pragma omp parallel for private(i)
		for(i = 0; i < numElems; i++) distNN[i] = -1;

		for(i = 0; i < numElems; i++)
		{
#pragma omp parallel for private(j, distance)
			for(j = i+1; j < numElems; j++)
			{
				if(((*meshElements)[i]->nodes()[0] == (*meshElements)[j]->nodes()[0]) || 
					((*meshElements)[i]->nodes()[0] == (*meshElements)[j]->nodes()[1]) || 
					((*meshElements)[i]->nodes()[0] == (*meshElements)[j]->nodes()[2]) || 
					((*meshElements)[i]->nodes()[1] == (*meshElements)[j]->nodes()[0]) || 
					((*meshElements)[i]->nodes()[1] == (*meshElements)[j]->nodes()[1]) || 
					((*meshElements)[i]->nodes()[1] == (*meshElements)[j]->nodes()[2]) || 
					((*meshElements)[i]->nodes()[2] == (*meshElements)[j]->nodes()[0]) || 
					((*meshElements)[i]->nodes()[2] == (*meshElements)[j]->nodes()[1]) || 
					((*meshElements)[i]->nodes()[2] == (*meshElements)[j]->nodes()[2]))
					continue;
				distance = minimumDistance((*(*meshElements)[i]), (*(*meshElements)[j]));
#pragma omp critical
				if((distNN[i] == -1) || (distance < distNN[i]))
				{
					distNN[i] = distance;
					nearestNeighbor[i] = j;
				}
#pragma omp critical
				if((distNN[j] == -1) || (distance < distNN[j]))
				{
					distNN[j] = distance;
					nearestNeighbor[j] = i;
				}
			}
		}
#if MEASURE_TIME
		cout<<"Time to execute MinimumDistance::updateNN() "<<(int)time(0)-t<<" s"<<endl;
#endif
	}
}

void MinimumDistance::updateNNKdSpherical()
{
	if(meshElements == NULL) return;
	else
	{
#if MEASURE_TIME
		int t = (int)time(0);
#endif
		int i, j, k, numElems, numNodes;
		double distance, radius; 

		nearestNeighbor.clear();
		distNN.clear();
		numElems = meshElements->numberOfItems();
		numNodes = nodes->numberOfItems();
		nearestNeighbor.resize(numElems);
		distNN.resize(numElems);
#pragma omp parallel for private(i)
		for(i = 0; i < numElems; i++) nearestNeighbor[i] = -1;
#pragma omp parallel for private(i)
		for(i = 0; i < numElems; i++) distNN[i] = -1;

		//Find the longest edge length:
		for(i = 0; i < numElems; i++)
		{
			g_Element * elem = (* meshElements)[i];
			for(j = 0; j < 3; j++)
			{
				distance = elem->nodes()[j]->coordinate().DistanceTo(elem->nodes()[(j+1)%3]->coordinate());
				if(((i == 0)&&(j == 0)) || (distance > radius))
					radius = distance;
			}
		}
		radius = 3.0/2.0 * radius;
		//Build a kd-tree for the nodes
#if USE_ANN
		ANNpointArray		dataPts;
		ANNpoint			queryPt;
		ANNidxArray			nnIdx;	
		ANNkd_tree*			kdTree;	
		queryPt = annAllocPt(3);				
		dataPts = annAllocPts(numNodes, 3);		
		nnIdx = new ANNidx[numNodes];	
		for(i = 0; i < numNodes; i++)
		{
			dataPts[i][0] = (*nodes)[i]->coordinate().x();
			dataPts[i][1] = (*nodes)[i]->coordinate().y();
			dataPts[i][2] = (*nodes)[i]->coordinate().z();
		}
		kdTree = new ANNkd_tree(dataPts, numNodes, 3);
#else
		kdTree * kd = kdTree_create();
		for(i = 0; i < numNodes; i++)
		{
			NODE * pNode = (NODE*)malloc(sizeof(NODE));
			pNode->coordinate.x = (*nodes)[i]->coordinate().x();
			pNode->coordinate.y = (*nodes)[i]->coordinate().y();
			pNode->coordinate.z = (*nodes)[i]->coordinate().z();
			pNode->local_id = i; 
			kdTree_insert(kd, pNode); 			
		}
#endif

		for(i = 0; i < numNodes; i++)
		{
			//Do spherical range searching and insert triangles into neighbors
#if USE_ANN
			queryPt[0] = (*nodes)[i]->coordinate().x();
			queryPt[1] = (*nodes)[i]->coordinate().y();
			queryPt[2] = (*nodes)[i]->coordinate().z();
			int numRange = kdTree->annkFRSearch(queryPt, radius*radius, numNodes, nnIdx);
			g_ElementContainer contToCheck;
			g_PtrAdapter<g_Element *> adapt(contToCheck);
			for(j = 0; j < numRange; j++)
			{
				contToCheck.append((*nodes)[nnIdx[j]]->elements());
			}
#else
			COORDINATE p;
			p.x = (*nodes)[i]->coordinate().x();
			p.y = (*nodes)[i]->coordinate().y();
			p.z = (*nodes)[i]->coordinate().z();
			nodeQueue * neighbors = kdTree_spherical_range_search(kd, p, radius);
			g_ElementContainer contToCheck;
			g_PtrAdapter<g_Element *> adapt(contToCheck);
			for(j = 0; j < nodeQueue_numberOfItems(neighbors); j++)
				contToCheck.append((*nodes)[nodeQueue_remove(neighbors)->local_id]->elements());
#endif
			adapt.removeDuplicates();
			int contNum = contToCheck.numberOfItems();
			//Find nearest neighbor
			int numElemNeighbors = (*nodes)[i]->elements().numberOfItems();
			for(j = 0; j < numElemNeighbors; j++)
			{
				int id1 = (*nodes)[i]->elements()[j]->id()-1;
#pragma omp parallel for private(k, distance)
				for(k = 0; k < contNum; k++)
				{
					int id2 = contToCheck[k]->id()-1;
					if(id2 <= id1) continue;
					if(((*meshElements)[id1]->nodes()[0] == (*meshElements)[id2]->nodes()[0]) || 
						((*meshElements)[id1]->nodes()[0] == (*meshElements)[id2]->nodes()[1]) || 
						((*meshElements)[id1]->nodes()[0] == (*meshElements)[id2]->nodes()[2]) || 
						((*meshElements)[id1]->nodes()[1] == (*meshElements)[id2]->nodes()[0]) || 
						((*meshElements)[id1]->nodes()[1] == (*meshElements)[id2]->nodes()[1]) || 
						((*meshElements)[id1]->nodes()[1] == (*meshElements)[id2]->nodes()[2]) || 
						((*meshElements)[id1]->nodes()[2] == (*meshElements)[id2]->nodes()[0]) || 
						((*meshElements)[id1]->nodes()[2] == (*meshElements)[id2]->nodes()[1]) || 
						((*meshElements)[id1]->nodes()[2] == (*meshElements)[id2]->nodes()[2]))
						continue;
					distance = minimumDistance((*(*meshElements)[id1]), (*(*meshElements)[id2]));
#pragma omp critical
					if((distNN[id1] == -1) || (distance < distNN[id1]))
					{
						distNN[id1] = distance;
						nearestNeighbor[id1] = id2;
					}
#pragma omp critical
					if((distNN[id2] == -1) || (distance < distNN[id2]))
					{
						distNN[id2] = distance;
						nearestNeighbor[id2] = id1;
					}
				}
			}
		}
		//Delete kd-tree
#if USE_ANN	
		delete [] nnIdx;
		delete kdTree;
		annClose();	
#else		
		kdTree_destroy(kd);
#endif

#if MEASURE_TIME
		cout<<"Time to execute MinimumDistance::updateNN() "<<(int)time(0)-t<<" s"<<endl;
#endif
	}
}

void MinimumDistance::updateNNKd()
{
	if(meshElements == NULL) return;
	else
	{
#if MEASURE_TIME
		int t = (int)time(0);
#endif
		int i, j, k, numElems, numNodes;
		double distance;

		nearestNeighbor.clear();
		distNN.clear();
		numElems = meshElements->numberOfItems();
		numNodes = nodes->numberOfItems();
		nearestNeighbor.resize(numElems);
		distNN.resize(numElems);
		int * kValueElement = new int[numElems];
		int * kValueNode = new int[numNodes];
#pragma omp parallel for private(i)
		for(i = 0; i < numElems; i++) nearestNeighbor[i] = -1;
#pragma omp parallel for private(i)
		for(i = 0; i < numElems; i++) distNN[i] = -1;
#pragma omp parallel for private(i)
		for(i = 0; i < numElems; i++) kValueElement[i] = 0;
#pragma omp parallel for private(i)
		for(i = 0; i < numNodes; i++) kValueNode[i] = 0;
		//Find the k-values (number of required nearest neighbors)
		for(i = 0; i < numElems; i++)
		{
			if((*meshElements)[i]->nodes().numberOfItems() == 3)
				kValueElement[i] = (*meshElements)[i]->nodes()[0]->elements().numberOfItems() + 
									(*meshElements)[i]->nodes()[1]->elements().numberOfItems() + 
									(*meshElements)[i]->nodes()[2]->elements().numberOfItems() - 3 + 1;
		}
		for(i = 0; i < numNodes; i++)
		{
			for(j = 0; j < (*nodes)[i]->elements().numberOfItems(); j++)
				kValueNode[i] = max(kValueNode[i], kValueElement[(*nodes)[i]->elements()[j]->id()-1]);
		}
		//Build a kd-tree for the nodes
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
			dataPts[i][0] = (*nodes)[i]->coordinate().x();
			dataPts[i][1] = (*nodes)[i]->coordinate().y();
			dataPts[i][2] = (*nodes)[i]->coordinate().z();
		}
		kdTree = new ANNkd_tree(dataPts, numNodes, 3);

		for(i = 0; i < numNodes; i++)
		{
			//Do k-nearest neighbor search and insert triangles into neighbors
			queryPt[0] = (*nodes)[i]->coordinate().x();
			queryPt[1] = (*nodes)[i]->coordinate().y();
			queryPt[2] = (*nodes)[i]->coordinate().z();
			kdTree->annkPriSearch(queryPt, kValueNode[i], nnIdx, dists);
			g_ElementContainer contToCheck;
			g_PtrAdapter<g_Element *> adapt(contToCheck);
			for(j = 0; j < kValueNode[i]; j++)
			{
				contToCheck.append((*nodes)[nnIdx[j]]->elements());
			}
			adapt.removeDuplicates();
			int contNum = contToCheck.numberOfItems();
			//Find nearest neighbor
			int numElemNeighbors = (*nodes)[i]->elements().numberOfItems();
			for(j = 0; j < numElemNeighbors; j++)
			{
				int id1 = (*nodes)[i]->elements()[j]->id()-1;
#pragma omp parallel for private(k, distance)
				for(k = 0; k < contNum; k++)
				{
					int id2 = contToCheck[k]->id()-1;
					if(id2 <= id1) continue;
					if(((*meshElements)[id1]->nodes()[0] == (*meshElements)[id2]->nodes()[0]) || 
						((*meshElements)[id1]->nodes()[0] == (*meshElements)[id2]->nodes()[1]) || 
						((*meshElements)[id1]->nodes()[0] == (*meshElements)[id2]->nodes()[2]) || 
						((*meshElements)[id1]->nodes()[1] == (*meshElements)[id2]->nodes()[0]) || 
						((*meshElements)[id1]->nodes()[1] == (*meshElements)[id2]->nodes()[1]) || 
						((*meshElements)[id1]->nodes()[1] == (*meshElements)[id2]->nodes()[2]) || 
						((*meshElements)[id1]->nodes()[2] == (*meshElements)[id2]->nodes()[0]) || 
						((*meshElements)[id1]->nodes()[2] == (*meshElements)[id2]->nodes()[1]) || 
						((*meshElements)[id1]->nodes()[2] == (*meshElements)[id2]->nodes()[2]))
						continue;
					distance = minimumDistance((*(*meshElements)[id1]), (*(*meshElements)[id2]));
#pragma omp critical
					if((distNN[id1] == -1) || (distance < distNN[id1]))
					{
						distNN[id1] = distance;
						nearestNeighbor[id1] = id2;
					}
#pragma omp critical
					if((distNN[id2] == -1) || (distance < distNN[id2]))
					{
						distNN[id2] = distance;
						nearestNeighbor[id2] = id1;
					}
				}
			}
		}
		//Delete kd-tree
		delete [] kValueElement;
		delete [] kValueNode;
		delete [] nnIdx;
		delete [] dists;
		annDeallocPt(queryPt);
		annDeallocPts(dataPts);
//		delete kdTree;
		annClose();	

#if MEASURE_TIME
		cout<<"Time to execute MinimumDistance::updateNN() "<<(int)time(0)-t<<" s"<<endl;
#endif
	}
}

double MinimumDistance::evaluateTrianglesOneNN()
{
	if(meshElements == NULL) return 0;
	else
	{
#if MEASURE_TIME
		int t = (int)time(0);
#endif
		int i, numElems;
		double stress;

		numElems = meshElements->numberOfItems();
		stress = 0.0;
		updateNNKd();

		for(i = 0; i < numElems; i++)
		{
			if(nearestNeighbor[i] == -1) continue;
			stress += 1.0/max(pow(distNN[i], 2), MD_DEFAULT_EPSILON);
		}
#if MEASURE_TIME
		cout<<"Time to execute MinimumDistance::evaluateTrianglesOneNN() "<<(int)time(0)-t<<" s"<<endl;
#endif
		return stress;
	}
}

void MinimumDistance::gradientTrianglesOneNN(double *& gradient)
{
	if(meshElements == NULL) return;
	else
	{
#if MEASURE_TIME
		int t = (int)time(0);
#endif
		int i, numElems, nNodes, index;
		nNodes = (int) nodes->numberOfItems();
		numElems = meshElements->numberOfItems();

#pragma omp parallel for private(i)
		for(i = 0; i < nNodes*3; i++) gradient[i] = 0.0;
		for(i = 0; i < numElems; i++)
		{
			if(nearestNeighbor[i] == -1) continue;

			minimumDistance((*(*meshElements)[i]), (*(*meshElements)[nearestNeighbor[i]]), &index);
			switch(index)
			{
			case 0:
				addSegmentGrad((*meshElements)[i]->nodes()[0], (*meshElements)[i]->nodes()[1], 
					(*meshElements)[nearestNeighbor[i]]->nodes()[0], (*meshElements)[nearestNeighbor[i]]->nodes()[1], gradient);
				break;
			case 1:
				addSegmentGrad((*meshElements)[i]->nodes()[0], (*meshElements)[i]->nodes()[1], 
					(*meshElements)[nearestNeighbor[i]]->nodes()[1], (*meshElements)[nearestNeighbor[i]]->nodes()[2], gradient);
				break;
			case 2:
				addSegmentGrad((*meshElements)[i]->nodes()[0], (*meshElements)[i]->nodes()[1], 
					(*meshElements)[nearestNeighbor[i]]->nodes()[2], (*meshElements)[nearestNeighbor[i]]->nodes()[0], gradient);
				break;
			case 3:
				addSegmentGrad((*meshElements)[i]->nodes()[1], (*meshElements)[i]->nodes()[2], 
					(*meshElements)[nearestNeighbor[i]]->nodes()[0], (*meshElements)[nearestNeighbor[i]]->nodes()[1], gradient);
				break;
			case 4:
				addSegmentGrad((*meshElements)[i]->nodes()[1], (*meshElements)[i]->nodes()[2], 
					(*meshElements)[nearestNeighbor[i]]->nodes()[1], (*meshElements)[nearestNeighbor[i]]->nodes()[2], gradient);
				break;
			case 5:
				addSegmentGrad((*meshElements)[i]->nodes()[1], (*meshElements)[i]->nodes()[2], 
					(*meshElements)[nearestNeighbor[i]]->nodes()[2], (*meshElements)[nearestNeighbor[i]]->nodes()[0], gradient);
				break;
			case 6:
				addSegmentGrad((*meshElements)[i]->nodes()[2], (*meshElements)[i]->nodes()[0], 
					(*meshElements)[nearestNeighbor[i]]->nodes()[0], (*meshElements)[nearestNeighbor[i]]->nodes()[1], gradient);
				break;
			case 7:
				addSegmentGrad((*meshElements)[i]->nodes()[2], (*meshElements)[i]->nodes()[0], 
					(*meshElements)[nearestNeighbor[i]]->nodes()[1], (*meshElements)[nearestNeighbor[i]]->nodes()[2], gradient);
				break;
			case 8:
				addSegmentGrad((*meshElements)[i]->nodes()[2], (*meshElements)[i]->nodes()[0], 
					(*meshElements)[nearestNeighbor[i]]->nodes()[2], (*meshElements)[nearestNeighbor[i]]->nodes()[0], gradient);
				break;
			case 9:
				addPointPlaneGrad((*meshElements)[i], (*meshElements)[nearestNeighbor[i]]->nodes()[0], gradient);
				break;
			case 10:
				addPointPlaneGrad((*meshElements)[i], (*meshElements)[nearestNeighbor[i]]->nodes()[1], gradient);
				break;
			case 11:
				addPointPlaneGrad((*meshElements)[i], (*meshElements)[nearestNeighbor[i]]->nodes()[2], gradient);
				break;
			case 12:
				addPointPlaneGrad((*meshElements)[nearestNeighbor[i]], (*meshElements)[i]->nodes()[0], gradient);
				break;
			case 13:
				addPointPlaneGrad((*meshElements)[nearestNeighbor[i]], (*meshElements)[i]->nodes()[1], gradient);
				break;
			case 14:
				addPointPlaneGrad((*meshElements)[nearestNeighbor[i]], (*meshElements)[i]->nodes()[2], gradient);
				break;
			}
		}
#if MEASURE_TIME
		cout<<"Time to execute MinimumDistance::gradientTrianglesOneNN(double *& gradient) "<<(int)time(0)-t<<" s"<<endl;
#endif
	}
}

void MinimumDistance::updateNNKd(int id)
{
	updateNNKd();
}

double MinimumDistance::evaluateTrianglesOneNN(int id)
{
	if(meshElements == NULL) return 0;
	else
	{
#if MEASURE_TIME
		int t = (int)time(0);
#endif
		int i, numElems;
		double stress;

		numElems = (*nodes)[id]->elements().numberOfItems();
		stress = 0.0;
		updateNNKd();

		for(i = 0; i < numElems; i++)
		{
			if(nearestNeighbor[(*nodes)[id]->elements()[i]->id()-1] == -1) continue;
			stress += 1.0/max(pow(distNN[(*nodes)[id]->elements()[i]->id()-1], 2), MD_DEFAULT_EPSILON);
		}
#if MEASURE_TIME
		cout<<"Time to execute MinimumDistance::evaluateTrianglesOneNN() "<<(int)time(0)-t<<" s"<<endl;
#endif
		return stress;
	}
}

void MinimumDistance::gradientTrianglesOneNN(int id, double *& gradient)
{
	if(meshElements == NULL) return;
	else
	{
		int i, index, localIdx, helpIdx;
		gradient[0] = 0.0; gradient[1] = 0.0; gradient[2] = 0.0;
		g_Element * nearestEl;
		g_ElementContainer localElems = (*nodes)[id]->elements();
		for(i = 0; i < (int)localElems.numberOfItems(); i++)
		{
			if(localElems[i]->nodes()[0] == (*nodes)[id]) localIdx = 0;
			else if(localElems[i]->nodes()[1] == (*nodes)[id]) localIdx = 1;
			else localIdx = 2;

			nearestEl = (*meshElements)[nearestNeighbor[localElems[i]->id()-1]];

			minimumDistance(*(localElems[i]), *nearestEl, &index);
			switch(index)
			{
			case 0:
				if(localIdx != 2)
				{
					if(localIdx == 0) helpIdx = 2;
					if(localIdx == 1) helpIdx = 3;
					addSegmentGrad(nearestEl->nodes()[0], nearestEl->nodes()[1], 
						localElems[i]->nodes()[0], localElems[i]->nodes()[1], gradient, helpIdx);
				}
				break;
			case 1:
				if(localIdx != 0)
				{
					if(localIdx == 1) helpIdx = 2;
					if(localIdx == 2) helpIdx = 3;
					addSegmentGrad(nearestEl->nodes()[0], nearestEl->nodes()[1], 
						localElems[i]->nodes()[1], localElems[i]->nodes()[2], gradient, helpIdx);
				}
				break;
			case 2:
				if(localIdx != 1)
				{
					if(localIdx == 2) helpIdx = 2;
					if(localIdx == 0) helpIdx = 3;
					addSegmentGrad(nearestEl->nodes()[0], nearestEl->nodes()[1], 
						localElems[i]->nodes()[2], localElems[i]->nodes()[0], gradient, helpIdx);
				}
				break;
			case 3:
				if(localIdx != 2)
				{
					if(localIdx == 0) helpIdx = 2;
					if(localIdx == 1) helpIdx = 3;
					addSegmentGrad(nearestEl->nodes()[1], nearestEl->nodes()[2], 
						localElems[i]->nodes()[0], localElems[i]->nodes()[1], gradient, helpIdx);
				}
				break;
			case 4:
				if(localIdx != 0)
				{
					if(localIdx == 1) helpIdx = 2;
					if(localIdx == 2) helpIdx = 3;
					addSegmentGrad(nearestEl->nodes()[1], nearestEl->nodes()[2], 
						localElems[i]->nodes()[1], localElems[i]->nodes()[2], gradient, helpIdx);
				}
				break;
			case 5:
				if(localIdx != 1)
				{
					if(localIdx == 2) helpIdx = 2;
					if(localIdx == 0) helpIdx = 3;
					addSegmentGrad(nearestEl->nodes()[1], nearestEl->nodes()[2], 
						localElems[i]->nodes()[2], localElems[i]->nodes()[0], gradient, helpIdx);
				}
				break;
			case 6:
				if(localIdx != 2)
				{
					if(localIdx == 0) helpIdx = 2;
					if(localIdx == 1) helpIdx = 3;
					addSegmentGrad(nearestEl->nodes()[2], nearestEl->nodes()[0], 
						localElems[i]->nodes()[0], localElems[i]->nodes()[1], gradient, helpIdx);
				}
				break;
			case 7:
				if(localIdx != 0)
				{
					if(localIdx == 1) helpIdx = 2;
					if(localIdx == 2) helpIdx = 3;
					addSegmentGrad(nearestEl->nodes()[2], nearestEl->nodes()[0], 
						localElems[i]->nodes()[1], localElems[i]->nodes()[2], gradient, helpIdx);
				}
				break;
			case 8:
				if(localIdx != 1)
				{
					if(localIdx == 2) helpIdx = 2;
					if(localIdx == 0) helpIdx = 3;
					addSegmentGrad(nearestEl->nodes()[2], nearestEl->nodes()[0], 
						localElems[i]->nodes()[2], localElems[i]->nodes()[0], gradient, helpIdx);
				}
				break;
			case 9:
				if(localIdx == 0)
					addPointPlaneGrad(nearestEl, localElems[i]->nodes()[0], gradient, 0);
				break;
			case 10:
				if(localIdx == 1)
					addPointPlaneGrad(nearestEl, localElems[i]->nodes()[1], gradient, 0);
				break;
			case 11:
				if(localIdx == 2)
					addPointPlaneGrad(nearestEl, localElems[i]->nodes()[2], gradient, 0);
				break;
			case 12:
				addPointPlaneGrad(localElems[i], nearestEl->nodes()[0], gradient, localIdx+1);
				break;
			case 13:
				addPointPlaneGrad(localElems[i], nearestEl->nodes()[1], gradient, localIdx+1);
				break;
			case 14:
				addPointPlaneGrad(localElems[i], nearestEl->nodes()[2], gradient, localIdx+1);
				break;
			}
		}
	}
}
