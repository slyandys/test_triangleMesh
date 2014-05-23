#include "RBFDeform.h"

RBFDeform::RBFDeform()
{
	centers = NULL;
	fctValues = NULL;
	rbf_coeff = NULL;
}

RBFDeform::RBFDeform(int numberCenters, int dimension, double * centers, double * functionValues, double standDev,
					 bool approximateRBF)
{
	initialize(numberCenters, dimension, centers, functionValues, standDev);
	approximation = approximateRBF;
}

RBFDeform::~RBFDeform()
{
	if(centers != NULL)
	{
		delete [] centers;
		centers = NULL;
	}
	if(fctValues != NULL)
	{
		delete [] fctValues;
		fctValues = NULL;
	}
	if(rbf_coeff != NULL)
	{
		delete [] rbf_coeff;
		rbf_coeff = NULL;
	}
}

void RBFDeform::initialize(int numberCenters, int dimension, double * centers, double * functionValues, double standDev)
{
   int i;
   rbf_coeff = NULL;
   this->numberCenters = numberCenters;
   this->dimension = dimension;
   srbfcoeff= numberCenters + 4;
   this->standDev = standDev;
   variance = standDev * standDev;
   
   this->centers = new double [numberCenters * dimension];
   for(i = 0; i < numberCenters * dimension; i++)
   {
	   this->centers[i] = centers[i];
   }
   fctValues = new double [numberCenters * dimension];
   for(i = 0; i < numberCenters * dimension; i++)
   {
	   fctValues[i] = functionValues[i];
   }
   
   compute_RBFCoeff();
}

void RBFDeform::compute_RBFCoeff()
{	
	int i, j;
	double *A = new double[srbfcoeff*srbfcoeff];
	//initialize everything to 0:
	for(i = 0; i < srbfcoeff*srbfcoeff; i++)
		A[i] = 0.0;
	rbf_coeff = new double[srbfcoeff * dimension];
	for(i = 0; i < srbfcoeff * dimension; i++)
		rbf_coeff[i] = 0.0;

	//set values for computation:
	for ( i = 0 ; i < numberCenters; i++ )
	{ 
		for ( j = 0; j < numberCenters; j++ )
			A[j * srbfcoeff + i] = RBFunction(i, j);
		for ( j = 0 ; j < dimension; j ++ )
		{
			// Added if rigid tranformation + RBF
			A[(j + numberCenters + 1) * srbfcoeff + i] = centers[i * dimension + j];
			A[i * srbfcoeff + (j + numberCenters + 1)] = centers[i * dimension + j];

			rbf_coeff[j * srbfcoeff + i] = fctValues[i * dimension + j];
		}
		// Added if rigid tranformation + RBF
		A[numberCenters * srbfcoeff + i] = 1;
		A[i * srbfcoeff + numberCenters] = 1;
	}

	//Only use this for approximate RBF
	if(approximation)
	{
		//make it an approximate RBF by setting rho: maybe use cross-validation later
		//use mean distance for now: like Jain, Zhang and van Kaick:
		double rho = 0;
		double squaredDist = 0;
		for (i = 0; i < numberCenters; i++ )
		{
			squaredDist = 0;
			for(j = 0; j < dimension; j++)
				squaredDist += pow((centers[i * dimension + j]-fctValues[i * dimension + j]), 2);
			rho += sqrt(squaredDist);
		}
		rho = rho / numberCenters;
		for ( i = 0 ; i < numberCenters; i++ )A[i * srbfcoeff + i]+=rho;
	}

	long int info;
	char uplo = 'U';
	long int * ipiv = new long int[srbfcoeff];
	long int n = (long int) srbfcoeff;
	long int m = (long int) dimension;
	char fact = 'E';
	char trans = 'N';
	char equed;
	double * af = new double[srbfcoeff*srbfcoeff];
	double * r = new double[srbfcoeff];
	double * c = new double[srbfcoeff];
	double * x = new double[srbfcoeff * dimension];
	double rcond;
	double * ferr = new double[dimension];
	double * berr = new double[dimension];
	double * work = new double [4*srbfcoeff];
	long int * iwork = new long int[srbfcoeff];
	clapack::dgesvx_(&fact, &trans, &n, &m, A, &n, af, &n, ipiv, &equed, r, c, rbf_coeff, &n, x, &n, &rcond, ferr, berr, work, 
		iwork, &info);
	//copy x to rbf_coeff:
	for(i = 0; i < srbfcoeff * dimension; i++)
		rbf_coeff[i] = x[i];
	//delete stuff:
	delete [] af;
	delete [] r;
	delete [] c;
	delete [] x;
	delete [] ferr;
	delete [] berr;
	delete [] work;
	delete [] iwork;

	delete [] A;
	delete [] ipiv;
}

inline double RBFDeform::RBFunction( unsigned int id1, unsigned int id2) 
{
	double squaredDist = 0;
	for(int i = 0; i < dimension; i++)
		squaredDist += pow((centers[id1 * dimension + i]-centers[id2 * dimension + i]), 2);
	
	//Thin plate spline:
	if(squaredDist < RBF_EPSILON) return 0;
	else return squaredDist * log(sqrt(squaredDist));
	
	//Gauss:
	//return exp( - squaredDist / (2*variance)) / (standDev * sqrt(2.0 * RBF_PI));
}

inline double RBFDeform::RBFunction(unsigned int id1, double * pointToEvaluate)
{
	double squaredDist = 0;
	for(int i = 0; i < dimension; i++)
		squaredDist += pow((pointToEvaluate[i]-centers[id1 * dimension + i]), 2);
	
	//Thin plate spline:
	if(squaredDist < RBF_EPSILON) return 0;
	else return squaredDist * log(sqrt(squaredDist));
	
	//Gauss:
	//return exp( - squaredDist / (2*variance)) / (standDev * sqrt(2.0 * RBF_PI));
}

void RBFDeform::evaluate(double * point, double * result)
{
	int i, j;
	double rbfunction;
	//initialize the result to 0:
	for(i = 0; i < dimension; i++) result[i] = 0;
	//evaluate the function using previously computed rbf_coeff
	for(j = 0; j < numberCenters; j++)
	{
		rbfunction = RBFunction(j, point);
		for(i = 0; i < dimension; i++)
		{
			result[i] += rbf_coeff[i * srbfcoeff + j] * rbfunction;
		}
	}
	for(i = 0; i < dimension; i++)
	{
		result[i] += rbf_coeff[i * srbfcoeff + numberCenters];
		for(j = 0; j < dimension; j++)
			result[i] += rbf_coeff[i * srbfcoeff + (numberCenters + j + 1)] * point[j];
	}
}