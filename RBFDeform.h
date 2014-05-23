// Deformation using a set of Radial Basis Functions
// Stefanie Wuhrer

#ifndef RBFDeform_H
#define RBFDeform_H

#include <iostream>
#include <math.h>

/*
//includes for CLAPACK
namespace clapack
{
	extern "C"
	{
		#include "blaswrap.h"
		#include "f2c.h"
		extern int dgemm_(char *transa, char *transb, integer *m, integer *
			n, integer *k, doublereal *alpha, doublereal *a, integer *lda, 
			doublereal *b, integer *ldb, doublereal *beta, doublereal *c, 
			integer *ldc);
		extern int dgemv_(char *trans, integer *m, integer *n, doublereal *
			alpha, doublereal *a, integer *lda, doublereal *x, integer *incx, 
			doublereal *beta, doublereal *y, integer *incy);
		extern int dsyevx_(char *jobz, char *range, char *uplo, integer *n, 
			doublereal *a, integer *lda, doublereal *vl, doublereal *vu, integer *
			il, integer *iu, doublereal *abstol, integer *m, doublereal *w, 
			doublereal *z__, integer *ldz, doublereal *work, integer *lwork, 
			integer *iwork, integer *ifail, integer *info);
		extern double dlamch_(char *cmach);
		extern int dgetrf_(integer *m, integer *n, doublereal *a, integer *
			lda, integer *ipiv, integer *info);
		extern int dgetrs_(char *trans, integer *n, integer *nrhs, 
			doublereal *a, integer *lda, integer *ipiv, doublereal *b, integer *
			ldb, integer *info);
		extern int dgetri_(integer *n, doublereal *a, integer *lda, integer 
			*ipiv, doublereal *work, integer *lwork, integer *info);
		extern int dgesvd_(char *jobu, char *jobvt, integer *m, integer *n, 
			doublereal *a, integer *lda, doublereal *s, doublereal *u, integer *
			ldu, doublereal *vt, integer *ldvt, doublereal *work, integer *lwork, 
			integer *info);
		extern int dgels_(char *trans, integer *m, integer *n, integer *
			nrhs, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
			doublereal *work, integer *lwork, integer *info);
		extern int dgesvx_(char *fact, char *trans, integer *n, integer *
			nrhs, doublereal *a, integer *lda, doublereal *af, integer *ldaf, 
			integer *ipiv, char *equed, doublereal *r__, doublereal *c__, 
			doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
			rcond, doublereal *ferr, doublereal *berr, doublereal *work, integer *
			iwork, integer *info);
	}
}
*/

#include "clapack_defines.h"

using namespace std; 

#define RBF_PI 3.14159265
#define RBF_EPSILON 0.000000001

class RBFDeform
{

public:
	RBFDeform();
	//Constructor: number of centers, dimensionality of RBF (always same dimensionality of input and output for us),
	//centers, and function of centers.
	RBFDeform(int numberCenters, int dimension, double * centers, double * functionValues, double standDev = 0, 
		bool approximateRBF = true);
	~RBFDeform();
	//Return the result in result pointer (memory must be allocated outside the function)
	void evaluate(double * point, double * result);

private:
	int numberCenters, dimension, srbfcoeff;
	double * centers, * fctValues;
	double * rbf_coeff;
	double variance, standDev;
	bool approximation;

	void initialize(int numberCenters, int dimension, double * centers, double * functionValues, double standDev);
	inline double RBFunction(unsigned int id1, unsigned int id2);
	inline double RBFunction(unsigned int id1, double * pointToEvaluate); 
	void compute_RBFCoeff();
};

#endif
