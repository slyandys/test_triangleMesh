//includes for CLAPACK
#ifndef NO_F2C
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
	  /* // Not used
		extern int dgemv_(char *trans, integer *m, integer *n, doublereal *
			alpha, doublereal *a, integer *lda, doublereal *x, integer *incx, 
			doublereal *beta, doublereal *y, integer *incy);
		extern int dsyevx_(char *jobz, char *range, char *uplo, integer *n, 
			doublereal *a, integer *lda, doublereal *vl, doublereal *vu, integer *
			il, integer *iu, doublereal *abstol, integer *m, doublereal *w, 
			doublereal *z__, integer *ldz, doublereal *work, integer *lwork, 
			integer *iwork, integer *ifail, integer *info);
		extern double dlamch_(char *cmach);
	  */
		extern int dgetrf_(integer *m, integer *n, doublereal *a, integer *
			lda, integer *ipiv, integer *info);
	  /* // Not used
		extern int dgetrs_(char *trans, integer *n, integer *nrhs, 
			doublereal *a, integer *lda, integer *ipiv, doublereal *b, integer *
			ldb, integer *info);
	  */
		extern int dgetri_(integer *n, doublereal *a, integer *lda, integer 
			*ipiv, doublereal *work, integer *lwork, integer *info);
	  /* // Not used
		extern int dgesvd_(char *jobu, char *jobvt, integer *m, integer *n, 
			doublereal *a, integer *lda, doublereal *s, doublereal *u, integer *
			ldu, doublereal *vt, integer *ldvt, doublereal *work, integer *lwork, 
			integer *info);
	  */
		extern int dgels_(char *trans, integer *m, integer *n, integer *
			nrhs, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
			doublereal *work, integer *lwork, integer *info);

		extern void dgesv_(integer *n, integer *nrhs, doublereal *a, integer *lda, 
			    integer *ipiv, doublereal *b, integer *ldb, integer *info);

		extern int dgesvx_(char *fact, char *trans, integer *n, integer *
			nrhs, doublereal *a, integer *lda, doublereal *af, integer *ldaf, 
			integer *ipiv, char *equed, doublereal *r__, doublereal *c__, 
			doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
			rcond, doublereal *ferr, doublereal *berr, doublereal *work, integer *
			iwork, integer *info);
	}
}
#else
namespace clapack
{
	extern "C"
	{
		extern void dgemm_(char *transa, char *transb, long *m, long *
			n, long *k, double *alpha, double *a, long *lda, 
			double *b, long *ldb, double *beta, double *c, 
			long *ldc);
	  /* // Not used
		extern void dgemv_(char *trans, long *m, long *n, double *
			alpha, double *a, long *lda, double *x, long *incx, 
			double *beta, double *y, long *incy);
		extern void dsyevx_(char *jobz, char *range, char *uplo, long *n, 
			double *a, long *lda, double *vl, double *vu, long *
			il, long *iu, double *abstol, long *m, double *w, 
			double *z__, long *ldz, double *work, long *lwork, 
			long *iwork, long *ifail, long *info);
		extern double dlamch_(char *cmach);
	  */
		extern void dgetrf_(long *m, long *n, double *a, long *
			lda, long *ipiv, long *info);
	  /* // Not used
		extern void dgetrs_(char *trans, long *n, long *nrhs, 
			double *a, long *lda, long *ipiv, double *b, long *
			ldb, long *info);
	  */
		extern void dgetri_(long *n, double *a, long *lda, long 
			*ipiv, double *work, long *lwork, long *info);
	  /* // Not used
	        extern void dgesvd_(char *jobu, char *jobvt, long *m, long *n, 
			double *a, long *lda, double *s, double *u, long *
			ldu, double *vt, long *ldvt, double *work, long *lwork, 
			long *info);
	  */
		extern void dgels_(char *trans, long *m, long *n, long *
			nrhs, double *a, long *lda, double *b, long *ldb, 
			double *work, long *lwork, long *info);
		extern void dgesv_(long *n, long *nrhs, double *a, long *lda, 
			    long *ipiv, double *b, long *ldb, long *info);

		extern int dgesvx_(char *fact, char *trans, long *n, long *
			nrhs, double *a, long *lda, double *af, long *ldaf, 
			long *ipiv, char *equed, double *r__, double *c__, 
			double *b, long *ldb, double *x, long *ldx, double *
			rcond, double *ferr, double *berr, double *work, long *
			iwork, long *info);

	}
}
#endif

