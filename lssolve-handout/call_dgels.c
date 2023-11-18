#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "msptools.h"

/* C prototype for LAPACK routine DGELS */
void dgels_(
	const char *trans, /* 'N' or 'T'             */
	const int *m,	   /* rows in A              */
	const int *n,	   /* cols in A              */
	const int *nrhs,   /* cols in B              */
	double *A,		   /* array A                */
	const int *lda,	   /* leading dimension of A */
	double *B,		   /* array B                */
	const int *ldb,	   /* leading dimension of B */
	double *work,	   /* workspace array        */
	int *lwork,		   /* workspace size         */
	int *info		   /* status code            */
);

/* call_dgels : wrapper for LAPACK's DGELS routine

Purpose:
Solves the least-squares problem

   minimize  || A*x-b ||_2^2

using LAPACK's DGELS routine. Upon exit, the input vector b is overwriten, i.e.,
the leading n elements of b contain the solution x.

Return value:
The function returns the output `info` from DGELS with the
following exceptions: the return value is

	-12 if the input A is NULL and/or the input b is NULL
	-13 if A->shape[0] < A->shape[1]
	-14 if the dimensions of the arrays that A and b point to are incompatible
	-15 in case of memory allocation errors.
*/


double ltwo_norm(array_t *x) {
    if (x == NULL || x->val == NULL) {
        return -1;
    }

    double sum = 0;
    for (size_t i = 0; i < x->len; i++) {
        sum += x->val[i] * x->val[i];
    }
    return sqrt(sum);
}

// Define max and min functions
static inline int max(int a, int b) {
    return (a > b) ? a : b;
}

static inline int min(int a, int b) {
    return (a < b) ? a : b;
}

int call_dgels(array2d_t * A, array_t * b, double * resnorm, double * rsquared)
{
	printf("doing stuff \n");
	if (A == NULL || b == NULL) {
		return -12;
	}
	int m = A->shape[0];
	int n = A->shape[1];
	int nrhs = 1;

	if (m < n) {
		return -13;
	}
	if (b->len != m) {
		printf("m is %d \n", m);
		printf("n is %d \n", n);
		printf("b length: %zu \n", b->len);
		return -14;
	}

	if (rsquared != NULL) {
		double mean = 0;
		*rsquared = 0;
			for (size_t i = 0; i < b->len; i++) {
				mean += b->val[i]/b->len;
			}
			for (size_t i = 0; i < b->len; i++) {
				*rsquared += (b->val[i]-mean)*(b->val[i]-mean);
			}
	}
	
	int lwork = max(1, min(m, n) + max(min(m, n), nrhs));

	//MEMORY ALLOCATION AND COPYING
	double* work_array = (double*)malloc(lwork * sizeof(double));


	if (!work_array) {
        return -15;
    }

	if (A->order == RowMajor) {
		char trans = 'T';
		int status_code = 1;
		int lda = n;
		int ldb = m;

		dgels_(&trans, &n, &m, &nrhs, A->val, &lda, b->val, &ldb, work_array, &lwork, &status_code);
		free(work_array);
		if (status_code != 0) {
			printf("failed status code");
			return status_code;
		}
	} else {
		char trans = 'N';
		int status_code = 1;
		int lda = m;
		int ldb = m;
		dgels_(&trans, &m, &n, &nrhs, A->val, &lda, b->val, &ldb, work_array, &lwork, &status_code);
		free(work_array);
		if (status_code != 0) {
			printf("failed status code");
			return status_code;
		}
	}
	b->len = n;

		if (resnorm != NULL) {
			array_t *residuals = array_alloc(b->capacity);
			double *dummy = b->val + n;
			memcpy(residuals->val, dummy, (m-n) * sizeof(double)); // this could be done better
			residuals->len = m-n;

			*resnorm = ltwo_norm(residuals);
			free(residuals);
		}
		if (rsquared != NULL) {
			*rsquared = 1 - (*resnorm * *resnorm) / *rsquared;
		}
	return 0;
}