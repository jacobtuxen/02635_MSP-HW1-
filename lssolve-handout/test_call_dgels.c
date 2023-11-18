#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "msptools.h"

#define MAX(a, b) ((a) > (b) ? (a) : (b));

int call_dgels(array2d_t *A, array_t *b, double *resnorm, double *rsquared);
int isclose(double a, double b, double rel_tol, double abs_tol);

int main(void)
{

    double arval[] = {-1.006, -1.634, -1.337,
                      0.791, 0.761, -1.474,
                      -0.117, 1.193, -0.042,
                      0.553, 1.632, -0.616,
                      -0.961, -1.532, 1.314};
    double acval[] = {-1.006, 0.791, -0.117, 0.553, -0.9610,
                      -1.634, 0.761, 1.193, 1.632, -1.5320,
                      -1.337, -1.474, -0.042, -0.616, 1.3140};
    double bval1[] = {4.97, -3.41, -2.90, -3.89, 3.06};
    double bval2[] = {4.97, -3.41, -2.90, -3.89, 3.06};
    double xval[] = {-1.057600048308403, -2.144472685786455, -0.145287705300387};

    array2d_t Ar = {
        .shape = {5, 3},
        .order = RowMajor,
        .val = arval,
    };
    array2d_t Ac = {
        .shape = {5, 3},
        .order = ColMajor,
        .val = acval,
    };
    array_t b1 = {.len = 5, .val = bval1, .capacity = 5};
    array_t b2 = {.len = 5, .val = bval2, .capacity = 5};
    array_t xref = {.len = 3, .val = xval, .capacity = 3};
    double rsquared_ref = 0.960227735451821;
    double resnorm_ref = 1.648078906183691;

    double resnorm = 0, rsquared = 0;

    printf("Testing with invalid inputs...\n");
    if ((call_dgels(NULL, &b1, &resnorm, &rsquared) != -12) || (call_dgels(&Ar, NULL, &resnorm, &rsquared) != -12))
    {
        fprintf(stderr, "  ***Test failed. Check if A or B is NULL.\n");
        return EXIT_FAILURE;
    }
    printf("  Passed.\n");

    printf("Testing with RowMajor input...\n");
    int ret = call_dgels(&Ar, &b1, &resnorm, &rsquared);
    if (ret != 0) {
        fprintf(stderr, "  ***Test failed. Expected return value 0, got %d.\n", ret);
        return EXIT_FAILURE;
    }
    if (!isclose(resnorm, resnorm_ref, 1e-8, 1e-8))
    {
        fprintf(stderr, "  ***Test failed. Unexpected resnorm.\n");
        return EXIT_FAILURE;
    }
    if (!isclose(rsquared, rsquared_ref, 1e-8, 1e-8))
    {
        fprintf(stderr, "  ***Test failed. Unexpected rsquared.\n");
        return EXIT_FAILURE;
    }
    if (b1.len != xref.len)
    {
        fprintf(stderr, "  ***Test failed. Unexpected length of b.\n");
        return EXIT_FAILURE;
    }
    for (size_t i = 0; i < b1.len; i++)
    {
        if (!isclose(b1.val[i], xref.val[i], 1e-8, 1e-8))
        {
            fprintf(stderr, "  ***Test failed. Unexpected result.\n");
            return EXIT_FAILURE;
        }
    }
    printf("  Passed.\n");

    printf("Testing with ColMajor input...\n");
    if (call_dgels(&Ac, &b2, &resnorm, &rsquared) != 0)
    {
        fprintf(stderr, "  ***Test failed. Expected return value 0.\n");
        return EXIT_FAILURE;
    }
    if (!isclose(resnorm, resnorm_ref, 1e-8, 1e-8))
    {
        fprintf(stderr, "  ***Test failed. Unexpected resnorm.\n");
        return EXIT_FAILURE;
    }
    if (!isclose(rsquared, rsquared_ref, 1e-8, 1e-8))
    {
        fprintf(stderr, "  ***Test failed. Unexpected rsquared.\n");
        return EXIT_FAILURE;
    }
    if (b2.len != xref.len)
    {
        fprintf(stderr, "  ***Test failed. Unexpected length of b.\n");
        return EXIT_FAILURE;
    }
    for (size_t i = 0; i < b2.len; i++)
    {
        if (!isclose(b2.val[i], xref.val[i], 1e-8, 1e-8))
        {
            fprintf(stderr, "  ***Test failed. Unexpected result.\n");
            return EXIT_FAILURE;
        }
    }
    printf("  Passed.\n");


    return EXIT_SUCCESS;
}

int isclose(double a, double b, double rel_tol, double abs_tol)
{
    if (isfinite(a) && isfinite(b))
    {
        // abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)
        double abs_a = fabs(a);
        double abs_b = fabs(b);
        double abs_ab_max = MAX(abs_a, abs_b);
        return fabs(a - b) <= MAX(rel_tol * abs_ab_max, abs_tol);
    }
    else if (isinf(a) && isinf(b))
    {
        return a == b; // a and b have the same sign
    }
    else if (isnan(a) && isnan(b))
    {
        return 1; // a and b are both NaN
    }
    else
    {
        return 0;
    }
}
