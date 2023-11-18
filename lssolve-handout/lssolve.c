#include <stdlib.h>
#include <stdio.h>
#include "msptools.h"

int call_dgels(array2d_t *A, array_t *b, double *resnorm, double *rsquared);

int main(int argc, char *argv[])
{

  if (argc != 4)
  {
    fprintf(stderr, "Usage: %s A b x\n", argv[0]);
    return EXIT_FAILURE;
  }

  array2d_t *A = array2d_from_file(argv[1]);
  array_t *b = array_from_file(argv[2]);
  if (A == NULL || b == NULL || (A->shape[0] < A->shape[1])) {
    return EXIT_FAILURE;
  }

  double resnorm = 0.0;
  double rsquared = 0.0;
  int status_code = 0;

  status_code = call_dgels(A, b, &resnorm, &rsquared);

  if (status_code != 0) {
    return 1;
  }
  
  
  array_to_file(argv[3], b);


  return EXIT_SUCCESS;
}
