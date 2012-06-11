/* -*- Mode: C; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

#include "defines.h"
#include "tests.h"
#include "dense_matrix.h"

static BOOL test_dense_matrix()
{
  BOOL result = TRUE;
  int i,j;
  /* input data */
  real A[3][3] = {{1, 2, 0}, {2, 0, 3}, {0, 2, 3}};
  real B[3][3] = {{0, 2, 1}, {1, 1, 1}, {3, 2, -1}};
  /* expected results */
  real result_matmul[3][3] = {{2, 4, 3}, {9, 10, -1}, {11, 8, -1}};
  real result_matmul_transp[3][3] = {{2, 4, 3}, {6, 8, 0}, {12, 9, 0}};
  real result_matmul_transp2[3][3] = {{4, 3, 7}, {3, 5, 3}, {7, 5, 1}};
  /* results array */
  real R[3][3];

  /* test A x B */
  matrix_mul3x3(A,B,R);
  for (i = 0; i < 3; ++ i)
    for (j = 0; j < 3; ++ j)
      result &= EQL(R[i][j],result_matmul[i][j]);

  if (result)
  {
    /* test A x B' */
    matrix_transpose_mul3x3(A,B,R);
    for (i = 0; i < 3; ++ i)
      for (j = 0; j < 3; ++ j)
        result &= EQL(R[i][j],result_matmul_transp[i][j]);
  }

  if (result)
  {
    /* test A' x B */
    matrix_transpose2_mul3x3(A,B,R);
    for (i = 0; i < 3; ++ i)
      for (j = 0; j < 3; ++ j)
        result &= EQL(R[i][j],result_matmul_transp2[i][j]);
  }
  printf("test_matrix result: *%s*\n",result ? "pass" : "fail");
  return result;
}

BOOL do_tests()
{
  return test_dense_matrix();
}
