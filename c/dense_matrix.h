/* -*- Mode: C; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#ifndef __DENSE_MATRIX_H__
#define __DENSE_MATRIX_H__

#include "defines.h"


/* A storage to keep 2nd rank tensor components  */
typedef struct tensor_tag {
  real components[MAX_DOF][MAX_DOF];
} tensor;
typedef tensor* tensor_ptr;

/*
 * Vector norm
 */
real vector_norm(real* vector, int size);
/* Scalar multiplication of 2 vectors of the same size */
real cdot(real* vector1, real* vector2, int size);

/*************************************************************/
/* Functions for operating on matrix 3x3                     */

/* Calculates the determinant of the matrix 3x3 */
real det3x3(real (*matrix3x3)[3]);

/*
 * Calculates in-place inverse of the matrix 3x3
 * Returns FALSE if the matrix is ill-formed
 * det shall store a determinant of the matrix
 */
BOOL inv3x3(real (*matrix3x3)[3], real* det);

/*
 * Performs matrix multiplication A x B writing result to R
 */
void matrix_mul3x3 (real (*A)[3],real (*B)[3],real (*R)[3]);

/*
 * Performs matrix multiplication A' x B writing result to R
 */
void matrix_transpose_mul3x3 (real (*A)[3],real (*B)[3],real (*R)[3]);

/*
 * Performs matrix multiplication A x B' writing result to R
 */
void matrix_transpose2_mul3x3 (real (*A)[3],real (*B)[3],real (*R)[3]);



#endif /* __DENSE_MATRIX_H__ */
