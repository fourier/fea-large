/* -*- Mode: C; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include <math.h>
#include "defines.h"
#include "dense_matrix.h"


real vector_norm(real* vector, int size)
{
  real norm = 0.0;
  int i = 0;
  for ( ; i < size; ++ i)
    norm += vector[i]*vector[i];
  return sqrt(norm);
}

real cdot(real* vector1, real* vector2, int size)
{
  real result = 0;
  int i = 0;
  for ( ; i < size; ++ i)
    result += vector1[i]*vector2[i];
  return result;
}

real det3x3(real (*m)[3])
{
  real result;
  result = m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1]) - 
    m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0]) + 
    m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]);
  return result;
}

BOOL inv3x3(real (*m)[3],real* det)
{
  real m00,m01,m02,m10,m11,m12,m20,m21,m22;
  *det = det3x3(m);
  if (EQL(*det,0.0))
    return FALSE;
	/* calculate components */
	/* first row */
  m00 = (m[1][1]*m[2][2]-m[1][2]*m[2][1])/(*det);
	m01 = (m[0][2]*m[2][1]-m[0][1]*m[2][2])/(*det);
	m02 = (m[0][1]*m[1][2]-m[0][2]*m[1][1])/(*det);
	/* second row */
	m10 = (m[1][2]*m[2][0]-m[1][0]*m[2][2])/(*det);
	m11 = (m[0][0]*m[2][2]-m[0][2]*m[2][0])/(*det);
	m12 = (m[0][2]*m[1][0]-m[0][0]*m[1][2])/(*det);
	/* third row */
	m20 = (m[1][0]*m[2][1]-m[1][1]*m[2][0])/(*det);
	m21 = (m[0][1]*m[2][0]-m[0][0]*m[2][1])/(*det);
	m22 = (m[0][0]*m[1][1]-m[0][1]*m[1][0])/(*det);
  
  /* perform in-place substitution of the result */
	m[0][0] = m00; 	m[0][1] = m01; 	m[0][2] = m02;
	m[1][0] = m10; 	m[1][1] = m11; 	m[1][2] = m12;
	m[2][0] = m20;	m[2][1] = m21;	m[2][2] = m22;
  
  return TRUE;
}

void matrix_mul3x3 (real (*A)[3],real (*B)[3],real (*R)[3])
{
  int i,j,k;
  real sum;
  for (i = 0; i < 3; ++ i)
  {
    for (j = 0; j < 3; ++ j)
    {
      sum = 0;
      for (k = 0; k < 3; ++ k)
        sum += A[i][k]*B[k][j];
      R[i][j] = sum;
    }
  }
}


void matrix_transpose_mul3x3 (real (*A)[3],real (*B)[3],real (*R)[3])
{
  int i,j,k;
  real sum;
  for (i = 0; i < 3; ++ i)
  {
    for (j = 0; j < 3; ++ j)
    {
      sum = 0;
      for (k = 0; k < 3; ++ k)
        sum += A[k][i]*B[k][j];
      R[i][j] = sum;
    }
  }
}


void matrix_transpose2_mul3x3 (real (*A)[3],real (*B)[3],real (*R)[3])
{
  int i,j,k;
  real sum;
  for (i = 0; i < 3; ++ i)
  {
    for (j = 0; j < 3; ++ j)
    {
      sum = 0;
      for (k = 0; k < 3; ++ k)
        sum += A[i][k]*B[j][k];
      R[i][j] = sum;
    }
  }
}
