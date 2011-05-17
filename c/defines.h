/* -*- Mode: C; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#ifndef __DEFINES_H__
#define __DEFINES_H__

/*#include <float.h>*/

/*************************************************************/
/* Type and constants definitions                            */

#ifndef DBL_MIN
#define DBL_MIN 2.2250738585072014e-308
#endif
#ifndef FLT_MIN
#define FLT_MIN 1.17549435e-38F
#endif

/* BOOL type as usual */
typedef int BOOL;
#define FALSE 0
#define TRUE 1

#define MAX_DOF 3
#define MAX_MATERIAL_PARAMETERS 10

/* Redefine type of the floating point values */
#ifdef SINGLE
typedef float real;
#else /* not SINGLE */
typedef double real;
#endif /* not SINGLE */

/* Equals macro for real values */
#ifdef SINGLE
#define EQL(x,y) ((fabs((x)-(y))<= (FLT_MIN)) ? TRUE:FALSE)
#else
#define EQL(x,y) ((fabs((x)-(y))<= (DBL_MIN)) ? TRUE:FALSE)
#endif

/* Kroneker delta */
#define DELTA(i,j) ((i)==(j) ? 1 : 0)

/* shortcut for adding of the matrix elements */
#define MTX(m,i,j,v) sp_matrix_element_add((m),(i),(j),(v));

#endif /* __DEFINES_H__ */
