/* -*- Mode: C; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#ifndef __DEFINES_H__
#define __DEFINES_H__

/*************************************************************/
/* Type and constants definitions                            */


#ifndef FLT_EPSILON
#define FLT_EPSILON 1.19209290e-07F
#endif
#ifndef DBL_EPSILON
#define DBL_EPSILON 2.2204460492503131e-16
#endif

/* BOOL type as usual */
typedef int BOOL;
#define FALSE 0
#define TRUE 1

#define MAX_DOF 3
#define MAX_MATERIAL_PARAMETERS 10

/* define specific macros used by GCC compiler */
#ifdef __GNUC__
/*
 * Use all these macros AFTER function declaration, like:
 * void print_me() DEPRECATED;
 */
#define DEPRECATED __attribute__ ((deprecated))
#define UNUSED __attribute__ ((unused))
#else
#define DEPRECATED
#define UNUSED
#endif /* __GNUC__ */

/* Redefine type of the floating point values */
#ifdef SINGLE
typedef float real;
#define REAL_EPSILON FLT_EPSILON
#else /* not SINGLE */
typedef double real;
#define REAL_EPSILON DBL_EPSILON
#endif /* SINGLE */

/*
 * Equals macro for real values
 * See http://www.rsdn.ru/forum/cpp/2640596.1.aspx for explanations
 */
#define EQL(x,y) ((fabs((x)-(y))<= fmax(fabs((x)),fabs((y)))*(REAL_EPSILON)) ? TRUE:FALSE)

/* Kroneker delta */
#define DELTA(i,j) ((i)==(j) ? 1 : 0)

/* shortcut for adding of the matrix elements */
#define MTX(m,i,j,v) sp_matrix_element_add((m),(i),(j),(v));

#endif /* __DEFINES_H__ */
