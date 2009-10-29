#include <stdio.h>
#include <string.h>
/* #include <pcre.h> */

/*************************************************************/
/* Type definitions                                          */

/* BOOL type as usual */
typedef int BOOL;
#define FALSE 0
#define TRUE 1

/* Redefine type of the floating point values */
typedef double real;

enum task_type { /* PLANE_STRESS, PLANE_STRAIN, AXISYMMETRIC,  */CARTESIAN3D};
enum element_type { /* TRIANGLE3, TRIANGLE6,TETRAHEDRA4, */TETRAHEDRA10 };
enum prescribed_boundary_type {
  FREE = 0,                    /* free */
  PRESCRIBEDX = 1,             /* x prescribed */
  PRESCRIBEDY = 2,             /* y prescribed */
  PRESCRIBEDXY = 3,            /* x, y prescribed */
  PRESCRIBEDZ = 4,             /* z prescribed*/
  PRESCRIBEDXZ = 5,            /* x, z prescribed*/
  PRESCRIBEDYZ = 6,            /* y, z prescribed*/
  PRESCRIBEDXYZ = 7            /* x, y, z prescribed.*/
}
/*
 * Task type declaration.
 * Defines an input parameters for the task, independent of 
 * the input geometry and loads
 */
typedef struct fea_task_tag {
  task_type type;               /* type of the task to solve */
  unsigned char dof;            /* number of degree of freedom */
  element_type ele_type;        /* type of the element */
  int load_increments_no;       /* number of load increments */
  real desired_tolerance;       /* desired energy tolerance */
  int linesearch_max;           /* maximum number of line searches */
  BOOL modified_newton;         /* use modified Newton's method or not */

} fea_task;


/* Calculated solution parameters */
typedef struct fea_solution_params_tag {
  int msize;                    /* size of the global stiffness matrix */
  int nodes_per_element;        /* number of nodes defined in element 
                                   based on fea_task::ele_type */
};

/*************************************************************/
/* Input geometry parameters                                 */

/* An array of nodes. */
typedef struct nodes_array_tag {
  int nodes_count;              /* number of input nodes */
  real **nodes;                 /* nodes array,sized as nodes_count x dof
                                 * so access is like nodes[node_number][dof]*/
};

/* An array of elements */
typedef struct elements_array_tag {
  int elements_count;           /* number of elements */
  int **elements;               /* elements array, each line represents an
                                 * element. Element is an array of node
                                 * indexes
                                 */
};

/* Particular prescribed boundary node */
typedef struct prescibed_boundary_node_tag {
  int node_number;
  real values[3];
  prescribed_boundary_type type;
} prescibed_boundary_node;
/* An array of prescribed boundary conditions */
typedef struct prescribed_boundary_array_tag {
  int prescribed_number;
  prescibed_boundary_node *prescribed_nodes;
} prescribed_boundary_array;


/*
 * Solution data parameters.
 * Most of them derived from fea_task
 */
typedef struct fea_solution_params_tag {

} fea_solution_params;


int main(int argc, char **argv)
{
  
  

  return 0;
}
