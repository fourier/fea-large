#include <stdio.h>
#include <string.h>
#include <pcre.h>

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
};


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
