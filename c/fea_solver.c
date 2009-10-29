#include <stdio.h>
#include <string.h>
/* #include <pcre.h> */

/*************************************************************/
/* Type definitions                                          */

/* BOOL type as usual */
typedef int BOOL;
#define FALSE 0
#define TRUE 1

#define MAX_DOF 3

/* Redefine type of the floating point values */
typedef double real;

typedef enum task_type_enum {
		/* PLANE_STRESS, PLANE_STRAIN, AXISYMMETRIC,  */
		CARTESIAN3D } task_type;

typedef enum element_type_enum {
		/* TRIANGLE3, TRIANGLE6,TETRAHEDRA4, */
		TETRAHEDRA10 } element_type;

typedef enum prescribed_boundary_type_enum {
		FREE = 0,                    /* free */
		PRESCRIBEDX = 1,             /* x prescribed */
		PRESCRIBEDY = 2,             /* y prescribed */
		PRESCRIBEDXY = 3,            /* x, y prescribed */
		PRESCRIBEDZ = 4,             /* z prescribed*/
		PRESCRIBEDXZ = 5,            /* x, z prescribed*/
		PRESCRIBEDYZ = 6,            /* y, z prescribed*/
		PRESCRIBEDXYZ = 7            /* x, y, z prescribed.*/
} prescribed_boundary_type;

		
/*
 * Task type declaration.
 * Defines an input parameters for the task, independent of 
 * the input geometry and loads
 */
typedef struct fea_task_tag {
		task_type type;               /* type of the task to solve */
		unsigned char dof;            /* number of degree of freedom */
		element_type ele_type;        /* type of the element */
		int load_increments_count;    /* number of load increments */
		real desired_tolerance;       /* desired energy tolerance */
		int linesearch_max;           /* maximum number of line searches */
		BOOL modified_newton;         /* use modified Newton's method or not */

} fea_task;


/* Calculated solution parameters */
typedef struct fea_solution_params_tag {
		int msize;                    /* size of the global stiffness matrix */
		int nodes_per_element;        /* number of nodes defined in element 
																		 based on fea_task::ele_type */
		int gauss_nodes_count;			  /* number of gauss nodes per element */
} fea_solution_params;

/*************************************************************/
/* Input geometry parameters                                 */

/* An array of nodes. */
typedef struct nodes_array_tag {
		int nodes_count;              /* number of input nodes */
		real **nodes;                 /* nodes array,sized as nodes_count x dof
																	 * so access is  nodes[node_number][dof] */
} nodes_array;

/* An array of elements */
typedef struct elements_array_tag {
		int elements_count;           /* number of elements */
		int **elements;               /* elements array, each line represents an
																	 * element. Element is an array of node
																	 * indexes
																	 */
} elements_array;

/* Particular prescribed boundary node */
typedef struct prescibed_boundary_node_tag {
		int node_number;
		real values[MAX_DOF];
		prescribed_boundary_type type;
} prescibed_boundary_node;

/*
 * An array of prescribed boundary conditions
 * either fixed and with prescribed displacements
 */
typedef struct prescribed_boundary_array_tag {
		int prescribed_number;
		prescibed_boundary_node *prescribed_nodes;
} prescribed_boundary_array;


/*************************************************************/
/* Application-specific structures                           */

/* Structure describing information for the gauss node
 * depending on number of shape functions N per element
 * TODO: add tables with layouts in comments */
typedef struct gauss_node_tag {
		real weigth; 								/* weight for the integration */
		real *forms; 								/* shape function values for gauss node, N */
		real **dforms;							/* derivatives of shape functions with
																 * respect to d.o.f.
																 * Rows represent d.o.f, columns represent
																 * derivatives in nodes */
} gauss_node;

/*
 * database of elements
 * Contains all gauss nodes for elements together with
 * derivatives
 */
typedef struct elements_database_tag {
		gauss_node **gauss_nodes;		/* Gauss nodes 2d array
																 * rows represent elements,
																 * columns are particular gauss nodes
																 * per element */
		/* For every gauss node */
} elements_database;


/*************************************************************/
/* Functions declarations declarations                            */



/* function for calculation value of shape function for 10-noded
 * tetrahedra 
 * by given element el, node number i,local coordinates r,s,t,  
 * where r,s,t from [0;1] 
 * all functions are taken from the book: 
 * "The Finite Element Method for 3D Thermomechanical Applications"
 * by - Guido Dhond p.72
 */
real isoform(int i,real r,real s,real t);


/*
 * function for calculation derivatives of shape 
 * function of 10noded tetrahedra element
 * with respect to local coordinate system
 * shape - number of node(and corresponding shape function)
 * dof - degree of freedom, dof = 1 is r, dof = 2 is s, dof = 3 is t
 * r,s,t is [0;1] - local coordinates
 */
real disoform(int shape,int dof,real r,real s,real t);
/* Particular derivatives */
real df_dr(int i, real r, real s, real t);
real df_ds(int i, real r, real s, real t);
real df_dt(int i, real r, real s, real t);


int main(int argc, char **argv)
{
		
		

		return 0;
}


real isoform(int i,real r,real s,real t)
{
		switch(i)
		{
		case 0: return (2*(1-r-s-t)-1)*(1-r-s-t);
		case 1: return (2*r-1)*r;
		case 2: return (2*s-1)*s;
		case 3: return (2*t-1)*t;
		case 4: return 4*r*(1-r-s-t);
		case 5: return 4*r*s;
		case 6: return 4*s*(1-r-s-t);
		case 7: return 4*t*(1-r-s-t);
		case 8: return 4*r*t;
		case 9: return 4*s*t;
		}
		return 0;
}


real disoform(int shape,int dof,real r,real s,real t)
{
		switch(dof)
		{
		case 0: return df_dr(shape,r,s,t);
		case 1: return df_ds(shape,r,s,t);
		case 2: return df_dt(shape,r,s,t);
		}
		return 0;
}


real df_dr(int i,real r,real s,real t)
{
		switch(i)
		{
		case 0: return 4*t+4*s+4*r-3;
		case 1: return 4*r-1;
		case 2: return 0;
		case 3: return 0;
		case 4: return -4*t-4*s-8*r+4;
		case 5: return 4*s;
		case 6: return -4*s;
		case 7: return -4*t;
		case 8: return 4*t;
		case 9:	return 0;
		}
		return 0;
}


real df_ds(int i, real r, real s, real t)
{
		switch(i)
		{
		case 0: return 4*t+4*s+4*r-3;
		case 1: return 0;
		case 2: return 4*s-1;
		case 3: return 0;
		case 4: return -4*r;
		case 5: return 4*r;
		case 6: return -4*t-8*s-4*r+4;
		case 7: return -4*t;
		case 8: return 0;
		case 9:	return 4*t;
		}
		return 0;
}


real df_dt(int i, real r, real s, real t)
{
		switch(i)
		{
		case 0: return 4*t+4*s+4*r-3;
		case 1: return 0;
		case 2: return 0;
		case 3: return 4*t-1;
		case 4: return -4*r;
		case 5: return 0;
		case 6: return -4*s;
		case 7: return -8*t-4*s-4*r+4;
		case 8: return 4*r;
		case 9:	return 4*s;
		}
		return 0;
}

