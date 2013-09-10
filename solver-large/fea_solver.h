/* -*- Mode: C; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#ifndef __FEA_SOLVER_H__
#define __FEA_SOLVER_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>
#include "defines.h"
#include "sp_matrix.h"
#include "sp_direct.h"
#include "dense_matrix.h"

/* default value of the tolerance for the iterative solvers */
#define MAX_ITERATIVE_TOLERANCE 1e-14
/* default value of the max number of iterations for the iterative solvers */
#define MAX_ITERATIVE_ITERATIONS 20000

/*************************************************************/
/* Global variables                                          */

typedef struct fea_solver_tag* fea_solver_ptr;

/*************************************************************/
/* Function pointers declarations                            */

/*
 * A pointer to the the isoparametric shape
 * function
 */
typedef real (*isoform_t)(int i,real r,real s,real t);

/*
 * A pointer to the derivative of the isoparametric shape
 * function
 */
typedef real (*disoform_t)(int shape,int dof,real r,real s,real t);

/*
 * A pointer to the function for exporting data from solver
 */ 
typedef void (*export_solution_t) (fea_solver_ptr, const char *filename);

/*
 * A pointer to the function for appling single BC for global
 * DOF[index], using argument arg
 */
typedef void (*apply_bc_t) (fea_solver_ptr self, int index, real arg);

/*
 * A pointer to the function for calculating Cauchy stresses by given
 * deformation gradient
 */
typedef void (*element_gauss_stress_t)(fea_solver_ptr self,
                                       real graddef_tensor[MAX_DOF][MAX_DOF],
                                       real stress_tensor[MAX_DOF][MAX_DOF]);
/*
 * A pointer to the function for calculating the C elasticity tensor
 * by given deformation gradient
 */
typedef void (*solver_ctensor_t)(fea_solver_ptr self,
                                 real (*graddef)[MAX_DOF],
                                 real (*ctensor)[MAX_DOF][MAX_DOF][MAX_DOF]);


/*************************************************************/
/* Enumerations declarations                                 */

typedef enum  {
  /* PLANE_STRESS, PLANE_STRAIN, AXISYMMETRIC,  */
  CARTESIAN3D
} task_type;

typedef enum {
  MODEL_A5,
  MODEL_COMPRESSIBLE_NEOHOOKEAN
} model_type;

typedef enum {
  CG,
  PCG_ILU,
  CHOLESKY
} slae_solver_type;
  
typedef enum  {
  /* TRIANGLE3, TRIANGLE6,TETRAHEDRA4, */
  TETRAHEDRA10
} element_type;


typedef enum  {
  FREE = 0,                    /* free */
  PRESCRIBEDX = 1,             /* x prescribed */
  PRESCRIBEDY = 2,             /* y prescribed */
  PRESCRIBEDXY = 3,            /* x, y prescribed */
  PRESCRIBEDZ = 4,             /* z prescribed*/
  PRESCRIBEDXZ = 5,            /* x, z prescribed*/
  PRESCRIBEDYZ = 6,            /* y, z prescribed*/
  PRESCRIBEDXYZ = 7            /* x, y, z prescribed.*/
} presc_boundary_type;


/*************************************************************/
/* Data structures                                           */

typedef struct {
  model_type model;                         /* model type */
  real parameters[MAX_MATERIAL_PARAMETERS]; /* model material parameters */
  int parameters_count;                     /* number of material params */
} fea_model;
typedef fea_model* fea_model_ptr;

/*
 * Task type declaration.
 * Defines an input parameters for the task, independent of 
 * the input geometry and loads
 */
typedef struct {
  task_type type;               /* type of the task to solve */
  fea_model model;              /* material model */
  slae_solver_type solver_type; /* SLAE solver */
  real solver_tolerance;        /* tolerance in case of iterative solver */
  int solver_max_iter;          /* max number of iters for iterative solver */
  int dof;                      /* number of degree of freedom */
  element_type ele_type;        /* type of the element */
  int load_increments_count;    /* number of load increments */
  real desired_tolerance;       /* desired energy tolerance */
  int max_newton_count;         /* maximum number of Newton's iterations */
  int linesearch_max;           /* maximum number of line searches */
  int arclength_max;            /* maximum number of arc lenght searches */
  BOOL modified_newton;         /* use modified Newton's method or not */
  const char* export_file;      /* export file name - guessing from input */
} fea_task;
typedef fea_task* fea_task_ptr;

/* Calculated solution parameters */
typedef struct {
  int nodes_per_element;        /* number of nodes defined in element 
                                   based on fea_task::ele_type */
  int gauss_nodes_count;        /* number of gauss nodes per element */
} fea_solution_params;
typedef fea_solution_params* fea_solution_params_ptr;

/*************************************************************/
/* Input geometry parameters                                 */

/* An array of nodes. */
typedef struct {
  int nodes_count;      /* number of input nodes */
  real **nodes;         /* nodes array,sized as nodes_count x MAX_DOF
                         * so access is  nodes[node_number][dof] */
} nodes_array;
typedef nodes_array* nodes_array_ptr;

/* An array of elements */
typedef struct {
  int elements_count;           /* number of elements */
  int **elements;               /* elements array, each line represents an
                                 * element. Element is an array of node
                                 * indexes
                                 */
} elements_array;
typedef elements_array* elements_array_ptr;

/* Particular prescribed boundary node */
typedef struct {
  int node_number;
  real values[MAX_DOF];
  presc_boundary_type type;
} prescribed_bnd_node;
typedef prescribed_bnd_node* prescribed_bnd_node_ptr;

/*
 * An array of prescribed boundary conditions
 * either fixed and with prescribed displacements
 */
typedef struct {
  int prescribed_nodes_count;
  prescribed_bnd_node* prescribed_nodes;
} presc_bnd_array;
typedef presc_bnd_array* presc_bnd_array_ptr;

/*************************************************************/
/* Application-specific structures                           */

/* Structure describing information for the gauss node
 * depending on number of shape functions N per element
 * TODO: add tables with layouts in comments */
typedef struct {
  real weight;                /* weight for the integration */
  real *forms;                /* shape function values for gauss node, N */
  real **dforms;              /* derivatives of shape functions with
                               * respect to d.o.f.
                               * Rows represent d.o.f, columns represent
                               * derivatives in nodes */
} gauss_node;
typedef gauss_node* gauss_node_ptr;

/*
 * database of elements
 * Contains all gauss nodes for elements together with
 * derivatives
 */
typedef struct {
  real (*gauss_nodes_data)[4];  /* pointer to an array of gauss
                                 * coefficients and weights */  
  gauss_node_ptr *gauss_nodes;   /* Gauss nodes array */
} elements_database;
typedef elements_database* elements_database_ptr;


/*
 * Matrix of gradients of the shape functions with respect to
 * global coordinates
 * Calculated in one node(gauss node) of an element
 * Layout: [d.o.f x nodes_count]
 *              dN_j(r,s,t)
 * grad[i][j] = -----------
 *                  dX_i
 * 
 * X_1 = x, X_2 = y, X_3 = z coordinates
 * 
 */
typedef struct {
  real **grads;                  /* array of derivatives of shape functions
                                  * [dof x nodes_per_element] */
  real detJ;                     /* determinant of Jacobi matrix */
} shape_gradients;
typedef shape_gradients* shape_gradients_ptr;

/*
 * Load increment step structure.
 * This structure holds all information needed on the current
 * load increment step
 */
typedef struct {
  int step_number;
  nodes_array_ptr nodes_p;      /* nodes in current configuration for step */
  tensor **graddefs;            /* Components of Deformation gradient tensor
                                 * in gauss nodes
                                 * array [number of elems] x [gauss nodes]
                                 */
  tensor **stresses;            /* Components of Cauchy stress tensor
                                 * in gauss nodes
                                 * array [number of elems] x [gauss nodes]
                                 */
} load_step;
typedef load_step* load_step_ptr;

/*
 * A main application structure which shall contain all
 * data necessary for solution
 */
typedef struct fea_solver_tag {
  disoform_t dshape;              /* a function pointer to derivative of the
                                   * shape function */
  isoform_t shape;                /* a function pointer to the shape
                                   * function */
  element_gauss_stress_t element_gauss_stress; /* a pointer to the function
                                                * for calculation of the
                                                * model-specific Cauchy
                                                * in gauss nodes per element */
  solver_ctensor_t ctensor; /* a function pointer to the C elasticity
                             * tensor
                             */
  export_solution_t export_function; /* a pointer to the export function */

  fea_task_ptr task_p;               
  fea_solution_params_ptr fea_params_p; 
  nodes_array_ptr nodes0_p;
  nodes_array_ptr nodes_p;              
  elements_array_ptr elements_p;
  presc_bnd_array_ptr presc_boundary_p;
  elements_database elements_db;  /* array of pre-constructed
                                   * values of derivatives of the
                                   * isoparametric shape functions
                                   * in gauss nodes */
  shape_gradients_ptr **shape_gradients0; /* an array of gradients of
                                           * shape functios per element
                                           * per gauss node in initial
                                           * configuration
                                           * [number of elems] x [gauss nodes]  
                                           */
  shape_gradients_ptr **shape_gradients; /* an array of gradients of
                                          * shape functios per element
                                          * per gauss node in current
                                          * configuration
                                          * [number of elems] x [gauss nodes]
                                          * shall be calculated after update of
                                          * nodes
                                          */

  tensor **graddefs;            /* Components of Deformation gradient tensor
                                 * in gauss nodes
                                 * array [number of elems] x [gauss nodes]
                                 */
  tensor **stresses;            /* Components of Cauchy stress tensor
                                 * in gauss nodes
                                 * array [number of elems] x [gauss nodes]
                                 */
  int current_load_step;
  load_step_ptr load_steps_p;   /* an array of stored load steps data
                                 * array size is task_p->load_increments_count
                                 * load_steps_p[0..current_load_step] shall be
                                 * filled during load steps iterations
                                 */
  sp_matrix global_mtx;         /* global stiffness matrix */
  sp_chol_symbolic_ptr symb_chol; /* symbolic Cholesky decomposition
                                   * of the global stiffness matrix
                                   */
  real* global_forces_vct;      /* external forces vector */
  real* global_reactions_vct;   /* reactions in fixed dofs */
  real* global_solution_vct;    /* vector of global solution */
} fea_solver;


/*************************************************************/
/* Functions declarations                                    */

/* error exit from the application with message */
void error(char* msg);

/*
 * Load initial data from file
 */
BOOL initial_data_load(char *filename,
                       fea_task_ptr *task,
                       fea_solution_params_ptr *fea_params,
                       nodes_array_ptr *nodes,
                       elements_array_ptr *elements,
                       presc_bnd_array_ptr *presc_boundary);


/*************************************************************/
/* Allocators for internal data structures                   */

/*************************************************************/
/* Allocators for structures with a data from file           */

/* Initializa fea task structure and fill with default values */
fea_task_ptr fea_task_alloc();
/* Initializes fea solution params with default values */
fea_solution_params_ptr fea_solution_params_alloc();
/* Initialize nodes array but not initialize particular arrays  */
nodes_array_ptr nodes_array_alloc();
/* create a copy of nodes array */
nodes_array_ptr nodes_array_copy_alloc(nodes_array_ptr nodes);
/* Initialize elements array but not initialize particular elements */
elements_array_ptr elements_array_alloc();
/* Initialize boundary nodes array but not initialize particular nodes */
presc_bnd_array_ptr presc_bnd_array_alloc();

/*
 * Constructor for the main application structure
 * all parameters shall be properly constructed and initialized
 * with data from file
 */
fea_solver_ptr fea_solver_alloc(fea_task_ptr task,
                                fea_solution_params_ptr fea_params,
                                nodes_array_ptr nodes,
                                elements_array_ptr elements,
                                presc_bnd_array_ptr presc);


/*************************************************************/
/* Deallocators for internal data structures                 */

fea_solution_params_ptr fea_solution_params_free(fea_solution_params_ptr ptr);
fea_task_ptr fea_task_free(fea_task_ptr task);
nodes_array_ptr nodes_array_free(nodes_array_ptr nodes);
elements_array_ptr elements_array_free(elements_array_ptr elements);
presc_bnd_array_ptr presc_bnd_array_free(presc_bnd_array_ptr presc);

/*
 * Destructor for the main solver
 * Will also clear all aggregated structures
 */
fea_solver_ptr fea_solver_free(fea_solver_ptr solver);

/*
 * Constructor for the shape functions gradients array
 * element - index of the element to calculate in
 * gauss - index of gauss node
 */
shape_gradients_ptr solver_shape_gradients_alloc(fea_solver_ptr self,
                                                 nodes_array_ptr nodes,
                                                 int element,
                                                 int gauss);
/* Destructor for the shape gradients array */
shape_gradients_ptr solver_shape_gradients_free(fea_solver_ptr self,
                                                shape_gradients_ptr grads);

/*
 * fills the self->shape_gradients or self->shape_gradients0 array
 * depending on flag current
 */
void solver_create_shape_gradients(fea_solver_ptr self,BOOL current);

/*
 * Returns a particular component of a node with local index 'node'
 * in element with index 'element' for the d.o.f. 'dof'
 */
real solver_node_dof(fea_solver_ptr self,
                     nodes_array_ptr nodes,
                     int element,
                     int node,
                     int dof);



/*
 * Constructor for the load step structure
 * It doesn't allocate a memory for a step itself,
 * just initializes the internal structures
 */
void solver_load_step_init(fea_solver_ptr self,
                           load_step_ptr step,
                           int step_number);
            
/*
 * Desctructor for the load step structure.
 * it doesn't deallocate a memory for a step itself,
 * just frees all step internal structures
 */
void solver_load_step_free(fea_solver_ptr self, load_step_ptr step);


/*************************************************************/
/* Functions describing concrete element types               */

/* function for calculation value of shape function for 10-noded
 * tetrahedra 
 * by given element el, node number i,local coordinates r,s,t,  
 * where r,s,t from [0;1] 
 * all functions are taken from the book: 
 * "The Finite Element Method for 3D Thermomechanical Applications"
 * by - Guido Dhond p.72
 */
real tetrahedra10_isoform(int i,real r,real s,real t);


/*
 * Element TETRAHEDRA10, isoparametric shape function derivative
 * with respece to the 1st variable r
 * i - node number
 */ 
real tetrahedra10_df_dr(int i,real r,real s,real t);

/*
 * Element TETRAHEDRA10, isoparametric shape function derivative
 * with respece to the 2nd variable s
 * i - node number
 */ 
real tetrahedra10_df_ds(int i, real r, real s, real t);

/*
 * Element TETRAHEDRA10, isoparametric shape function derivative
 * with respece to the 3rd variable t
 * i - node number
 */ 
real tetrahedra10_df_dt(int i, real r, real s, real t);

/*
 * function for calculation derivatives of shape 
 * function of 10noded tetrahedra element
 * with respect to local coordinate system
 * shape - number of node(and corresponding shape function)
 * dof - degree of freedom, dof = 1 is r, dof = 2 is s, dof = 3 is t
 * r,s,t is [0;1] - local coordinates
 */
real tetrahedra10_disoform(int shape,int dof,real r,real s,real t);


/*************************************************************/
/* Functions particular material models                      */

/*
 * Calculate stress tensor of the material model A5 by given
 * deformation gradient
 */
void solver_element_gauss_stress_A5(fea_solver_ptr self,
                                    real F[MAX_DOF][MAX_DOF],
                                    real stress_tensor[MAX_DOF][MAX_DOF]);

/*
 * Calculate stress tensor of the Neo-hookean compressible model by given
 * deformation gradient
 */
void solver_element_gauss_stress_compr_neohookean(fea_solver_ptr self,
                                                  real F[MAX_DOF][MAX_DOF],
                                                  real S[MAX_DOF][MAX_DOF]);
/*
 * Calculate 4th rank tensor C of the elastic material model A5
 * T = C(4)**S
 * S - deformation tensor
 * by given deformation gradient
 */
void solver_ctensor_A5(fea_solver_ptr self,
                       real (*graddef)[MAX_DOF],
                       real (*ctensor)[MAX_DOF][MAX_DOF][MAX_DOF]);
/*
 * Calculate 4th rank tensor C of the Neo-hookean compressible material model
 * T = C(4)**S
 * S - deformation tensor 
 * by given deformation gradient
 */
void solver_ctensor_compr_neohookean(fea_solver_ptr self,
                                     real (*graddef)[MAX_DOF],
                                     real (*ctensor)[MAX_DOF][MAX_DOF][MAX_DOF]);

/*************************************************************/
/* Functions for exporting data in different formats         */
void solver_export_tetrahedra10_gmsh(fea_solver_ptr solver,
                                     const char *filename);


/*************************************************************/
/* General functions                                         */

/*
 * Parse command line parameters and return input file name
 * into the filename variable
 */
int parse_cmdargs(int argc, char **argv,char **filename);

/*
 * Real main function with input filename as a parameter
 */
int do_main(char* filename);


/*
 * Solver function which shall be called
 * when all data read to an appropriate structures
 */
void solve(fea_task_ptr task,
           fea_solution_params_ptr fea_params,
           nodes_array_ptr nodes,
           elements_array_ptr elements,
           presc_bnd_array_ptr presc_boundary);

/*
 * Solver wrapper function to solve SLAE
 */
BOOL solver_solve_slae(fea_solver_ptr solver);

#ifdef DUMP_DATA
/* Dump input data to check if parser works correctly */
void dump_input_data( char* filename,
                      fea_task_ptr task,
                      fea_solution_params_ptr fea_params,
                      nodes_array_ptr nodes,
                      elements_array_ptr elements,
                      presc_bnd_array_ptr presc_boundary);
#endif

/*************************************************************/
/* Functions for fea_solver structure                        */

/*
 * Fills solver with pointers to functions and pointers to arrays
 * of gauss nodes for particular element type.
 * This function is used on a construction phase
 */
void solver_create_element_params(fea_solver_ptr self);

/*
 * Initialize material model functions
 */
void solver_create_model_params(fea_solver_ptr self);

/* Allocates memory and construct elements database for solver */
void solver_create_element_database(fea_solver_ptr self);
/* Destructor for the element database */
void solver_free_element_database(fea_solver_ptr self);


/*
 * Create an array of shape functions gradients
 * in gauss nodes per elements in initial configuration
 */
void solver_create_initial_shape_gradients(fea_solver_ptr self);

/*
 * Create an array of shape functions gradients
 * in gauss nodes per elements in current configuration
 */
void solver_create_current_shape_gradients(fea_solver_ptr self);

/* Create components of the Cauchy stresses in gauss nodes  */
void solver_create_stresses(fea_solver_ptr self);

/* Create global residual forces vector */
void solver_create_residual_forces(fea_solver_ptr self);

/* Create global stiffness matrix */
void solver_create_stiffness(fea_solver_ptr self);


/* Update global forces vector with residual forces for the element */
void solver_local_residual_forces(fea_solver_ptr self,int element);


/* Create constitutive component of the stiffness matrix */
void solver_local_constitutive_part(fea_solver_ptr self,int element);

/* Create initial stress component of the stiffness matrix */
void solver_local_initial_stess_part(fea_solver_ptr self,int element);


/*
 * Calculate Deformation gradient in gauss node
 * element - element number
 * gauss - gauss node number in element
 * graddef - MAX_DOF x MAX_DOF array of components of Deformation
 * gradient tensor
 */
void solver_element_gauss_graddef(fea_solver_ptr self,
                                  int element,
                                  int gauss,
                                  real graddef[MAX_DOF][MAX_DOF]);


/*
 * Calculate stress tensor in gauss node
 * element - element number
 * gauss - gauss node number in element
 * graddef - MAX_DOF x MAX_DOF array of components of Deformation gradient
 * stress - MAX_DOF x MAX_DOF array of components of Cauchy stress tensor
 */
void solver_element_gauss_stress(fea_solver_ptr self,
                                 int element,
                                 int gauss,
                                 real graddef_tensor[MAX_DOF][MAX_DOF],
                                 real stress_tensor[MAX_DOF][MAX_DOF]);



/* Fill greate global forces vector */
void solver_create_forces_bc(fea_solver_ptr self);

/*
 * Apply BC in form of prescribed displacements
 * lambda - multiplier for the prescribed displacements
 */
void solver_apply_prescribed_bc(fea_solver_ptr self,real lambda);

/*
 * Call the 'apply' function for every prescribed displacement
 * with an argument lambda
 */
void solver_apply_bc_general(fea_solver_ptr self,
                             apply_bc_t apply,
                             real lambda);

/* Apply BC in form of prescribed displacements to a single specified
 * global d.o.f.
 * This function is called from solver_apply_prescribed_bc
 */
void solver_apply_single_bc(fea_solver_ptr self,
                            int index, real value);

/* Add BC in form of prescribed displacements to a single specified
 * global d.o.f. of a global nodes vector
 * This function is called from solver_update_nodes_with_bc */
void solver_update_node_with_bc(fea_solver_ptr self,
                                int index,
                                real value);

/*
 * Update array solver->nodes with displacements from vector x
 * This function may be used for:
 * 1) calculation of the deformation gradients
 * 2) create deformed configuration by applying prescribed displacements
 */
void solver_update_nodes_with_solution(fea_solver_ptr self,
                                       real* x);

/*
 * Update solver->nodes array with prescribed displacements
 * lambda - multiplier for prescribed displacements
 */
void solver_update_nodes_with_bc(fea_solver_ptr self, real lambda);

/*
 * Creates a particular gauss node for the element
 * with index element_index and gauss node number gauss_node_index
 */
gauss_node_ptr solver_gauss_node_alloc(fea_solver_ptr self,
                                       int gauss_node_index);

/* Deallocate gauss node */
gauss_node_ptr solver_gauss_node_free(fea_solver_ptr self,
                                      gauss_node_ptr node);




#endif /* __FEA_SOLVER_H__ */
