/* -*- Mode: C; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#ifdef USE_EXPAT
#include <expat.h>
#endif

/*************************************************************/
/* Type and constants definitions                            */

/* BOOL type as usual */
typedef int BOOL;
#define FALSE 0
#define TRUE 1

#define MAX_DOF 3
#define MAX_MATERIAL_PARAMETERS 10

/*
 * SLAE solver constants
 * possibly shall go to the initial data in future
 */
const int MAX_ITER = 10000;
const double TOLERANCE = 1e-10;

/* Redefine type of the floating point values */
#ifdef SINGLE
typedef float real;
#else /* not SINGLE */
typedef double real;
#endif /* not SINGLE */

/* Equals macro for real values */
#ifdef SINGLE
#define EQL(x,y) (fabs((x)-(y))<= (FLT_MIN))
BOOL eql(float x,float y);
inline BOOL eql(float x,float y) {return (fabs((x)-(y))<= (FLT_MIN))?TRUE:FALSE;}
#else
#define EQL(x,y) (fabs((x)-(y))<= (DBL_MIN))
BOOL eql(double x,double y);
BOOL eql(double x,double y) {return (fabs((x)-(y))<= (DBL_MIN))?TRUE:FALSE;}
#endif

/* Kroneker delta */
#define DELTA(i,j) ((i)==(j) ? 1 : 0)

/* shortcut for adding of the matrix elements */
#define MTX(m,i,j,v) sp_matrix_element_add((m),(i),(j),(v));

/*************************************************************/
/* Global variables                                          */

extern int errno;
struct fea_solver_tag* global_solver;
typedef struct fea_solver_tag* fea_solver_ptr;

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
typedef void (*export_solution_t) (fea_solver_ptr, char *filename);

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


/*
 * arrays of gauss nodes with coefficients                   
 * layout: [number_of_nodes x 4], with values:               
 * {weight, r,s,t}                                           
 * per gauss node. For 2d cases t = 0
 * Note what divisor 6 for tetraheadras and 2 for triangles
 * shall be already taken into account in weights.
 * See below.
 */

/* Element: TETRAHEDRA10, 4 nodes */
real gauss_nodes4_tetr10[4][4] = { {(1/4.)/6.,   /* weight */
                                    0.58541020,  /* a */
                                    0.13819660,  /* b */
                                    0.13819660}, /* b */
                                   {(1/4.)/6.,
                                    0.13819660,  /* b */
                                    0.58541020,  /* a */
                                    0.13819660}, /* b */
                                   {(1/4.)/6.,
                                    0.13819660,  /* b */
                                    0.13819660,  /* b */
                                    0.58541020}, /* a */
                                   {(1/4.)/6.,
                                    0.13819660,  /* b */
                                    0.13819660,  /* b */
                                    0.13819660}  /* b */
};
/* Element: TETRAHEDRA10, 5 nodes */
real gauss_nodes5_tetr10[5][4] = { {(-4/5.)/6., 1/4., 1/4., 1/4.},
                                   {(9/20.)/6., 1/2., 1/6., 1/6.},
                                   {(9/20.)/6., 1/6., 1/2., 1/6.},
                                   {(9/20.)/6., 1/6., 1/6., 1/2.},
                                   {(9/20.)/6., 1/6., 1/6., 1/6.} };



typedef enum task_type_enum {
  /* PLANE_STRESS, PLANE_STRAIN, AXISYMMETRIC,  */
  CARTESIAN3D
} task_type;

typedef enum model_type_enum {
  MODEL_A5,
  MODEL_COMPRESSIBLE_NEOHOOKEAN
} model_type;
  
typedef enum element_type_enum {
  /* TRIANGLE3, TRIANGLE6,TETRAHEDRA4, */
  TETRAHEDRA10
} element_type;


typedef enum presc_boundary_type_enum {
  FREE = 0,                    /* free */
  PRESCRIBEDX = 1,             /* x prescribed */
  PRESCRIBEDY = 2,             /* y prescribed */
  PRESCRIBEDXY = 3,            /* x, y prescribed */
  PRESCRIBEDZ = 4,             /* z prescribed*/
  PRESCRIBEDXZ = 5,            /* x, z prescribed*/
  PRESCRIBEDYZ = 6,            /* y, z prescribed*/
  PRESCRIBEDXYZ = 7            /* x, y, z prescribed.*/
} presc_boundary_type;

typedef enum sparse_storage_type_enum {
  CRS = 0,                      /* Compressed Row Storage */
  CCS = 1                       /* Compressed Column Storage */
} sparse_storage_type;

typedef struct fea_material_model_tag {
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
typedef struct fea_task_tag {
  task_type type;               /* type of the task to solve */
  fea_model model;              /* material model */
  unsigned char dof;            /* number of degree of freedom */
  element_type ele_type;        /* type of the element */
  int load_increments_count;    /* number of load increments */
  real desired_tolerance;       /* desired energy tolerance */
  int linesearch_max;           /* maximum number of line searches */
  int arclength_max;            /* maximum number of arc lenght searches */
  BOOL modified_newton;         /* use modified Newton's method or not */

} fea_task;
typedef fea_task* fea_task_ptr;

/* Calculated solution parameters */
typedef struct fea_solution_params_tag {
  int nodes_per_element;        /* number of nodes defined in element 
                                   based on fea_task::ele_type */
  int gauss_nodes_count;        /* number of gauss nodes per element */
} fea_solution_params;
typedef fea_solution_params* fea_solution_params_ptr;

/*************************************************************/
/* Input geometry parameters                                 */

/* An array of nodes. */
typedef struct nodes_array_tag {
  int nodes_count;      /* number of input nodes */
  real **nodes;         /* nodes array,sized as nodes_count x MAX_DOF
                         * so access is  nodes[node_number][dof] */
} nodes_array;
typedef nodes_array* nodes_array_ptr;

/* An array of elements */
typedef struct elements_array_tag {
  int elements_count;           /* number of elements */
  int **elements;               /* elements array, each line represents an
                                 * element. Element is an array of node
                                 * indexes
                                 */
} elements_array;
typedef elements_array* elements_array_ptr;

/* Particular prescribed boundary node */
typedef struct prescibed_boundary_node_tag {
  int node_number;
  real values[MAX_DOF];
  presc_boundary_type type;
} prescibed_boundary_node;
typedef prescibed_boundary_node* prescibed_boundary_node_ptr;

/*
 * An array of prescribed boundary conditions
 * either fixed and with prescribed displacements
 */
typedef struct presc_boundary_array_tag {
  int prescribed_nodes_count;
  prescibed_boundary_node* prescribed_nodes;
} presc_boundary_array;
typedef presc_boundary_array* presc_boundary_array_ptr;

/*************************************************************/
/* Application-specific structures                           */

/* Structure describing information for the gauss node
 * depending on number of shape functions N per element
 * TODO: add tables with layouts in comments */
typedef struct gauss_node_tag {
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
typedef struct elements_database_tag {
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
typedef struct shape_gradients_tag {
  real **grads;                  /* array of derivatives of shape functions
                                  * [dof x nodes_per_element] */
  real detJ;                     /* determinant of Jacobi matrix */
} shape_gradients;
typedef shape_gradients* shape_gradients_ptr;

/* A storage to keep 2nd rank tensor components  */
typedef struct tensor_tag {
  real components[MAX_DOF][MAX_DOF];
} tensor;
typedef tensor* tensor_ptr;

/*
 * Load increment step structure.
 * This structure holds all information needed on the current
 * load increment step
 */
typedef struct load_step_tag {
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
 * Sparse matrix row/column storage array
 */
typedef struct indexed_array_tag {
  int width;                    /* size of an array */
  int last_index;               /* last stored index, i.e. if width = 20
                                 * it will be 9 if only 10 nonzero elements
                                 * stored */
  int  *indexes;                /* array of column/row indexes */
  real *values;                 /* array of values */
} indexed_array;
typedef indexed_array* indexed_array_ptr;

/*
 * Sparse matrix row storage 
 * Internal format based on CRS or CCS
 */
typedef struct sp_matrix_tag {
  int rows_count;
  int cols_count;
  indexed_array* storage;
  BOOL ordered;                              /* if matrix was finalized */
  sparse_storage_type storage_type;          /* Storage type */
} sp_matrix;
typedef sp_matrix* sp_matrix_ptr;

/*
 * Sparse matrix CSLR(Skyline) format
 * used in sparse iterative solvers
 * Constructed from Sparse Matrix in assumption of the symmetric
 * matrix portrait
 */
typedef struct sp_matrix_skyline_tag {
  int rows_count;
  int cols_count;
  int nonzeros;                 /* number of nonzero elements in matrix */
  int tr_nonzeros;  /* number of nonzero elements in
                     * upper or lower triangles */
  real *diag;                   /* rows_count elements in matrix diagonal */
  real *lower_triangle;         /* nonzero elements of the lower triangle */
  real *upper_triangle;         /* nonzero elements of the upper triangle */
  int *jptr;                    /* array of column/row indexes of the
                                 * lower/upper triangles */
  int *iptr;                    /* array of row/column offsets in jptr
                                 * for lower or upper triangles */
} sp_matrix_skyline;
typedef sp_matrix_skyline* sp_matrix_skyline_ptr;

/*
 * ILU decomposition of the sparse matrix in Skyline (CSLR) format
 * ILU decomposition keeps the symmetric portrait of the sparse matrix
 */
typedef struct sp_matrix_skyline_ilu_tag {
  sp_matrix_skyline parent;
  real *ilu_diag;              /* U matrix diagonal */
  real *ilu_lowertr;           /* nonzero elements of the lower(L) matrix */
  real *ilu_uppertr;           /* nonzero elements of the upper(U) matrix */
} sp_matrix_skyline_ilu;
typedef sp_matrix_skyline_ilu* sp_matrix_skyline_ilu_ptr;

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
  export_solution_t export;     /* a pointer to the export function */

  fea_task_ptr task_p;               
  fea_solution_params_ptr fea_params_p; 
  nodes_array_ptr nodes0_p;
  nodes_array_ptr nodes_p;              
  elements_array_ptr elements_p;
  presc_boundary_array_ptr presc_boundary_p;
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
  sp_matrix global_mtx;     /* global stiffness matrix */
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
                       presc_boundary_array_ptr *presc_boundary);


/*************************************************************/
/* Allocators for internal data structures                   */

/*************************************************************/
/* Allocators for structures with a data from file           */

/* Initializa fea task structure and fill with default values */
fea_task_ptr new_fea_task();
/* Initializes fea solution params with default values */
fea_solution_params_ptr new_fea_solution_params();
/* Initialize nodes array but not initialize particular arrays  */
nodes_array_ptr new_nodes_array();
/* create a copy of nodes array */
nodes_array_ptr new_copy_nodes_array(nodes_array_ptr nodes);
/* Initialize elements array but not initialize particular elements */
elements_array_ptr new_elements_array();
/* Initialize boundary nodes array but not initialize particular nodes */
presc_boundary_array_ptr new_presc_boundary_array();

/*
 * Constructor for the main application structure
 * all parameters shall be properly constructed and initialized
 * with data from file
 */
fea_solver_ptr new_fea_solver(fea_task_ptr task,
                              fea_solution_params_ptr fea_params,
                              nodes_array_ptr nodes,
                              elements_array_ptr elements,
                              presc_boundary_array_ptr presc);


/*************************************************************/
/* Deallocators for internal data structures                 */

void free_fea_solution_params(fea_solution_params_ptr params);
void free_fea_task(fea_task_ptr task);
void free_nodes_array(nodes_array_ptr nodes);
void free_elements_array(elements_array_ptr elements);
void free_presc_boundary_array(presc_boundary_array_ptr presc);

/*
 * Destructor for the main solver
 * Will also clear all aggregated structures
 */
void free_fea_solver(fea_solver_ptr solver);

/*
 * Constructor for the shape functions gradients array
 * element - index of the element to calculate in
 * gauss - index of gauss node
 */
shape_gradients_ptr solver_new_shape_gradients(fea_solver_ptr self,
                                               nodes_array_ptr nodes,
                                               int element,
                                               int gauss);
/* Destructor for the shape gradients array */
void solver_free_shape_gradients(fea_solver_ptr self,
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
void init_load_step(fea_solver_ptr self,
                    load_step_ptr step,
                    int step_number);
            
/*
 * Desctructor for the load step structure.
 * it doesn't deallocate a memory for a step itself,
 * just frees all step internal structures
 */
void free_load_step(fea_solver_ptr self, load_step_ptr step);


/*************************************************************/
/* Sparse matrix operations                                  */

/* indexed_arrays operations */
/* Swap i and j elements in the indexed array.
 * Used in indexed_array_sort*/
void indexed_array_swap(indexed_array_ptr self,int i, int j);
/* Performs in-place sort of the indexed array */
void indexed_array_sort(indexed_array_ptr self, int l, int r);
/* Print contents of the indexed array to the stdout  */
void indexed_array_printf(indexed_array_ptr self);

/*
 * Initializer for a sparse matrix with specified rows and columns
 * number.
 * This function doesn't allocate the memory for the matrix itself;
 * only for its structures. Matrix mtx shall be already allocated
 * bandwdith - is a start bandwidth of a matrix row
 * type - CRS or CCS sparse matrix storage types
 */
void init_sp_matrix(sp_matrix_ptr mtx,
                    int rows,
                    int cols,
                    int bandwidth,
                    sparse_storage_type type);
/*
 * Destructor for a sparse matrix
 * This function doesn't deallocate memory for the matrix itself,
 * only for its structures.
 */
sp_matrix_ptr free_sp_matrix(sp_matrix_ptr mtx);

/*
 * Clear the sparse matrix.
 * Set the element values to zero keeping sparsity portrait
 */
void clear_sp_matrix(sp_matrix_ptr mtx);

/*
 * Copy sparse matrix from mtx_from to mtx_to
 * This function assumes what mtx_to is already cleared by free_sp_matrix
 * or mtx_to is a pointer to uninitialized sp_matrix structure
 */
void copy_sp_matrix(sp_matrix_ptr mtx_from,
                    sp_matrix_ptr mtx_to);

/*
 * Converts matrix storage format CRS <=> CCS
 * mtx_to shall be uninitialized sp_matrix structure
 */
void sp_matrix_convert(sp_matrix_ptr mtx_from,
                       sp_matrix_ptr mtx_to,
                       sparse_storage_type type);

/*
 * Creates ILU decomposition of the sparse matrix 
 */
void sp_matrix_create_ilu(sp_matrix_ptr self,sp_matrix_skyline_ilu_ptr ilu);

/*
 * Construct CSLR sparse matrix based on sp_matrix format
 * mtx - is the (reordered) sparse matrix to take data from
 * Acts as a copy-constructor
 */
void init_sp_matrix_skyline(sp_matrix_skyline_ptr self,
                            sp_matrix_ptr mtx);
/*
 * Destructor for a sparse matrix in CSLR format
 * This function doesn't deallocate memory for the matrix itself,
 * only for its structures.
 */
void free_sp_matrix_skyline(sp_matrix_skyline_ptr self);

/* getters/setters for a sparse matrix */

/* returns a pointer to the specific element
 * zero pointer if not found */
real* sp_matrix_element(sp_matrix_ptr self,int i, int j);
/* adds an element value to the matrix node (i,j) and return (i,j) */
real sp_matrix_element_add(sp_matrix_ptr self,
                           int i, int j, real value);

/* rearrange columns of a matrix to prepare for solving SLAE */
void sp_matrix_compress(sp_matrix_ptr self);

/*
 * Implements BLAS level 2 function SAXPY: y = A*x+b
 * All vectors shall be already allocated
 void sp_matrix_saxpy((sp_matrix_ptr self,real* b,real* x, real* y);
*/

/* Matrix-vector multiplication
 * y = A*x*/
void sp_matrix_mv(sp_matrix_ptr self,real* x, real* y);

/*
 * Solves SLAE L*x = b
 * by given L sparse matrix 
 * n - is the size of the x vector, and therefore
 * the matrix L will be used up to nth row & column.
 */
void sp_matrix_lower_solve(sp_matrix_ptr self,
                           int n,
                           real* b,
                           real* x);


/*
 * Solve SLAE for a matrix self with right-part b
 * Store results to the vector x. It shall be already allocated
 */
void sp_matrix_solve(sp_matrix_ptr self,real* b,real* x);
/*
 * Conjugate Grade solver
 * self - matrix
 * b - right-part vector
 * x0 - first approximation of the solution
 * max_iter - pointer to maximum number of iterations, MAX_ITER if zero;
 * will contain a number of iterations passed
 * tolerance - pointer to desired tolerance value, TOLERANCE if zero;
 * will contain norm of the residual at the end of iteration
 * x - output vector
 */
void sp_matrix_solve_cg(sp_matrix_ptr self,
                        real* b,
                        real* x0,
                        int* max_iter,
                        real* tolerance,
                        real* x);

/*
 * Preconditioned Conjugate Grade solver
 * Preconditioner in form of the ILU decomposition
 * self - matrix
 * b - right-part vector
 * x0 - first approximation of the solution
 * max_iter - pointer to maximum number of iterations, MAX_ITER if zero;
 * will contain a number of iterations passed
 * tolerance - pointer to desired tolerance value, TOLERANCE if zero;
 * will contain norm of the residual at the end of iteration
 * x - output vector
 */
void sp_matrix_solve_pcg_ilu(sp_matrix_ptr self,
                             sp_matrix_skyline_ilu_ptr ilu,
                             real* b,
                             real* x0,
                             int* max_iter,
                             real* tolerance,
                             real* x);


/*
 * Create ILU decomposition of the sparse matrix in skyline format
 * lu_diag - ILU decomposition diagonal
 * lu_lowertr - lower triangle of the ILU decomposition
 * lu_uppertr - upper triangle of the ILU decomposition
 */
void init_copy_sp_matrix_skyline_ilu(sp_matrix_skyline_ilu_ptr self,
                                     sp_matrix_skyline_ptr parent);

/* Free the sparse matrix skyline & ilu decomposition structure */
void free_sp_matrix_skyline_ilu(sp_matrix_skyline_ilu_ptr self);

/*
 * by given L,U - ILU decomposition of the matrix A
 * calculates L*x = y
 */
void sp_matrix_skyline_ilu_lower_mv(sp_matrix_skyline_ilu_ptr self,
                                    real* x,
                                    real* y);
/*
 * by given L,U - ILU decomposition of the matrix A
 * calculates U*x = y
 */
void sp_matrix_skyline_ilu_upper_mv(sp_matrix_skyline_ilu_ptr self,
                                    real* x,
                                    real* y);

/*
 * by given L,U - ILU decomposition of the matrix A
 * Solves SLAE L*x = b
 * Warning! Side-Effect: modifies b
 */
void sp_matrix_skyline_ilu_lower_solve(sp_matrix_skyline_ilu_ptr self,
                                       real* b,
                                       real* x);

/*
 * by given L,U - ILU decomposition of the matrix A
 * Solves SLAE U*x = b
 * Warning! Side-Effect: modifies b 
 */
void sp_matrix_skyline_ilu_upper_solve(sp_matrix_skyline_ilu_ptr self,
                                       real* b,
                                       real* x);

/* Print contens of the matrix in index form to the stdout */
void sp_matrix_printf(sp_matrix_ptr self);
#ifdef DUMP_DATA
void sp_matrix_dump(sp_matrix_ptr self, char* filename);
void sp_matrix_skyline_dump(sp_matrix_skyline_ptr self, char* filename);
#endif 

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
void solver_export_tetrahedra10_gmsh(fea_solver_ptr solver, char *filename);


/*************************************************************/
/* General functions                                         */

/*
 * A function which will be called in case of error to
 * clear all memory occupied by internal structures.
 * Will clear global variable global_solver in a proper way
 * by calling free_fea_solver.
 * Also shall clear other allocated resources
 */
void application_done(void);


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
 * Vector norm
 */
real vector_norm(real* vector, int size);
/* Scalar multiplication of 2 vectors of the same size */
real cdot(real* vector1, real* vector2, int size);


/*
 * Test function to prove what all matrix manipulations are correct
 * returns FALSE if fail
 */
BOOL do_tests();
/* test dense matrix operations */
BOOL test_matrix();
/* test sparse matrix operations */
BOOL test_sp_matrix();
/* test set of triangle solvers */
BOOL test_triangle_solver();
/* test matrix/vector Conjugate Gradient SLAE solver */
BOOL test_cg_solver();
/*
 * test matrix/vector Preconditioned Conjugate Gradient SLAE solver
 * where preconditioner M = ILU decomposition
 */
BOOL test_pcg_ilu_solver();
/* test ILU decomposition */
BOOL test_ilu();
/* test Cholesky decomposition */
BOOL test_cholesky();

/*
 * Solver function which shall be called
 * when all data read to an appropriate structures
 */
void solve(fea_task_ptr ask,
           fea_solution_params_ptr fea_params,
           nodes_array_ptr nodes,
           elements_array_ptr elements,
           presc_boundary_array_ptr presc_boundary);

#ifdef DUMP_DATA
/* Dump input data to check if parser works correctly */
void dump_input_data( char* filename,
                      fea_task_ptr task,
                      fea_solution_params_ptr fea_params,
                      nodes_array_ptr nodes,
                      elements_array_ptr elements,
                      presc_boundary_array_ptr presc_boundary);
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
gauss_node_ptr solver_new_gauss_node(fea_solver_ptr self,
                                     int gauss_node_index);

/* Deallocate gauss node */
void solver_free_gauss_node(fea_solver_ptr self,
                            gauss_node_ptr node);

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
