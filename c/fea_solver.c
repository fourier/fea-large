#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef USE_EXPAT
#include <expat.h>
#endif

/*************************************************************/
/* Constants                                                 */
#define READ_BUF_SIZE 10000

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
  int arclength_max;            /* maximum number of arc lenght searches */
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
  int prescribed_nodes_count;
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
/* Functions declarations                                    */


/*
 * Load initial data from file
 */
BOOL initial_data_load(char *filename,
                       fea_task **task,
                       fea_solution_params **fea_params,
                       nodes_array **nodes,
                       elements_array **elements,
                       prescribed_boundary_array **presc_boundary);

 

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


/*************************************************************/
/* Global variables                                          */


int parse_cmdargs(int argc, char **argv,char **filename)
{
  if (argc < 2)
  {
    printf("Usage: fea_solve input_data.xml\n");
    return 1;
  }
  *filename = argv[1];
  return 0;
}

int do_main(char* filename)
{
  int result = 0;
  fea_task *task;
  fea_solution_params *fea_params;
  nodes_array *nodes;
  elements_array *elements;
  prescribed_boundary_array *presc_boundary;
  
  if(!initial_data_load(filename,&task,&fea_params,&nodes,&elements,&presc_boundary))
  {
    printf("Error. Unable to load %s.\n",filename);
    result = 1;
  }
  return result;
}

int main(int argc, char **argv)
{
  char* filename = 0;
  int result = 0;
  
  do
  {
    if ( result = parse_cmdargs(argc, argv,&filename))
      break;
    if ( result = do_main(filename))
      break;
  } while(0);

  return result;
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

/* Initializa fea task structure and fill with default values */
static fea_task* initialize_fea_task()
{
  /* allocate memory */
  fea_task *task = (fea_task *)malloc(sizeof(fea_task));
  /* set default values */
  task->desired_tolerance = 1e-8;
  task->dof = 3;
  task->ele_type = TETRAHEDRA10;
  task->linesearch_max = 0;
  task->arclength_max = 0;
  task->load_increments_count = 0;
  task->type = CARTESIAN3D;
  task->modified_newton = TRUE;
  return task;
}

static void free_fea_task(fea_task* task)
{
  /* TODO: implement this */

}

/* Initializes fea solution params with default values */
static fea_solution_params* initialize_fea_solution_params()
{
  /* allocate memory */
  fea_solution_params *fea_params = (fea_solution_params *)
    malloc(sizeof(fea_solution_params));
  /* set default values */
  fea_params->gauss_nodes_count = 5;
  fea_params->nodes_per_element = 10;
  fea_params->msize = 0;
  return fea_params;
}


/* Initialize nodes array but not initialize particular arrays  */
static nodes_array* initialize_nodes_array()
{
  /* allocate memory */
  nodes_array *nodes = (nodes_array*)malloc(sizeof(nodes_array));
  /* set zero values */
  nodes->nodes = (void*)0;
  nodes->nodes_count = 0;
  return nodes;
}

static void free_nodes_array(nodes_array* array)
{
  /* TODO: implement this */
}


/* Initialize elements array but not initialize particular elements */
static elements_array* initialize_elements_array()
{
  /* allocate memory */
  elements_array *elements = (elements_array*)malloc(sizeof(elements_array));
  /* set zero values */
  elements->elements = (void*)0;
  elements->elements_count = 0;
  return elements;
}

static free_elements_array(elements_array *elements)
{
  /* TODO: implement this */
}

/* Initialize boundary nodes array but not initialize particular nodes */
static prescribed_boundary_array* initialize_prescribed_boundary_array()
{
  /* allocate memory */
  prescribed_boundary_array *presc_boundary = (prescribed_boundary_array*)
    malloc(sizeof(prescribed_boundary_array));
  /* set zero values */
  presc_boundary->prescribed_nodes = (void*)0;
  presc_boundary->prescribed_nodes_count = 0;
  return presc_boundary;
}

static free_prescribed_boundary_array(prescribed_boundary_array* presc)
{
  /* TODO: implement this */
}


#ifdef USE_EXPAT

typedef enum xml_format_tags_enum {
  UNKNOWN_TAG,
  TASK,
  MODEL,
  MODEL_TYPE,
  MODEL_PARAMETERS,
  LAMBDA,
  MU,
  SOLUTION,
  MODIFIED_NEWTON,
  TASK_TYPE,
  ELEMENT_TYPE,
  GAUSS_NODES_COUNT,
  LOAD_INCREMENTS_COUNT,
  DESIRED_TOLERANCE,
  LINE_SEARCH,
  LINE_SEARCH_MAX,
  ARC_LENGTH,
  ARC_LENGTH_MAX,
  INPUT_DATA,
  GEOMETRY,
  NODES,
  NODE,
  X,
  Y,
  Z,
  ELEMENTS,
  ELEMENT,
  NODE_ID,
  BOUNDARY_CONDITIONS,
  PRESCRIBED_DISPLACEMENTS,
  PRESC_NODE,
  PRESC_ID,
  PRESC_X,
  PRESC_Y,
  PRESC_Z,
  PRESC_TYPE,
} xml_format_tags;

enum xml_format_tags tagname_to_enum(const XML_Char* name)
{
  if (!strcmp(name,"task")) return TASK;
  if (!strcmp(name,"model")) return MODEL;
  if (!strcmp(name,"model_type")) return MODEL_TYPE;
  if (!strcmp(name,"model_parameters")) return MODEL_PARAMETERS;
  if (!strcmp(name,"lambda")) return LAMBDA;
  if (!strcmp(name,"mu")) return MU;
  if (!strcmp(name,"solution")) return SOLUTION;
  if (!strcmp(name,"modified_newton")) return MODIFIED_NEWTON;
  if (!strcmp(name,"task_type")) return TASK_TYPE;
  if (!strcmp(name,"element_type")) return ELEMENT_TYPE;
  if (!strcmp(name,"gauss_nodes_count")) return GAUSS_NODES_COUNT;
  if (!strcmp(name,"load_increments_count")) return LOAD_INCREMENTS_COUNT;
  if (!strcmp(name,"desired_tolerance")) return DESIRED_TOLERANCE;
  if (!strcmp(name,"line_search")) return LINE_SEARCH;
  if (!strcmp(name,"line_search_max")) return LINE_SEARCH_MAX;
  if (!strcmp(name,"arc_length")) return ARC_LENGTH;
  if (!strcmp(name,"arc_length_max")) return ARC_LENGTH_MAX;
  if (!strcmp(name,"input_data")) return INPUT_DATA;
  if (!strcmp(name,"geometry")) return GEOMETRY;
  if (!strcmp(name,"nodes")) return NODES;
  if (!strcmp(name,"node")) return NODE;
  if (!strcmp(name,"x")) return X;
  if (!strcmp(name,"y")) return Y;
  if (!strcmp(name,"z")) return Z;
  if (!strcmp(name,"elements")) return ELEMENTS;
  if (!strcmp(name,"element")) return ELEMENT;
  if (!strcmp(name,"node_id")) return NODE_ID;
  if (!strcmp(name,"boundary_conditions")) return BOUNDARY_CONDITIONS;
  if (!strcmp(name,"prescribed_displacements")) return PRESCRIBED_DISPLACEMENTS;
  if (!strcmp(name,"presc_node")) return PRESC_NODE;
  if (!strcmp(name,"presc_id")) return PRESC_ID;
  if (!strcmp(name,"presc_x")) return PRESC_X;
  if (!strcmp(name,"presc_y")) return PRESC_Y;
  if (!strcmp(name,"presc_z")) return PRESC_Z;
  if (!strcmp(name,"presc_type")) return PRESC_TYPE;
  return UNKNOWN_TAG;
}

/* An input data structure used in parser */
typedef struct parse_data_tag {
  fea_task *task;
  fea_solution_params *fea_params;
  nodes_array *nodes;
  elements_array *elements;
  prescribed_boundary_array *presc_boundary;
  xml_format_tags last_tag;
  char* current_text;
  int current_size;
} parse_data;


void expat_start_tag_handler(void *userData,
                      const XML_Char *name,
                      const XML_Char **atts)
{
  parse_data* data = (parse_data*)userData;
  data->last_tag = tagname_to_enum(name);
  printf("open: %s\n",name);
}

void expat_end_tag_handler(void *userData,
                   const XML_Char *name)
{
  parse_data* data = (parse_data*)userData;
  printf("total text len = %d\n",data->current_size);
  data->current_text = (char*)0;
  data->current_size = 0;
  printf("close %s\n",name);
}

void expat_text_handler(void *userData,
                        const XML_Char *s,
                        int len)
{
  parse_data* data = (parse_data*)userData;
  char *ptr;
  if (!data->current_text)
  {
    data->current_text = (char*)malloc(len);
    ptr = data->current_text;
  }
  else
  {
    /* FIXME: after every reallocation data->current_text is another pointer!!!*/
    ptr = data->current_text;
    data->current_text = (char*)realloc(data->current_text,data->current_size+len);
    ptr += data->current_size;  
  }
  memcpy(ptr,s,len);
  data->current_size += len;
  printf("fragment text len = %d\n",len);
}


static BOOL expat_data_load(char *filename,
                            fea_task **task,
                            fea_solution_params **fea_params,
                            nodes_array **nodes,
                            elements_array **elements,
                            prescribed_boundary_array **presc_boundary)
{
  BOOL result = FALSE;
  FILE* xml_document_file;
  XML_Parser parser;
  size_t file_size = 0;
  size_t read_bytes = 0;
  parse_data parse;
  enum XML_Status status;
  char *file_contents = (char*)0;
  char *file_current_position = file_contents;
  
  /* Try to open file */
  if (!(xml_document_file = fopen(filename,"rt")))
  {
    printf("Error, could not open file %s\n",filename);
    return FALSE;
  }
  /* Determine file size */
  if (fseek(xml_document_file,0,SEEK_END))
  {
    printf("Error reading file %s\n",filename);
    return FALSE;
  }
  file_size = ftell(xml_document_file);
  /* rewind to the begin of file */
  fseek(xml_document_file,0,SEEK_SET);

  /* Create parser */
  parser = XML_ParserCreate(NULL);
  /* Set handlers */
  XML_SetElementHandler(parser, &expat_start_tag_handler, &expat_end_tag_handler);
  XML_SetCharacterDataHandler(parser,expat_text_handler);

  /* initialize data */
  parse.current_text = (char*)0;
  parse.current_size = 0;
  /* TODO: add initialization here */
  XML_SetUserData(parser,&parse);

  /* read whole file */
  file_contents = (char*)malloc(file_size);
  read_bytes = fread(file_contents,1,file_size,xml_document_file);
  if (errno)
  {
    free(file_contents);
    /* deallocate parse data */
    return FALSE;
  }
  /* call parser */
  status = XML_Parse(parser,file_contents,(int)read_bytes,1);
  free(file_contents);
  
  return result;
}



#endif



BOOL initial_data_load(char *filename,
                       fea_task **task,
                       fea_solution_params **fea_params,
                       nodes_array **nodes,
                       elements_array **elements,
                       prescribed_boundary_array **presc_boundary)
{
  BOOL result = FALSE;
#ifdef USE_EXPAT
  result = expat_data_load(filename,task,fea_params,nodes,elements,presc_boundary);
#endif
  
  return result;
}

