/* -*- Mode: C; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include "fea_solver.h"
#include "sp_matrix.h"
#include "dense_matrix.h"
#include "tests.h"
#include "xml_loader.h"


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


void error(char* msg)
{
  fprintf(stderr,"feasolve error encountered: %s",msg);
  exit(1);
}


int main(int argc, char **argv)
{
  char* filename = 0;
  int result = 0;

  /* Perform tests before start */
  if (!do_tests())
  {
    printf("Error! Tests failed!\n");
    return 1;
  }
  
  do
  {
    if ( TRUE == (result = parse_cmdargs(argc, argv,&filename)))
      break;
    if ( TRUE == (result = do_main(filename)))
      break;
  } while(0);

  return result;
}

int do_main(char* filename)
{
  /* initialize variables */
  int result = 0;
  fea_task_ptr task = (fea_task_ptr)0;
  fea_solution_params_ptr fea_params = (fea_solution_params_ptr)0;
  nodes_array_ptr nodes = (nodes_array_ptr)0;
  elements_array_ptr elements = (elements_array_ptr)0;
  presc_boundary_array_ptr presc_boundary = (presc_boundary_array_ptr)0;
  
  /* load geometry and solution details */
  if(!initial_data_load(filename,
                        &task,
                        &fea_params,
                        &nodes,
                        &elements,
                        &presc_boundary))
  {
    printf("Error. Unable to load %s.\n",filename);
    result = 1;
  }
  
  /* solve task */
  solve(task, fea_params, nodes, elements, presc_boundary);
  
  return result;
}

void solve( fea_task_ptr task,
            fea_solution_params_ptr fea_params,
            nodes_array_ptr nodes,
            elements_array_ptr elements,
            presc_boundary_array_ptr presc_boundary)
{
  /* initialize variables */
  fea_solver_ptr solver = (fea_solver_ptr)0;
  int it = 0;
  real tolerance;
  sp_matrix stiffness;
#ifdef DUMP_DATA
  /* Dump all data in debug version */
  dump_input_data("input.txt",task,fea_params,nodes,elements,presc_boundary);
#endif
  /* Prepare solver instance */
  solver = new_fea_solver(task,
                          fea_params,
                          nodes,
                          elements,
                          presc_boundary);
#ifdef DUMP_DATA
  /* solver_update_nodes_with_bc(solver, 1); */
  /* dump_input_data("input1.txt",task,fea_params,solver->nodes_p,elements,
     presc_boundary); */
#endif
  
  /* Create elements database */
  solver_create_element_database(solver);
  /* create an array of shape functions gradients in initial configuration */
  solver_create_initial_shape_gradients(solver);

  /* Increment loop starts here */
  for (; solver->current_load_step < solver->task_p->load_increments_count;
       ++ solver->current_load_step)
  {
    it = 0;
    /* apply prescribed displacements */
    solver_update_nodes_with_bc(solver, 1);

    /* Create an array of shape functions gradients in current configuration */
    solver_create_current_shape_gradients(solver);
    /* create stresses in order to use them in residual forces and in
     * initial stress component of the stiffness matrix */
    solver_create_stresses(solver);

    /* create global stiffness matrix K */
    solver_create_stiffness(solver);
    /* store global stiffness matrix for modified Newton method */
    copy_sp_matrix(&solver->global_mtx,&stiffness);
    do 
    {
      it ++;

      /* create right-side vector of residual forces (-R) */
      solver_create_residual_forces(solver);

      /* create global stiffness matrix K */
      if (solver->task_p->modified_newton) 
      {
        /*
         * use stored  stiffness matrix in
         * modified Newton method
         */
        free_sp_matrix(&solver->global_mtx);
        copy_sp_matrix(&stiffness,&solver->global_mtx);
      }
      else                        
      {
        /* create global stiffness otherwise */
        solver_create_stiffness(solver);
      }
      /* apply prescribed boundary conditions */
      solver_apply_prescribed_bc(solver,0);

      /* solve global equation system K*u=-R */
      sp_matrix_solve(&solver->global_mtx,
                      solver->global_forces_vct,
                      solver->global_solution_vct);
      /* check for convergence */

      tolerance = cdot(solver->global_forces_vct,
                       solver->global_solution_vct,
                       solver->global_mtx.rows_count);
    
      printf("Tolerance <X,R> = %e\n",tolerance);
      printf("Newton iteration %d finished\n",it);
    
      /* update nodes array with solution */
      solver_update_nodes_with_solution(solver,solver->global_solution_vct);
      solver_create_current_shape_gradients(solver);
      solver_create_stresses(solver);

    } while ( fabs(tolerance) > solver->task_p->desired_tolerance);
    /* store current load step */
    init_load_step(solver,
                   &solver->load_steps_p[solver->current_load_step],
                   solver->current_load_step);

    /* clear stored stiffness matrix */
    free_sp_matrix(&stiffness);
    printf("Load increment %d finished\n",solver->current_load_step);
  }
  /* export solution */
  printf("Exporting data...\n");
  solver->export_function(solver,"deformed.msh");
  
  free_fea_solver(solver);
}


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

#ifdef DUMP_DATA
void dump_input_data( char* filename,
                      fea_task_ptr task,
                      fea_solution_params_ptr fea_params,
                      nodes_array_ptr nodes,
                      elements_array_ptr elements,
                      presc_boundary_array_ptr presc_boundary)
{
  int i,j;
  FILE *f;
  if ((f = fopen(filename,"w+")))
  {
    fprintf(f,"nodes\n");
    for ( i = 0; i < nodes->nodes_count; ++ i)
    {
      for ( j = 0; j < MAX_DOF; ++ j)
        fprintf(f,"%f ", nodes->nodes[i][j]);
      fprintf(f,"\n");
    }
    fprintf(f,"elements\n");
    for ( i = 0; i < elements->elements_count; ++ i)
    {
      for ( j = 0; j < fea_params->nodes_per_element; ++ j)
        fprintf(f,"%d ", elements->elements[i][j]);
      fprintf(f,"\n");
    }
    fprintf(f,"boundary\n");
    for ( i = 0; i < presc_boundary->prescribed_nodes_count; ++ i)
    {
      fprintf(f,"%d %f %f %f %d\n",
              presc_boundary->prescribed_nodes[i].node_number,
              presc_boundary->prescribed_nodes[i].values[0],
              presc_boundary->prescribed_nodes[i].values[1],
              presc_boundary->prescribed_nodes[i].values[2],
              presc_boundary->prescribed_nodes[i].type);
    }
    fclose(f);
  }
}
#endif


real solver_node_dof(fea_solver_ptr self,
                     nodes_array_ptr nodes,
                     int element,
                     int node,
                     int dof)
{
  return nodes->nodes[self->elements_p->elements[element][node]][dof];
}


fea_solver* new_fea_solver(fea_task_ptr task,
                           fea_solution_params_ptr fea_params,
                           nodes_array_ptr nodes,
                           elements_array_ptr elements,
                           presc_boundary_array_ptr prs_boundary)
{
  int msize,bandwidth,elnum,gauss_count,i,j,k,l;
  /* Allocate structure */
  fea_solver_ptr solver = (fea_solver_ptr)malloc(sizeof(fea_solver));
  /* Copy pointers to the solver structure */
  solver->task_p = task;
  solver->fea_params_p = fea_params;
  solver->nodes0_p = nodes;
  solver->nodes_p = new_copy_nodes_array(nodes);
  solver->elements_p = elements;
  solver->presc_boundary_p = prs_boundary;

  solver->elements_db.gauss_nodes = (gauss_node**)0;
  solver_create_element_params(solver);
  solver_create_model_params(solver);
  
  /* initialize an array of gradients of shape functions per
   * element/gauss node and arrays of deformation gradients/stresses */
  elnum = elements->elements_count;
  gauss_count = solver->fea_params_p->gauss_nodes_count;
  solver->shape_gradients0 =
    (shape_gradients***)malloc(sizeof(shape_gradients_ptr*)*elnum);
  solver->shape_gradients  =
    (shape_gradients***)malloc(sizeof(shape_gradients_ptr*)*elnum);
  solver->stresses = (tensor**)malloc(sizeof(tensor*)*elnum);
  solver->graddefs = (tensor**)malloc(sizeof(tensor*)*elnum);
  for (i = 0; i < elnum; ++ i)
  {
    solver->stresses[i] = (tensor*)malloc(sizeof(tensor)*gauss_count);
    solver->graddefs[i] = (tensor*)malloc(sizeof(tensor)*gauss_count);
    solver->shape_gradients0[i] =
      (shape_gradients**)malloc(sizeof(shape_gradients_ptr)*gauss_count);
    solver->shape_gradients[i] =
      (shape_gradients**)malloc(sizeof(shape_gradients_ptr)*gauss_count);
    
    for (j = 0; j < gauss_count; ++ j)
    {
      solver->shape_gradients0[i][j] = (shape_gradients_ptr)0;
      solver->shape_gradients[i][j] = (shape_gradients_ptr)0;
      for ( k = 0; k < MAX_DOF; ++ k)
        for (l = 0; l < MAX_DOF; ++ l)
        {
          solver->stresses[i][j].components[k][l] = 0.0;
          solver->graddefs[i][j].components[k][l] = 0.0;
        }
    }
  }
  solver->current_load_step = 0;
  solver->load_steps_p = (load_step_ptr)malloc(sizeof(load_step)*
                                               task->load_increments_count);
  /* allocate resources initialize global stiffness matrix */
  /* global matrix size */
  msize = nodes->nodes_count*solver->task_p->dof;
  /* approximate bandwidth of a global matrix
   * usually sqrt(msize)*2*/
  bandwidth = (int)sqrt(msize)*2;
  init_sp_matrix(&solver->global_mtx,msize,msize,bandwidth,CRS);
  /* allocate memory for global forces and solution vectors */
  solver->global_forces_vct = (real*)malloc(sizeof(real)*msize);
  solver->global_solution_vct = (real*)malloc(sizeof(real)*msize);
  memset(solver->global_forces_vct,0,sizeof(real)*msize);
  memset(solver->global_solution_vct,0,sizeof(real)*msize);
  return solver;
}


void free_fea_solver(fea_solver_ptr solver)
{
  /* deallocate resources */
  /* free shape gradients, graddefs and stresses */
  int elnum = solver->elements_p->elements_count;
  int gauss_count = solver->fea_params_p->gauss_nodes_count;
  int i,j;
  for (i = 0; i < elnum; ++ i)
  {
    for (j = 0; j < gauss_count; ++ j)
    {
      if (solver->shape_gradients0[i][j])
        solver_free_shape_gradients(solver,solver->shape_gradients0[i][j]);
      if (solver->shape_gradients[i][j])
        solver_free_shape_gradients(solver,solver->shape_gradients[i][j]);
    }
    free(solver->shape_gradients0[i]);
    free(solver->shape_gradients[i]);
    free(solver->stresses[i]);
    free(solver->graddefs[i]);
  }
  free(solver->shape_gradients0);
  free(solver->shape_gradients);  
  free(solver->stresses);
  free(solver->graddefs);
  /* free stored load steps */
  for ( i = 0; i < solver->current_load_step; ++ i)
    free_load_step(solver, &solver->load_steps_p[i]);
  free(solver->load_steps_p);
  /* deallocate all other resources */
  solver_free_element_database(solver);
  free_fea_task(solver->task_p);
  free_fea_solution_params(solver->fea_params_p);
  free_nodes_array(solver->nodes0_p);
  free_nodes_array(solver->nodes_p);
  free_elements_array(solver->elements_p);
  free_presc_boundary_array(solver->presc_boundary_p);
  free_sp_matrix(&solver->global_mtx);
  free(solver->global_forces_vct);
  free(solver->global_solution_vct);
  free(solver);
}

gauss_node_ptr solver_new_gauss_node(fea_solver_ptr self,
                                     int gauss_node_index)
{
  gauss_node_ptr node = (gauss_node_ptr)0;
  int i,j;
  real r,s,t;
  /* Check for array bounds*/
  if (gauss_node_index >= 0 &&
      gauss_node_index < self->fea_params_p->gauss_nodes_count)
  {
    node = (gauss_node_ptr)malloc(sizeof(gauss_node));
    /* set the weight for this gauss node */
    node->weight = self->elements_db.gauss_nodes_data[gauss_node_index][0];
    /* set shape function values and their derivatives for this node */
    node->forms =
      (real*)malloc(sizeof(real)*(self->fea_params_p->nodes_per_element));
    node->dforms = (real**)malloc(sizeof(real*)*(self->task_p->dof));
    for ( i = 0; i < self->task_p->dof; ++ i)
      node->dforms[i] =
        (real*)malloc(sizeof(real)*(self->fea_params_p->nodes_per_element));
    for ( i = 0; i < self->fea_params_p->nodes_per_element; ++ i)
    {
      r = self->elements_db.gauss_nodes_data[gauss_node_index][1];
      s = self->elements_db.gauss_nodes_data[gauss_node_index][2];
      t = self->elements_db.gauss_nodes_data[gauss_node_index][3];
      node->forms[i] = self->shape(i,r,s,t);
      for ( j = 0; j < self->task_p->dof; ++ j)
        node->dforms[j][i] = self->dshape(i,j,r,s,t);
    }

  }
  return node;
}

/* Deallocate gauss node */
void solver_free_gauss_node(fea_solver_ptr self,
                            gauss_node_ptr node)
{
  int i;
  if (node)
  {
    /* clear forms and dforms arrays */
    free(node->forms);
    for ( i = 0; i < self->task_p->dof; ++ i)
      free(node->dforms[i]);
    free(node->dforms);
    /* free the node itself */
    free(node);
  }
}
  

void solver_create_element_database(fea_solver_ptr self)
{
  int gauss;
  int gauss_count = self->fea_params_p->gauss_nodes_count;
  /* Create database only if not created yet */
  if (!self->elements_db.gauss_nodes)
  {
    /* allocate memory for gauss nodes array */
    self->elements_db.gauss_nodes =
      (gauss_node_ptr*)malloc(sizeof(gauss_node*)*gauss_count);
    for (gauss = 0; gauss < gauss_count; ++ gauss)
      self->elements_db.gauss_nodes[gauss] =
        solver_new_gauss_node(self,gauss);
    
  }
}

void solver_free_element_database(fea_solver_ptr solver)
{
  int gauss;
  if (solver->elements_db.gauss_nodes)
  {
    for (gauss = 0; gauss < solver->fea_params_p->gauss_nodes_count; ++gauss)
      solver_free_gauss_node(solver,solver->elements_db.gauss_nodes[gauss]);
    
    free(solver->elements_db.gauss_nodes);
  }
}

void solver_create_element_params_tetrahedra10(fea_solver_ptr solver);

/*
 * Creates particular element-dependent data in fea_solver
 * All new element types shall be added here 
 */
void solver_create_element_params(fea_solver_ptr solver)
{
  switch (solver->task_p->ele_type)
  {
  case TETRAHEDRA10:
    solver_create_element_params_tetrahedra10(solver);
    break;
  default:
    /* TODO: add error handling here */
    printf("Error: unknown element type");
    exit(1);
  };
  
}

void init_load_step(fea_solver_ptr self,
                    load_step_ptr step,
                    int step_number)
{
  int elnum = self->elements_p->elements_count;
  int gauss_count = self->fea_params_p->gauss_nodes_count;
  int i,j,k,l;
  if (step)
  {
    step->step_number = step_number;
    step->nodes_p = new_copy_nodes_array(self->nodes_p);
    step->stresses = (tensor**)malloc(sizeof(tensor*)*elnum);
    step->graddefs = (tensor**)malloc(sizeof(tensor*)*elnum);

    for (i = 0; i < elnum; ++ i)
    {
      step->stresses[i] = (tensor*)malloc(sizeof(tensor)*gauss_count);
      step->graddefs[i] = (tensor*)malloc(sizeof(tensor)*gauss_count);    
      for (j = 0; j < gauss_count; ++ j)
      {
        for ( k = 0; k < MAX_DOF; ++ k)
          for (l = 0; l < MAX_DOF; ++ l)
          {
            step->stresses[i][j].components[k][l] =
              self->stresses[i][j].components[k][l];
            step->graddefs[i][j].components[k][l] =
              self->graddefs[i][j].components[k][l];
          }
      }
    }
  }
}

void free_load_step(fea_solver_ptr self, load_step_ptr step)
{
  int elnum = self->elements_p->elements_count;
  int i;
  if (step)
  {
    for (i = 0; i < elnum; ++ i)
    {
      free(step->stresses[i]);
      free(step->graddefs[i]);
    }
    free(step->stresses);
    free(step->graddefs);
    free_nodes_array(step->nodes_p);
  }
}


shape_gradients_ptr solver_new_shape_gradients(fea_solver_ptr self,
                                               nodes_array_ptr nodes,
                                               int element,
                                               int gauss)
{
  int i,j,k;
  int row_size;
  real detJ;
  /* J is a Jacobi matrix of transformation btw local and global */
  /* coordinate systems */
  real J[MAX_DOF][MAX_DOF];
  shape_gradients_ptr grads = (shape_gradients_ptr)0;

  /* Fill an array using Bonet & Wood 7.6(a,b) p.198, 1st edition */
  /* also see Zienkiewitz v1, 6th edition, p.146-147 */

  for (i = 0; i < MAX_DOF; ++ i)
    memset(&J[i],0,sizeof(real)*MAX_DOF);
  /* First, fill the Jacobi matrix (3x3) */
  /* I = 1..n, n - number of nodes per element
   * x_I, y_I, z_I - nodal coordinates for the element
   *           dN_1(r,s,t)            dN_n(r,s,t)     
   * J(1,1) =  ---------- * x_1 + ... ---------- * x_n
   *               dr                      dr
   *
   *           dN_1(r,s,t)            dN_n(r,s,t)     
   * J(1,2) =  ---------- * y_1 + ... ---------- * y_n
   *               dr                      dr
   *
   *           dN_1(r,s,t)            dN_n(r,s,t)     
   * J(2,1) =  ---------- * x_1 + ... ---------- * x_n
   *               ds                      ds
   * ...
   */
  for (i = 0; i < MAX_DOF; ++ i)
    for (j = 0; j < MAX_DOF; ++ j)
    {
      for (k = 0; k < self->fea_params_p->nodes_per_element; ++ k)
        J[i][j] += self->elements_db.gauss_nodes[gauss]->dforms[i][k]* \
          solver_node_dof(self,nodes,element,k,j);
    }
  if (inv3x3(J,&detJ))                /* inverse exists */
  {
    /* Allocate memory for shape gradients */
    grads = (shape_gradients_ptr)malloc(sizeof(shape_gradients));
    grads->grads = (real**)malloc(sizeof(real*)*(self->task_p->dof));
    row_size = sizeof(real)*(self->fea_params_p->nodes_per_element);
    for (i = 0; i < self->task_p->dof; ++ i)
    {
      grads->grads[i] = (real*)malloc(row_size);
      memset(grads->grads[i],0,row_size);
    }
    /* Store determinant of the Jacobi matrix */
    grads->detJ = detJ;
    
    /* [ dN/dx ]           [ dN/dr ] */
    /* [ dN/dy ]  = J^-1 * [ dN/ds ] */
    /* [ dN/dz ]           [ dN/dt ] */
    for ( i = 0; i < MAX_DOF; ++ i)
      for ( j = 0; j < self->fea_params_p->nodes_per_element; ++ j)
        for ( k = 0; k < MAX_DOF; ++ k)
          grads->grads[i][j] += J[i][k]* \
            self->elements_db.gauss_nodes[gauss]->dforms[k][j];
  }

  return grads;
}

/* Destructor for the shape gradients array */
void solver_free_shape_gradients(fea_solver_ptr self,
                                 shape_gradients_ptr grads)
{
  int i;
  /* for (i = 0; i < self->fea_params->nodes_per_element; ++ i */
  for (i = 0; i < self->task_p->dof; ++ i)
    free(grads->grads[i]);
  free(grads->grads);
  grads->grads = (real**)0;
  free(grads);
}

#ifdef DUMP_DATA
void solver_dump_local_stiffness(fea_solver* self,real **stiff,int el)
{
  int i,j;
  FILE* f;
  char fname[50];
  int size = self->fea_params_p->nodes_per_element*self->task_p->dof;
  sprintf(fname,"elements/K%d.txt",el);
  if ((f = fopen(fname,"w+")))
  {
    for ( i = 0; i < size; ++ i)
    {
      for ( j = 0; j < size; ++ j)
        fprintf(f,"%e ",stiff[i][j]);
      fprintf(f,"\n");
    }
    fclose(f);
  }
}
#endif

#ifdef DUMP_DATA
void matrix_tensor_mapping(int I, int* i, int* j)
{
  switch (I)
  {
  case 0:
    *i = 0; *j = 0;
    break;
  case 1:
    *i = 1; *j = 1;
    break;
  case 2:
    *i = 2; *j = 2;
    break;
  case 3:
    *i = 0; *j = 1;
    break;
  case 4:
    *i = 1; *j = 2;
    break;
  case 5:
    *i = 0; *j = 2;
  }
}

#endif


void solver_create_shape_gradients(fea_solver_ptr self,BOOL current)
{
  /* prepare an array of shape functions gradients in
   * gauss nodes per element */
  shape_gradients_ptr grads = (shape_gradients_ptr)0;
  int gauss,element;
  /* loop by elements */
  for ( element = 0;
        element < self->elements_p->elements_count;
        ++ element)
  {
    /* loop by gauss nodes per element */
    for (gauss = 0;
         gauss < self->fea_params_p->gauss_nodes_count;
         ++ gauss)
    {
      /* create shape gradients either in initial or current configuration */
      grads = solver_new_shape_gradients(self,
                                         current ?
                                         self->nodes_p : self->nodes0_p,
                                         element,gauss);
      if (grads)
      {
        /* free previous shape gradients array */
        if (current)            /* current configuration */
        {
          if ( self->shape_gradients[element][gauss] )
            solver_free_shape_gradients(self,
                                        self->shape_gradients[element][gauss]);
          self->shape_gradients[element][gauss] = grads;
        }
        else                    /* initial configuration */
        {
          if ( self->shape_gradients0[element][gauss] )
            solver_free_shape_gradients(self,
                                        self->shape_gradients0[element][gauss]);
          self->shape_gradients0[element][gauss] = grads;
        }
      }
    }
  }
}


void solver_create_current_shape_gradients(fea_solver_ptr self)
{
  solver_create_shape_gradients(self,TRUE);
}

void solver_create_initial_shape_gradients(fea_solver_ptr self)
{
  solver_create_shape_gradients(self,FALSE);
}



void solver_create_stresses(fea_solver_ptr self)
{
  int gauss,el;
  /* loop by elements */
  for ( el = 0;
        el < self->elements_p->elements_count;
        ++ el)
  {
    /* loop by gauss nodes per element */
    for (gauss = 0;
         gauss < self->fea_params_p->gauss_nodes_count;
         ++ gauss)
    {
      solver_element_gauss_stress(self, el, gauss,
                                  self->graddefs[el][gauss].components,
                                  self->stresses[el][gauss].components);
    }
  }
}

void solver_create_residual_forces(fea_solver_ptr self)
{
  int el = 0;
  memset(self->global_forces_vct,0,sizeof(real)*self->global_mtx.rows_count);

  for (; el < self->elements_p->elements_count; ++ el)
    solver_local_residual_forces(self, el);
}

/* Create global stiffness matrix */
void solver_create_stiffness(fea_solver_ptr self)
{
  int el = 0;
  /* clear global stiffness matrix before constructing a new one */
  clear_sp_matrix(&self->global_mtx);
  
  for (; el < self->elements_p->elements_count; ++ el)
  {
    solver_local_constitutive_part(self,el);
    solver_local_initial_stess_part(self,el);
  }
}



void solver_local_constitutive_part(fea_solver_ptr self,int element)
{
  /* matrix of gradients of shape functions */
  shape_gradients_ptr grads = (shape_gradients_ptr)0;
  int gauss,a,b,i,j,k,l,I,J,globalI,globalJ;
  real sum;
  /* size of a local stiffness matrix */
  int size;
  /* number of nodes per element */
  int nelem;
  /* current number of d.o.f */
  int dof;
  /* local stiffness matrix */
  real **stiff = (real**)0;
  /* C tensor depending on material model */
  real ctens[MAX_DOF][MAX_DOF][MAX_DOF][MAX_DOF];
  
  real cikjl = 0;
  /* allocate memory for a local stiffness matrix */
  size = self->fea_params_p->nodes_per_element*self->task_p->dof;
  stiff = (real**)malloc(sizeof(real*)*size);
  for (i = 0; i < size; ++ i)
  {
    stiff[i] = (real*)malloc(sizeof(real)*size);
    memset(stiff[i],0,sizeof(real)*size);
  }
    
  dof = self->task_p->dof;
  nelem = self->fea_params_p->nodes_per_element;
  
  /* loop by gauss nodes - numerical integration */
  for (gauss = 0; gauss < self->fea_params_p->gauss_nodes_count ; ++ gauss)
  {
    /* obtain a C tensor */
    self->ctensor(self,self->graddefs[element][gauss].components,ctens);

    grads = self->shape_gradients[element][gauss];
    if (grads)
    {
      /* Construct components of stiffness matrix in
       * indical form using Bonet & Wood 7.35 p.207, 1st edition */
      
      /* loop by nodes */
      for ( a = 0; a < nelem; ++ a)
        for (b = 0; b < nelem; ++ b)
        {
          /* loop by d.o.f in a stiffness matrix block [K_{ab}]ij, 3x3 */
          for (i = 0; i < dof; ++ i)
            for (j = 0; j < dof; ++ j)
            {
              /* indicies in a local stiffness matrix */
              I = a*dof + i;
              J = b*dof + j;
              sum = 0.0;
              /* sum of particular derivatives and components of C tensor */
              for (k = 0; k < dof; ++ k)
                for (l = 0; l < dof; ++ l)
                {
                  /* cikjl = ctens[i][k][j][l]; */
                  cikjl = (ctens[i][k][j][l]+ctens[i][k][l][j]+
                           ctens[k][i][j][l]+ctens[k][i][l][j])/4.;
                  sum += 
                    grads->grads[k][a]*cikjl*grads->grads[l][b];
                }
              /*
               * multiply by volume of an element = det(J)
               * where divider 6 or 2 or others already accounted in
               * weights of gauss nodes
               */
              sum *= fabs(grads->detJ);
              /* ... and weight of the gauss nodes for  */
              sum *= self->elements_db.gauss_nodes[gauss]->weight;
              /* append to the local stiffness */
              stiff[I][J] += sum;
              /* finally distribute to the global matrix */
              globalI = self->elements_p->elements[element][a]*dof + i;
              globalJ = self->elements_p->elements[element][b]*dof + j;
              sp_matrix_element_add(&self->global_mtx,
                                    globalI,
                                    globalJ,
                                    sum);
            }
        }
    }
  }
  
#ifdef DUMP_DATA
  solver_dump_local_stiffness(self,stiff,element);
#endif
  
  /* clear local stiffness */
  for ( i = 0; i < size; ++ i )
    free(stiff[i]);
  free(stiff);
}

/* Create initial stress component of the stiffness matrix */
void solver_local_initial_stess_part(fea_solver_ptr self,int element)
{
  /* matrix of gradients of shape functions */
  shape_gradients_ptr grads = (shape_gradients_ptr)0;
  int gauss,a,b,i,j,k,l,I,J,globalI,globalJ;
  real sum;
  /* size of a local stiffness matrix */
  int size;
  /* number of nodes per element */
  int nelem;
  /* current number of d.o.f */
  int dof;
  /* local stiffness matrix */
  real **stiff = (real**)0;
  
  /* allocate memory for a local stiffness matrix */
  size = self->fea_params_p->nodes_per_element*self->task_p->dof;
  stiff = (real**)malloc(sizeof(real*)*size);
  for (i = 0; i < size; ++ i)
  {
    stiff[i] = (real*)malloc(sizeof(real)*size);
    memset(stiff[i],0,sizeof(real)*size);
  }
  
  dof = self->task_p->dof;
  nelem = self->fea_params_p->nodes_per_element;
  
  /* loop by gauss nodes - numerical integration */
  for (gauss = 0; gauss < self->fea_params_p->gauss_nodes_count ; ++ gauss)
  {
    grads = self->shape_gradients[element][gauss];
    if (grads)
    {
      /* Construct components of stiffness matrix in
       * indical form using Bonet & Wood 7.35 p.207, 1st edition */
      
      /* loop by nodes */
      for ( a = 0; a < nelem; ++ a)
        for (b = 0; b < nelem; ++ b)
        {
          /* loop by d.o.f in a stiffness matrix block [K_{ab}]ij, 3x3 */
          for (i = 0; i < dof; ++ i)
            for (j = 0; j < dof; ++ j)
            {
              /* indicies in a local stiffness matrix */
              I = a*dof + i;
              J = b*dof + j;
              sum = 0.0;
              /* sum of particular derivatives and components of C tensor */
              for (k = 0; k < dof; ++ k)
                for (l = 0; l < dof; ++ l)
                  sum += 
                    grads->grads[k][a] *
                    self->stresses[element][gauss].components[k][l] *
                    grads->grads[l][b] *
                    DELTA(i,j);
              /*
               * multiply by volume of an element = det(J)
               * where divider 6 or 2 or others already accounted in
               * weights of gauss nodes
               */
              sum *= fabs(grads->detJ);
              /* ... and weight of the gauss nodes for  */
              sum *= self->elements_db.gauss_nodes[gauss]->weight;
              /* append to the local stiffness */
              stiff[I][J] += sum;
              /* finally distribute to the global matrix */
              globalI = self->elements_p->elements[element][a]*dof + i;
              globalJ = self->elements_p->elements[element][b]*dof + j;
              sp_matrix_element_add(&self->global_mtx,
                                    globalI,
                                    globalJ,
                                    sum);
            }
        }
    }
  }
    
  /* clear local stiffness */
  for ( i = 0; i < size; ++ i )
    free(stiff[i]);
  free(stiff);
}



void solver_local_residual_forces(fea_solver_ptr self,int element)
{
  /*
   * Calculate residual force vector using formula
   * Bonet & Wood, 1st edition, 7.15
   */
  int a,i,j,I,gauss;
  real sum;
  shape_gradients_ptr grads = (shape_gradients_ptr)0;
  int nelem = self->fea_params_p->nodes_per_element;
  int dof = self->task_p->dof;

  for (gauss = 0; gauss < self->fea_params_p->gauss_nodes_count ; ++ gauss)
  {
    grads = self->shape_gradients[element][gauss];
    if (grads)
    {
      /* loop by nodes */
      for ( a = 0; a < nelem; ++ a)
        /* loop by d.o.f in a residual vector matrix block T_{ai}, 3x1 */
        for (i = 0; i < dof; ++ i)
        {
          sum = 0.0;
          /* sum of particular derivatives and components of C tensor */
          for (j = 0; j < dof; ++ j)
            sum += self->stresses[element][gauss].components[i][j] *
              grads->grads[j][a];
          /*
           * multiply by volume of an element = det(J)
           * where divider 6 or 2 or others already accounted in
           * weights of gauss nodes
           */
          sum *= fabs(grads->detJ);
          /* ... and weight of the gauss node */
          sum *= self->elements_db.gauss_nodes[gauss]->weight;
          /* finally distribute to the global residual forces vector */
          I = self->elements_p->elements[element][a]*dof + i;
          self->global_forces_vct[I] += -sum;
        }
    }
  }

}



void solver_element_gauss_graddef(fea_solver_ptr self,
                                  int element,
                                  int gauss,
                                  real graddef[MAX_DOF][MAX_DOF])
{
  int i,j,k;
  /*
   * There are 2 ways to calculate Deformation gradient
   * First, by using macro CURRENT_SHAPE_GRADIENTS, calculate
   * using gradients of shape functions in current configuration,
   * therefore they shall be obtained using solver_new_shape_gradients
   * function
   */
#ifdef CURRENT_SHAPE_GRADIENTS
  /*
   * Deformation gradient using formula:
   *                              dX_I 
   * F^-1 = \sum\limits_{I,i=1}^3 ----      E_I \otimes e_i
   *                              dx_i
   * See Bonet & Wood 7.6(a,b), 7.7 p.198, 1st edition
   */
  real detF = 0;
  shape_gradients_ptr grads = self->shape_gradients[element][gauss];
  for (i = 0; i < MAX_DOF; ++ i)
  {
    for (j = 0; j < MAX_DOF; ++ j)
    {
      graddef[i][j] = 0;
      for (k = 0; k < self->fea_params_p->nodes_per_element; ++ k)
        graddef[i][j] +=
          grads->grads[j][k] * 
          self->nodes0_p->nodes[self->elements_p->elements[element][k]][i];
    }
  }
  inv3x3(graddef,&detF);
#else /* Second way is to use gradients of shapes in initial configuration */
  /*
   * Deformation gradient could be calculated using the following
   * formula:
   *                                        N_k
   * F_{ij} = \sum\limits_{k=1}^{n} x_{k,i}------
   *                                        dX_j
   * See Bonet & Wood 7.6(a,b), 7.7 p.198, 1st edition
   */

  for (i = 0; i < MAX_DOF; ++ i)
  {
    for (j = 0; j < MAX_DOF; ++ j)
    {
      graddef[i][j] = 0;
      for (k = 0; k < self->fea_params_p->nodes_per_element; ++ k)
        graddef[i][j] +=
          self->shape_gradients0[element][gauss]->grads[j][k] *
          self->nodes_p->nodes[self->elements_p->elements[element][k]][i];
    }
  }
#endif /* CURRENT_SHAPE_GRADIENTS */
}

void solver_element_gauss_stress_A5(fea_solver_ptr self,
                                    real F[MAX_DOF][MAX_DOF],
                                    real stress_tensor[MAX_DOF][MAX_DOF])
{
  int i,j,k;
  real C[MAX_DOF][MAX_DOF];
  real G[MAX_DOF][MAX_DOF];
  real Sn[MAX_DOF][MAX_DOF];
  real lambda,mu;
  real detF = 0;
  real I1 = 0;

  lambda = self->task_p->model.parameters[0];
  mu = self->task_p->model.parameters[1];
  
  /* G = F'*F */
  for (i = 0; i < MAX_DOF; ++ i)
    for (j = 0; j < MAX_DOF; ++ j)
    {
      G[i][j] = 0;
      for ( k = 0; k < MAX_DOF; ++ k)
        G[i][j] += F[k][i]*F[k][j];
    }

  /* C = 0.5(G - E) */
  for (i = 0; i < MAX_DOF; ++ i)
    for (j = 0; j < MAX_DOF; ++ j)
    {
      C[i][j] = 0.5*(G[i][j] - DELTA(i,j));
    }

  /* I1 = 1st invariant of C */
  for (i = 0; i < MAX_DOF; ++ i)
    I1 += C[i][i];
  
  detF = det3x3(F);
  /* Sigma = ( lambda*I1(C)+2mu*C )/det(F)*/
  for (i = 0; i < MAX_DOF; ++ i)
    for (j = 0; j < MAX_DOF; ++ j)
    {
      Sn[i][j] =
        (lambda*I1*DELTA(i,j)+ 2*mu*C[i][j])/detF;
    }
  /* S = F*S*F'; for model A5 */
  /*
   * Calculate in 2 steps, using C as for a temporary results 
   * 1) C = F*Sn
   */
  matrix_mul3x3(F,Sn,C);
  /* 2) S = C*F' */
  matrix_transpose2_mul3x3(C,F,stress_tensor);
}

void solver_element_gauss_stress_compr_neohookean(fea_solver_ptr self,
                                                  real F[MAX_DOF][MAX_DOF],
                                                  real S[MAX_DOF][MAX_DOF])
{
  int i,j,k;
  real B[MAX_DOF][MAX_DOF];
  real lambda,mu;
  real J = det3x3(F);

  lambda = self->task_p->model.parameters[0];
  mu = self->task_p->model.parameters[1];
  
  /* B = F*F' */
  for (i = 0; i < MAX_DOF; ++ i)
    for (j = 0; j < MAX_DOF; ++ j)
    {
      B[i][j] = 0;
      for ( k = 0; k < MAX_DOF; ++ k)
        B[i][j] += F[i][k]*F[j][k];
    }

  /* S = mu/J*(B-E)+lambda/J*log(J)*E; */
  for (i = 0; i < MAX_DOF; ++ i)
    for (j = 0; j < MAX_DOF; ++ j)
    {
      S[i][j] = mu*(B[i][j]-DELTA(i,j))/J + 
        lambda*log(J)*DELTA(i,j)/J;
    }
}


void solver_element_gauss_stress(fea_solver_ptr self,
                                 int element,
                                 int gauss,
                                 real graddef_tensor[MAX_DOF][MAX_DOF],
                                 real stress_tensor[MAX_DOF][MAX_DOF])
  
{
  /* get deformation gradient */
  solver_element_gauss_graddef(self,element,gauss,graddef_tensor);
  self->element_gauss_stress(self,graddef_tensor,stress_tensor);
}


void solver_ctensor_A5(fea_solver_ptr self,
                       real (*graddef)[MAX_DOF],
                       real (*ctensor)[MAX_DOF][MAX_DOF][MAX_DOF])
{
  int i,j,k,l;
  real lambda,mu;
  real detF = det3x3(graddef);
  lambda = self->task_p->model.parameters[0];
  mu = self->task_p->model.parameters[1];
  for ( i = 0; i < MAX_DOF; ++ i )
    for ( j = 0; j < MAX_DOF; ++ j )
      for ( k = 0; k < MAX_DOF; ++ k )
        for ( l = 0; l < MAX_DOF; ++ l )
          ctensor[i][j][k][l] = 
            ( lambda * DELTA (i, j) * DELTA (k, l)  \
              + mu * DELTA (i, k) * DELTA (j, l)    \
              + mu * DELTA (i, l) * DELTA (j, k) ) / detF;
}

void solver_ctensor_compr_neohookean(fea_solver_ptr self,
                                     real (*graddef)[MAX_DOF],
                                     real (*ctensor)[MAX_DOF][MAX_DOF][MAX_DOF])
{
  int i,j,k,l;
  real lambda,mu,lambda1,mu1;
  real J = det3x3(graddef);
  lambda = self->task_p->model.parameters[0];
  mu = self->task_p->model.parameters[1];
  lambda1 = lambda/J;
  mu1 = (mu - lambda*log(J))/J;
  for ( i = 0; i < MAX_DOF; ++ i )
    for ( j = 0; j < MAX_DOF; ++ j )
      for ( k = 0; k < MAX_DOF; ++ k )
        for ( l = 0; l < MAX_DOF; ++ l )
          ctensor[i][j][k][l] = 
            lambda1 * DELTA (i, j) * DELTA (k, l)    \
            + 2*mu1 * DELTA (i, k) * DELTA (j, l);
  
}


void solver_create_forces_bc(fea_solver_ptr self)
{
  /* TODO: implement this */
  if (self)
  {
    /* solver->global_forces_vct */
  }
}

void solver_apply_prescribed_bc(fea_solver_ptr self, real lambda)
{
  solver_apply_bc_general(self,solver_apply_single_bc,lambda);
}

void solver_apply_bc_general(fea_solver_ptr self,apply_bc_t apply,real lambda)
{
  int i,j;
  int type,index,offset,node_number;
  real presc[MAX_DOF];
  for ( i =0; i < self->presc_boundary_p->prescribed_nodes_count; ++ i)
  {
    node_number = self->presc_boundary_p->prescribed_nodes[i].node_number;
    for ( j = 0; j < MAX_DOF; ++ j)
      presc[j] = self->presc_boundary_p->prescribed_nodes[i].values[j]*lambda;
    
    type = self->presc_boundary_p->prescribed_nodes[i].type;

    /* set the index offset depending on condition type */
    if ( type == PRESCRIBEDX || type == PRESCRIBEDXY || 
         type == PRESCRIBEDXZ || type == PRESCRIBEDXYZ )
    {
      offset = 0;
      index = node_number*self->task_p->dof+offset;
      apply(self, index, presc[offset]);
    }
    if ( type == PRESCRIBEDY || type == PRESCRIBEDXY || 
         type == PRESCRIBEDYZ || type == PRESCRIBEDXYZ )
    {
      offset = 1;
      index = node_number*self->task_p->dof+offset;
      apply(self, index, presc[offset]);
    }
    if ( type == PRESCRIBEDZ || type == PRESCRIBEDXZ || 
         type == PRESCRIBEDYZ || type == PRESCRIBEDXYZ )
    {
      offset = 2;
      index = node_number*self->task_p->dof+offset;
      apply(self, index, presc[offset]);
    }
  }

}

void solver_apply_single_bc(fea_solver_ptr self, int index, real presc)
{
  real *pvalue,*pvalue1,value;
  int size = self->global_mtx.rows_count;
  int j;
  pvalue = sp_matrix_element(&self->global_mtx,index,index);
  /* global matrix always shall have
   * nonzero diagonal elements */
  assert(pvalue);
  pvalue1 = pvalue;
  value = *pvalue;
  
  for (j = 0; j < size; ++ j)
  {
    pvalue = sp_matrix_element(&self->global_mtx,j,index);
    if (pvalue)
    {
      self->global_forces_vct[j] -= *pvalue*presc;
      *pvalue = 0;
    }
    pvalue = sp_matrix_element(&self->global_mtx,index,j);
    if (pvalue)
      *pvalue = 0;
  }
  
  *pvalue1 = value;
  self->global_forces_vct[index] = value*presc;
}

void solver_update_node_with_bc(fea_solver_ptr self,
                                int index,
                                real value)
{
  int i = index / self->task_p->dof;
  int j = index % self->task_p->dof;
  self->nodes_p->nodes[i][j] += value;
}



void solver_update_nodes_with_solution(fea_solver_ptr self,
                                       real* x)
{
  int i,j;
  for ( i = 0; i < self->nodes_p->nodes_count; ++ i)
  {
    for ( j = 0; j < self->task_p->dof; ++ j)
      self->nodes_p->nodes[i][j] += x[i*self->task_p->dof + j];
  }
}

void solver_update_nodes_with_bc(fea_solver_ptr self, real lambda)
{
  solver_apply_bc_general(self,solver_update_node_with_bc,lambda) ;
}


real tetrahedra10_isoform(int i,real r,real s,real t)
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
  default: error("tetrahedra10_isoform: wrong index");
  }
  return 0;
}

real tetrahedra10_df_dr(int i,real r,real s,real t)
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
  case 9: return 0;
  default: error("tetrahedra10_df_dr: wrong index");
  }
  return 0;
}

real tetrahedra10_df_ds(int i, real r, real s, real t)
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
  case 9: return 4*t;
  default: error("tetrahedra10_df_ds: wrong index");
  }
  return 0;
}

real tetrahedra10_df_dt(int i, real r, real s, real t)
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
  case 9: return 4*s;
  default: error("tetrahedra10_df_dt: wrong index");    
  }
  return 0;
}

real tetrahedra10_disoform(int shape,int dof,real r,real s,real t)
{
  switch(dof)
  {
  case 0: return tetrahedra10_df_dr(shape,r,s,t);
  case 1: return tetrahedra10_df_ds(shape,r,s,t);
  case 2: return tetrahedra10_df_dt(shape,r,s,t);
  default: error("tetrahedra10_disoform: wrong dof");
  }
  return 0;
}

void solver_export_tetrahedra10_gmsh(fea_solver_ptr solver, char *filename)
{
  /* Our(left) and Gmsh(Right) nodal ordering.
   * 
   * see http://geuz.org/gmsh/doc/texinfo/#Node-ordering
   *
   *    Our Tetrahedron10:                       Gmsh Tetrahedron10:
   * 
   *                    v
   *                  .
   *                ,/
   *               /
   *             2                                     2
   *           ,/|`\                                 ,/|`\
   *         ,/  |  `\                             ,/  |  `\
   *       ,6    '.   `5                         ,6    '.   `5
   *     ,/       9     `\                     ,/       8     `\
   *   ,/         |       `\                 ,/         |       `\
   *  0--------4--'.--------1 --> u         0--------4--'.--------1
   *   `\.         |      ,/                 `\.         |      ,/
   *      `\.      |    ,8                      `\.      |    ,9
   *         `7.   '. ,/                           `7.   '. ,/
   *            `\. |/                                `\. |/
   *               `3                                    `3
   *                `\.
   *                    ` w
   *
   *                    Difference in nodes 8 <=> 9
   */
  FILE* f;
  int i,j,k;
  int load;
 
  f = fopen(filename,"w+");
  if ( f )
  {
    /* Header */
    fprintf(f,"$MeshFormat\n");
    fprintf(f,"2.0 0 8\n");
    fprintf(f,"$EndMeshFormat\n");
    /* Geometry */
    /* Nodes section */
    fprintf(f,"$Nodes\n");
    fprintf(f,"%d\n",solver->nodes_p->nodes_count);
    for (i = 0; i < solver->nodes_p->nodes_count; ++ i)
      fprintf(f,"%d %f %f %f\n",i+1,
              solver->nodes0_p->nodes[i][0],
              solver->nodes0_p->nodes[i][1],
              solver->nodes0_p->nodes[i][2]);
    fprintf(f,"$EndNodes\n");
    /* Elements section */
    fprintf(f,"$Elements\n");
    fprintf(f,"%d\n", solver->elements_p->elements_count);
    for (i = 0; i < solver->elements_p->elements_count; ++ i)
    {
      fprintf(f,"%d 11 3 1 1 1 ",i+1);
      for (j = 0; j < 8; ++ j)
        fprintf(f,"%d ",solver->elements_p->elements[i][j]+1);
      fprintf(f,"%d ",solver->elements_p->elements[i][9]+1);
      fprintf(f,"%d ",solver->elements_p->elements[i][8]+1);
      fprintf(f,"\n");
    }
    fprintf(f,"$EndElements\n");
    
    /* loop by load increments */
    for ( load = 0; load <= solver->current_load_step; ++ load)
    {
      /* Export displacements */
      fprintf(f,"$NodeData\n");
      fprintf(f,"1\n");
      fprintf(f,"\"Displacements\"\n");
      fprintf(f,"1\n");           /* number-of-real-tags */
      fprintf(f,"%f\n", load*0.83333333);    /* timestamp */
      fprintf(f,"3\n");           /* number-of-integer-tags */
      fprintf(f,"%d\n", load);           /* step index (starting at 0) */
      fprintf(f,"3\n");           /* number of field components (1, 3 or 9)*/
      /* number of entities */
      fprintf(f,"%d\n",solver->nodes_p->nodes_count);
      for (i = 0; i < solver->nodes_p->nodes_count; ++ i)
        fprintf(f,"%d %f %f %f\n",i+1,
                load ? solver->load_steps_p[load-1].nodes_p->nodes[i][0] -
                solver->nodes0_p->nodes[i][0] : 0.0,
                load ? solver->load_steps_p[load-1].nodes_p->nodes[i][1] -
                solver->nodes0_p->nodes[i][1] : 0.0,
                load ? solver->load_steps_p[load-1].nodes_p->nodes[i][2] -
                solver->nodes0_p->nodes[i][2] : 0.0);
      fprintf(f,"$EndNodeData\n");
    
      /* Export stresses */
      fprintf(f,"$ElementData\n");
      fprintf(f,"1\n");           /* number-of-string-tags */
      fprintf(f,"\"Stress tensor\"\n"); /* string tag */
      fprintf(f,"1\n");           /* number-of-real-tags */
      fprintf(f,"%f\n",load*0.83333333);         /* timestamp */
      fprintf(f,"3\n");           /* number-of-integer-tags */
      fprintf(f,"%d\n",load);           /* step index (starting at 0) */
      fprintf(f,"9\n");           /* number of field components (1, 3 or 9) */
      /* number of entities */
      fprintf(f,"%d\n",solver->elements_p->elements_count);
      for (i = 0; i < solver->elements_p->elements_count; ++ i)
      {
        fprintf(f,"%d ",i+1);     /* element index */
        for ( j = 0; j < MAX_DOF; ++ j)
          for ( k = 0; k < MAX_DOF; ++ k)
            fprintf(f,"%f ", load ?
                    solver->load_steps_p[load-1].stresses[i][0].components[j][k]
                    : 0.0); 
        fprintf(f,"\n");
      }
      fprintf(f,"$EndElementData\n");
    }
    fclose(f);
  }
}
                
void solver_create_model_params(fea_solver_ptr self)
{
  switch(self->task_p->model.model)
  {
  case MODEL_A5:
    self->element_gauss_stress = solver_element_gauss_stress_A5;
    self->ctensor = solver_ctensor_A5;
    break;
  case MODEL_COMPRESSIBLE_NEOHOOKEAN:
    self->element_gauss_stress =
      solver_element_gauss_stress_compr_neohookean;
    self->ctensor = solver_ctensor_compr_neohookean;
    break;
  default:
    assert(FALSE);
  };

}



void solver_create_element_params_tetrahedra10(fea_solver* solver)
{
  solver->shape = tetrahedra10_isoform;
  solver->dshape = tetrahedra10_disoform;
  switch (solver->fea_params_p->gauss_nodes_count)
  {
  case 4:
    solver->elements_db.gauss_nodes_data = gauss_nodes4_tetr10;
    break;
  case 5:
    solver->elements_db.gauss_nodes_data = gauss_nodes5_tetr10;
    break;
  default: error("solver_create_element_params_tetrahedra10: gauss nodes");
  }
  solver->export_function = solver_export_tetrahedra10_gmsh;
}


fea_task_ptr new_fea_task()
{
  /* allocate memory */
  fea_task_ptr task = (fea_task_ptr)malloc(sizeof(fea_task));
  /* set default values */
  task->desired_tolerance = 1e-8;
  task->dof = 3;
  task->ele_type = TETRAHEDRA10;
  task->linesearch_max = 0;
  task->arclength_max = 0;
  task->load_increments_count = 0;
  task->type = CARTESIAN3D;
  task->modified_newton = TRUE;
  task->model.model = MODEL_A5;
  task->model.parameters_count = 2;
  task->model.parameters[0] = 100;
  task->model.parameters[1] = 100;
  return task;
}

void free_fea_task(fea_task_ptr task)
{
  free(task);
}

/* Initializes fea solution params with default values */
fea_solution_params_ptr new_fea_solution_params()
{
  /* allocate memory */
  fea_solution_params_ptr fea_params = (fea_solution_params_ptr)
    malloc(sizeof(fea_solution_params));
  /* set default values */
  fea_params->gauss_nodes_count = 5;
  fea_params->nodes_per_element = 10;
  return fea_params;
}

/* clear fea solution params */
void free_fea_solution_params(fea_solution_params_ptr params)
{
  free(params);
}

/* Initialize nodes array but not initialize particular arrays  */
nodes_array_ptr new_nodes_array()
{
  /* allocate memory */
  nodes_array_ptr nodes = (nodes_array_ptr)malloc(sizeof(nodes_array));
  /* set zero values */
  nodes->nodes = (real**)0;
  nodes->nodes_count = 0;
  return nodes;
}

nodes_array_ptr new_copy_nodes_array(nodes_array_ptr nodes)
{
  int i;
  /* allocate memory */
  nodes_array_ptr copy = (nodes_array_ptr)malloc(sizeof(nodes_array));
  /* set zero values */
  copy->nodes = (real**)0;
  copy->nodes_count = nodes->nodes_count;
  /* copy nodes */
  if ( nodes->nodes_count && nodes->nodes)
  {
    copy->nodes = (real**)malloc(sizeof(real*)*copy->nodes_count);
    for ( i = 0; i < copy->nodes_count; ++ i)
    {
      copy->nodes[i] = (real*)malloc(sizeof(real)*MAX_DOF);
      memcpy(copy->nodes[i],nodes->nodes[i],sizeof(real)*MAX_DOF);
    }
  }
  return copy;
}

/* carefully deallocate nodes array */
void free_nodes_array(nodes_array_ptr nodes)
{
  int counter = 0;
  if (nodes)
  {
    if (nodes->nodes_count && nodes->nodes)
    {
      for (; counter < nodes->nodes_count; ++ counter)
        free(nodes->nodes[counter]);
      free(nodes->nodes);
    }
    free(nodes);
  }
}


/* Initialize elements array but not initialize particular elements */
elements_array_ptr new_elements_array()
{
  /* allocate memory */
  elements_array_ptr elements = (elements_array_ptr)
    malloc(sizeof(elements_array));
  /* set zero values */
  elements->elements = (int**)0;
  elements->elements_count = 0;
  return elements;
}

void free_elements_array(elements_array_ptr elements)
{
  if(elements)
  {
    int counter = 0;
    if (elements->elements_count && elements->elements)
    {
      for (; counter < elements->elements_count; ++ counter)
        free(elements->elements[counter]);
      free(elements->elements);
    }
    free(elements);
  }
}

/* Initialize boundary nodes array but not initialize particular nodes */
presc_boundary_array_ptr new_presc_boundary_array()
{
  /* allocate memory */
  presc_boundary_array_ptr presc_boundary =
    (presc_boundary_array_ptr)malloc(sizeof(presc_boundary_array));
  /* set zero values */
  presc_boundary->prescribed_nodes = (prescibed_boundary_node*)0;
  presc_boundary->prescribed_nodes_count = 0;
  return presc_boundary;
}

void free_presc_boundary_array(presc_boundary_array_ptr presc)
{
  if (presc)
  {
    if (presc->prescribed_nodes_count && presc->prescribed_nodes)
    {
      free(presc->prescribed_nodes);
    }
    free(presc);
  }
}


BOOL initial_data_load(char *filename,
                       fea_task_ptr *task,
                       fea_solution_params_ptr *fea_params,
                       nodes_array_ptr *nodes,
                       elements_array_ptr *elements,
                       presc_boundary_array_ptr *presc_boundary)
{
  BOOL result = FALSE;
#ifdef USE_EXPAT
  result = expat_data_load(filename,task,fea_params,nodes,elements,
                           presc_boundary);
#endif
  return result;
}

