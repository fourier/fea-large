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
#include "fea_solver.h"

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

  /* Set the application exit handler */
  /* atexit(application_done); */
  
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
  
  /* backup solver to the global_solver for the case of emergency exit */
  global_solver = solver;

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
  solver->export(solver,"deformed.msh");
  
  free_fea_solver(solver);
  global_solver = (fea_solver*)0;
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

void application_done(void)
{
  /* TODO: add additional finalization routines here */
  if (global_solver)
  {
    free_fea_solver(global_solver);
  }
  
}


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

void init_sp_matrix(sp_matrix_ptr mtx,
                    int rows,
                    int cols,
                    int bandwidth,
                    sparse_storage_type type)
{
  int i,n;
  if (mtx)
  {
    mtx->rows_count = rows;
    mtx->cols_count = cols;
    mtx->ordered = FALSE;
    mtx->storage_type = type;
    n = type == CRS ? rows : cols;
    mtx->storage = (indexed_array*)malloc(sizeof(indexed_array)*n);
    /* create rows or cols with fixed bandwidth */
    for (i = 0; i < n; ++ i)
    {
      mtx->storage[i].width = bandwidth;
      mtx->storage[i].last_index = -1;
      mtx->storage[i].indexes = (int*)malloc(sizeof(int)*bandwidth);
      mtx->storage[i].values = (real*)malloc(sizeof(real)*bandwidth);
      memset(mtx->storage[i].indexes,0,sizeof(int)*bandwidth);
      memset(mtx->storage[i].values,0,sizeof(real)*bandwidth);
    }
  }
}


sp_matrix_ptr free_sp_matrix(sp_matrix_ptr mtx)
{
  int i,n;
  if (mtx)
  {
    n = mtx->storage_type  == CRS ? mtx->rows_count : mtx->cols_count;
    for (i = 0; i < n; ++ i)
    {
      free(mtx->storage[i].indexes);
      free(mtx->storage[i].values);
    }
    free(mtx->storage);
    mtx->storage = (indexed_array*)0;
    mtx->cols_count = 0;
    mtx->rows_count = 0;
  }
  return (sp_matrix_ptr)0;
}

void clear_sp_matrix(sp_matrix_ptr mtx)
{
  int i,n;
  if (mtx)
  {
    n = mtx->storage_type  == CRS ? mtx->rows_count : mtx->cols_count;
    for (i = 0; i < n; ++ i)
    {
      memset(mtx->storage[i].values,0,sizeof(real)*(mtx->storage[i].width));
    }
  }
}

void copy_sp_matrix(sp_matrix_ptr mtx_from, sp_matrix_ptr mtx_to)
{
  int i,n;
  assert(mtx_from && mtx_to);
  n = mtx_from->storage_type  == CRS ? mtx_from->rows_count :
    mtx_from->cols_count;
  mtx_to->rows_count = mtx_from->rows_count;
  mtx_to->cols_count = mtx_from->cols_count;
  mtx_to->ordered = mtx_from->ordered;
  mtx_to->storage_type = mtx_from->storage_type;
  mtx_to->storage =
    (indexed_array*)malloc(sizeof(indexed_array)*n);
  /* copy rows */
  for (i = 0; i < n; ++ i)
  {
    memset(&mtx_to->storage[i],0,sizeof(indexed_array));
    mtx_to->storage[i].width = mtx_from->storage[i].width;
    mtx_to->storage[i].last_index = mtx_from->storage[i].last_index;
    mtx_to->storage[i].indexes =
      (int*)malloc(sizeof(int)*mtx_from->storage[i].width);
    mtx_to->storage[i].values =
      (real*)malloc(sizeof(real)*mtx_from->storage[i].width);
    memcpy(mtx_to->storage[i].indexes, mtx_from->storage[i].indexes,
           sizeof(int)*mtx_from->storage[i].width);
    memcpy(mtx_to->storage[i].values, mtx_from->storage[i].values,
           sizeof(real)*mtx_from->storage[i].width);
  }
}

void sp_matrix_convert(sp_matrix_ptr mtx_from,
                       sp_matrix_ptr mtx_to,
                       sparse_storage_type type)
{
  int i,j;
  if (type == mtx_from->storage_type)
    return;
  init_sp_matrix(mtx_to,
                 mtx_from->rows_count,
                 mtx_from->cols_count,
                 mtx_from->storage[0].width,
                 type);
  if (type == CCS)              /* CRS -> CCS */
  {
    for (i = 0; i < mtx_from->rows_count; ++ i)
    {
      for (j = 0; j <= mtx_from->storage[i].last_index; ++ j)
        MTX(mtx_to,i,mtx_from->storage[i].indexes[j],
            mtx_from->storage[i].values[j]);
    }
  }
  else                          /* CCS -> CRS*/
  {
    for (i = 0; i < mtx_from->cols_count; ++ i)
    {
      for (j = 0; j <= mtx_from->storage[i].last_index; ++ j)
        MTX(mtx_to,mtx_from->storage[i].indexes[j],i,
            mtx_from->storage[i].values[j]);
    }
  }
  
}


void sp_matrix_create_ilu(sp_matrix_ptr self,sp_matrix_skyline_ilu_ptr ilu)
{
  sp_matrix_skyline A;
  /* reorder self if not already reordered */
  if (!self->ordered)
    sp_matrix_compress(self);
  /* initialize skyline matrix for ILU decomposition */
  init_sp_matrix_skyline(&A,self);
  /*
   * create ILU decomposition of the sparse matrix in skyline format
   * taking ownership of the skyline A matrix
   */
  init_copy_sp_matrix_skyline_ilu(ilu,&A);
  /*
   * since init_copy_sp_matrix_skyline_ilu takes the ownership
   * of the A matrix it is not needed to free A matrix
   */
}

void init_sp_matrix_skyline(sp_matrix_skyline_ptr self,sp_matrix_ptr mtx)
{
  /*
   * Construct CSLR matrix from the sp_matrix
   * with symmetric portrait
   */
  int i,j,k,iptr,l_count,u_count,column;
  real* pvalue = 0;

  /* assert what mtx is already reordered */
  assert(mtx->ordered == TRUE);
  /* currenty implemented conversion only from CRS format */
  assert(mtx->storage_type == CRS);
  
  self->rows_count = mtx->rows_count;
  self->cols_count = mtx->cols_count;
  /*
   * get an information about number of nonzero elements
   */
  self->nonzeros = 0;
  for (i = 0; i < mtx->rows_count; ++ i)
    self->nonzeros += mtx->storage[i].last_index + 1;
  
  /* calculate number of upper-triangle elements */
  l_count = 0;
  u_count = 0;
  for (i = 0; i < mtx->rows_count; ++ i)
    for (j = 0; j <= mtx->storage[i].last_index; ++ j)
      if ( mtx->storage[i].indexes[j] > i)
        u_count ++;
      else if (mtx->storage[i].indexes[j] < i)
        l_count ++;
  /*
   * check if the number of upper triangle nonzero elements
   * is the same as number of lower triangle nonzero elements
   */
  assert(l_count == u_count);
  self->tr_nonzeros = l_count;
  
  /* allocate memory for arrays */
  self->diag = (real*)malloc(sizeof(real)*mtx->rows_count);
  self->lower_triangle = l_count ? (real*)malloc(sizeof(real)*l_count) : 0;
  self->upper_triangle = u_count ? (real*)malloc(sizeof(real)*u_count) : 0;
  self->jptr = l_count ? (int*)malloc(sizeof(int)*l_count)  : 0;
  self->iptr = (int*)malloc(sizeof(int)*(mtx->rows_count+1));

  /* fill diagonal */
  for (i = 0; i < mtx->rows_count; ++ i)
  {
    pvalue = sp_matrix_element(mtx,i,i);
    self->diag[i] = pvalue ? *pvalue : 0;
  }
  /* now fill arrays with proper values */
  u_count = 0,l_count = 0;
  for (i = 0; i < mtx->rows_count; ++ i)
  {
    iptr = -1;
    self->iptr[i] = 0;
    for (j = 0; j <= mtx->storage[i].last_index; ++ j)
    {
      if ( mtx->storage[i].indexes[j] < i)
      {
        /*
         * set a flag what we found the first nonzero element in
         * current row in lower triangle
         */
        if (iptr == -1)
          iptr  = l_count;
        /* fill lower triangle values */
        column = mtx->storage[i].indexes[j];
        self->jptr[l_count] = column;
        self->lower_triangle[l_count] = mtx->storage[i].values[j];
        /* fill upper triangle values - column-wise */
        for ( k = 0; k <= mtx->storage[column].last_index; ++ k)
          if (mtx->storage[column].indexes[k] == i)
          {
            self->upper_triangle[l_count] =
              mtx->storage[column].values[k];
            break;
          }
        l_count ++;
      }
    }
    self->iptr[i] = iptr == -1 ? l_count : iptr;
  }
  /* finalize iptr array */
  self->iptr[i] = self->tr_nonzeros;
}

void free_sp_matrix_skyline(sp_matrix_skyline_ptr self)
{
  if (self)
  {
    self->rows_count = 0;
    self->cols_count = 0;
    self->nonzeros = 0;
    self->tr_nonzeros = 0;
    free(self->diag);
    free(self->lower_triangle);
    free(self->upper_triangle);
    free(self->jptr);
    free(self->iptr);
  }
}


real* sp_matrix_element(sp_matrix_ptr self,int i, int j)
{
  int index;
  /* check for matrix and if i,j are proper indicies */
  assert(self && 
         (i >= 0 && i < self->rows_count ) &&
         (j >= 0 && j < self->cols_count ));
  {
    if (self->storage_type == CRS)
    {
      /* loop by nonzero columns in row i */
      for (index = 0; index <= self->storage[i].last_index; ++ index)
        if (self->storage[i].indexes[index] == j)
          return &self->storage[i].values[index];
    }
    else                        /* CCS */
    {
      /* loop by nonzero rows in column i */
      for (index = 0; index <= self->storage[j].last_index; ++ index)
        if (self->storage[j].indexes[index] == i)
          return &self->storage[j].values[index];
    }
  }
  return (real*)0;
}

real sp_matrix_element_add(sp_matrix_ptr self,int i, int j, real value)
{
  int index,new_width,I,J;
  int* indexes = (int*)0;
  real* values = (real*)0;
  /* check for matrix and if i,j are proper indicies */
  assert (self && 
          (i >= 0 && i < self->rows_count ) &&
          (j >= 0 && j < self->cols_count ));
  /* set I and J to be i and j in case of CRS or j and i otherwise */
  I = self->storage_type == CRS ? i : j;
  J = self->storage_type == CRS ? j : i;
  /* loop by nonzero columns in row/col i */
  for (index = 0; index <= self->storage[I].last_index; ++ index)
    if (self->storage[I].indexes[index] == J)
    {
      /* nonzerod element found, add to it */
      self->storage[I].values[index] += value;
      return self->storage[I].values[index];
    }
  /* needed to add a new element to the row/col */
    
  /*
   * check if bandwidth is not exceed and reallocate memory
   * if necessary
   */
  if (self->storage[I].last_index == self->storage[I].width - 1)
  {
    new_width = self->storage[I].width*2;
    indexes = (int*)realloc(self->storage[I].indexes,new_width*sizeof(int));
    assert(indexes);
    self->storage[I].indexes = indexes;
    values = (real*)realloc(self->storage[I].values,new_width*sizeof(real));
    assert(values);
    self->storage[I].values = values;
    self->storage[I].width = new_width;
    self->ordered = FALSE;
  }
  /* add an element to the row/col */
  self->storage[I].last_index++;
  self->storage[I].values[self->storage[I].last_index] = value;
  self->storage[I].indexes[self->storage[I].last_index] = J;
  return value;
}

/* Swap 2 elements of the indexed array */
void indexed_array_swap(indexed_array_ptr self,int i, int j)
{
  int tmp_idx;
  real tmp_val;
  tmp_idx = self->indexes[i];
  self->indexes[i] = self->indexes[j];
  self->indexes[j] = tmp_idx;
  tmp_val = self->values[i];
  self->values[i] = self->values[j];
  self->values[j] = tmp_val;
}

void indexed_array_sort(indexed_array_ptr self, int l, int r)
{
  /*
   * Quick sort procedure for indexed(compressed) arrays
   * for example rows for CRS sparse matrix or columns for CSC
   * sparse matrix
   */
  int pivot,i;
  int tmp_idx;

  /* boundary checks */
  if (l < r)
  {
    if ( r - l == 1)
    {
      if (self->indexes[l] > self->indexes[r])
        indexed_array_swap(self,r,l);
      return;
    }
    /* choose the pivoting element */
    pivot = (int)((r+l)/2.);
    /* in-place partition procedure - move all elements
     * lower than pivoting to the left, greater to the right */
    tmp_idx  = self->indexes[pivot];
    indexed_array_swap(self,pivot,r);
    pivot = l;
    for ( i = l; i < r; ++ i)
    {
      if (self->indexes[i] <= tmp_idx )
      {
        indexed_array_swap(self,i,pivot);
        pivot++;
      }
    }
    indexed_array_swap(self,r,pivot);
    /* repeat procedure for the left and right parts of an array */
    indexed_array_sort(self,l,pivot-1);
    indexed_array_sort(self,pivot+1,r);
  }
}

void indexed_array_printf(indexed_array_ptr self)
{
  int i;
  if (self)
  {
    printf("indexes = [");
    for (i = 0; i <= self->last_index; ++ i)
      printf("%d,\t",self->indexes[i]);
    printf("%d]\n",self->indexes[i]);
    printf("values  = [");
    for (i = 0; i <= self->last_index; ++ i)
      printf("%f,\t",self->values[i]);
    printf("%f]\n",self->values[i]);
  }
}


void sp_matrix_compress(sp_matrix_ptr self)
{
  int i,j,n;
  n = self->storage_type == CRS ? self->rows_count : self->cols_count;
  for (i = 0; i < n; ++ i)
    indexed_array_sort(&self->storage[i],0,self->storage[i].last_index);
  
  for (i = 0; i < n; ++ i)
    for (j = 0; j < self->storage[i].last_index; ++ j)
      /* assert(self->storage[i].indexes[j] < self->storage[i].indexes[j+1]); */
      if (self->storage[i].indexes[j] >= self->storage[i].indexes[j+1])
      {
        printf("row %d\n",i);
        indexed_array_printf(&self->storage[i]);
        assert(FALSE);
      }
  self->ordered = TRUE;
}

void sp_matrix_mv(sp_matrix_ptr self,real* x, real* y)
{
  int i,j;
  memset(y,0,sizeof(real)*self->rows_count);
  if (self->storage_type == CRS)
  {
    for ( i = 0; i < self->rows_count; ++ i)
    {
      for ( j = 0; j <= self->storage[i].last_index; ++ j)
        y[i] += self->storage[i].values[j]*x[self->storage[i].indexes[j]];
    }
  }
  else                          /* CCS */
  {
    for ( j = 0; j < self->cols_count; ++ j)
    {
      for ( i = 0; i<= self->storage[j].last_index; ++ i)
        y[self->storage[j].indexes[i]] += self->storage[j].values[i]*x[j];
    }
  }
}

void sp_matrix_lower_solve(sp_matrix_ptr self,
                           int n,
                           real* b,
                           real* x)
{
  int i,j;
  assert(self);
  assert(b);
  assert(x);
  assert(n>0 && n <= self->rows_count);

  memset(x,0,sizeof(real)*n);
  
  if (!self->ordered)
    sp_matrix_compress(self);
  if (self->storage_type == CCS)
  {
    for ( j = 0; j < n; ++ j)
      x[j] = b[j];
    for ( j = 0; j < n; ++ j)
    {
      x[j] /= self->storage[j].values[0]; 
      for (i = 1; i <= self->storage[j].last_index; ++ i)
        x[self->storage[j].indexes[i]] -= x[j]*self->storage[j].values[i];
    }
  }
  else                          /* CRS */
  {
    for ( i = 0; i < n; ++ i)
    {
      x[i] = b[i];
      for (j = 0; j <= self->storage[i].last_index &&
             self->storage[i].indexes[j] <= i-1; ++ j)
        x[i] -= x[self->storage[i].indexes[j]]*self->storage[i].values[j];
      x[i] /= self->storage[i].values[self->storage[i].last_index];
    }
  }
}

void sp_matrix_solve(sp_matrix_ptr self,real* b,real* x)
{
  real tolerance = 1e-15;
  int max_iter = 20000;
  int i;
  real tol = 0;
  real* r = malloc(self->rows_count*sizeof(real));
  memset(r,0,self->rows_count*sizeof(real));
  /* reorder columns for to prepare to solve SLAE */
  sp_matrix_compress(self);
  /* sp_matrix_solve_pcg(self,b,b,&max_iter,&tolerance,x); */
  sp_matrix_solve_cg(self,b,b,&max_iter,&tolerance,x);
  /* Calculare residual r = A*x-b */
  sp_matrix_mv(self,x,r);
  for ( i = 0; i < self->rows_count; ++ i)
    r[i] -= b[i];
  /* Find a norm of residual vector */
  for ( i = 0; i < self->rows_count; ++ i)
    tol += r[i]*r[i];
  tol = sqrt(tol);
  /* TODO: move iter, tolerance1 and tolerance2 to the output parameters */
  printf("iter = %d, tolerance1 = %e, tolerance2 = %e\n",
         max_iter,tolerance,tol);
  free(r);
}

void sp_matrix_solve_cg(sp_matrix_ptr self,
                        real* b,
                        real* x0,
                        int* max_iter,
                        real* tolerance,
                        real* x)
{
  /* Conjugate Gradient Algorithm */
  /*
   * Taken from the book:
   * Saad Y. Iterative methods for sparse linear systems (2ed., 2000)
   * page 178
   */
   
  /* variables */
  int i,j;
  real alpha, beta,a1,a2;
  real residn = 0;
  int size = sizeof(real)*self->rows_count;
  int msize = self->rows_count;
  int max_iterations = max_iter ? *max_iter : MAX_ITER;
  real tol = tolerance ? *tolerance : TOLERANCE;
  real* r;              /* residual */
  real* p;              /* search direction */
  real* temp;

  /* check if matrix is reordered */
  if (!self->ordered)
    sp_matrix_compress(self);
  
  /* allocate memory for vectors */
  r = malloc(size);
  p = malloc(size);
  temp = malloc(size);
  /* clear vectors */
  memset(r,0,size);
  memset(p,0,size);
  memset(temp,0,size);

  /* x = x_0 */
  for ( i = 0; i < msize; ++ i)
    x[i] = x0[i];

  /* r_0 = b - A*x_0 */
  sp_matrix_mv(self,b,r);
  for ( i = 0; i < msize; ++ i)
    r[i] = b[i] - r[i];

  /* p_0 = r_0 */
  memcpy(p,r,size);
  
  /* CG loop */
  for ( j = 0; j < max_iterations; j ++ )
  {
    /* temp = A*p_j */
    sp_matrix_mv(self,p,temp);
    /* compute (r_j,r_j) and (A*p_j,p_j) */
    a1 = 0; a2 = 0;
    for (i = 0; i < msize; ++ i)
    {
      a1 += r[i]*r[i]; /* (r_j,r_j) */
      a2 += p[i]*temp[i];      /* (A*p_j,p_j) */
    }

    /*            (r_j,r_j) 
     * alpha_j = -----------
     *           (A*p_j,p_j)
     */                     
    alpha = a1/a2;              
                                
    /* x_{j+1} = x_j+alpha_j*p_j */
    for (i = 0; i < msize; ++ i)
      x[i] += alpha*p[i];
    
    /* r_{j+1} = r_j-alpha_j*A*p_j */
    for (i = 0; i < msize; ++ i)
      r[i] -= alpha*temp[i]; 

    /* check for convergence */
    residn = fabs(r[0]);
    for (i = 1; i < msize; ++ i )
      if (fabs(r[i]) > residn) residn = fabs(r[i]);
    if (residn < tol )
      break;

    /* compute (r_{j+1},r_{j+1}) */
    a2 = 0;
    for (i = 0; i < msize; ++ i)
      a2 += r[i]*r[i];

    /* b_j = (r_{j+1},r_{j+1})/(r_j,r_j) */
    beta = a2/a1;
    
    /* d_{j+1} = r_{j+1} + beta_j*d_j */
    for (i = 0; i < msize; ++ i)
      p[i] = r[i] + beta*p[i];
  }
  *max_iter = j;
  *tolerance = residn;
  
  free(r);
  free(p);
  free(temp);
}

void sp_matrix_solve_pcg_ilu(sp_matrix_ptr self,
                             sp_matrix_skyline_ilu_ptr ILU,                         
                             real* b,
                             real* x0,
                             int* max_iter,
                             real* tolerance,
                             real* x)
{
  /* Preconditioned Conjugate Gradient Algorithm */
  /*
   * Taken from the book:
   * Saad Y. Iterative methods for sparse linear systems (2ed., 2000)
   * page 246
   *
   * Preconditioner: Incomplete LU decomposition (ILU)
   * M = L*U, A = M-R
   */

  /* variables */
  int i,j;
  real alpha, beta,a1,a2;
  real residn = 0;
  int size = sizeof(real)*self->rows_count;
  int msize = self->rows_count;
  int max_iterations = max_iter ? *max_iter : MAX_ITER;
  real tol = tolerance ? *tolerance : TOLERANCE;
  
  real* r;              /* residual */
  real* r1;             /* backup of the residual */
  real* p;              /* search direction */
  real* z;              /* z = M^{-1}*r */
  real* temp;
  
  /* allocate memory for vectors */
  r = malloc(size);
  r1 = malloc(size);
  p = malloc(size);
  z = malloc(size);
  temp = malloc(size);
  
  /* clear vectors */
  memset(r,0,size);
  memset(r1,0,size);
  memset(p,0,size);
  memset(z,0,size);
  memset(temp,0,size);

  /* x = x_0 */
  for ( i = 0; i < msize; ++ i)
    x[i] = x0[i];

  /* r_0 = b - A*x_0 */
  sp_matrix_mv(self,b,r);
  for ( i = 0; i < msize; ++ i)
    r[i] = b[i] - r[i];
  
  /* backup residual */
  memcpy(r1,r,size);
  /* z_0 = M^{-1}*r_0 */
  /*
   * to solve system L*U*x = b
   * y = U*x, => L*y = b
   * U*x = y => x
   */ 
  sp_matrix_skyline_ilu_lower_solve(ILU,r1,temp); /* temp = L^{-1}*r */
  /* r1 now changed, temp contains solution */
  sp_matrix_skyline_ilu_upper_solve(ILU,temp,z); /* z = U^{-1}*temp */
  /* temp now changed, z contains solution*/
  
  /* p_0 = z_0 */
  memcpy(p,z,size);
  
  /* CG loop */
  for ( j = 0; j < max_iterations; j ++ )
  {
    /* temp = A*p_j */
    memset(temp,0,size);
    sp_matrix_mv(self,p,temp);
    /* compute (r_j,z_j) and (A*p_j,p_j) */
    a1 = 0; a2 = 0;
    for (i = 0; i < msize; ++ i)
    {
      a1 += r[i]*z[i]; /* (r_j,z_j) */
      a2 += p[i]*temp[i];      /* (A*p_j,p_j) */
    }

    /*            (r_j,z_j) 
     * alpha_j = -----------
     *           (A*p_j,p_j)
     */                     
    alpha = a1/a2;              
                                
    /* x_{j+1} = x_j+alpha_j*p_j */
    for (i = 0; i < msize; ++ i)
      x[i] += alpha*p[i];
    
    /* r_{j+1} = r_j-alpha_j*A*p_j */
    for (i = 0; i < msize; ++ i)
      r[i] -= alpha*temp[i]; 

    /* check for convergence */
    residn = fabs(r[0]);
    for (i = 1; i < msize; ++ i )
      if (fabs(r[i]) > residn) residn = fabs(r[i]);
    if (residn < tol )
      break;

    /* z_{j+1} = M^{-1}*r_{j+1} */
    memcpy(r1,r,size);
    memset(temp,0,size);
    sp_matrix_skyline_ilu_lower_solve(ILU,r1,temp); /* temp = L^{-1}*r */
    sp_matrix_skyline_ilu_upper_solve(ILU,temp,z); /* z = U^{-1}*temp */

    
    /* compute (r_{j+1},z_{j+1}) */
    a2 = 0;
    for (i = 0; i < msize; ++ i)
      a2 += r[i]*z[i];

    /* b_j = (r_{j+1},z_{j+1})/(r_j,z_j) */
    beta = a2/a1;
    
    /* d_{j+1} = r_{j+1} + beta_j*d_j */
    for (i = 0; i < msize; ++ i)
      p[i] = z[i] + beta*p[i];
  }
  *max_iter = j;
  *tolerance = residn;
  
  /* free vectors */
  free(r);
  free(r1);
  free(z);
  free(p);
  free(temp);
}

void init_copy_sp_matrix_skyline_ilu(sp_matrix_skyline_ilu_ptr self,
                                     sp_matrix_skyline_ptr parent)
{
  int i,j,k,l,q;
  real sum;
  
  /* copy parent member-wise */
  self->parent = *parent;
  /* allocate memory for ILU decomposition arrays */
  self->ilu_diag = (real*)malloc(sizeof(real)*parent->rows_count);
  self->ilu_lowertr = (real*)malloc(sizeof(real)*parent->tr_nonzeros);
  self->ilu_uppertr = (real*)malloc(sizeof(real)*parent->tr_nonzeros);
  /* clear arrays before construction of the ILU decomposition */
  memset(self->ilu_diag,0,sizeof(real)*parent->rows_count);
  memset(self->ilu_lowertr,0,sizeof(real)*parent->tr_nonzeros);
  memset(self->ilu_uppertr,0,sizeof(real)*parent->tr_nonzeros);

  for (k = 0; k < parent->rows_count; ++ k)
  {
    for ( j = parent->iptr[k]; j < parent->iptr[k+1]; ++ j)
    {
      /*
       * L_{kj} = (A_{kj} - \sum\limits_{i=1}^{j-1}L_{ki}U_{ij}/U_{jj}
       * calculate using L_{k,jptr[j]}
       */
      sum = 0;
      q = parent->jptr[j];        /* column index */
      for ( i = parent->iptr[k]; i < parent->iptr[k+1]; ++ i)
      {
        for ( l = parent->iptr[q]; l < parent->iptr[q+1]; ++ l)
        {
          /* if row and column indicies are the same */
          if ( parent->jptr[i] == parent->jptr[l] )
            sum += self->ilu_lowertr[i]*self->ilu_uppertr[l];
        }
      }
      self->ilu_lowertr[j] =
        (parent->lower_triangle[j] - sum)/self->ilu_diag[q];
    }

    /*
     * U_{kk} = A_{kk} -
     * \sum\limits_{i=1}^{k-1} L_{ki}U_{ik}
     */
    sum = 0;
    for ( i = parent->iptr[k]; i < parent->iptr[k+1]; ++ i)
      sum += self->ilu_lowertr[i]*self->ilu_uppertr[i];
    self->ilu_diag[k] = parent->diag[k] - sum;

    for (j = k; j < parent->rows_count; ++ j)
    {
      for ( q = parent->iptr[j]; q < parent->iptr[j+1]; ++ q)
        if (k == parent->jptr[q])
        {
          /*
           * U_{kj} = A_{kj} -
           * \sum\limits_{i=1}^{k-1}L_{ki}U_{ij}
           */
          sum = 0;
          /*
           * i = iptr[k]:iptr[k+1]-1 are coordinates of the
           * k-th row in lower matrix array (self->ilu_lowertr)
           * l = iptr[j]:iptr[j+1]-1 are coordinates of the
           * j-th column in upper matrix array (self->ilu_uppertr)
           */
          for ( i = parent->iptr[k]; i < parent->iptr[k+1]; ++ i)
            for ( l = parent->iptr[j]; l < parent->iptr[j+1]; ++ l)
            {
              /* if row and column indicies are the same */
              if ( parent->jptr[i] == parent->jptr[l] )
                sum += self->ilu_lowertr[i]*self->ilu_uppertr[l];
            }
          self->ilu_uppertr[q] = parent->upper_triangle[q] - sum;
        }
    }
  }
}

void free_sp_matrix_skyline_ilu(sp_matrix_skyline_ilu_ptr self)
{
  free(self->ilu_diag);
  free(self->ilu_lowertr);
  free(self->ilu_uppertr);
  free_sp_matrix_skyline(&self->parent);
}


void sp_matrix_skyline_ilu_lower_mv(sp_matrix_skyline_ilu_ptr self,
                                    real* x,
                                    real* y)
{
  int i,j;
  memset(y,0,sizeof(real)*self->parent.rows_count);
  for (i = 0; i < self->parent.rows_count; ++ i)
  {
    y[i] = x[i];
    for (j = self->parent.iptr[i]; j < self->parent.iptr[i+1]; ++ j)
      y[i] += x[self->parent.jptr[j]]*self->ilu_lowertr[j];
  }
}

void sp_matrix_skyline_ilu_upper_mv(sp_matrix_skyline_ilu_ptr self,
                                    real* x,
                                    real* y)
{
  int i,j;
  memset(y,0,sizeof(real)*self->parent.rows_count);

  for (i = 0; i < self->parent.rows_count; ++ i)
    y[i] = x[i]*self->ilu_diag[i];
  for (i = 0; i < self->parent.rows_count; ++ i)
  {
    for ( j = self->parent.iptr[i]; j < self->parent.iptr[i+1]; ++j )
      y[self->parent.jptr[j]] += x[i]*self->ilu_uppertr[j];
  }
}

void sp_matrix_skyline_ilu_lower_solve(sp_matrix_skyline_ilu_ptr self,
                                       real* b,
                                       real* x)
{
  int i,j;
  memset(x,0,sizeof(real)*self->parent.rows_count);
  
  for ( i = 0; i < self->parent.rows_count; ++ i)
  {
    for (j = self->parent.iptr[i]; j < self->parent.iptr[i+1]; ++ j)
      b[i] -= x[self->parent.jptr[j]]*self->ilu_lowertr[j];
    x[i] = b[i];
  }
}

void sp_matrix_skyline_ilu_upper_solve(sp_matrix_skyline_ilu_ptr self,
                                       real* b,
                                       real* x)
{
  int i,j;
  memset(x,0,sizeof(real)*self->parent.rows_count);

  for ( i = self->parent.rows_count-1; i >= 0; -- i)
  {
    x[i] = b[i]/self->ilu_diag[i];
    for (j = self->parent.iptr[i]; j < self->parent.iptr[i+1]; ++ j)
      b[self->parent.jptr[j]] -= x[i]*self->ilu_uppertr[j];
  }
}

void sp_matrix_printf(sp_matrix_ptr self)
{
  int i,n;
  n = self->storage_type == CRS ? self->rows_count : self->cols_count;
  for (i = 0; i < n; ++ i)
    indexed_array_printf(&self->storage[i]);
}

#ifdef DUMP_DATA
void sp_matrix_dump(sp_matrix_ptr self,char* filename)
{
  FILE* f;
  int i,j;
  real *pvalue,value;

  if ((f = fopen(filename,"w+")))
  {
    for (i = 0; i < self->rows_count; ++ i)
    {
      for (j = 0; j < self->rows_count; ++ j)
      {
        pvalue = sp_matrix_element(self,i,j);
        value = pvalue ? *pvalue : 0;
        fprintf(f,"%e ",value);
      }
      fprintf(f,"\n");
    }
    fflush(f);
    fclose(f);
  }
}


void sp_matrix_skyline_dump(sp_matrix_skyline_ptr self,char* filename)
{
  int i;
  FILE* f;
  if ((f = fopen("global_matrix_skyline.txt","w+")))
  {
    fprintf(f,"adiag = [");
    for ( i = 0; i < self->rows_count; ++ i )
      fprintf(f,"%f ",self->diag[i]);
    fprintf(f,"]\n");

    fprintf(f,"altr = [");
    for ( i = 0; i < self->tr_nonzeros; ++ i )
      fprintf(f,"%f ",self->lower_triangle[i]);
    fprintf(f,"]\n");

    fprintf(f,"autr = [");
    for ( i = 0; i < self->tr_nonzeros; ++ i )
      fprintf(f,"%f ",self->upper_triangle[i]);
    fprintf(f,"]\n");
  
    fprintf(f,"jptr = [");
    for ( i = 0; i < self->tr_nonzeros; ++ i )
      fprintf(f,"%d ",self->jptr[i]+1);
    fprintf(f,"]\n");

    fprintf(f,"iptr = [");
    for ( i = 0; i < self->rows_count; ++ i )
      fprintf(f,"%d ",self->iptr[i]+1);
    fprintf(f,"]\n");
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
  fea_solver_ptr solver = malloc(sizeof(fea_solver));
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
  solver->shape_gradients0 = malloc(sizeof(shape_gradients_ptr*)*elnum);
  solver->shape_gradients  = malloc(sizeof(shape_gradients_ptr*)*elnum);
  solver->stresses = malloc(sizeof(tensor*)*elnum);
  solver->graddefs = malloc(sizeof(tensor*)*elnum);
  for (i = 0; i < elnum; ++ i)
  {
    solver->stresses[i] = malloc(sizeof(tensor)*gauss_count);
    solver->graddefs[i] = malloc(sizeof(tensor)*gauss_count);
    solver->shape_gradients0[i] =
      malloc(sizeof(shape_gradients_ptr)*gauss_count);
    solver->shape_gradients[i] =
      malloc(sizeof(shape_gradients_ptr)*gauss_count);
    
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
    node = malloc(sizeof(gauss_node));
    /* set the weight for this gauss node */
    node->weight = self->elements_db.gauss_nodes_data[gauss_node_index][0];
    /* set shape function values and their derivatives for this node */
    node->forms =
      malloc(sizeof(real)*(self->fea_params_p->nodes_per_element));
    node->dforms = malloc(sizeof(real*)*(self->task_p->dof));
    for ( i = 0; i < self->task_p->dof; ++ i)
      node->dforms[i] =
        malloc(sizeof(real)*(self->fea_params_p->nodes_per_element));
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
      malloc(sizeof(gauss_node*)*gauss_count);
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
    step->stresses = malloc(sizeof(tensor*)*elnum);
    step->graddefs = malloc(sizeof(tensor*)*elnum);

    for (i = 0; i < elnum; ++ i)
    {
      step->stresses[i] = malloc(sizeof(tensor)*gauss_count);
      step->graddefs[i] = malloc(sizeof(tensor)*gauss_count);    
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
  solver->export = solver_export_tetrahedra10_gmsh;
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
    copy->nodes = malloc(sizeof(real*)*copy->nodes_count);
    for ( i = 0; i < copy->nodes_count; ++ i)
    {
      copy->nodes[i] = malloc(sizeof(real)*MAX_DOF);
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
    malloc(sizeof(presc_boundary_array));
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



#ifdef USE_EXPAT

/* All known XML tags */
typedef enum xml_format_tags_enum {
  UNKNOWN_TAG,
  TASK,
  MODEL,
  MODEL_PARAMETERS,
  SOLUTION,
  ELEMENT_TYPE,
  LINE_SEARCH,
  ARC_LENGTH,
  INPUT_DATA,
  GEOMETRY,
  NODES,
  NODE,
  ELEMENTS,
  ELEMENT,
  BOUNDARY_CONDITIONS,
  PRESCRIBED_DISPLACEMENTS,
  PRESC_NODE
} xml_format_tags;

/* An input data structure used in parser */
typedef struct parse_data_tag {
  fea_task *task;
  fea_solution_params *fea_params;
  nodes_array *nodes;
  elements_array *elements;
  presc_boundary_array *presc_boundary;
  xml_format_tags parent_tag;
  char* current_text;
  int current_size;
} parse_data;


/*************************************************************/
/* Declarations of functions used                            */

/* Case-insensitive string comparsion procedure */
int istrcmp(char *s1,char *s2);

/* Convert particular string to the XML tag enum */
xml_format_tags tagname_to_enum(const XML_Char* name);

/*
 * Remove leading and trailing whitespaces from the string,
 * allocating null-terminated string as a result
 */
char *trim_whitespaces(const char* string,size_t size);

/* Expat start tag handler */
void expat_start_tag_handler(void *userData,
                             const XML_Char *name,
                             const XML_Char **atts);

/* Expat End tag handler */
void expat_end_tag_handler(void *userData,
                           const XML_Char *name);


/*    
 * Expat handler for the text between tags
 * Since this function could be called several times for the current tag,
 * it is necessary to store text somewhere. We use parse_data->current_text
 * pointer and parse_data->current_size for these purposes
 */
void expat_text_handler(void *userData,
                        const XML_Char *s,
                        int len);

/* test if an attribute name is what expected(attribute_name)
 * and increase pointer to the next attribute if yes*/
BOOL check_attribute(const char* attribute_name, const XML_Char ***atts);

/* model tag handler */
void process_model_type(parse_data* data, const XML_Char **atts);

/* model-parameters tag handler */
void process_model_params(parse_data* data, const XML_Char **atts);

/* solution tag handler */
void process_solution(parse_data* data, const XML_Char **atts);

/* element-type tag handler */
void process_element_type(parse_data* data, const XML_Char **atts);

/* line-search tag handler */
void process_line_search(parse_data* data, const XML_Char **atts);

/* arc-length tag handler */
void process_arc_length(parse_data* data, const XML_Char **atts);

/* nodes tag handler */
void process_nodes(parse_data* data, const XML_Char **atts);

/* node tag handler */
void process_node(parse_data* data, const XML_Char **atts);

/* elements tag handler */
void process_elements(parse_data* data, const XML_Char **atts);

/* take the node id/position from the element attributes
 * like 'node1' or 'node10'
 * returns -1 in case of wrong attribute name
 * but not skip it in this case!
 */
int node_position_from_attr(const XML_Char ***atts);

/* element tag handler */
void process_element(parse_data* data, const XML_Char **atts);

/* prescribed-displacements tag handler */
void process_prescribed_displacements(parse_data* data, const XML_Char **atts);

/* node tag handler */
void process_prescribed_node(parse_data* data, const XML_Char **atts);

/* main function for loading data using expat XML processor */
BOOL expat_data_load(char *filename,
                     fea_task **task,
                     fea_solution_params **fea_params,
                     nodes_array **nodes,
                     elements_array **elements,
                     presc_boundary_array **presc_boundary);

/*************************************************************/
/* Definition of functions used                              */

int istrcmp(char *s1,char *s2)
{
  /* case insensitive comparison */
  int d;
  for (;;) {
#ifdef ASCII_CTYPE
    if (!isascii(*s1) || !isascii(*s2))
      d = *s1 - *s2;
    else
#endif
      d = (tolower((unsigned char) *s1) - tolower((unsigned char)*s2));
    if ( d != 0 || *s1 == '\0' || *s2 == '\0' )
      return d;
    ++s1;
    ++s2;
  }
  /*NOTREACHED*/
}

xml_format_tags tagname_to_enum(const XML_Char* name)
{
  if (!istrcmp((char*)name,"TASK")) return TASK;
  if (!istrcmp((char*)name,"MODEL")) return MODEL;
  if (!istrcmp((char*)name,"MODEL-PARAMETERS")) return MODEL_PARAMETERS;
  if (!istrcmp((char*)name,"SOLUTION")) return SOLUTION;
  if (!istrcmp((char*)name,"ELEMENT-TYPE")) return ELEMENT_TYPE;
  if (!istrcmp((char*)name,"LINE-SEARCH")) return LINE_SEARCH;
  if (!istrcmp((char*)name,"ARC-LENGTH")) return ARC_LENGTH;
  if (!istrcmp((char*)name,"INPUT-DATA")) return INPUT_DATA;
  if (!istrcmp((char*)name,"GEOMETRY")) return GEOMETRY;
  if (!istrcmp((char*)name,"NODES")) return NODES;
  if (!istrcmp((char*)name,"NODE")) return NODE;
  if (!istrcmp((char*)name,"ELEMENTS")) return ELEMENTS;
  if (!istrcmp((char*)name,"ELEMENT")) return ELEMENT;
  if (!istrcmp((char*)name,"BOUNDARY-CONDITIONS"))
    return BOUNDARY_CONDITIONS;
  if (!istrcmp((char*)name,"PRESCRIBED-DISPLACEMENTS"))
    return PRESCRIBED_DISPLACEMENTS;
  if (!istrcmp((char*)name,"PRESC-NODE")) return PRESC_NODE;
  return UNKNOWN_TAG;
}

char *trim_whitespaces(const char* string,size_t size)
{
  const char* end = string+size;
  char* result = (char*)0;
  int not_ws_start = 0;
  int not_ws_end = 0;
  const char* ptr = string;
  /* find starting non-whitespace character */
  while( isspace(*ptr++) && size-- ) not_ws_start++;
  if (size != 0 || not_ws_start == 0)
  {
    ptr--;
    /* find trailing non-whitespace character */
    while(isspace(*--end) && end != ptr) not_ws_end++;
    size = end-ptr+1;
    result = (char*)malloc(size+1);
    memcpy(result,ptr,size);
    result[size] = '\0';
  }
  return result;
}

/*
 * Functions called from expat_start/end_tag_handler
 * when the tag is known
 */
void process_begin_tag(parse_data* data, int tag,const XML_Char **atts);
void process_end_tag(parse_data* data, int tag);

void expat_start_tag_handler(void *userData,
                             const XML_Char *name,
                             const XML_Char **atts)
{
  parse_data* data = (parse_data*)userData;
  int tag = tagname_to_enum(name);
  if(tag != UNKNOWN_TAG)
    process_begin_tag(data,tag,atts);
}

void expat_end_tag_handler(void *userData,
                           const XML_Char *name)
{
  parse_data* data = (parse_data*)userData;
  int tag = tagname_to_enum(name);

  if (tag != UNKNOWN_TAG)
    process_end_tag(data,tag);
  /* clear tag text data at tag close */
  if (data->current_text)
    free(data->current_text);
  data->current_text = (char*)0;
  data->current_size = 0;
}

void expat_text_handler(void *userData,
                        const XML_Char *s,
                        int len)
{
  parse_data* data = (parse_data*)userData;
  char *ptr;
  if (len)
  {
    if (!data->current_text)    /* allocate memory for the text in tag */
    {
      data->current_text = (char*)malloc(len);
      ptr = data->current_text;
    }
    else                        /* reallocate/widen memory alread allocated */
    {
      data->current_text = (char*)realloc(data->current_text,
                                          data->current_size+len);
      ptr = data->current_text;
      ptr += data->current_size;  
    }
    /* append text to the end of allocated/reallocad buffer */
    /* and increase size variable */
    memcpy(ptr,s,len);
    data->current_size += len;
  }
}


BOOL expat_data_load(char *filename,
                     fea_task **task,
                     fea_solution_params **fea_params,
                     nodes_array **nodes,
                     elements_array **elements,
                     presc_boundary_array **presc_boundary)
{
  BOOL result = FALSE;
  FILE* xml_document_file;
  XML_Parser parser;
  size_t file_size = 0;
  size_t read_bytes = 0;
  parse_data parse;
  enum XML_Status status;
  char *file_contents = (char*)0;

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
  XML_SetElementHandler(parser, &expat_start_tag_handler,
                        &expat_end_tag_handler);
  XML_SetCharacterDataHandler(parser,expat_text_handler);

  /* initialize data */
  parse.current_text = (char*)0;
  parse.current_size = 0;

  /* read whole file */
  file_contents = (char*)malloc(file_size);
  read_bytes = fread(file_contents,1,file_size,xml_document_file);
  if (errno)
  {
    free(file_contents);
    return FALSE;
  }
  fclose(xml_document_file);

  /* allocate parse data */
  parse.task = new_fea_task();
  parse.fea_params = new_fea_solution_params();
  parse.nodes = new_nodes_array();
  parse.elements = new_elements_array();
  parse.presc_boundary = new_presc_boundary_array();
  parse.current_size = 0;
  parse.current_text = (char*)0;
  /* set user data */
  XML_SetUserData(parser,&parse);

  /* call parser */
  status = XML_Parse(parser,file_contents,(int)read_bytes,1);
  free(file_contents);
  XML_ParserFree(parser);
  
  *task = parse.task;
  *fea_params = parse.fea_params;
  *nodes = parse.nodes;
  *elements = parse.elements;
  *presc_boundary = parse.presc_boundary;

  
  result = TRUE;
  return result;
}


BOOL check_attribute(const char* attribute_name, const XML_Char ***atts)
{
  BOOL result = FALSE;
  if(!istrcmp((char*)(**atts),(char*)attribute_name))
  {
    (*atts)++;
    result = TRUE;
  }
  return result;
}

void process_model_type(parse_data* data, const XML_Char **atts)
{
  char* text = (char*)0;
  for (; *atts; atts++ )
  {
    if (check_attribute("name",&atts))
    {
      text = trim_whitespaces(*atts,strlen(*atts));
      if (!istrcmp(text,"A5")) 
      {
        data->task->model.model = MODEL_A5;
        data->task->model.parameters_count = 2;
      }
      else if (!istrcmp(text,"COMPRESSIBLE_NEOHOOKEAN"))
      {
        data->task->model.model = MODEL_COMPRESSIBLE_NEOHOOKEAN;
        data->task->model.parameters_count = 2;
      }
      else
      {
        
        printf("unknown model type %s\n",text);
      }
      if(text)
        free(text);
    }
  }
}

/* model-parameters tag handler */
void process_model_params(parse_data* data, const XML_Char **atts)
{
  char* text = (char*)0;
  int count = 0;
  for (; *atts && count < data->task->model.parameters_count; atts++ )
  {
    atts++;
    text = trim_whitespaces(*atts,strlen(*atts));
    data->task->model.parameters[count] = atof(text);
    if(text)
      free(text);
    count++;
  }
}

/* solution tag handler */
void process_solution(parse_data* data, const XML_Char **atts)
{
  char* text = (char*)0;
  for (; *atts; atts++ )
  {
    if (check_attribute("modified-newton",&atts))
    {
      text = trim_whitespaces(*atts,strlen(*atts));
      data->task->modified_newton =
        (!istrcmp(text,"yes") || !istrcmp(text,"true"))? TRUE: FALSE;
      if (text)
        free(text);
    }
    else if (check_attribute("task-type",&atts))
    {
      text = trim_whitespaces(*atts,strlen(*atts));
      if (!istrcmp(text,"CARTESIAN3D"))
        data->task->type = CARTESIAN3D;
      if (text)
        free(text);
    }
    else if (check_attribute("load-increments-count",&atts))
    {
      text = trim_whitespaces(*atts,strlen(*atts));
      data->task->load_increments_count = atoi(text);
      if (text)
        free(text);
    }
    else if (check_attribute("desired-tolerance",&atts))
    {
      text = trim_whitespaces(*atts,strlen(*atts));
      data->task->desired_tolerance = atof(text);
      if (text)
        free(text);
    }
  }
  data->parent_tag = SOLUTION;
}

/* element-type tag handler */
void process_element_type(parse_data* data, const XML_Char **atts)
{
  char* text = (char*)0;
  for (; *atts; atts++ )
  {
    if (check_attribute("name",&atts))
    {
      text = trim_whitespaces(*atts,strlen(*atts));
      if (!istrcmp(text,"TETRAHEDRA10"))
        data->task->ele_type = TETRAHEDRA10;
      if (text)
        free(text);
    }
    else if (check_attribute("nodes-count",&atts))
    {
      text = trim_whitespaces(*atts,strlen(*atts));
      data->fea_params->nodes_per_element = atoi(text);
      if (text)
        free(text);
    }
    else if (check_attribute("nodes-count",&atts))
    {
      text = trim_whitespaces(*atts,strlen(*atts));
      data->fea_params->nodes_per_element = atoi(text);
      if (text)
        free(text);
    }
    else if (check_attribute("gauss-nodes-count",&atts))
    {
      text = trim_whitespaces(*atts,strlen(*atts));
      data->fea_params->gauss_nodes_count = atoi(text);
      if (text)
        free(text);
    }
  }
}

/* line-search tag handler */
void process_line_search(parse_data* data, const XML_Char **atts)
{
  char* text = (char*)0;
  for (; *atts; atts++ )
  {
    if (check_attribute("max",&atts))
    {
      text = trim_whitespaces(*atts,strlen(*atts));
      data->task->linesearch_max = atoi(text);
      if (text)
        free(text);
    }
  }
}

/* arc-length tag handler */
void process_arc_length(parse_data* data, const XML_Char **atts)
{
  char* text = (char*)0;
  for (; *atts; atts++ )
  {
    if (check_attribute("max",&atts))
    {
      text = trim_whitespaces(*atts,strlen(*atts));
      data->task->arclength_max = atoi(text);
      if(text)
        free(text);
    }
  }
}

/* nodes tag handler */
void process_nodes(parse_data* data, const XML_Char **atts)
{
  char* text = (char*)0;
  int i = 0;
  /* set parameters only when 'nodes' tag is a child of the 'geometry' tag */
  if ( data->parent_tag == GEOMETRY )
  {
    for (; *atts; atts++ )
    {
      if (check_attribute("count",&atts))
      {
        text = trim_whitespaces(*atts,strlen(*atts));
        data->nodes->nodes_count = atoi(text);
        /* allocate storage for nodes */
        data->nodes->nodes =
          (real**)malloc(data->nodes->nodes_count*sizeof(real*));
        for (; i < data->nodes->nodes_count; ++ i)
          data->nodes->nodes[i] = (real*)malloc(MAX_DOF*sizeof(real));
        if (text)
          free(text);
      }
    }
    /* set parent tag to 'nodes' to recoginze an appropriate 'node' tag */
    data->parent_tag = NODES;
  }
}

/* node tag handler */
void process_node(parse_data* data, const XML_Char **atts)
{
  char* text = (char*)0;
  real dofs[MAX_DOF];
  int id = -1;
  if ( data->parent_tag == NODES )
  {
    for (; *atts; atts++ )
    {
      if (check_attribute("id",&atts))
      {
        text = trim_whitespaces(*atts,strlen(*atts));
        id = atoi(text);
        if (text) free(text);
      }
      else if (check_attribute("x",&atts))
      {
        text = trim_whitespaces(*atts,strlen(*atts));
        dofs[0] = atof(text);
        if (text) free(text);
      }
      else if (check_attribute("y",&atts))
      {
        text = trim_whitespaces(*atts,strlen(*atts));
        dofs[1] = atof(text);
        if (text) free(text);
      }
      else if (check_attribute("z",&atts))
      {
        text = trim_whitespaces(*atts,strlen(*atts));
        dofs[2] = atof(text);
        if (text) free(text);
      }
    }
    if (id != -1)
      memcpy(data->nodes->nodes[id],dofs,sizeof(dofs));
  }
}

/* elements tag handler */
void process_elements(parse_data* data, const XML_Char **atts)
{
  char* text = (char*)0;
  int i = 0;
  /* set parameters only when 'nodes' tag is a child of the 'geometry' tag */
  if ( data->parent_tag == GEOMETRY )
  {
    for (; *atts; atts++ )
    {
      if (check_attribute("count",&atts))
      {
        text = trim_whitespaces(*atts,strlen(*atts));
        data->elements->elements_count = atoi(text);
        /* allocate storage for elements */
        data->elements->elements = 
          (int**)malloc(data->elements->elements_count*sizeof(int*));
        for (; i < data->elements->elements_count; ++ i)
          data->elements->elements[i] =
            (int*)malloc(data->fea_params->nodes_per_element*sizeof(int));
        if (text) free(text);
      }
    }
    /* set parent tag to 'ELEMENTS' to recoginze an appropriate 'ELEMENT' tag */
    data->parent_tag = ELEMENTS;
  }
}

int node_position_from_attr(const XML_Char ***atts)
{
  int result = -1;
  char* pos = strstr(**atts,"node");
  if (pos == **atts)
  {
    result = atoi(pos+strlen("node"))-1;
    (*atts)++;
  }
  return result;
}

void process_element(parse_data* data, const XML_Char **atts)
{
  char* text = (char*)0;
  int pos = -1;
  int element_size_bytes = data->fea_params->nodes_per_element*sizeof(int);
  int* element = (int*)malloc(element_size_bytes);
  int id = -1;
  if ( data->parent_tag == ELEMENTS )
  {
    for (; *atts; atts++ )
    {
      if (check_attribute("id",&atts))
      {
        text = trim_whitespaces(*atts,strlen(*atts));
        id = atoi(text);
        if (text) free(text);
      }
      if (-1 != (pos = node_position_from_attr(&atts)))
      {
        text = trim_whitespaces(*atts,strlen(*atts));
        element[pos] = atoi(text); 
        if (text) free(text);
      }
    }
    if (id != -1)
      memcpy(data->elements->elements[id],element,element_size_bytes);
  }
  free(element);
}


/* prescribed-displacements tag handler */
void process_prescribed_displacements(parse_data* data, const XML_Char **atts)
{
  char* text = (char*)0;
  int size;
  /* set parameters only when 'nodes' tag is a child of the 'geometry' tag */
  if ( data->parent_tag == BOUNDARY_CONDITIONS )
  {
    for (; *atts; atts++ )
    {
      if (check_attribute("count",&atts))
      {
        text = trim_whitespaces(*atts,strlen(*atts));
        data->presc_boundary->prescribed_nodes_count = atoi(text);
        size = data->presc_boundary->prescribed_nodes_count;
        size = size*sizeof(prescibed_boundary_node);
        /* allocate storage for prescribed nodes */
        data->presc_boundary->prescribed_nodes =
          (prescibed_boundary_node*)malloc(size);
        if (text) free(text);
      }
    }
    /* set parent tag to 'nodes' to recoginze an appropriate 'node' tag */
    data->parent_tag = PRESCRIBED_DISPLACEMENTS;
  }
}

/* node tag handler */
void process_prescribed_node(parse_data* data, const XML_Char **atts)
{
  /* <presc-node id="1" node-id="10" x="0" y="0" z="0" type="7"/> */
  char* text = (char*)0;
  prescibed_boundary_node node;
  int id = -1;
  if ( data->parent_tag == PRESCRIBED_DISPLACEMENTS )
  {
    for (; *atts; atts++ )
    {
      if (check_attribute("id",&atts))
      {
        text = trim_whitespaces(*atts,strlen(*atts));
        id = atoi(text);
        if (text) free(text);
      }
      else if (check_attribute("node-id",&atts))
      {
        text = trim_whitespaces(*atts,strlen(*atts));
        node.node_number = atoi(text);
        if (text) free(text);
      }
      else if (check_attribute("x",&atts))
      {
        text = trim_whitespaces(*atts,strlen(*atts));
        node.values[0] = atof(text);
        if (text) free(text);
      }
      else if (check_attribute("y",&atts))
      {
        text = trim_whitespaces(*atts,strlen(*atts));
        node.values[1] = atof(text);
        if (text) free(text);
      }
      else if (check_attribute("z",&atts))
      {
        text = trim_whitespaces(*atts,strlen(*atts));
        node.values[2] = atof(text);
        if (text) free(text);
      }
      else if (check_attribute("type",&atts))
      {
        text = trim_whitespaces(*atts,strlen(*atts));
        /* TODO: add proper conversion */
        node.type= (presc_boundary_type)atoi(text);
        if (text) free(text);
      }
    }
    if (id != -1)
      memcpy(&data->presc_boundary->prescribed_nodes[id],&node,sizeof(node));
  }
}

void process_begin_tag(parse_data* data, int tag,const XML_Char **atts)
{
  switch(tag)
  {
  case TASK:
    break;
  case MODEL:
    process_model_type(data,atts);
    break;
  case MODEL_PARAMETERS:
    process_model_params(data,atts);
    break;
  case SOLUTION:
    process_solution(data,atts);
    break;
  case ELEMENT_TYPE:
    process_element_type(data,atts);
    break;
  case LINE_SEARCH:
    process_line_search(data,atts);
    break;
  case ARC_LENGTH:
    process_arc_length(data,atts);
    break;
  case INPUT_DATA:
    data->parent_tag = INPUT_DATA;
    break;
  case GEOMETRY:
    data->parent_tag = GEOMETRY;
    break;
  case NODES:
    process_nodes(data,atts);
    break;
  case NODE:
    process_node(data,atts);
    break;
  case ELEMENTS:
    process_elements(data,atts);
    break;
  case ELEMENT:
    process_element(data,atts);
    break;
  case BOUNDARY_CONDITIONS:
    data->parent_tag = BOUNDARY_CONDITIONS;
    break;
  case PRESCRIBED_DISPLACEMENTS:
    process_prescribed_displacements(data,atts);
    data->parent_tag = PRESCRIBED_DISPLACEMENTS;
    break;
  case PRESC_NODE:
    process_prescribed_node(data,atts);
    break;
  default:
    break;
  };
}

void process_end_tag(parse_data* data, int tag)
{
  switch(tag)
  {
  case NODE:
    break;
  case ELEMENT:
    break;
  case PRESC_NODE:
    break;
  case MODEL:
  case SOLUTION:
  case INPUT_DATA:    
    data->parent_tag = TASK;
    break;
  case MODEL_PARAMETERS:
    data->parent_tag = MODEL;
    break;
  case ELEMENT_TYPE:
  case LINE_SEARCH:
  case ARC_LENGTH:
    data->parent_tag = SOLUTION;
    break;
  case GEOMETRY:
  case BOUNDARY_CONDITIONS:
    data->parent_tag = INPUT_DATA;
    break;
  case NODES:
  case ELEMENTS:
    data->parent_tag = GEOMETRY;
    break;
  case TASK:
  default:
    data->parent_tag = UNKNOWN_TAG;
    break;
  };
}


#endif


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

BOOL test_matrix()
{
  BOOL result = TRUE;
  int i,j;
  /* input data */
  real A[3][3] = {{1, 2, 0}, {2, 0, 3}, {0, 2, 3}};
  real B[3][3] = {{0, 2, 1}, {1, 1, 1}, {3, 2, -1}};
  /* expected results */
  real result_matmul[3][3] = {{2, 4, 3}, {9, 10, -1}, {11, 8, -1}};
  real result_matmul_transp[3][3] = {{2, 4, 3}, {6, 8, 0}, {12, 9, 0}};
  real result_matmul_transp2[3][3] = {{4, 3, 7}, {3, 5, 3}, {7, 5, 1}};
  /* results array */
  real R[3][3];

  /* test A x B */
  matrix_mul3x3(A,B,R);
  for (i = 0; i < 3; ++ i)
    for (j = 0; j < 3; ++ j)
      result &= EQL(R[i][j],result_matmul[i][j]);

  if (result)
  {
    /* test A x B' */
    matrix_transpose_mul3x3(A,B,R);
    for (i = 0; i < 3; ++ i)
      for (j = 0; j < 3; ++ j)
        result &= EQL(R[i][j],result_matmul_transp[i][j]);
  }

  if (result)
  {
    /* test A' x B */
    matrix_transpose2_mul3x3(A,B,R);
    for (i = 0; i < 3; ++ i)
      for (j = 0; j < 3; ++ j)
        result &= EQL(R[i][j],result_matmul_transp2[i][j]);
  }
  printf("test_matrix result: *%s*\n",result ? "pass" : "fail");
  return result;
}

BOOL test_sp_matrix()
{
  BOOL result = TRUE;
  sp_matrix mtx,mtx2,mtx3;
  real b[] = {1, 2, 3, 4, 3, 2, 1};
  real expected[] = {25, 34, 40, 45, 42, 16, 23};
  int size = sizeof(expected)/sizeof(real);
  int i;
  real x[] = {0,0,0,0,0,0,0};
  /*
   * Sparse matrix
   * 9  0  0  3  1  0  1
   * 0  11 2  1  0  0  2
   * 0  1  10 2  0  0  0
   * 0  0  2  9  1  0  0
   * 1  0  0  1  12 0  1
   * 0  0  0  0  0  8  0
   * 2  2  0  0  3  0  8
   */

  init_sp_matrix(&mtx,7,7,5,CRS);

  MTX(&mtx,0,0,9);MTX(&mtx,0,3,3);MTX(&mtx,0,4,1);MTX(&mtx,0,6,1);
  MTX(&mtx,1,1,11);MTX(&mtx,1,2,2);MTX(&mtx,1,3,1);MTX(&mtx,1,6,2);
  MTX(&mtx,2,1,1);MTX(&mtx,2,2,10);MTX(&mtx,2,3,2);
  MTX(&mtx,3,2,2);MTX(&mtx,3,3,9);MTX(&mtx,3,4,1);
  MTX(&mtx,4,0,1);MTX(&mtx,4,3,1);MTX(&mtx,4,4,12);MTX(&mtx,4,6,1);
  MTX(&mtx,5,5,8);
  MTX(&mtx,6,0,2);MTX(&mtx,6,1,2);MTX(&mtx,6,4,3);MTX(&mtx,6,6,8);

  sp_matrix_compress(&mtx);

  /* 1st test - matrix-vector multiplication */
  sp_matrix_mv(&mtx,b,x);
  for (i = 0; i < size; ++ i)
    result &= eql(x[i],expected[i]);
  
  /* 2nd test - conversion btw different storage types */
  if (result)
  {
    sp_matrix_convert(&mtx,&mtx2,CCS);
    sp_matrix_mv(&mtx2,b,x);
    for (i = 0; i < size; ++ i)
    {
      result &= eql(x[i],expected[i]);
    }
  }
  if (result)
  {
    sp_matrix_convert(&mtx2,&mtx3,CRS);
    sp_matrix_mv(&mtx3,b,x);
    for (i = 0; i < size; ++ i)
    {
      result &= eql(x[i],expected[i]);
    }
    free_sp_matrix(&mtx2);
    free_sp_matrix(&mtx3);
  }
  
  free_sp_matrix(&mtx);
  printf("test_sp_matrix result: *%s*\n",result ? "pass" : "fail");
  return result;
}

BOOL test_triangle_solver()
{
  BOOL result = TRUE;
  int i;
  sp_matrix mtx,mtx2;
  real x[5] = {0};
  real x_expected[] = {1,2,-3,5,-7};
  real b[] = {-1, 5, -10, 40, -71};

  /*
   * |-1  0  0  0  0 |   | 1 |   |-1 |
   * | 1  2  0  0  0 |   | 2 |   | 5 |
   * |-1  0  3  0  0 | x |-3 | = |-10|
   * | 0  5  0  6  0 |   | 5 |   | 40|
   * | 0  0 -2  0 11 |   |-7 |   |-71|
   */
  init_sp_matrix(&mtx,5,5,3,CRS);
  MTX(&mtx,0,0,-1);
  MTX(&mtx,1,0,1);MTX(&mtx,1,1,2);
  MTX(&mtx,2,0,-1);MTX(&mtx,2,2,3);
  MTX(&mtx,3,1,5);MTX(&mtx,3,3,6);
  MTX(&mtx,4,2,-2);MTX(&mtx,4,4,11);
  sp_matrix_lower_solve(&mtx,5,b,x);
  for (i = 0; i < 5; ++ i)
    result &= EQL(x_expected[i],x[i]);
  
  if (result)
  {
    sp_matrix_convert(&mtx,&mtx2,CCS);
    memset(x,0,sizeof(real)*5);
    sp_matrix_lower_solve(&mtx,5,b,x);
    for (i = 0; i < 5; ++ i)
      result &= EQL(x_expected[i],x[i]);
    free_sp_matrix(&mtx2);
  }
  
  free_sp_matrix(&mtx);
  
  printf("test_triangle_solver result: *%s*\n",result ? "pass" : "fail");
  return result;
}

BOOL test_cg_solver()
{
  BOOL result = TRUE;
  sp_matrix mtx;
  real v[3] = {0}, x[3] = {0};
  int max_iter = 20000;
  real tolerance = 1e-15;

  /* matrix solver test  */

  /* Test 1: */
  /*
   * | 1 0 -2 |   | 1 |   |-5 |
   * | 0 1  0 | x | 2 | = | 2 | 
   * |-2 0  5 |   | 3 |   |13 |
   */
  memset(x,0,3);
  v[0] = -5;
  v[1] = 2;
  v[2] = 13;
  init_sp_matrix(&mtx,3,3,2,CRS);
  
  MTX(&mtx,0,0,1);MTX(&mtx,0,2,-2);
  MTX(&mtx,1,1,1);
  MTX(&mtx,2,0,-2);MTX(&mtx,2,2,5);

  sp_matrix_compress(&mtx);
  sp_matrix_solve_cg(&mtx,v,v,&max_iter,&tolerance,x);

  result = !( fabs(x[0]-1) > TOLERANCE ||
              fabs(x[1]-2) > TOLERANCE ||
              fabs(x[2]-3) > TOLERANCE);
  free_sp_matrix(&mtx);
  
  printf("test_cg_solver result: *%s*\n",result ? "pass" : "fail");
  return result;
}

BOOL test_ilu()
{
  BOOL result = TRUE;
  sp_matrix mtx;
  sp_matrix_skyline m;
  sp_matrix_skyline_ilu ILU;
  real x_exact[] = {1,2,3,0,3,2,1};
  real x[7];
  real b[7];
  int i;
  /* test data for ILU decomposition test */
  real lu_diag_expected[] = {9.000000,
                             11.000000,
                             9.818182,
                             7.888889,
                             11.823161,
                             8.000000,
                             7.205303};
  real lu_lowertr_expected[] = {0.090909,
                                0.222222,
                                0.090909,
                                0.185185,
                                0.111111,
                                0.084507,
                                0.222222,
                                0.181818,
                                0.234944};
  real lu_uppertr_expected[] = {2.000000,
                                3.000000,
                                1.000000,
                                1.909091,
                                1.000000,
                                0.777778,
                                1.000000,
                                2.000000,
                                0.888889};

  memset(b,0,sizeof(b));
  memset(x,0,sizeof(x));
  
  /* Sparse matrix from Balandin
   * 9  0  0  3  1  0  1
   * 0  11 2  1  0  0  2
   * 0  1  10 2  0  0  0
   * 2  1  2  9  1  0  0
   * 1  0  0  1  12 0  1
   * 0  0  0  0  0  8  0
   * 2  2  0  0  3  0  8
   *
   * Test for
   * 1) skyline format
   * 2) ILU decomposition
   * 3) LU - solvers for ILU decomposition
   */

  init_sp_matrix(&mtx,7,7,5,CRS);

  MTX(&mtx,0,0,9);MTX(&mtx,0,3,3);MTX(&mtx,0,4,1);MTX(&mtx,0,6,1);
  MTX(&mtx,1,1,11);MTX(&mtx,1,2,2);MTX(&mtx,1,3,1);MTX(&mtx,1,6,2);
  MTX(&mtx,2,1,1);MTX(&mtx,2,2,10);MTX(&mtx,2,3,2);
  MTX(&mtx,3,0,2);MTX(&mtx,3,1,1);MTX(&mtx,3,2,2);MTX(&mtx,3,3,9);
  MTX(&mtx,3,4,1);
  MTX(&mtx,4,0,1);MTX(&mtx,4,3,1);MTX(&mtx,4,4,12);MTX(&mtx,4,6,1);
  MTX(&mtx,5,5,8);
  MTX(&mtx,6,0,2);MTX(&mtx,6,1,2);MTX(&mtx,6,4,3);MTX(&mtx,6,6,8);

  sp_matrix_compress(&mtx);
  init_sp_matrix_skyline(&m,&mtx);
  init_copy_sp_matrix_skyline_ilu(&ILU,&m);

  for (i = 0; i <  m.rows_count; ++ i)
    result &= fabs(ILU.ilu_diag[i] - lu_diag_expected[i]) < 1e-5;
  
  if (result)
  {
    for (i = 0; i <  m.tr_nonzeros; ++ i)
      result &= fabs(ILU.ilu_lowertr[i] - lu_lowertr_expected[i]) < 1e-5;
  }
  
  if (result)
  {       
    for (i = 0; i <  m.tr_nonzeros; ++ i)
      result &= fabs(ILU.ilu_uppertr[i] - lu_uppertr_expected[i]) < 1e-5;
  }

  /*
   * test for solving Lx=b
   */
  if (result)
  {
    /* prepare a right-part vector */
    sp_matrix_skyline_ilu_lower_mv(&ILU,x_exact,b);
    
    /* solve for x */
    sp_matrix_skyline_ilu_lower_solve(&ILU,b,x);
    /* test result */
    for ( i = 0; i < m.rows_count; ++ i)
      result &= EQL(x[i],x_exact[i]);
  }
  
  /*
   * test for solving Ux=b
   */
  if (result)
  {
    memset(b,0,sizeof(b));
    memset(x,0,sizeof(x));
    /* prepare a right-part vector */
    sp_matrix_skyline_ilu_upper_mv(&ILU,x_exact,b);
    
    /* solve for x */
    sp_matrix_skyline_ilu_upper_solve(&ILU,b,x);
    /* test result */
    for ( i = 0; i < m.rows_count; ++ i)
      result &= EQL(x[i],x_exact[i]);

    free_sp_matrix(&mtx);
    free_sp_matrix_skyline_ilu(&ILU);
  }
  
  printf("test_ilu result: *%s*\n",result ? "pass" : "fail");
  return result;
}

BOOL test_pcg_ilu_solver()
{
  BOOL result = TRUE;
  sp_matrix mtx;
  sp_matrix_skyline_ilu ilu;
  real v[3] = {0}, x[3] = {0};
  int max_iter = 20000;
  real tolerance = 1e-15;

  /* matrix solver test  */

  /* Test 1: */
  /*
   * | 1 0 -2 |   | 1 |   |-5 |
   * | 0 1  0 | x | 2 | = | 2 | 
   * |-2 0  5 |   | 3 |   |13 |
   */
  memset(x,0,3);
  v[0] = -5;
  v[1] = 2;
  v[2] = 13;
  init_sp_matrix(&mtx,3,3,2,CRS);
  
  MTX(&mtx,0,0,1);MTX(&mtx,0,2,-2);
  MTX(&mtx,1,1,1);
  MTX(&mtx,2,0,-2);MTX(&mtx,2,2,5);

  sp_matrix_compress(&mtx);
  sp_matrix_create_ilu(&mtx,&ilu);

  sp_matrix_solve_pcg_ilu(&mtx,&ilu,v,v,&max_iter,&tolerance,x);

  result = !( fabs(x[0]-1) > TOLERANCE ||
              fabs(x[1]-2) > TOLERANCE ||
              fabs(x[2]-3) > TOLERANCE);

  free_sp_matrix_skyline_ilu(&ilu);
  free_sp_matrix(&mtx);
  
  printf("test_pcg_ilu_solver result: *%s*\n",result ? "pass" : "fail");
  return result;
}

BOOL test_cholesky()
{
  BOOL result = TRUE;
  /* initial matrix */
  /* {90, 6, 4, 46, 29, 0, 26}, */
  /* {6, 127, 34, 22, 7, 0, 38}, */
  /* {4, 34, 108, 40, 2, 0, 4}, */
  /* {46, 22, 40, 96, 24, 0, 6}, */
  /* {29, 7, 2, 24, 155, 0, 37}, */
  /* {0, 0, 0, 0, 0, 64, 0}, */
  /* {26, 38, 4, 6, 37, 0, 70} */
  sp_matrix mtx;
  sp_matrix_skyline m;

  /* expected decomposition */
  real cholesky_expected[7][7] = 
    {{9.48683, 0.632456, 0.421637, 4.84883, 3.05687, 0., 2.74064},
     {0.,11.2517, 2.99807, 1.68271, 0.450304, 0., 3.22323},
     {0., 0., 9.94152,3.31043, -0.0642691, 0., -0.685914},
     {0., 0., 0., 7.66149, 1.12678,0., -1.36292},
     {0., 0., 0., 0., 12.0075, 0., 2.38705},
     {0., 0., 0.,0., 0., 8., 0.},
     {0., 0., 0., 0., 0., 0., 6.6388}};

  /* fill initial matrix */
  init_sp_matrix(&mtx,7,7,5,CRS);

/* {90, 6, 4, 46, 29, 0, 26}, */
  MTX(&mtx,0,0,90);MTX(&mtx,0,1,6);MTX(&mtx,0,2,4);MTX(&mtx,0,3,46);
  MTX(&mtx,0,4,29);MTX(&mtx,0,6,26);
/* {6, 127, 34, 22, 7, 0, 38}, */
  MTX(&mtx,1,0,6);MTX(&mtx,1,1,127);MTX(&mtx,1,2,34);MTX(&mtx,1,3,22);
  MTX(&mtx,1,4,7);MTX(&mtx,1,6,38);
/* {4, 34, 108, 40, 2, 0, 4}, */
  MTX(&mtx,2,0,4);MTX(&mtx,2,1,34);MTX(&mtx,2,2,108);MTX(&mtx,2,3,40);
  MTX(&mtx,2,4,2);MTX(&mtx,2,6,4);
/* {46, 22, 40, 96, 24, 0, 6}, */
  MTX(&mtx,3,0,46);MTX(&mtx,3,1,22);MTX(&mtx,3,2,40);MTX(&mtx,3,3,96);
  MTX(&mtx,3,4,24);MTX(&mtx,3,6,6);
/* {29, 7, 2, 24, 155, 0, 37}, */
  MTX(&mtx,4,0,29);MTX(&mtx,4,1,7);MTX(&mtx,4,2,2);MTX(&mtx,4,3,24);
  MTX(&mtx,4,4,155);MTX(&mtx,4,6,37);
/* {0, 0, 0, 0, 0, 64, 0}, */
  MTX(&mtx,5,5,64);
/* {26, 38, 4, 6, 37, 0, 70} */
  MTX(&mtx,6,0,26);MTX(&mtx,6,1,38);MTX(&mtx,6,2,4);MTX(&mtx,6,3,6);
  MTX(&mtx,6,4,37);MTX(&mtx,6,6,70);

  /* prepare initial matrix for conversion to Skyline format */
  sp_matrix_compress(&mtx);
  /* initialize skyline format from given CRS format */
  init_sp_matrix_skyline(&m,&mtx);


  /* clear matrix */
  free_sp_matrix(&mtx);
  free_sp_matrix_skyline(&m);
  
  printf("test_cholesky result: *%s*\n",result ? "pass" : "fail");
  return result;
}


BOOL do_tests()
{
  return test_matrix() && 
    test_sp_matrix() &&
    test_triangle_solver() &&
    test_cg_solver() &&
    test_ilu() &&
    test_pcg_ilu_solver() &&
    test_cholesky();
}
