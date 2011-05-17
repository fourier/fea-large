/* -*- Mode: C; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>

#include "sp_matrix.h"
#include "defines.h"

/*
 * SLAE solver constants
 * possibly shall go to the initial data in future
 */
const int MAX_ITER = 10000;
const real TOLERANCE = 1e-10;


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
    memset(mtx_to->storage[i].indexes,0,sizeof(int)*mtx_from->storage[i].width);
    mtx_to->storage[i].values =
      (real*)malloc(sizeof(real)*mtx_from->storage[i].width);
    memset(mtx_to->storage[i].values,0,sizeof(real)*mtx_from->storage[i].width);
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
  real* r = (real*)malloc(self->rows_count*sizeof(real));
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
  r = (real*)malloc(size);
  p = (real*)malloc(size);
  temp = (real*)malloc(size);
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
  r = (real*)malloc(size);
  r1 = (real*)malloc(size);
  p = (real*)malloc(size);
  z = (real*)malloc(size);
  temp = (real*)malloc(size);
  
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
