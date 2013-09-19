/* -*- Mode: C; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include <assert.h>
#include <math.h>
#include "fea_model.h"


void fea_model_init(fea_model_ptr self, model_type type)
{
  switch(type)
  {
  case MODEL_A5:
    self->stress = fea_model_stress_A5;
    self->ctensor = fea_model_ctensor_A5;
    break;
  case MODEL_COMPRESSIBLE_NEOHOOKEAN:
    self->stress = fea_model_stress_compr_neohookean;
    self->ctensor = fea_model_ctensor_compr_neohookean;
    break;
  default:
    assert(FALSE);
  };
   
}


void fea_model_stress_A5(fea_model_ptr self,
                         real (*F)[MAX_DOF],
                         real (*stress_tensor)[MAX_DOF])
{
  int i,j,k;
  real C[MAX_DOF][MAX_DOF];
  real G[MAX_DOF][MAX_DOF];
  real Sn[MAX_DOF][MAX_DOF];
  real lambda,mu;
  real detF = 0;
  real I1 = 0;

  lambda = self->parameters[0];
  mu = self->parameters[1];
  
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

void fea_model_stress_compr_neohookean(fea_model_ptr self,
                                       real (*F)[MAX_DOF],
                                       real (*S)[MAX_DOF])
{
  int i,j,k;
  real B[MAX_DOF][MAX_DOF];
  real lambda,mu;
  real J = det3x3(F);

  lambda = self->parameters[0];
  mu = self->parameters[1];
  
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


void fea_model_ctensor_A5(fea_model_ptr self,
                          real (*graddef)[MAX_DOF],
                          real (*ctensor)[MAX_DOF][MAX_DOF][MAX_DOF])
{
  int i,j,k,l;
  real lambda,mu;
  real detF = det3x3(graddef);
  lambda = self->parameters[0];
  mu = self->parameters[1];
  for ( i = 0; i < MAX_DOF; ++ i )
    for ( j = 0; j < MAX_DOF; ++ j )
      for ( k = 0; k < MAX_DOF; ++ k )
        for ( l = 0; l < MAX_DOF; ++ l )
          ctensor[i][j][k][l] = 
            ( lambda * DELTA (i, j) * DELTA (k, l)  \
              + mu * DELTA (i, k) * DELTA (j, l)    \
              + mu * DELTA (i, l) * DELTA (j, k) ) / detF;
}

void fea_model_ctensor_compr_neohookean(fea_model_ptr self,
                                        real (*graddef)[MAX_DOF],
                                        real (*ctensor)[MAX_DOF][MAX_DOF][MAX_DOF])
{
  int i,j,k,l;
  real lambda,mu,lambda1,mu1;
  real J = det3x3(graddef);
  lambda = self->parameters[0];
  mu = self->parameters[1];
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
