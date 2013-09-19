/* -*- Mode: C; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#ifndef _FEA_MODEL_H_
#define _FEA_MODEL_H_

#include "dense_matrix.h"


/*************************************************************/
/* Forward declarations                                      */

typedef struct fea_model fea_model;
typedef fea_model* fea_model_ptr;


/*************************************************************/
/* Function pointers declarations                            */

/*
 * A pointer to the function for calculating Cauchy stresses by given
 * deformation gradient
 */
typedef void (*stress_func_t)(fea_model_ptr self,
                              real (*graddef_tensor)[MAX_DOF],
                              real (*stress_tensor)[MAX_DOF]);
/*
 * A pointer to the function for calculating the C elasticity tensor
 * by given deformation gradient
 */
typedef void (*ctensor_func_t)(fea_model_ptr self,
                               real (*graddef)[MAX_DOF],
                               real (*ctensor)[MAX_DOF][MAX_DOF][MAX_DOF]);


/*************************************************************/
/* Enumerations declarations                                 */

typedef enum {
  MODEL_A5,
  MODEL_COMPRESSIBLE_NEOHOOKEAN
} model_type;


/*************************************************************/
/* Data structures                                           */

struct fea_model {
  model_type model;                         /* model type */
  real parameters[MAX_MATERIAL_PARAMETERS]; /* model material parameters */
  int parameters_count;                     /* number of material params */
  stress_func_t stress; /* a pointer to the function
                         * for calculation of the
                         * model-specific Cauchy
                         * in gauss nodes per element */
  ctensor_func_t ctensor; /* a function pointer to the C elasticity
                           * tensor
                           */
};


/*************************************************************/
/* C'tor/D'tor of model structure                            */
void fea_model_init(fea_model_ptr self, model_type type);


/*************************************************************/
/* Functions particular material models                      */

/*
 * Calculate stress tensor of the material model A5 by given
 * deformation gradient
 */
void fea_model_stress_A5(fea_model_ptr self,
                         real (*F)[MAX_DOF],
                         real (*stress_tensor)[MAX_DOF]);

/*
 * Calculate stress tensor of the Neo-hookean compressible model by given
 * deformation gradient
 */
void fea_model_stress_compr_neohookean(fea_model_ptr self,
                                       real (*F)[MAX_DOF],
                                       real (*S)[MAX_DOF]);
/*
 * Calculate 4th rank tensor C of the elastic material model A5
 * T = C(4)**S
 * S - deformation tensor
 * by given deformation gradient
 */
void fea_model_ctensor_A5(fea_model_ptr self,
                          real (*graddef)[MAX_DOF],
                          real (*ctensor)[MAX_DOF][MAX_DOF][MAX_DOF]);
/*
 * Calculate 4th rank tensor C of the Neo-hookean compressible material model
 * T = C(4)**S
 * S - deformation tensor 
 * by given deformation gradient
 */
void fea_model_ctensor_compr_neohookean(fea_model_ptr self,
                                        real (*graddef)[MAX_DOF],
                                        real (*ctensor)[MAX_DOF][MAX_DOF][MAX_DOF]);




#endif /* _FEA_MODEL_H_ */
