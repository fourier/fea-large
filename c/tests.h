/* -*- Mode: C; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#ifndef __TESTS_H__
#define __TESTS_H__

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


#endif /* __TESTS_H__ */
