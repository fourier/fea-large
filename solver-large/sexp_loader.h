/* -*- Mode: C; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */

#include "defines.h"
#include "fea_solver.h"

/* loader for the data from the .sexp file */
BOOL sexp_data_load(char *filename,
                    fea_task **task,
                    fea_solution_params **fea_params,
                    nodes_array **nodes,
                    elements_array **elements,
                    presc_bnd_array **presc_boundary);
