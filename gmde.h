/******************************************************************************

    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#ifndef GMDE_H
#define GMDE_H

#include "generics.h"
#include "mat.h"
#include "fmpz_poly.h"
#include "fmpq_mat.h"

void gmde_convert_matrix(mat_t C, fmpz_poly_t u, const ctx_t ctxC, 
                         const mat_t M, const ctx_t ctxM);

void gmde_solve_series(fmpq_mat_struct *C, long N, 
                       const mat_t M, const ctx_t ctxM);

void gmde_solve_convert(mat_t A, const ctx_t ctxA, 
                        const fmpq_mat_struct *C, long N);

#endif

