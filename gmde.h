/******************************************************************************

    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#ifndef GMDE_H
#define GMDE_H

#include "generics.h"
#include "mat.h"
#include "fmpz_poly.h"
#include "fmpq_mat.h"

void gmde_convert_gmc(fmpq_mat_struct **B, long *lenB, fmpz_poly_t denB, 
                      const mat_t M, const ctx_t ctxM);

void gmde_solve_fmpq(fmpq_mat_struct *C, long N, 
                     const mat_t M, const ctx_t ctxM);

void gmde_solve_inv_fmpq(fmpq_mat_struct *C, long N, 
                         const mat_t M, const ctx_t ctxM);

void gmde_convert_soln_fmpq(mat_t A, const ctx_t ctxA, 
                            const fmpq_mat_struct *C, long N);

#endif

