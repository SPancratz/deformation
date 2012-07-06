/******************************************************************************

    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#ifndef GMDE_H
#define GMDE_H

#include "generics.h"
#include "mat.h"
#include "fmpz_poly.h"
#include "fmpq_mat.h"
#include "padic_mat.h"

void gmde_convert_gmc_fmpq(fmpq_mat_struct **numM, long *lenM, fmpz_poly_t denM, 
                           const mat_t M, const ctx_t ctxM);

void gmde_convert_gmc(padic_mat_struct **B, long *lenB, fmpz_poly_t denB, 
                      const padic_ctx_t pctx, 
                      const mat_t M, const ctx_t ctxM);

void gmde_solve_fmpq(fmpq_mat_struct *C, long N, 
                     const mat_t M, const ctx_t ctxM);

void gmde_solve(padic_mat_struct *C, long N, const padic_ctx_t pctx, 
                const mat_t M, const ctx_t ctxM);

void gmde_convert_soln_fmpq(mat_t A, const ctx_t ctxA, 
                            const fmpq_mat_struct *C, long N);

void gmde_convert_soln(mat_t A, const ctx_t ctxA, 
                       const padic_mat_struct *C, long N);


/* Checks ********************************************************************/

void gmde_check_soln(const mat_t C, const ctx_t Zpt, long K, 
                     const mat_t M, const ctx_t FracZt);

#endif

