/******************************************************************************

    Copyright (C) 2010, 2012 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include <assert.h>
#include <limits.h>

#include "gmde.h"

void gmde_convert_soln_fmpq(mat_t A, const ctx_t ctxA, 
                            const fmpq_mat_struct *C, long N)
{
    long i, j, k;

    assert(N > 0);
    assert(A->m == C->r && A->n == C->c);

    for (i = 0; i < A->m; i++)
        for (j = 0; j < A->n; j++)
        {
            ctxA->zero(ctxA, mat_entry(A, i, j, ctxA));

            for (k = N - 1; k >= 0; k--)
            {
                if (!fmpq_is_zero(fmpq_mat_entry(C + k, i, j)))
                {
                    fmpq_poly_set_coeff_fmpq(
                        (fmpq_poly_struct *) mat_entry(A, i, j, ctxA), 
                        k, fmpq_mat_entry(C + k, i, j));
                }
            }
        }
}

void gmde_convert_soln(fmpz_poly_mat_t A, long *vA, 
                       const padic_mat_struct *C, long N, const fmpz_t p)
{
    long i, j, k;
    fmpz_t s, t;

    assert(N > 0);
    assert(A->r == padic_mat(C)->r && A->c == padic_mat(C)->c);

    fmpz_init(s);
    fmpz_init(t);

    /* Find valuation */
    *vA = LONG_MAX;
    for (k = 0; k < N; k++)
        *vA = FLINT_MIN(*vA, padic_mat_val(C + k));

    fmpz_poly_mat_zero(A);

    for (k = N - 1; k >= 0; k--)
    {
        if (padic_mat_val(C + k) == *vA)
        {
            for (i = 0; i < A->r; i++)
                for (j = 0; j < A->c; j++)
                    if (!fmpz_is_zero(padic_mat_entry(C + k, i, j)))
                        fmpz_poly_set_coeff_fmpz(fmpz_poly_mat_entry(A, i, j), k, 
                                                 padic_mat_entry(C + k, i, j));
        }
        else
        {
            fmpz_pow_ui(s, p, padic_mat_val(C + k) - *vA);

            for (i = 0; i < A->r; i++)
                for (j = 0; j < A->c; j++)
                    if (!fmpz_is_zero(padic_mat_entry(C + k, i, j)))
                    {
                        fmpz_mul(t, s, padic_mat_entry(C + k, i, j));
                        fmpz_poly_set_coeff_fmpz(fmpz_poly_mat_entry(A, i, j), 
                                                 k, t);
                    }
        }
    }

    fmpz_clear(s);
    fmpz_clear(t);
}

