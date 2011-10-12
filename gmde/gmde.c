/******************************************************************************

    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include <assert.h>

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

void gmde_convert_soln(mat_t A, const ctx_t ctxA, 
                       const padic_mat_struct *C, long N)
{
    long i, j, k;
    padic_t t;

    assert(N > 0);
    assert(A->m == padic_mat(C)->r && A->n == padic_mat(C)->c);

    _padic_init(t);

    for (i = 0; i < A->m; i++)
        for (j = 0; j < A->n; j++)
        {
            ctxA->zero(ctxA, mat_entry(A, i, j, ctxA));

            for (k = N - 1; k >= 0; k--)
            {
                if (!fmpz_is_zero(padic_mat_unit(C + k, i, j)))
                {
                    _padic_mat_get_entry_padic(t, C + k, i, j, ctxA->pctx);
                    padic_poly_set_coeff_padic(
                        (padic_poly_struct *) mat_entry(A, i, j, ctxA), k, t, 
                        ctxA->pctx);
                }
            }
        }

    _padic_clear(t);
}

