/******************************************************************************

    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include <assert.h>

#include "gmde.h"

void gmde_convert_gmc(fmpq_mat_struct **B, long *lenB, fmpz_poly_t denB, 
                      const mat_t M, const ctx_t ctxM)
{
    const long n = M->m;

    mat_t Mnum;
    ctx_t ctxNum;
    fmpz_poly_t t;
    long i, j, k;

    assert(M->m == M->n);

    ctx_init_fmpz_poly(ctxNum);
    mat_init(Mnum, n, n, ctxNum);
    fmpz_poly_init(t);

    fmpz_poly_set_ui(denB, 1);
    for (i = 0; i < M->m; i++)
        for (j = 0; j < M->n; j++)
        {
            fmpz_poly_lcm(t, denB, fmpz_poly_q_denref(
                (fmpz_poly_q_struct *) mat_entry(M, i, j, ctxM)));
            fmpz_poly_swap(denB, t);
        }

    for (i = 0; i < M->m; i++)
        for (j = 0; j < M->n; j++)
        {
            fmpz_poly_divides(t, denB, fmpz_poly_q_denref(
                (fmpz_poly_q_struct *) mat_entry(M, i, j, ctxM)));
            fmpz_poly_mul((fmpz_poly_struct *) mat_entry(Mnum, i, j, ctxNum), 
                          fmpz_poly_q_numref(
                (fmpz_poly_q_struct *) mat_entry(M, i, j, ctxM)), t);
        }

    /* Let lenB be the maximum length of a polynomial in Mnum */
    *lenB = 0;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            *lenB = FLINT_MAX(*lenB, fmpz_poly_length(
                (fmpz_poly_struct *) mat_entry(Mnum, i, j, ctxNum)));
    
    /* Rewrite Mnum (matrix of polys) as b (list of matrices / Q) */
    *B = malloc(*lenB * sizeof(fmpq_mat_struct));

    if (*B == NULL)
    {
        printf("ERROR (gmde_convert_gmc).  Cannot allocate memory.\n");
        abort();
    }

    for (i = 0; i < *lenB; i++)
        fmpq_mat_init(*B + i, n, n);
    
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
        {
            const fmpz_poly_struct *poly = 
                (fmpz_poly_struct *) mat_entry(Mnum, i, j, ctxNum);

            for (k = 0; k < fmpz_poly_length(poly); k++)
            {
                fmpz_set(fmpq_mat_entry_num(*B + k, i, j), 
                         fmpz_poly_get_coeff_ptr(poly, k));
            }
        }

    mat_clear(Mnum, ctxNum);
    ctx_clear(ctxNum);
    fmpz_poly_clear(t);
}

