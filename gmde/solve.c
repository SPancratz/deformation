/******************************************************************************

    Copyright (C) 2011, 2012 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>

#include "gmde.h"

void gmde_solve(padic_mat_struct *C, long N, const padic_ctx_t pctx, 
                const mat_t M, const ctx_t ctxM)
{
    const long n = M->m;

    padic_mat_struct *B;
    long lenB;

    fmpz_poly_t r;
    fmpz * r0;
    long lenR;

    long i, j;

    /* Initialisation */
    fmpz_poly_init(r);

    /* Express M as B / r */
    gmde_convert_gmc(&B, &lenB, r, pctx, M, ctxM);

    r0   = fmpz_poly_get_coeff_ptr(r, 0);
    lenR = fmpz_poly_length(r);

    /* Solve the differential system iteratively */
    {
        padic_mat_t mat;
        fmpz_t coeff;

        padic_mat_init(mat, n, n);
        fmpz_init(coeff);

        _padic_mat_one(C + 0);

        for (i = 0; i < N - 1; i++)
        {
            padic_mat_zero(C + (i + 1));

            j = FLINT_MAX(0, i - lenB + 1);
            for ( ; j <= i; j++)
            {
                /* C[i+1] = C[i+1] + b[i-j] * C[j]; */
                padic_mat_mul(mat, B + (i - j), C + j, pctx);
                padic_mat_reduce(mat, pctx);
                padic_mat_add(C + (i + 1), C + (i + 1), mat, pctx);
            }

            j = FLINT_MAX(0, i - lenR + 1) + 1;
            for ( ; j <= i; j++)
            {
                /* C[i+1] = C[i+1] + r[i-j+1] * j * C[j]; */
                fmpz_mul_ui(coeff, fmpz_poly_get_coeff_ptr(r, i - j + 1), j);
                padic_mat_scalar_mul_fmpz(mat, C + j, coeff, pctx);
                padic_mat_add(C + (i + 1), C + (i + 1), mat, pctx);
            }

            fmpz_mul_si(coeff, r0, -(i + 1));
            padic_mat_scalar_div_fmpz(C + (i + 1), C + (i + 1), coeff, pctx);
        }

        padic_mat_clear(mat);
        fmpz_clear(coeff);
    }

    /* Clean-up */
    for (i = 0; i < lenB; i++)
        padic_mat_clear(B + i);
    free(B);

    fmpz_poly_clear(r);
}

