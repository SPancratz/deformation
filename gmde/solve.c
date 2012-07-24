/******************************************************************************

    Copyright (C) 2011, 2012 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>

#include "gmde.h"

void gmde_solve(padic_mat_struct *C, long K, const fmpz_t p, long N, long Nw, 
                const mat_t M, const ctx_t ctxM)
{
    const long n = M->m;

    padic_ctx_t pctx, pctxW;

    padic_mat_struct *B;
    long lenB;

    fmpz_poly_t r;
    fmpz * r0;
    long lenR;

    long i, j;

    /* Initialisation */
    fmpz_poly_init(r);
    padic_ctx_init(pctx,  p,  N, PADIC_SERIES);
    padic_ctx_init(pctxW, p, Nw, PADIC_SERIES);

    /* Express M as B / r */
    gmde_convert_gmc(&B, &lenB, r, pctxW, M, ctxM);

    r0   = fmpz_poly_get_coeff_ptr(r, 0);
    lenR = fmpz_poly_length(r);

    /* Solve the differential system iteratively */
    {
        padic_mat_t mat;
        fmpz_t coeff;

        padic_mat_init(mat, n, n);
        fmpz_init(coeff);

        _padic_mat_one(C + 0);

        for (i = 0; i < K - 1; i++)
        {
            padic_mat_zero(C + (i + 1));

            j = FLINT_MAX(0, i - lenB + 1);
            for ( ; j <= i; j++)
            {
                /* C[i+1] = C[i+1] + b[i-j] * C[j]; */
                _padic_mat_mul(mat, B + (i - j), C + j, pctxW);
                _padic_mat_add(C + (i + 1), C + (i + 1), mat, pctxW);
            }

            j = FLINT_MAX(0, i - lenR + 1) + 1;
            for ( ; j <= i; j++)
            {
                /* C[i+1] = C[i+1] + r[i-j+1] * j * C[j]; */
                fmpz_mul_ui(coeff, fmpz_poly_get_coeff_ptr(r, i - j + 1), j);
                _padic_mat_scalar_mul_fmpz(mat, C + j, coeff, pctxW);
                _padic_mat_add(C + (i + 1), C + (i + 1), mat, pctxW);
            }

            fmpz_mul_si(coeff, r0, -(i + 1));
            padic_mat_scalar_div_fmpz(C + (i + 1), C + (i + 1), coeff, pctxW);
        }

        padic_mat_clear(mat);
        fmpz_clear(coeff);
    }

    for (i = 0; i < K; i++)
    {
        padic_mat_reduce(C + i, pctx);
    }

    /* Clean-up */
    for (i = 0; i < lenB; i++)
        padic_mat_clear(B + i);
    free(B);

    fmpz_poly_clear(r);
    padic_ctx_clear(pctx);
    padic_ctx_clear(pctxW);
}

