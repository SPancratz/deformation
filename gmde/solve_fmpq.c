/******************************************************************************

    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include <assert.h>

#include "gmde.h"

void gmde_solve_fmpq(fmpq_mat_struct *C, long N, 
                     const mat_t M, const ctx_t ctxM)
{
    const long n = M->m;

    fmpq_mat_struct *B;
    long lenB;

    fmpz_poly_t r;
    fmpz * r0;
    long lenR;

    long i, j;

    /* Initialisation */
    fmpz_poly_init(r);

    /* Express M as B / r */
    gmde_convert_gmc(&B, &lenB, r, M, ctxM);

    r0   = fmpz_poly_get_coeff_ptr(r, 0);
    lenR = fmpz_poly_length(r);

    /* Solve the differential system iteratively */
    {
        fmpq_mat_t mat;
        fmpz_t coeff;

        fmpq_mat_init(mat, n, n);
        fmpz_init(coeff);

        fmpq_mat_one(C + 0);

        for (i = 0; i < N - 1; i++)
        {
            fmpq_mat_zero(C + (i + 1));

            j = FLINT_MAX(0, i - lenB + 1);
            for ( ; j <= i; j++)
            {
                /* C[i+1] = C[i+1] + b[i-j] * C[j]; */
                fmpq_mat_mul(mat, B + (i - j), C + j);
                fmpq_mat_add(C + (i + 1), C + (i + 1), mat);
            }

            j = FLINT_MAX(0, i - lenR + 1) + 1;
            for ( ; j <= i; j++)
            {
                /* C[i+1] = C[i+1] + r[i-j+1] * j * C[j]; */
                fmpz_mul_ui(coeff, fmpz_poly_get_coeff_ptr(r, i - j + 1), j);
                fmpq_mat_scalar_mul_fmpz(mat, C + j, coeff);
                fmpq_mat_add(C + (i + 1), C + (i + 1), mat);
            }

            fmpz_mul_si(coeff, r0, -(i + 1));
            fmpq_mat_scalar_div_fmpz(C + (i + 1), C + (i + 1), coeff);
        }

        fmpq_mat_clear(mat);
        fmpz_clear(coeff);
    }

    /* Clean-up */
    for (i = 0; i < lenB; i++)
        fmpq_mat_clear(B + i);
    free(B);

    fmpz_poly_clear(r);
}

