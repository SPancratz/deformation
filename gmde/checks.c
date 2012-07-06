/******************************************************************************

    Copyright (C) 2012 Sebastian Pancratz

******************************************************************************/

#include "flint.h"
#include "fmpq_poly.h"

#include "gmde.h"

/*
    Checks whether the matrix C over Zp[[t]] satisfies (d/dt + M) C = 0.
 */

void gmde_check_soln(const mat_t C, const ctx_t Zpt, long K, 
                     const mat_t M, const ctx_t FracZt)
{
    long b = M->m, i, j, k, v;
    mat_t T;
    padic_ctx_struct *pctx = Zpt->pctx;

    mat_init(T, b, b, Zpt);

    for (i = 0; i < b; i++)
        for (j = 0; j < b; j++)
        {
            padic_poly_derivative((padic_poly_struct *) mat_entry(T, i, j, Zpt), 
                                  (padic_poly_struct *) mat_entry(C, i, j, Zpt), pctx);

            for (k = 0; k < b; k++)
            {
                fmpq_poly_t t1, t2;
                padic_poly_t t;

                fmpq_poly_init(t1);
                fmpq_poly_init(t2);
                padic_poly_init(t);

                fmpq_poly_set_fmpz_poly(t1, fmpz_poly_q_numref(
                    (fmpz_poly_q_struct *) mat_entry(M, i, k, FracZt)));
                fmpq_poly_set_fmpz_poly(t2, fmpz_poly_q_denref(
                    (fmpz_poly_q_struct *) mat_entry(M, i, k, FracZt)));

                fmpq_poly_inv_series(t2, t2, K);
                fmpq_poly_mul(t1, t1, t2);

                padic_poly_set_fmpq_poly(t, t1, pctx);
                padic_poly_mul(t, t, 
                    (padic_poly_struct *) mat_entry(C, k, j, Zpt), pctx);
                padic_poly_add((padic_poly_struct *) mat_entry(T, i, j, Zpt), 
                               (padic_poly_struct *) mat_entry(T, i, j, Zpt), 
                               t, pctx);

                padic_poly_truncate(
                    (padic_poly_struct *) mat_entry(T, i, j, Zpt), K - 1, 
                    pctx->p);

                fmpq_poly_clear(t1);
                fmpq_poly_clear(t2);
                padic_poly_clear(t);
            }
        }

    v = LONG_MAX;
    for (i = 0; i < b; i++)
        for (j = 0; j < b; j++)
        {
            const padic_poly_struct *poly = 
                (padic_poly_struct *) mat_entry(T, i, j, Zpt);
            if (!padic_poly_is_zero(poly))
                v = FLINT_MIN(padic_poly_val(poly), v);
        }

    printf("(d/dt + M) * C (mod t^%ld)\n", K - 1);
    mat_print(T, Zpt);
    printf("\n");
    printf("ord_p(...) = %ld\n", v);

    mat_clear(T, Zpt);
}

