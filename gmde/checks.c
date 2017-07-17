/******************************************************************************

    Copyright (C) 2012 Sebastian Pancratz

******************************************************************************/

#include "flint/flint.h"
#include "flint/fmpq_poly.h"

#include "gmde.h"

/*
    Checks whether the matrix C over Zp[[t]] satisfies the equation 
    (d/dt + M) C = 0 modulo t^K.
 */

void gmde_check_soln(const fmpz_poly_mat_t C, long vC, const fmpz_t p, long N, 
                     long K, const mat_t M, const ctx_t FracZt)
{
    long b = M->m, i, j, v;

    fmpz_poly_t r;
    fmpz_poly_mat_t Mn, T;

    fmpz_poly_init(r);
    fmpz_poly_mat_init(Mn, b, b);
    fmpz_poly_mat_init(T, b, b);

    /* Compute denominator r(t) */
    {
        fmpz_poly_t t;

        fmpz_poly_init(t);
        fmpz_poly_set_ui(r, 1);
        for (i = 0; i < M->m; i++)
            for (j = 0; j < M->n; j++)
            {
                fmpz_poly_lcm(t, r, fmpz_poly_q_denref(
                    (fmpz_poly_q_struct *) mat_entry(M, i, j, FracZt)));
                fmpz_poly_swap(r, t);
            }
        fmpz_poly_clear(t);
    }

    /* Compute Mn = r M */
    for (i = 0; i < b; i++)
        for (j = 0; j < b; j++)
        {
            fmpz_poly_div(fmpz_poly_mat_entry(Mn, i, j), r, fmpz_poly_q_denref((fmpz_poly_q_struct *) mat_entry(M, i, j, FracZt)));
            fmpz_poly_mul(fmpz_poly_mat_entry(Mn, i, j), fmpz_poly_mat_entry(Mn, i, j), fmpz_poly_q_numref((fmpz_poly_q_struct *) mat_entry(M, i, j, FracZt)));
        }

    fmpz_poly_mat_mul(T, Mn, C);

    for (i = 0; i < b; i++)
        for (j = 0; j < b; j++)
        {
            fmpz_poly_t t;

            fmpz_poly_init(t);
            fmpz_poly_derivative(t, fmpz_poly_mat_entry(C, i, j));
            fmpz_poly_mul(t, r, t);
            fmpz_poly_add(fmpz_poly_mat_entry(T, i, j), t, fmpz_poly_mat_entry(T, i, j));
            fmpz_poly_truncate(fmpz_poly_mat_entry(T, i, j), K - 1);
            fmpz_poly_clear(t);
        }

    v = LONG_MAX;
    for (i = 0; i < b; i++)
        for (j = 0; j < b; j++)
            if (!fmpz_poly_is_zero(fmpz_poly_mat_entry(T, i, j)))
            {
                fmpz_poly_struct *poly = fmpz_poly_mat_entry(T, i, j);
                long w = _fmpz_vec_ord_p(poly->coeffs, poly->length, p);
                v = FLINT_MIN(v, w);
            }

    printf("(d/dt + M) * C (mod t^%ld)\n", K - 1);
    fmpz_poly_mat_print(T, "t");
    printf("ord_p(.) = %ld\n", v);
    printf("ord_p(.) >= N - ord_p(C) ? %d\n", v >= N - vC);

    fmpz_poly_clear(r);
    fmpz_poly_mat_clear(Mn);
    fmpz_poly_mat_clear(T);
}

