/******************************************************************************

    Copyright (C) 2011, 2012 Sebastian Pancratz

******************************************************************************/

#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_mat.h"
#include "padic_mat.h"

/*
    TODO:  Currently only works with non-negative valuations.
 */

static 
void _deformation_revcharpoly(fmpz *rop, const padic_mat_t op, long n, 
                              const fmpz_t p, long a, long N0)
{
    const long b  = padic_mat(op)->r;
    const long v  = padic_mat_val(op);
    const long hi = (b % 2L == 0) ? (b + 2) / 2 : b;

    long i, j;
    fmpz_t t, q, pN;
    fmpz_mat_t mat;

    fmpz_init(t);
    fmpz_init(q);
    fmpz_init(pN);
    fmpz_mat_init(mat, b, b);

    fmpz_pow_ui(q, p, a);
    fmpz_pow_ui(pN, p, N0);

    /* Step 1.  Compute exact reverse charpoly *******************************/

    fmpz_pow_ui(t, p, v);
    for (i = 0; i < b; i++)
        for (j = 0; j < b; j++)
            fmpz_mul(fmpz_mat_entry(mat, i, j), padic_mat_unit(op, i, j), t);

    _fmpz_mat_charpoly(rop, mat);
    _fmpz_poly_reverse(rop, rop, b + 1, b + 1);

    /* Step 2.  Compute the lower half of the coefficients *******************/

    /* t = p^N - (p^N / 2) */
    fmpz_fdiv_q_ui(t, pN, 2);
    fmpz_sub(t, pN, t);

    for (i = 0; i <= hi; i++)
    {
        fmpz_mod(rop + i, rop + i, pN);
        if (fmpz_cmp(rop + i, t) >= 0)
            fmpz_sub(rop + i, rop + i, pN);
    }

    if (b % 2L == 0)
    {
        const int sgn = 1;

        /* Step 3.  Use a{i} = sgn * a{b-i} q^{(n-1 i - (n-1) b / 2} *********/

        for (i = (b + 2) / 2 + 1; i <= b; i++)
        {
            fmpz_pow_ui(t, q, (n - 1) * i - ((n - 1) * b) / 2);
            fmpz_mul(rop + i, rop + (b - i), t);
            if (sgn == -1)
                fmpz_neg(rop + i, rop + i);
        }
    }

    fmpz_clear(t);
    fmpz_clear(q);
    fmpz_clear(pN);
    fmpz_mat_clear(mat);
}

/*
    Assumes that the matrix is integral and that its entries lie 
    in the interval $[0,p^{N_0})$.
 */

void deformation_revcharpoly(fmpz_poly_t rop, const padic_mat_t op, long n, 
                             const fmpz_t p, long a, long N0, long r, long s)
{
    const long b = padic_mat(op)->r;
    const long v = padic_mat_val(op);

    if (v < - (r + s))
    {
        printf("Exception (deformation_revcharpoly).\n");
        printf("The valuation of the action of F_q is %ld, less than\n", v);
        printf("- (r + s) = %ld.\n", - (r + s));
        abort();
    }

    if (v < 0)
    {
        printf("Exception (deformation_revcharpoly).\n");
        printf("Case of negative valuation not implemented yet.\n");
        abort();
    }

    fmpz_poly_fit_length(rop, b + 1);
    _fmpz_poly_set_length(rop, b + 1);

    _deformation_revcharpoly(rop->coeffs, op, n, p, a, N0);
}

