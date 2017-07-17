/******************************************************************************

    Copyright (C) 2011, 2012 Sebastian Pancratz

******************************************************************************/

#include "flint/fmpz.h"
#include "flint/fmpz_mat.h"
#include "flint/fmpz_poly.h"
#include "flint/fmpz_poly_mat.h"
#include "flint/padic_mat.h"
#include "flint/qadic.h"

/*
    Computes the reverse characteristic polynomial (cp, n+1) 
    of the n-by-n matrix mat.
 */

static 
void _fmpz_poly_mat_revcharpoly(fmpz_poly_struct *cp, 
                                const fmpz_poly_mat_t mat)
{
    const long n = mat->r;

    if (n == 0)
    {
        fmpz_poly_one(cp + 0);
    }
    else if (n == 1)
    {
        fmpz_poly_neg(cp + 0, fmpz_poly_mat_entry(mat, 0, 0));
        fmpz_poly_one(cp + 1);
    }
    else
    {
        long i, j, k, t;
        fmpz_poly_struct *a, *A, *s;
        fmpz_poly_t w;

        a = flint_malloc(n * n * sizeof(fmpz_poly_struct));
        A = a + (n - 1) * n;

        for (i = 0; i < n * n; i++)
            fmpz_poly_init(a + i);
        fmpz_poly_init(w);

        for (i = 0; i <= n; i++)
            fmpz_poly_zero(cp + i);

        fmpz_poly_neg(cp + 0, fmpz_poly_mat_entry(mat, 0, 0));

        for (t = 1; t < n; t++)
        {
            for (i = 0; i <= t; i++)
            {
                fmpz_poly_set(a + 0 * n + i, fmpz_poly_mat_entry(mat, i, t));
            }

            fmpz_poly_set(A + 0, fmpz_poly_mat_entry(mat, t, t));

            for (k = 1; k < t; k++)
            {
                for (i = 0; i <= t; i++)
                {
                    s = a + k * n + i;
                    fmpz_poly_zero(s);
                    for (j = 0; j <= t; j++)
                    {
                        fmpz_poly_mul(w, fmpz_poly_mat_entry(mat, i, j), 
                                         a + (k - 1) * n + j);
                        fmpz_poly_add(s, s, w);
                    }
                }
                fmpz_poly_set(A + k, a + k * n + t);
            }

            fmpz_poly_zero(A + t);
            for (j = 0; j <= t; j++)
            {
                fmpz_poly_mul(w, fmpz_poly_mat_entry(mat, t, j), 
                                 a + (t - 1) * n + j);
                fmpz_poly_add(A + t, A + t, w);
            }

            for (k = 0; k <= t; k++)
            {
                for (j = 0; j < k; j++)
                {
                    fmpz_poly_mul(w, A + j, cp + (k - j - 1));
                    fmpz_poly_sub(cp + k, cp + k, w);
                }
                fmpz_poly_sub(cp + k, cp + k, A + k);
            }
        }

        /* Shift all coefficients up by one */
        for (i = n; i > 0; i--)
        {
            fmpz_poly_swap(cp + i, cp + (i - 1));
        }
        fmpz_poly_one(cp + 0);

        for (i = 0; i < n * n; i++)
            fmpz_poly_clear(a + i);
        flint_free(a);
        fmpz_poly_clear(w);
    }
}

static 
void _deformation_revcharpoly(fmpz *rop, const fmpz_poly_mat_t op, long v, long n, 
                              long N0, const qadic_ctx_t Qq)
{
    const long a  = qadic_ctx_degree(Qq);
    const long b  = op->r;
    const long hi = (n % 2L == 0) ? (b / 2) : b;
    const fmpz *p = (&Qq->pctx)->p;

    long i, j;
    fmpz_t t, q, pN;
    fmpz_poly_struct *cp;

    fmpz_init(t);
    fmpz_init(q);
    fmpz_init(pN);
    cp = flint_malloc((b + 1) * sizeof(fmpz_poly_struct));
    for (i = 0; i <= b; i++)
        fmpz_poly_init(cp + i);

    fmpz_pow_ui(q, p, a);
    fmpz_pow_ui(pN, p, N0);

    /* Step 1.  Compute exact reverse charpoly *******************************/

    if (v >= 0)
    {
        fmpz_poly_mat_t mat;

        fmpz_poly_mat_init(mat, b, b);
        fmpz_pow_ui(t, p, v);
        fmpz_poly_mat_scalar_mul_fmpz(mat, op, t);
        _fmpz_poly_mat_revcharpoly(cp, mat);

        for (i = 0; i <= hi; i++)
        {
            _fmpz_mod_poly_reduce((cp + i)->coeffs, (cp + i)->length, 
                                  Qq->a, Qq->j, Qq->len, pN);
            _fmpz_poly_set_length(cp + i, FLINT_MIN((cp + i)->length, a));
            _fmpz_poly_normalise(cp + i);
        }
        fmpz_poly_mat_clear(mat);
    }
    else
    {
        _fmpz_poly_mat_revcharpoly(cp, op);

        for (i = 0; i <= hi; i++)
        {
            _fmpz_poly_reduce((cp + i)->coeffs, (cp + i)->length, 
                              Qq->a, Qq->j, Qq->len);
            _fmpz_poly_set_length(cp + i, FLINT_MIN((cp + i)->length, a));
            _fmpz_poly_normalise(cp + i);

            fmpz_pow_ui(t, p, (-v) * i);
            fmpz_poly_scalar_divexact_fmpz(cp + i, cp + i, t);
            fmpz_poly_scalar_mod_fmpz(cp + i, cp + i, pN);
        }
    }

    /* Step 2.  Compute the lower half of the coefficients *******************/

    /* t = p^N - (p^N / 2) */
    fmpz_fdiv_q_ui(t, pN, 2);
    fmpz_sub(t, pN, t);

    for (i = 0; i <= hi; i++)
    {
        if (fmpz_poly_length(cp + i) > 1)
        {
            printf("Exception (deformation_revcharpoly).\n");
            printf("The coefficient of T^{%ld} is not constant modulo %ld^%ld, \n", i, *p, N0);
            printf("it is equal to {"), fmpz_poly_print_pretty(cp + i, "X"), printf(".\n");
            abort();
        }

        fmpz_poly_get_coeff_fmpz(rop + i, cp + i, 0);
        fmpz_mod(rop + i, rop + i, pN);
        if (fmpz_cmp(rop + i, t) >= 0)
        {
            fmpz_sub(rop + i, rop + i, pN);
        }
    }

    if (n % 2L == 0)
    {
        const int sgn = 1;

        /* Step 3.  Use a{i} = sgn * a{b-i} q^{(n-1 i - (n-1) b / 2} *********/

        for (i = hi + 1; i <= b; i++)
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
    for (i = 0; i <= b; i++)
        fmpz_poly_clear(cp + i);
    flint_free(cp);
}

static 
void _deformation_revcharpoly_surfaces(
    fmpz *rop, const fmpz_poly_mat_t op, long v, long d, long N0, const qadic_ctx_t Qq)
{
    const long n  = 3;
    const long a  = qadic_ctx_degree(Qq);
    const long b  = op->r;
    const fmpz *p = (&Qq->pctx)->p;

    long i, j;
    fmpz_t h02, s, t, q, pN;
    fmpz_poly_mat_t mat;
    fmpz_poly_struct *cp;

    if (fmpz_cmp_ui(p, 2) == 0)
    {
        printf("Exception (_deformation_revcharpoly_surfaces). \n");
        printf("Assumes that p is not 2, but found p = "), fmpz_print(p), printf(".\n");
        abort();
    }

    if (v < 0)
    {
        printf("Exception (_deformation_revcharpoly_surfaces). \n");
        printf("Assumes that v >= 0, but found v = %ld.\n", v);
        abort();
    }

    fmpz_init(h02);
    fmpz_init(s);
    fmpz_init(t);
    fmpz_init(q);
    fmpz_init(pN);
    fmpz_poly_mat_init(mat, b, b);
    cp = flint_malloc((b + 1) * sizeof(fmpz_poly_struct));
    for (i = 0; i <= b; i++)
        fmpz_poly_init(cp + i);

    fmpz_pow_ui(q, p, a);
    fmpz_pow_ui(pN, p, N0);

    /* Step 1.  Compute exact reverse charpoly *******************************/

    fmpz_pow_ui(t, p, v);
    fmpz_poly_mat_scalar_mul_fmpz(mat, op, t);
    _fmpz_poly_mat_revcharpoly(cp, mat);

    for (i = 0; i <= b; i++)
    {
        _fmpz_poly_reduce((cp + i)->coeffs, (cp + i)->length, 
                          Qq->a, Qq->j, Qq->len);
        _fmpz_poly_set_length(cp + i, FLINT_MIN((cp + i)->length, a));
        _fmpz_poly_normalise(cp + i);
    }

    /* Step 2.  Compute the coefficients of q^{h_{0,2}} f(T/q) ***************/

    /* t = p^N - (p^N / 2) */
    fmpz_fdiv_q_ui(t, pN, 2);
    fmpz_sub(t, pN, t);

    /* h02 = h_{0,2} */
    fmpz_bin_uiui(h02, d-1, 3);

    for (i = 0; i <= b; i++)
    {
        fmpz_poly_get_coeff_fmpz(rop + i, cp + i, 0);

        fmpz_pow_ui(s, q, FLINT_ABS(*h02 - i));
        if (*h02 > i)
        {
            fmpz_mul(rop + i, rop + i, s);
        }
        else if (*h02 < i)
        {
            if (!fmpz_divisible(rop + i, s))
            {
                printf("Exception (_deformation_revcharpoly_surfaces). \n");
                printf("Divisibility condition is not satisfied.\n");
                abort();
            }

            fmpz_divexact(rop + i, rop + i, s);
        }

        fmpz_mod(rop + i, rop + i, pN);
        if (fmpz_cmp(rop + i, t) >= 0)
        {
            fmpz_sub(rop + i, rop + i, pN);
        }

        if (*h02 > i)
        {
            fmpz_divexact(rop + i, rop + i, s);
        }
        else if (*h02 < i)
        {
            fmpz_mul(rop + i, rop + i, s);
        }
    }

    fmpz_clear(h02);
    fmpz_clear(s);
    fmpz_clear(q);
    fmpz_clear(pN);
    fmpz_poly_mat_clear(mat);
    for (i = 0; i <= b; i++)
        fmpz_poly_clear(cp + i);
    flint_free(cp);
}

/*
    Assumes that the matrix is integral and that its entries lie 
    in the interval $[0,p^{N_0})$.
 */

void deformation_revcharpoly(fmpz_poly_t rop, const fmpz_poly_mat_t op, long v, long n, long d, 
                             long N0, long r, long s, const qadic_ctx_t Qq)
{
    const long b  = op->r;
    const fmpz *p = (&Qq->pctx)->p;

    if (v < - (r + s))
    {
        printf("Exception (deformation_revcharpoly).\n");
        printf("The valuation of the action of q^{-1} F_q is %ld,\n", v);
        printf("less than - (r + s) = %ld.\n", - (r + s));
        abort();
    }

    fmpz_poly_fit_length(rop, b + 1);
    _fmpz_poly_set_length(rop, b + 1);

    if (n == 3 && fmpz_cmp_ui(p, 2) != 0)
        _deformation_revcharpoly_surfaces(rop->coeffs, op, v, d, N0, Qq);
    else
        _deformation_revcharpoly(rop->coeffs, op, v, n, N0, Qq);
}

