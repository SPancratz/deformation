/******************************************************************************

    Copyright (C) 2010, 2011, 2012, 2013 Sebastian Pancratz

******************************************************************************/

#ifndef DIAGFROB_H
#define DIAGFROB_H

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "flint.h"
#include "nmod_mat.h"
#include "padic.h"
#include "padic_mat.h"

#include "gmconnection.h"
#include "mon.h"

static __inline__
long diagfrob_delta(long n, const fmpz_t p)
{
    if (fmpz_cmp_ui(p, n) < 0)
    {
        long delta, i;

        delta = padic_val_fac_ui(n - 1, p);

        for (i = 1; i <= n - 1; i++)
            delta += n_flog(i, *p);

        return delta;
    }
    else
    {
        return 0;
    }
}

static __inline__ 
long diagfrob_k(const long *u, long n, long d)
{
    long i, ud;

    ud = n + 1;
    for (i = 0; i <= n; i++)
        ud += u[i];
    ud /= d;

    return ud;
}

/* 
    Observe that 
    \[
    \floor{\log_p (2 (b/i) \mathfrak{q}^{i (n-1)/2})} + 1
        = \floor{1/2 (\floor{\log_p (4 b^2/i^2)} + a i (n-1))} + 1.
    \]
 */

static __inline__
long diagfrob_prec_chi(long n, long b, const fmpz_t p, long a, long i)
{
    return (n_flog((4*b*b)/(i*i), fmpz_get_ui(p)) + a*i*(n-1)) / 2 + 1;
}

static long * diagfrob_gamma(long n, long d)
{
    long i, k, *gamma;
    mon_t *B;
    long *iB, l, u, lenB;

    gmc_basis_sets(&B, &iB, &lenB, &l, &u, n, d);

    gamma = malloc((lenB + 1) * sizeof(long));

    gamma[0] = 0;

    for (k = l; k <= u; k++)
        for (i = iB[k] + 1; i <= iB[k+1]; i++)
            gamma[i] = gamma[i-1] + (k - 1);

    free(B);
    free(iB);

    return gamma;
}

static __inline__ 
long diagfrob_prec_phi(long n, long d, const fmpz_t p, long a)
{
    const long b     = gmc_basis_size(n, d);
    const long delta = diagfrob_delta(n, p);

    long curr, i, prec = 0;

    if (fmpz_cmp_ui(p, n) < 0)
    {
        for (i = 1; i <= b; i++)
        {
            curr = diagfrob_prec_chi(n, b, p, a, i) + delta;
            prec = FLINT_MAX(prec, curr);
        }
    }
    else  /* p >= n */
    {
        long *gamma = diagfrob_gamma(n, d);

        for (i = 1; i <= b; i++)
        {
            curr = diagfrob_prec_chi(n, b, p, a, i) - a * gamma[i-1];
            prec = FLINT_MAX(prec, curr);
        }
        free(gamma);
    }
    return prec;
}

void diagfrob(padic_mat_t F, const fmpz *a, long n, long d, 
              const padic_ctx_t ctx, int verbose);

void diagfrob_revcharpoly(fmpz_poly_t chi, 
                          const padic_mat_t phi, const padic_ctx_t ctx);

void _diagfrob_zetafunction(fmpz *chi, long n, long d, const fmpz_t p, long a);

void diagfrob_zetafunction(fmpz_poly_t chi, 
                           long n, long d, const fmpz_t p, long a);

static 
int diagfrob_verify_functional_eq(const fmpz_poly_t chi, 
                                  long n, long d, const fmpz_t p, long a)
{
    const long b = gmc_basis_size(n, d);

    if (fmpz_poly_length(chi) != b + 1)
    {
        return 0;
    }
    else if (!fmpz_is_one(fmpz_poly_get_coeff_ptr(chi, 0)))
    {
        return 0;
    }
    else
    {
        const int sgn = fmpz_sgn(fmpz_poly_lead(chi));
        long i;
        fmpz_t q, t;

        fmpz_init(q);
        fmpz_init(t);
        fmpz_pow_ui(q, p, a);

        /* want i <= b - i, hence 2 i <= b hence i <= b/2 hence i <= floor(b/2) */
        for (i = 0; i <= b / 2; i++)
        {
            fmpz_pow_ui(t, q, ((n - 1) * b) / 2 - (n - 1) * i);
            fmpz_mul(t, t, fmpz_poly_get_coeff_ptr(chi, i));
            if (sgn < 0)
                fmpz_neg(t, t);
            if (!fmpz_equal(t, fmpz_poly_get_coeff_ptr(chi, b - i)))
                return 0;
        }

        fmpz_clear(q);
        fmpz_clear(t);

        return 1;
    }
}

void nmod_mat_hessenberg(nmod_mat_t rop, const nmod_mat_t op);

void nmod_mat_charpoly(nmod_poly_t rop, const nmod_mat_t op);

void fmpz_mat_charpoly_modular(fmpz_poly_t rop, const fmpz_mat_t op);

#endif

