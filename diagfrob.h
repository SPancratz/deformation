/******************************************************************************

    Copyright (C) 2010, 2011 Sebastian Pancratz

******************************************************************************/

#ifndef DIAGFROB_H
#define DIAGFROB_H

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "flint.h"
#include "padic.h"

#include "generics.h"
#include "gmconnection.h"
#include "mat.h"
#include "mon.h"

#define DIAGFROB_MOD(x, m)  (((x) % (m) >= 0) ? (x) % (m) : (x) % (m) + (m))

static __inline__ 
long fdiv_si(long n, long d)
{
    return ((n < 0) && (n % d) ? n / d - 1 : n / d);
}

void diagfrob_falling_fac_mpq(mpq_t rop, long u, long d, long r);

void diagfrob_coefficient_mpq(mpq_t rop, long m, long p);

void diagfrob_coefficient(padic_t rop, long m, const padic_ctx_t ctx);

void diagfrob_alpha_mpq(mpq_t rop, const fmpz *a, long n, long d, 
                        const mon_t u, const mon_t v, 
                        const padic_ctx_t ctx);

void diagfrob_alpha(padic_t rop, const fmpz *a, long n, long d, 
                    const mon_t u, const mon_t v, 
                    const padic_ctx_t ctx);

void diagfrob_entry_mpq(mpq_t rop, const fmpz *a, long n, long d, 
                        const mon_t u, const mon_t v, long bound, 
                        const padic_ctx_t ctx);

void diagfrob_entry(padic_t rop, const fmpz *a, long n, long d, 
                    const mon_t u, const mon_t v, long bound, 
                    const padic_ctx_t ctx);

void diagfrob(mat_t F, const fmpz *a, long n, long d, const ctx_t ctx);

void diagfrob_revcharpoly(fmpz *poly, const char *f, long n, long b, 
                          const ctx_t ctx);

static __inline__
long fmpz_clog(const fmpz_t op, const fmpz_t b)
{
    assert((fmpz_cmp_ui(op, 1) >= 0) && (fmpz_cmp_ui(b, 2) >= 0));

    if (fmpz_is_one(op))
    {
        return 0;
    }
    else
    {
        long n;
        fmpz_t x;

        fmpz_init(x);
        fmpz_set(x, b);
        for (n = 1; fmpz_cmp(x, op) < 0; n++)
        {
            fmpz_mul(x, x, b);
        }
        fmpz_clear(x);

        return n;
    }
}

static __inline__
long fmpz_flog(const fmpz_t op, const fmpz_t b)
{
    assert((fmpz_cmp_ui(op, 1) >= 0) && (fmpz_cmp_ui(b, 2) >= 0));

    if (fmpz_is_one(op))
    {
        return 0;
    }
    else
    {
        long n;
        fmpz_t x;

        fmpz_init(x);
        fmpz_set(x, b);
        for (n = 1; fmpz_cmp(x, op) <= 0; n++)
        {
            fmpz_mul(x, x, b);
        }
        fmpz_clear(x);

        return n - 1;
    }
}

/*
    Computes the precision to which one needs to compute the 
    bottom half of the coefficients.  We relax the bound from 
    [Ger07, Theorem 3.2] to 
    \begin{equation*}
    \ceil{\log_p \binom{b}{\floor{b/2}}} + b \ceil{(n-1)/4} a + 2.
    \end{equation*}

    Here, $n, d, p$ are as in the statement of the theorem, 
    and $a$ is such that $q = p^a$.

    [Ger07] - ``Relative rigid cohomology and deformation of 
    hypersurfaces``
  */
static __inline__ 
long diagfrob_charpoly_prec(long n, long d, const fmpz_t p, long a)
{
    fmpz_t x;
    long logx;
    long b = gmc_basis_size(n, d);

    fmpz_init(x);

    fmpz_bin_uiui(x, b, b / 2);
    logx = fmpz_clog(x, p);

    fmpz_clear(x);

    return logx + b * ((n - 1 + 3) / 4) * a + 2;
}

/*
    Sets $r$ and $s$ to non-negative integers suitable for the 
    application of [Ger07, Lemma 3.4].

    Thus, to determine the coefficients of the characteristic 
    polynomial to $p$-padic precision $N$ it suffices to 
    compute a $p$-adic approximation to the matrix $q^{-1} F$ 
    to precision $N + (r + s)$.
 */
static __inline__ 
void diagfrob_matrix_prec(long *r, long *s, long n, const fmpz_t p)
{
    fmpz_t x;

    fmpz_init(x);

    fmpz_fac_ui(x, n - 1);
    *r = fmpz_remove(x, x, p);

    fmpz_set_ui(x, n - 1);
    *s = fmpz_flog(x, p);
    *s = (n + 1) * *s;

    fmpz_clear(x);
}

#endif

