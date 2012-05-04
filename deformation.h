#ifndef DEFORMATION_H
#define DEFORMATION_H

#include <stdlib.h>
#include <mpir.h>

#include "generics.h"
#include "fmpz.h"
#include "mpoly.h"
#include "mat.h"

#include "gmconnection.h"

/*
    Extracts diagonal fibre from the multivariate polynomial $P$, 
    which is expected to be at $t = 0$.
 */
static __inline__ 
void mpoly_diagonal_fibre(fmpz *a, const mpoly_t P, const ctx_t ctx)
{
    mpoly_iter_t iter;
    mpoly_term t;

    const long d = mpoly_degree(P, -1, ctx);
    const long n = P->n;

    mpq_t x, y;

    mpq_init(x);
    mpq_init(y);

    mpoly_iter_init(iter, P);
    while ((t = mpoly_iter_next(iter)))
    {
        long diag = 1, e, i, var = -1;

        for (i = 0; diag && (i < n); i++)
        {
            e = mon_get_exp(t->key, i);
            if (e)
            {
                if (var >= 0)
                    diag = 0;
                else
                    var = i;
            }
        }

        if (diag)
        {
            fmpz_poly_q_evaluate(y, t->val, x);

            if (mpq_sgn(y))
            {
                fmpz_set_mpz(a + var, mpq_numref(y));
            }
            else
            {
                printf("ERROR (mpoly_diagonal_fibre).  Zero diagonal term.\n");
                abort();
            }
        }
    }
    mpoly_iter_clear(iter);

    mpq_clear(x);
    mpq_clear(y);
}

/*
    Returns a $p$-adic precision bound for the computation of 
    the zeta-function of a smooth projective hypersurface in 
    projective space $\mathbf{P}^n$ of degree $d$, defined over 
    $\mathbf{F}_q$ where $q = p^a$.

    If the reverse characteristic polynomial of $q^{-1} F_q$ is 
    computed to $p$-adic precision at least the return value, 
    then uniquely determines the correct integer polynomial.

    Let $b$ denote the dimension of the cohomology vector space, 
    or equivalently the degree of the characteristic polynomial. 
    Whenever $b$ is even, this function returns a bound such 
    that one can compute the bottom half of the coefficients 
    correctly and then use the functional equation to complete 
    the computation.  Whenever $b$ is odd, the sign in the 
    functional equation is generally unknown and this function 
    returns a bound sufficient to determine all coefficients 
    directly.
 */

static __inline__ 
long deformation_prec_zeta_function(const fmpz_t p, long a, 
    long n, long d)
{
    const long b = gmc_basis_size(n, d);
    long N = 0;

    if (b % 2L == 0)
    {
        if (n == 2 && a == 1)
        {
            if (fmpz_cmp_ui(p, 2) == 0 && d <= 5)
            {
                if (d == 3) 
                    N = 2 + 1;
                else if (d == 4) 
                    N = 4 + 1;
                else  /* d == 5 */
                    N = 5 + 1;
            }
            else if (fmpz_cmp_ui(p, 3) == 0 && d == 3)
            {
                N = 1 + 1;
            }
        }

        if (N == 0)
        {
            fmpz_t t;

            fmpz_init(t);
            fmpz_set_ui(t, 6);
            N = fmpz_flog(t, p) + 1 + ((b / 2) * (n - 1) * a + 1) / 2;
            fmpz_clear(t);
        }
    }
    else
    {
        if (n == 2 && a == 1)
        {
            if (fmpz_cmp_ui(p, 2) == 0)
            {
                if (d == 3)
                    N = 2 + 1;
                else if (d == 4)
                    N = 5 + 1;
            }
            else if (fmpz_cmp_ui(p, 3) && d == 3)
            {
                N = 1 + 1;
            }
        }

        if (N == 0)
        {
            N = 1 + (b * (n - 1) * a) / 2 + 1;
        }
    }

    return N;
}

/*
    Returns a number $N_1$ such that, in order to determine the 
    coefficients of the characteristic polynomial of $p^{-1} F_p$ 
    to $p$-adic precision $N$ it suffices to compute the matrix 
    to precision $N_1$.
 */
static __inline__ 
long deformation_prec_frob(const fmpz_t p, long n, long N0)
{
    long r, s;
    fmpz_t t;

    r = padic_val_fac_ui(n - 1, p);

    fmpz_init_set_ui(t, n - 1);
    s = (n + 1) * fmpz_flog(t, p);
    fmpz_clear(t);

    return N0 + r + s;
}

/*
    Currently returns the heuristic value $m = 1.10 \times p N_1$.

    XXX:

    Assumes that $p$ is small and that $p N_1$ 
    fits into a signed long.
 */

static __inline__ 
long deformation_prec_pole_order(const fmpz_t p, long N1)
{
    return (11 * (*p) * N1) / 10;
}

/*
    Returns the required $t$-adic precision.

    If one has bounds $m_{0}$ and $m_{\infty}$ for the 
    pole orders of the matrix $p^{-1} F_p(t)$ modulo $p^{N_1}$ 
    at finite points and infinity, then we could choose 
    $K_1 = \deg(r) m_{0} + m_{\infty}$.

    Without this more precise information, we simply choose 
    to return $K_1 = (\deg(r) + 1) N_1$.
 */

static __inline__
long deformation_prec_tadic(const fmpz_t p, long degR, long N1)
{
    return (degR + 1) * ((11 * (*p) * N1) / 10);
}

static __inline__ 
long deformation_prec_local_solution(const fmpz_t p, long n, long N1, long K1)
{
    return N1 + (n - 1) * n_flog(K1, *p);
}

static __inline__ 
long deformation_prec_diagfrob(const fmpz_t p, long n, long N1, long K1)
{
    return N1 + (n - 1) * (2 * n_flog(K1, *p) - 1);
}

void frob(const mpoly_t P, const ctx_t ctxFracQt, const fmpz_t p);

void frob_with_precisions(mat_t F, const ctx_t ctxF, 
                          const mpoly_t P, const ctx_t ctxFracQt, 
                          long NWork, long Kfinite, long Kinfinite);

void frob_with_precisions_fmpq(mat_t F, const ctx_t ctxF, 
                               const mpoly_t P, const ctx_t ctxFracQt, 
                               long NWork, long K);

#endif

