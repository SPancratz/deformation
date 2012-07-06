/******************************************************************************

    Copyright (C) 2012 Sebastian Pancratz

******************************************************************************/

#include "flint.h"
#include "fmpz.h"
#include "padic.h"

#include "deformation.h"
#include "gmconnection.h"

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
    Returns data $N_1 = N_0 + r + s$ such that, in order to determine the 
    coefficients of the characteristic polynomial of $p^{-1} F_p$ to $p$-adic 
    precision $N_0$ it suffices to compute the matrix to precision $N_1$.
 */
static __inline__ 
void deformation_prec_frob(long *r, long *s, const fmpz_t p, long n, long N0)
{
    fmpz_t t;

    *r = padic_val_fac_ui(n - 1, p);

    fmpz_init_set_ui(t, n - 1);
    *s = (n + 1) * fmpz_flog(t, p);
    fmpz_clear(t);
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
    to return $K_1 = 1.1 (\deg(r) + 1) p N_1$.
 */

static __inline__
long deformation_prec_tadic(const fmpz_t p, long degR, long N1)
{
    long r;

    r = (degR + 1) * (*p) * N1;

    return (11 * r + 9) / 10;
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

void deformation_precisions(prec_struct *prec, 
                            const fmpz_t p, long a, long n, long d, long degR)
{
    prec->N0 = deformation_prec_zeta_function(p, a, n, d);
    deformation_prec_frob(&(prec->r), &(prec->s), p, n, prec->N0);
    prec->N1 = prec->N0 + prec->r + prec->s;
    prec->m  = deformation_prec_pole_order(p, prec->N1);
    prec->K  = deformation_prec_tadic(p, degR, prec->N1);
    prec->N2 = deformation_prec_local_solution(p, n, prec->N1, prec->K);
    prec->N3 = deformation_prec_diagfrob(p, n, prec->N1, prec->K);
}

