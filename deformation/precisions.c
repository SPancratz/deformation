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
    then this uniquely determines the correct integer polynomial.

    Let $b$ denote the dimension of the cohomology vector space, 
    or equivalently the degree of the characteristic polynomial. 
    Whenever $b$ is even, this function returns a bound such 
    that one can compute the bottom half of the coefficients 
    correctly and then use the functional equation to complete 
    the computation.  Whenever $b$ is odd, the sign in the 
    functional equation is unknown in general and this function 
    returns a bound sufficient to determine all coefficients 
    directly.
 */

static __inline__ 
long _zeta_function(const fmpz_t p, long a, 
    long n, long d)
{
    const long b = gmc_basis_size(n, d);
    long N = 0;

    if (b % 2L == 0)
    {
        if (n == 2 && a == 1)
        {
            if (*p == 2L && d <= 5)
            {
                if (d == 3) 
                    N = 2 + 1;
                else if (d == 4) 
                    N = 4 + 1;
                else  /* d == 5 */
                    N = 5 + 1;
            }
            else if (*p == 3L && d == 3)
            {
                N = 1 + 1;
            }
        }

        if (N == 0)
        {
            fmpz_t t = {6L};  /* 2 b / (b // 2) <= 6 */

            N = fmpz_clog(t, p) + (a * (b / 2) * (n - 1) + 1) / 2 + 1;
        }
    }
    else
    {
        if (n == 2 && a == 1)
        {
            if (*p == 2L)
            {
                if (d == 3)
                    N = 2 + 1;
                else if (d == 4)
                    N = 5 + 1;
            }
            else if (*p == 3L && d == 3)
            {
                N = 1 + 1;
            }
        }

        if (N == 0)
        {
            N = 1 + (a * b * (n - 1) + 1) / 2 + 1;
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
long _frobq(long *r, long *s, const fmpz_t p, long n, long N0)
{
    *r = padic_val_fac_ui(n - 1, p);

    *s = (n + 1) * n_flog(n - 1, *p);

    return N0 + *r + *s;
}

static __inline__ 
long _frobp(long a, long N1, long r, long s)
{
    return N1 + (a - 1) * (r + s);
}

void deformation_precisions(prec_t *prec, 
                            const fmpz_t p, long a, long n, long d, long degR)
{
    long f;

    prec->N0 = _zeta_function(p, a, n, d);
    prec->N1 = _frobq(&(prec->r), &(prec->s), p, n, prec->N0);
    prec->N2 = _frobp(a, prec->N1, prec->r, prec->s);

    prec->m = (long) (1.1 * (*p * prec->N1));
    prec->K = (degR + 1) * prec->m;

    f = 2 * (prec->r + prec->s) + (n - 1);

    prec->N3   = prec->N2  + (prec->r + prec->s) + f * n_clog((prec->K + *p - 1) / *p, *p);
    prec->N3i  = prec->N2  + (prec->r + prec->s) + f * n_clog(prec->K, *p);
    prec->N3w  = prec->N3  + (f + 1) * n_clog(prec->K, *p);
    prec->N3iw = prec->N3i + (f + 1) * n_clog((prec->K + *p - 1) / *p, *p);
    prec->N4   = prec->N2  + f * (n_clog(prec->K, *p) + n_clog((prec->K + *p - 1) / *p, *p));
}

