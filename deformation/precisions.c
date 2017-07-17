/******************************************************************************

    Copyright (C) 2012 Sebastian Pancratz

******************************************************************************/

#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/padic.h"

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
    long i, N;

    fmpz_t f, g, max;

    fmpz_init(f);
    fmpz_init(g);
    fmpz_init(max);

    if (n == 3 && fmpz_cmp_ui(p, 2) != 0)
    {
        fmpz_bin_uiui(f, d-1, 3);
        fmpz_bin_uiui(g, b, b / 2);
        fmpz_mul_ui(g, g, 2);
        N = a * (*f) + fmpz_clog(g, p);
    }
    else if (n % 2L == 0)  /* n even implies b even */
    {
        fmpz_bin_uiui(f, b, b / 2);
        fmpz_pow_ui(g, p, (a * (b / 2) * (n - 1) + 1) / 2);
        fmpz_mul(f, f, g);
        fmpz_mul_ui(f, f, 2);

        N = fmpz_flog(f, p) + 1;
    }
    else
    {
        for (i = b / 2; i <= b; i++)
        {
            fmpz_bin_uiui(f, b, i);
            fmpz_pow_ui(g, p, (a * i * (n - 1) + 1) / 2);
            fmpz_mul(f, f, g);
            fmpz_mul_ui(f, f, 2);

            if (fmpz_cmp(max, f) < 0)
                fmpz_swap(max, f);
        }

        N = fmpz_flog(max, p) + 1;
    }

    fmpz_clear(f);
    fmpz_clear(g);
    fmpz_clear(max);

    return N;
}

/*
    Returns data $N_1 = N_0 + r + s$ such that, in order to determine the 
    coefficients of the characteristic polynomial of $p^{-1} F_p$ to $p$-adic 
    precision $N_0$ it suffices to compute the matrix to precision $N_1$.

    When $X$ is a smooth projective surface and $p > 2$, we can do better 
    by instead computing $q^{h_{0,2}} f(T/q)$, in which case the extraction 
    of the polynomial requires an extra $a$ digits of precision.
 */
static __inline__ 
long _frobq(long *r, long *s, const fmpz_t p, long a, long n, long d, long N0)
{
    if (n == 3 && fmpz_cmp_ui(p, 2) != 0)
    {
        *r = 0;
        *s = 0;

        return N0 + a;
    }
    else
    {
        *r = padic_val_fac_ui(n - 1, p);
        *s = (n + 1) * n_flog(n - 1, *p);

        return N0 + *r + *s;
    }
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
    prec->N1 = _frobq(&(prec->r), &(prec->s), p, a, n, d, prec->N0);
    prec->N2 = _frobp(a, prec->N1, prec->r, prec->s);

    prec->m = (long) (1.1 * (*p * prec->N1));
    prec->K = (degR + 1) * prec->m;

    f = 2 * (prec->r + prec->s) + (n - 1);

    prec->N3   = prec->N2  + (prec->r + prec->s) + f * n_clog((prec->K + *p - 1) / *p, *p);
    prec->N3i  = prec->N2  + (prec->r + prec->s) + f * n_clog(prec->K, *p);
    prec->N3w  = prec->N3  + (f + 1) * n_clog(prec->K, *p);
    prec->N3iw = prec->N3i + (f + 1) * n_clog((prec->K + *p - 1) / *p, *p);
    prec->N4   = prec->N2  + f * (n_clog(prec->K, *p) + n_clog((prec->K + *p - 1) / *p, *p));

    prec->denR = NULL;
}

