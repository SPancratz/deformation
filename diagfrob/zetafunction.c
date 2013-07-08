/******************************************************************************

    Copyright (C) 2011, 2012, 2013 Sebastian Pancratz

******************************************************************************/

#include "diagfrob.h"

#include <math.h>

#define DEBUG  1

void _diagfrob_zetafunction(fmpz *chi, long n, long d, const fmpz_t p, long a)
{
    const long b = gmc_basis_size(n, d);

    fmpz *s;
    fmpz_t f, g, pN, sum, lo, hi;
    long i, j, n_chi;

#if DEBUG
    double centre, radius;
#endif

#if DEBUG
    printf("Enter _diagfrob_zetafunction:\n");
    printf("  chi = "), _fmpz_vec_print(chi, b + 1), printf("\n");
    printf("  p   = %ld\n", *p);
    printf("  a   = %ld\n", a);
#endif

    assert(fmpz_is_one(chi + 0));

    s = _fmpz_vec_init(b + 1);
    fmpz_init(f);
    fmpz_init(pN);
    fmpz_init(sum);
    fmpz_init(lo);
    fmpz_init(hi);

    /* Compute s[j] and hence round chi[j] for j = 1,...,b */
    for (j = 1; j <= b; j++)
    {
        fmpz_zero(sum);
        for (i = 1; i <= j - 1; i++)
        {
            fmpz_submul(sum, s + (j - i), chi + i);
        }

        /* The exact value of chi[j] lies in the ball with centre sum/j 
           and radius (b/j) q^{j(n-1)/2}.  The intersection of this ball 
           with the integers is equal to [lo, hi].
         */
        fmpz_pow_ui(f, p, a * j * (n - 1));
        fmpz_mul_ui(f, f, b * b);
        fmpz_sqrt(f, f);
        fmpz_sub(lo, sum, f);
        fmpz_cdiv_q_ui(lo, lo, j);
        fmpz_add(hi, sum, f);
        fmpz_fdiv_q_ui(hi, hi, j);

        n_chi = diagfrob_prec_chi(n, b, p, a, j);
        fmpz_pow_ui(pN, p, n_chi);

        fmpz_sub(f, chi + j, lo);
        fmpz_mod(f, f, pN);
        fmpz_add(chi + j, f, lo);

        /* Now s[j] + j chi[j] = sum */
        fmpz_mul_ui(f, chi + j, j);
        fmpz_sub(s + j, sum, f);

#if DEBUG
    centre = fmpz_get_d(sum) / (double)j;
    radius = (double)b / (double)j * pow(*p, (a * j * (n - 1)) / 2.0);

    printf("    j             = %ld\n", j);
    printf("    centre of B/j = %f\n", centre);
    printf("    radius of B/j = %f\n", radius);
    printf("    lo (real)     = %f\n", centre - radius);
    printf("    hi (real)     = %f\n", centre + radius);
    printf("    lo (int)      = "), fmpz_print(lo), printf("\n");
    printf("    hi (int)      = "), fmpz_print(hi), printf("\n");
    printf("    hi - lo       = "), fmpz_sub(f, hi, lo), fmpz_print(f), printf("\n");
    printf("    n_chi         = %ld\n", n_chi);
    printf("    p^{n_chi}     = "), fmpz_print(pN), printf("\n");
    printf("    chi[j]        = "), fmpz_print(chi + j), printf("\n");
    printf("    s[j]          = "), fmpz_print(s + j), printf("\n");
#endif
    }

    _fmpz_vec_clear(s, b + 1);
    fmpz_clear(f);
    fmpz_clear(pN);
    fmpz_clear(sum);
    fmpz_clear(lo);
    fmpz_clear(hi);
}

void diagfrob_zetafunction(fmpz_poly_t chi, 
                           long n, long d, const fmpz_t p, long a)
{
    _diagfrob_zetafunction(chi->coeffs, n, d, p, a);
}

