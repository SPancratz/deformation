/******************************************************************************

    Copyright (C) 2011, 2012, 2013 Sebastian Pancratz

******************************************************************************/

#include "diagfrob.h"

void _diagfrob_zetafunction(fmpz *z, const fmpz *chi, long n, long d, 
                            const fmpz_t p, long a)
{
    const long b = gmc_basis_size(n, d);

    fmpz *s;
    fmpz_t f, pN, pN_2, sum;
    long i, j, n_chi;

    s = _fmpz_vec_init(b + 1);
    fmpz_init(f);
    fmpz_init(pN);
    fmpz_init(pN_2);
    fmpz_init(sum);

    fmpz_one(z + 0);

    /* Compute z[j] for j = 1,...,b */
    for (j = 1; j <= b; j++)
    {
        fmpz_zero(sum);
        for (i = 1; i <= j - 1; i++)
        {
            fmpz_mul(f, s + (j - i), z + i);
            fmpz_add(sum, sum, f);
        }
        fmpz_neg(sum, sum);

        /* Now s[j] + j z[j] = sum */
        fmpz_mul_ui(f, chi + j, j);
        fmpz_sub(f, sum, f);

        /* Now f is an approximation to s[j] */
        n_chi = diagfrob_prec_chi(n, b, p, a, j);

        fmpz_pow_ui(pN, p, n_chi);
        fmpz_fdiv_q_ui(pN_2, pN, 2);
        fmpz_sub(pN_2, pN, pN_2);

        fmpz_mod(s + j, f, pN);
        if (fmpz_cmp(s + j, pN_2) >= 0)
            fmpz_sub(s + j, s + j, pN);

        /* Now f is equal to s[j] and we can recover z[j] */
        fmpz_sub(f, sum, s + j);

        assert(fmpz_divisible_si(f, j));

        fmpz_divexact_si(z + j, f, j);
    }

    _fmpz_vec_clear(s, b + 1);
    fmpz_clear(f);
    fmpz_clear(pN);
    fmpz_clear(pN_2);
    fmpz_clear(sum);
}

void diagfrob_zetafunction(fmpz_poly_t z, 
                           const fmpz_poly_t chi, long n, long d, 
                           const fmpz_t p, long a)
{
    const long b = gmc_basis_size(n, d);

    fmpz_poly_fit_length(z, b + 1);

    _diagfrob_zetafunction(z->coeffs, chi->coeffs, n, d, p, a);

    _fmpz_poly_set_length(z, b + 1);
    _fmpz_poly_normalise(z);
}

