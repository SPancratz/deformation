/******************************************************************************

    Copyright (C) 2011, 2012, 2013 Sebastian Pancratz

******************************************************************************/

#include "diagfrob.h"

void _diagfrob_zetafunction(fmpz *chi, long n, long d, const fmpz_t p, long a)
{
    const long b = gmc_basis_size(n, d);

    fmpz *s;
    fmpz_t f, pN, pN_2, sum;
    long i, j, n_chi;

    assert(fmpz_is_one(chi + 0));

    s = _fmpz_vec_init(b + 1);
    fmpz_init(f);
    fmpz_init(pN);
    fmpz_init(pN_2);
    fmpz_init(sum);

    /* Compute s[j] and hence round chi[j] for j = 1,...,b */
    for (j = 1; j <= b; j++)
    {
        fmpz_zero(sum);
        for (i = 1; i <= j - 1; i++)
        {
            fmpz_mul(f, s + (j - i), chi + i);
            fmpz_sub(sum, sum, f);
        }

        /* Now s[j] + j chi[j] = sum, but we don't know s[j] yet */
        fmpz_mul_ui(f, chi + j, j);
        fmpz_sub(f, sum, f);

        /* Round f to s[j] in $[- \floor{p^N/2}, \ceil{p^N/2})$ */
        n_chi = diagfrob_prec_chi(n, b, p, a, j);

        fmpz_pow_ui(pN, p, n_chi);
        fmpz_fdiv_q_ui(pN_2, pN, 2);
        fmpz_sub(pN_2, pN, pN_2);

        fmpz_mod(s + j, f, pN);
        if (fmpz_cmp(s + j, pN_2) >= 0)
            fmpz_sub(s + j, s + j, pN);

        /* Now we have s[j], we can recover chi[j] exactly */
        fmpz_sub(sum, sum, s + j);

        assert(fmpz_divisible_si(sum, j));

        fmpz_divexact_si(chi + j, sum, j);
    }

    _fmpz_vec_clear(s, b + 1);
    fmpz_clear(f);
    fmpz_clear(pN);
    fmpz_clear(pN_2);
    fmpz_clear(sum);
}

void diagfrob_zetafunction(fmpz_poly_t chi, 
                           long n, long d, const fmpz_t p, long a)
{
    _diagfrob_zetafunction(chi->coeffs, n, d, p, a);
}

