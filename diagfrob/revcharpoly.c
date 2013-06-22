/******************************************************************************

    Copyright (C) 2011, 2012, 2013 Sebastian Pancratz

******************************************************************************/

#include "diagfrob.h"

void diagfrob_revcharpoly(fmpz_poly_t chi, 
                          const padic_mat_t phi, const padic_ctx_t ctx)
{
    const long b = padic_mat_nrows(phi);

    if (padic_mat_val(phi) >= 0)
    {
        fmpz_t t;
        fmpz_mat_t mat;

        fmpz_init(t);
        fmpz_mat_init(mat, b, b);

        fmpz_pow_ui(t, ctx->p, padic_mat_val(phi));
        fmpz_mat_scalar_mul_fmpz(mat, padic_mat(phi), t);
        fmpz_mat_charpoly(chi, mat);
        fmpz_poly_reverse(chi, chi, b + 1);

        fmpz_clear(t);
        fmpz_mat_clear(mat);
    }
    else  /* v < 0 */
    {
        fmpz_t t;
        long i;

        fmpz_init(t);

        fmpz_mat_charpoly(chi, padic_mat(phi));
        fmpz_poly_reverse(chi, chi, b + 1);
        for (i = 0; i <= b; i++)
        {
            fmpz_pow_ui(t, ctx->p, (-padic_mat_val(phi)) * i);

            assert(fmpz_divisible(chi->coeffs + i, t));

            fmpz_divexact(chi->coeffs + i, chi->coeffs + i, t);
        }

        fmpz_clear(t);
    }
}

