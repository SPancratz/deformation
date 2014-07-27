/******************************************************************************

    Copyright (C) 2012 Sebastian Pancratz
 
******************************************************************************/

#include "fmpz_mod_poly.h"
#include "ulong_extras.h"
#include "qadic.h"

/*
    Assumes that \code{len1} and \code{len2} are positive but at 
    most~$d$, and also that \code{len1} is at least $6$.

    The latter assumption guarantees that $\ceil{n/B} \geq 2$, 
    i.e.\ $n \geq 2B$ so $n \geq 2 \ceil{\sqrt{n}}$.
 */

static void 
_fmpz_mod_poly_compose_smod_rectangular(fmpz *rop, 
                           const fmpz *op1, slong len1, 
                           const fmpz *op2, slong len2, 
                           const fmpz *a, const slong *j, slong lena, 
                           const fmpz_t p)
{
    const slong d = j[lena - 1];

    if (len2 == 1)
    {
        _fmpz_mod_poly_evaluate_fmpz(rop, op1, len1, op2, p);
        _fmpz_vec_zero(rop + 1, d - 1);
    }
    else
    {
        const slong B = n_sqrt(len1);
        slong i, k;
        fmpz *pows, *t;

        pows = _fmpz_vec_init((B + 2) * d);
        t    = _fmpz_vec_init(2 * d - 1);

        fmpz_one(pows + 0 * d + 0);
        _fmpz_vec_set(pows + 1 * d, op2, len2);
        for (i = 2; i <= B; i++)
        {
            _fmpz_poly_mul(pows + i * d, pows + (i - 1) * d, d, op2, len2);
            _fmpz_poly_reduce(pows + i * d, d + len2 - 1, a, j, lena);
            _fmpz_vec_scalar_mod_fmpz(pows + i * d, pows + i * d, d, p);
        }

        _fmpz_vec_zero(rop, d);

        for (i = (len1 + B - 1) / B - 1; i >= 0; i--)
        {
            _fmpz_poly_mul(t, rop, d, pows + B * d, d);
            _fmpz_poly_reduce(t, 2 * d - 1, a, j, lena);

            _fmpz_vec_set(rop, t, d);
            fmpz_add(rop + 0, rop + 0, op1 + i*B);
            for (k = FLINT_MIN(B, len1 - i*B) - 1; k > 0; k--)
            {
                _fmpz_vec_scalar_addmul_fmpz(rop, pows + k * d, d, op1 + (i*B + k));
            }

            _fmpz_vec_scalar_mod_fmpz(rop, rop, d, p);
        }

        _fmpz_vec_clear(pows, (B + 2) * d);
        _fmpz_vec_clear(t, 2 * d - 1);
    }
}

static void 
_fmpz_mod_poly_compose_smod_horner(fmpz *rop, 
                           const fmpz *op1, slong len1, 
                           const fmpz *op2, slong len2, 
                           const fmpz *a, const slong *j, slong lena, 
                           const fmpz_t p)
{
    const slong d = j[lena - 1];

    if (len1 == 1)
    {
        fmpz_set(rop, op1);
        _fmpz_vec_zero(rop + 1, d - 1);
    }
    else if (len2 == 1)
    {
        _fmpz_mod_poly_evaluate_fmpz(rop, op1, len1, op2, p);
        _fmpz_vec_zero(rop + 1, d - 1);
    }
    else
    {
        slong i;
        fmpz *t;

        t = _fmpz_vec_init(2*d - 1);

        _fmpz_vec_zero(rop, d);

        for (i = len1 - 1; i >= 0; i--)
        {
            _fmpz_poly_mul(t, rop, d, op2, len2);
            _fmpz_poly_reduce(t, d + len2 - 1, a, j, lena);
            _fmpz_poly_add(rop, t, d, op1 + i, 1);
            _fmpz_vec_scalar_mod_fmpz(rop, rop, d, p);
        }

        _fmpz_vec_clear(t, 2*d - 1);
    }
}

/* 
    Computes the composition $f(g(X))$ modulo the sparse polynomial 
    given by the data \code{(a, j, lena)}, which is assumed to be 
    of degree~$d \geq 2$.

    Sets the vector \code{(rop, d)}.

    Assumes that \code{len1} and \code{len2} are positive but at 
    most~$d$.

    Does not support aliasing.
 */

void 
_fmpz_mod_poly_compose_smod(fmpz *rop, 
                           const fmpz *op1, slong len1, 
                           const fmpz *op2, slong len2, 
                           const fmpz *a, const slong *j, slong lena, 
                           const fmpz_t p)
{
    if (len1 < 6)
    {
        _fmpz_mod_poly_compose_smod_horner(rop, op1, len1, op2, len2, a, j, lena, p);
    }
    else
    {
        _fmpz_mod_poly_compose_smod_rectangular(rop, op1, len1, op2, len2, a, j, lena, p);
    }
}

