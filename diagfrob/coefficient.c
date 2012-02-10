#include "diagfrob.h"

void diagfrob_coefficient(padic_t rop, long m, const padic_ctx_t ctx)
{
    if (m == 0)
    {
        padic_one(rop, ctx);
    }
    else if (m == 1)
    {
        padic_set_si(rop, fmpz_cmp_ui(ctx->p, 2) ? 1 : -2, ctx);
    }
    else
    {
        const fmpz *prime = ctx->p;
        const long p      = *prime;
        const long N      = ctx->N;

        const long low   = (m + (p - 1)) / p;  /* Ceiling of m / p */
        const long Nplus = (m / p) / (p - 1);  /* Excess of the precision */

        long ell;     /* Index variable for updating factorials */
        long n, k;    /* Indices of the double sum */
        long e;       /* Counting powers of p */

        fmpz_t r;     /* Term that is added to the coefficient at each step */
        fmpz_t s, t;  /* Factor that changes the above term */
        fmpz_t ppow;  /* Power of p */

        fmpz_init(r);
        fmpz_init(s);
        fmpz_init(t);
        fmpz_init(ppow);
        fmpz_pow_ui(ppow, prime, N + Nplus);

        fmpz_fac_ui(r, m);
        e  = m / (p - 1);
        e -= fmpz_remove(r, r, prime);
        _padic_inv(r, r, prime, N + Nplus);
        fmpz_pow_ui(t, prime, e);
        fmpz_mul(r, r, t);
        fmpz_mod(r, r, ppow);
        fmpz_set(padic_unit(rop), r);

        for (n = m - (p - 1), k = 1; n >= low; n -= (p - 1), k++)
        {
            fmpz_set_si(s, 1);                   /* Update (n-k)! */
            for (ell = 1; ell <= p; ell++)
                fmpz_mul_si(s, s, n - k + ell);
            fmpz_divexact_ui(s, s, p);           /* Update p^{q'} */

            fmpz_set_si(t, k);                   /* Update k! */
            e = fmpz_remove(t, t, prime);
            _padic_inv(t, t, prime, N + Nplus);

            fmpz_mul(r, r, s);                   /* Update product */
            fmpz_pow_ui(s, prime, e);
            fmpz_divexact(r, r, s);
            fmpz_mod(r, r, ppow);
            fmpz_mul(r, r, t);
            fmpz_mod(r, r, ppow);
            
            fmpz_add(padic_unit(rop), padic_unit(rop), r);  /* Update sum */
            fmpz_mod(padic_unit(rop), padic_unit(rop), ppow);
        }

        padic_val(rop) = 0;
        padic_reduce(rop, ctx);

        if ((m / (p - 1)) & 1L)
            padic_neg(rop, rop, ctx);
        
        fmpz_clear(r);
        fmpz_clear(s);
        fmpz_clear(t);
        fmpz_clear(ppow);
    }
}

