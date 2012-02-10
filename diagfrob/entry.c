#include "fmpq.h"
#include "diagfrob.h"

/*
    Computes an entry in the matrix representing $p^{-1} F$.

    This is equal to 
    \begin{equation*}
    (-1)^{u'+v'} \frac{(v'-1)!}{(u'-1)!} p^n \alpha_{u+1,v+1}^{-1}.
    \end{equation*}

    Expects \code{bound} to be a lower bound for the valuations 
    of the entries in the matrix $p^{-1} F$.

    The entry is computed modulo $p^N$.
 */
void 
diagfrob_entry(padic_t rop, const fmpz *a, long n, long d, 
               const mon_t u, const mon_t v, long bound, 
               const padic_ctx_t ctx)
{
    mon_t up1, vp1;       /* u + 1, v + 1 */
    fmpq_t x;
    long i;
    long e, N1, N2;
    padic_ctx_t ctx2;

    const long p = *(ctx->p), N = ctx->N;

    mon_init(up1);
    mon_init(vp1);

    /* Check obvious condition for non-zero entry */
    for (i = 0; i <= n; i++)
    {
        mon_set_exp(up1, i, mon_get_exp(u, i) + 1);
        mon_set_exp(vp1, i, mon_get_exp(v, i) + 1);
        
        if ((p * mon_get_exp(up1, i) - mon_get_exp(vp1, i)) % d != 0)
        {
            padic_zero(rop);
            mon_clear(up1);
            mon_clear(vp1);
            return;
        }
    }

    fmpq_init(x);

    /* Compute (-1)^{u'+v'} (v'-1)! / (u'-1)! p^n */
    {
        long up = ((long) mon_degree(u) + (n + 1)) / d;
        long vp = ((long) mon_degree(v) + (n + 1)) / d;

        fmpz_fac_ui(fmpq_numref(x), (vp - 1 <= 0) ? 1 : vp - 1);
        fmpz_fac_ui(fmpq_denref(x), (up - 1 <= 0) ? 1 : up - 1);
        e  = _fmpz_remove(fmpq_numref(x), ctx->p, ctx->pinv);
        e -= _fmpz_remove(fmpq_denref(x), ctx->p, ctx->pinv);
        if ((up + vp) & 1L)
            fmpz_neg(fmpq_numref(x), fmpq_numref(x));
        e += n;
    }

    /* Update the precision */
    N1 = ctx->N - e;
    N2 = ctx->N + e - 2 * bound;

    padic_ctx_init(ctx2, ctx->p, N2, PADIC_SERIES);

    /* Set rop to alpha_{u+1,v+1}^{-1} */
    diagfrob_alpha(rop, a, n, d, up1, vp1, ctx2);

    _padic_inv(padic_unit(rop), 
               padic_unit(rop), ctx->p, N1 + padic_val(rop));
    e -= padic_val(rop);
    padic_val(rop) = 0;

    /* Compute final product */
    {
        _padic_inv(fmpq_denref(x), fmpq_denref(x), ctx->p, N - e);
        fmpz_mul(padic_unit(rop), padic_unit(rop), fmpq_numref(x));
        fmpz_mul(padic_unit(rop), padic_unit(rop), fmpq_denref(x));
        padic_shift(rop, rop, e, ctx);
    }

    mon_clear(up1);
    mon_clear(vp1);
    fmpq_clear(x);
    padic_ctx_clear(ctx2);
}

