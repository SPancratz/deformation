#include "diagfrob.h"

/*
    Computes an entry in the matrix representing $p^{-1} F$.

    This is equal to 
    \begin{equation*}
    (-1)^{u'+v'} \frac{(v'-1)!}{(u'-1)!} p^n \alpha_{u+1,v+1}^{-1}.
    \end{equation*}

    The entry is computed modulo $p^N$.
 */
void 
diagfrob_entry_mpq(mpq_t rop, const fmpz *a, long n, long d, 
                   const mon_t u, const mon_t v, 
                   const padic_ctx_t ctx)
{
    long i;
    long up, vp;
    mon_t up1, vp1;
    mpq_t x;
    mpz_t y;
    
    const long p = *(ctx->p);
    
    mon_init(up1);
    mon_init(vp1);
    
    /* Check obvious condition for non-zero entry */
    for (i = 0; i <= n; i++)
    {
        mon_set_exp(up1, i, mon_get_exp(u, i) + 1);
        mon_set_exp(vp1, i, mon_get_exp(v, i) + 1);
        
        if ((p * mon_get_exp(up1, i) - mon_get_exp(vp1, i)) % d != 0)
        {
            mpq_set_si(rop, 0, 1);
            mon_clear(up1);
            mon_clear(vp1);
            return;
        }
    }
    
    mpq_init(x);
    mpz_init(y);
    
    up = ((long) mon_degree(u) + (n + 1)) / d;
    vp = ((long) mon_degree(v) + (n + 1)) / d;
    
    /* Set rop to alpha_{u+1,v+1}^{-1} */
    diagfrob_alpha_mpq(rop, a, n, d, up1, vp1, ctx);
    if (mpq_sgn(rop) == 0)
    {
        printf("ERROR (diagfrob_entry).  Found a zero value of alpha.\n");
        abort();
    }
    mpq_inv(rop, rop);
    
    /* Multiply by p^n */
    mpz_set_si(y, p);
    mpz_pow_ui(y, y, n);
    mpz_mul(mpq_numref(rop), mpq_numref(rop), y);
    mpq_canonicalize(rop);
    
    /* Multiply by (-1)^{u'+v'} (v'-1)! / (u'-1)! */
    mpz_fac_ui(mpq_numref(x), vp - 1 <= 0 ? 1 : vp - 1);
    mpz_fac_ui(mpq_denref(x), up - 1 <= 0 ? 1 : up - 1);
    mpq_canonicalize(x);
    if ((up - vp) % 2 != 0)
    {
        mpz_neg(mpq_numref(x), mpq_numref(x));
    }
    mpq_mul(rop, x, rop);
    
    mon_clear(up1);
    mon_clear(vp1);
    mpq_clear(x);
    mpz_clear(y);
}

