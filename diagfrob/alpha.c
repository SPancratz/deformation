#include <limits.h>
#include "diagfrob.h"

/*
    Computes the entry $\alpha_{uv}$.

    This is defined by 
    \begin{equation*}
    \alpha_{uv} = \pi^{v'-u'} \prod_{i=0}^n \sum_{m,r \geq 0} \lambda_m \binom{u_i}{d}_r (-1)^r \pi^{-r} \bigl(\hat{a}_i\bigr)^{m-r}.
    \end{equation*}

    Although this is defined over $\mathbf{Q}_p(\pi)$, it turns out to 
    actually lie in $\mathbf{Q}_p$.  Since we are working with finite 
    precision, everything is happening over $\mathbf{Q}$.

    Result is computed modulo $p^N$.  Following [Lau], it suffices the 
    inner sum for non-negative values of $m, r$ such that 
    $m \leq 2 p^2 (N + (2n - 1)) / (p-1)$.

    \param[out]  rop  $\alpha_{u,v}$
    \param[in]   a    Array of coefficients, of length $n+1$
    \param[in]   n    One less than the number of variables
    \param[in]   d    Homogeneous degree
    \param[in]   p    Prime number
    \param[in]   u    Monomial, first index
    \param[in]   v    Monomial, second index

    \see  [Lau] - Lauder, <i>Counting solutions to equations in many variables 
                  over finite fields</i>, 2003
 */
void 
diagfrob_alpha(padic_t rop, const fmpz *a, long n, long d, 
               const mon_t u, const mon_t v, 
               const padic_ctx_t ctx0)
{
    long m, r;                 /* Index variables */
    long M;                    /* Bounds */
    long Nt;                   /* Temporarily increased precision */
    long i, j;                 /* Index variables */
    long ui, vi;               /* u[i], v[i] */
    long lhs;                  /* (p u[i] - v[i]) / d */
    long quot;
    
    fmpz_t ppow;                /* Power of p */
    fmpz_t dinv;                /* d^{-1} */
    padic_t factor;             /* ith factor */
    long e_product, e_factor;   /* Exponents of pi */
    padic_t summand;
    padic_t lambda;
    padic_t ai;                 /* Teichmuller lift */
    fmpz_t falling_fac1;        /* u[i] (u[i] + d) ... (u[i] + (r-1) d) */
    fmpz_t falling_fac2;        /* u[i] + (r-1) d */
    long e_falling_fac;

    const long p = *(ctx0->p), N = ctx0->N;
    padic_ctx_t ctx;

    /* Check the size of the input parameters */
    {
        mpz_t t;

        if (N > LONG_MAX - (2*n - 1))
        {
            printf("ERROR (diagfrob_alpha).  N too large\n");
            abort();
        }

        mpz_init_set_ui(t, 2);
        mpz_mul_ui(t, t, p);
        mpz_mul_ui(t, t, p);
        mpz_mul_ui(t, t, N + (2*n - 1));
        mpz_fdiv_q_ui(t, t, p - 1);

        if (!mpz_fits_slong_p(t))
        {
            printf("ERROR (diagfrob_alpha).  M too large\n");
            abort();
        }
        else
            M = mpz_get_si(t);

        mpz_clear(t);
    }

    /* Ensure that no factor in the product is (obviously) zero. */
    for (i = 0; i <= n; i++)
    {
        ui  = (long) mon_get_exp(u, i);
        vi  = (long) mon_get_exp(v, i);
        lhs = p * ui - vi;
        if (lhs % d != 0)
        {
            padic_zero(rop);
            return;
        }
    }

    /* Find the minimum of 0 and the (p u_i - v_i) / d, call it Nt for now */
    Nt = 0;
    for (i = 0; i <= n; i++)
    {
        ui  = (long) mon_get_exp(u, i);
        vi  = (long) mon_get_exp(v, i);
        lhs = (p * ui - vi) / d;
        Nt  = FLINT_MIN(Nt, lhs);
    }
    
    /* Set bounds */
    Nt = (Nt >= 0) ? N : (N + (n + 1) * (-Nt * (p - 1) / (p * p)) + 1);

    /* Initialisation */
    padic_ctx_init(ctx, ctx0->p, Nt, PADIC_SERIES);

    padic_init(factor, ctx);
    padic_init(summand, ctx);
    padic_init(lambda, ctx);
    padic_init(ai, ctx);

    fmpz_init(ppow);
    fmpz_init(falling_fac1);
    fmpz_init(falling_fac2);
    fmpz_init(dinv);

    fmpz_pow_ui(ppow, ctx->p, Nt);
    fmpz_set_ui(dinv, d);
    _padic_inv(dinv, dinv, ctx->p, ctx->N);

    /* Outer product.. */
    padic_one(rop, ctx);
    e_product = ((long) mon_degree(v) - (long) mon_degree(u)) / d;
    
    for (i = 0; i <= n; i++)
    {
        ui  = (long) mon_get_exp(u, i);
        vi  = (long) mon_get_exp(v, i);
        lhs = (p * ui - vi) / d;
        
        /*
           Inner sum over all non-negative m, r such that lhs = m - p r
           and m <= M
         */
        
        /* Initial values of m, r */
        r = (lhs >= 0 ? 0 : (lhs % p == 0 ? (-lhs) / p : (-lhs) / p + 1));
        m = p * r + lhs;
        if (m > M)
        {
            printf("ERROR (alpha_padic).  Bound on m is too low.");
            abort();
        }

        /* Power of Teichmuller lift */
        padic_set_fmpz(ai, a + i, ctx);
        padic_teichmuller(ai, ai, ctx);
        fmpz_powm_ui(padic_unit(ai), 
                     padic_unit(ai), DIAGFROB_MOD(m - r, p - 1), ppow);
        
        padic_zero(factor);
        e_factor = DIAGFROB_MOD(m - r, p - 1);
        
        /* Compute the pre-initial value of (u_i / d)_r */
        if (r == 0)
        {
            fmpz_set_ui(falling_fac1, 1);
            e_falling_fac = 0;
        }
        else
        {
            fmpz_powm_ui(falling_fac1, dinv, r - 1, ppow);
            for (j = 0; j < r - 1; j++)
            {
                fmpz_mul_ui(falling_fac1, falling_fac1, ui + j * d);
                fmpz_mod(falling_fac1, falling_fac1, ppow);
            }
            e_falling_fac = _fmpz_remove(falling_fac1, ctx->p, ctx->pinv);
        }

        for ( ; m <= M; r++, m += p)
        {
            diagfrob_coefficient(lambda, m, ctx);
            
            if (r > 0)
            {
                if (ui + (r - 1) * d == 0)
                    break;

                fmpz_set_ui(falling_fac2, ui + (r - 1) * d);
                e_falling_fac += _fmpz_remove(falling_fac2, ctx->p, ctx->pinv);
                fmpz_mul(falling_fac1, falling_fac1, falling_fac2);
                fmpz_mul(falling_fac1, falling_fac1, dinv);
                fmpz_mod(falling_fac1, falling_fac1, ppow);
            }
            
            /*
               The relevant monomial is pi^{(m mod (p-1)) - r}, we write this 
               exponent as quot * (p-1) + rem with 0 <= rem < p-1
             */
            quot = fdiv_si(m % (p - 1) - r, p - 1);

            fmpz_mul(padic_unit(summand), padic_unit(lambda), falling_fac1);
            fmpz_mod(padic_unit(summand), padic_unit(summand), ppow);
            padic_val(summand) = 0;

            if ((r + quot) % 2 != 0)
                padic_neg(summand, summand, ctx);
            
            padic_val(summand) = padic_val(lambda) + e_falling_fac + quot;
            padic_reduce(summand, ctx);
            
            /* Now summand is the p-adic coefficient of pi^{rem} */
            padic_add(factor, factor, summand, ctx);
        }

        fmpz_mul(padic_unit(factor), padic_unit(factor), padic_unit(ai));
        padic_reduce(factor, ctx);
        
        /* Set the product */
        padic_mul(rop, rop, factor, ctx);
        e_product += e_factor;
        
        if (padic_is_zero(rop, ctx))
            break;
    }

    if (!padic_is_zero(rop, ctx))
    {
        /* DEBUG.  This condition should never be met */
        if (e_product % (p - 1) != 0)
        {
            printf("ERROR (alpha).  Result not defined over Qp.\n");
            abort();
        }
        
        /* Multiply the product by (-p)^quot */
        quot = e_product / (p - 1);
        if (quot % 2 != 0)
            padic_neg(rop, rop, ctx);
        padic_shift(rop, rop, quot, ctx0);  /* XXX: N, not Nt */
    }
    
    padic_clear(factor, ctx);
    padic_clear(summand, ctx);
    padic_clear(lambda, ctx);
    padic_clear(ai, ctx);
    padic_ctx_clear(ctx);

    fmpz_clear(ppow);
    fmpz_clear(falling_fac1);
    fmpz_clear(falling_fac2);
    fmpz_clear(dinv);
}

