#include <limits.h>

#include "diagfrob.h"

#define DEBUG 0

/*
    Computes the entry $\alpha_{uv}$.

    This is defined by 
    \begin{equation*}
    \alpha_{uv} = \pi^{v'-u'} \prod_{i=0}^n \sum_{m,r \geq 0} \lambda_m \binom{u_i}{d}_r (-1)^r \pi^{-r} \bigl(\hat{a}_i\bigr)^{m-r}.
    \end{equation*}

    Although this is defined over $\mathbf{Q}_p(\pi)$, it turns out to 
    actually lie in $\mathbf{Q}_p$.  Since we are working with finite 
    precision, everything is happening over $\mathbf{Q}$.

    The length to which the infinite sum above is computed depends on 
    the precision bound $N$.  Following [Lau], it suffices to consider 
    non-negative values of $m, r$ such that 
    $m \leq 2 p^2 (N + (2n - 1)) / (p-1)$.

    Here, all computations are carried out over the rationals.

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
diagfrob_alpha_mpq(mpq_t rop, const fmpz *a, long n, long d, 
                   const mon_t u, const mon_t v, 
                   const padic_ctx_t ctx)
{
    long M;                    /* Bound on m */
    long i, m, r;              /* Index variables */
    long lhs;                  /* p u[i] - v[i] */
    long quot;
    long ui, vi;               /* u[i], v[i] */
    mpq_t s1, s2;              /* Temporary terms */
    mpz_t pmpz;                /* Integer p */
    mpq_t factor;              /* ith factor */
    long e_product, e_factor;  /* Exponents */
    mpz_t ai;                  /* Teichmuller lift */
    
    const long p = *(ctx->p), N = ctx->N;

    M = 2 * p * p * (N + 2 * n - 1) / (p - 1);

    /* Check the size of the input parameters */
    {
        mpz_t x;

        if (N > LONG_MAX - (2*n - 1))
        {
            printf("ERROR (diagfrob_alpha_mpq).  N is too large.\n");
            abort();
        }

        mpz_init_set_ui(x, 2);
        mpz_mul_ui(x, x, p);
        mpz_mul_ui(x, x, p);
        mpz_mul_ui(x, x, N + (2*n - 1));
        mpz_fdiv_q_ui(x, x, p - 1);

        if (!mpz_fits_slong_p(x))
        {
            printf("ERROR (diagfrob_alpha_mpq).  M is too large.\n");
            abort();
        }

        mpz_clear(x);
    }
 
    /* Ensure that no factor in the product is (obviously) zero */
    for (i = 0; i <= n; i++)
    {
        ui  = (long) mon_get_exp(u, i);
        vi  = (long) mon_get_exp(v, i);
        lhs = p * ui - vi;
        if (lhs % d != 0)
        {
            mpq_set_si(rop, 0, 1);
            return;
        }
    }
    
    if (DEBUG)
    {
        printf("diagfrob_alpha_mpq(");
        printf("u = %s, v = %s)\n", mon_get_str_pretty(u, n + 1, NULL), 
                                    mon_get_str_pretty(v, n + 1, NULL));
    }
    
    mpq_init(factor);
    mpq_init(s1);
    mpq_init(s2);
    mpz_init_set_si(pmpz, p);
    mpz_init(ai);
    
    /* Outer product */
    mpq_set_si(rop, 1, 1);
    e_product = ((long) mon_degree(v) - (long) mon_degree(u)) / d;
    
    if (DEBUG)
        printf("  First factor:  (pi)^%ld\n", e_product);
    
    for (i = 0; i <= n; i++)
    {
        if (DEBUG)
            printf("  Loop i = %ld\n", i);
        
        ui  = (long) mon_get_exp(u, i);
        vi  = (long) mon_get_exp(v, i);
        lhs = (p * ui - vi) / d;

        /*
           Inner sum over all non-negative m, r such that lhs = m - pr
           and m <= M
         */
        
        /* Initial values of m, r */
        r = (lhs >= 0 ? 0 : (lhs % p == 0 ? (-lhs)/p : (-lhs)/p + 1));
        m = p * r + lhs;
        if (m > M)
        {
            printf("ERROR (diagfrob_alpha_mpq).  Bound on m is too low.");
            abort();
        }
        
        mpq_set_si(factor, 0, 1);
        e_factor = (m - r) % (p - 1);
        if (e_factor < 0)
            e_factor += (p - 1);
        
        {
            padic_t x;

            padic_init(x, ctx);
            padic_set_fmpz(x, a + i, ctx);
            padic_teichmuller(x, x, ctx);
            padic_get_mpz(ai, x, ctx);
            padic_clear(x, ctx);
        }
        
        while (m <= M)
        {
            diagfrob_coefficient_mpq(s1, m, p);
            diagfrob_falling_fac_mpq(s2, ui, d, r);
            mpq_mul(s1, s1, s2);
            if (r & 1L)
                mpq_neg(s1, s1);
            if (m >= r)
            {
                mpz_pow_ui(mpq_numref(s2), ai, m - r);
                mpz_set_si(mpq_denref(s2), 1);
            }
            else
            {
                mpz_set_si(mpq_numref(s2), 1);
                mpz_pow_ui(mpq_denref(s2), ai, r - m);
                if (mpz_sgn(mpq_denref(s2)) < 0)
                {
                    mpz_neg(mpq_numref(s2), mpq_numref(s2));
                    mpz_neg(mpq_denref(s2), mpq_denref(s2));
                }
            }
            mpq_mul(s1, s1, s2);
            
            /* 
               Now s1 is the coefficient of a monomial term with monomial 
               pi^{(m mod (p-1)) - r}.  Write this exponent as 
               quot * (p-1) + rem
             */
            quot = fdiv_si(m % (p - 1) - r, p - 1);
            
            mpq_set_si(s2, -p, 1);
            if (quot >= 0)
            {
                mpz_pow_ui(mpq_numref(s2), mpq_numref(s2), quot);
            }
            else
            {
                mpz_pow_ui(mpq_numref(s2), mpq_numref(s2), -quot);
                mpz_swap(mpq_numref(s2), mpq_denref(s2));
                if (mpz_sgn(mpq_denref(s2)) < 0)
                {
                    mpz_neg(mpq_numref(s2), mpq_numref(s2));
                    mpz_neg(mpq_denref(s2), mpq_denref(s2));
                }
            }
            mpq_mul(s1, s1, s2);
            
            /* Now s1 is the coefficient of pi^{rem} */
            mpq_add(factor, factor, s1);
            
            /* Increment m, r */
            r++;
            m += p;
        }
        
        /* Set the product */
        mpq_mul(rop, rop, factor);
        e_product += e_factor;
        
        if (mpq_sgn(rop) == 0)
            break;
    }
    
    if (mpq_sgn(rop) != 0)
    {
        /* DEBUG.  This condition should never be met */
        if (DEBUG && (e_product % (p - 1) != 0))
        {
            printf("ERROR (diagfrob_alpha_mpq).  Result not defined over Qp.\n");
            abort();
        }
        
        /* Multiply the product by (-p)^quot */
        quot = e_product / (p - 1);
        mpq_set_si(factor, -p, 1);
        if (quot >= 0)
        {
            mpz_pow_ui(mpq_numref(factor), mpq_numref(factor), quot);
        }
        else
        {
            mpz_pow_ui(mpq_numref(factor), mpq_numref(factor), -quot);
            mpz_swap(mpq_numref(factor), mpq_denref(factor));
            if (mpz_sgn(mpq_denref(factor)) < 0)
            {
                mpz_neg(mpq_numref(factor), mpq_numref(factor));
                mpz_neg(mpq_denref(factor), mpq_denref(factor));
            }
        }
        mpq_mul(rop, rop, factor);
    }
    
    mpq_clear(factor);
    mpq_clear(s1);
    mpq_clear(s2);
    mpz_clear(pmpz);
    mpz_clear(ai);
}

