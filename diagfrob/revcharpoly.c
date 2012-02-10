#include "diagfrob.h"

/*
    Given the lower $\ceil{(b + 1)/2}$ coefficients of the reverse 
    characteristic polynomial $f$ of $F$ to some sufficiently 
    high $p$-adic precision $N$, computes the reverse 
    characteristic polynomial $a$ over the integers.

    The number of variables of the diagonal hypersurfaces is $n + 1$ 
    and we have $q = p^r$.
 */

void diagfrob_revcharpoly(fmpz *a, const char *f, long n, long b, 
                          const ctx_t ctx)
{
    const fmpz *p = ctx->pctx->p;
    const long N  = ctx->pctx->N;
    const long r  = 1;

    /* Compute the lower half of the coefficients */
    {
        fmpz_t B, l, u;
        long i;

        fmpz_init(B);
        fmpz_init(l);
        fmpz_init(u);

        /* B = p^N, u = p^N - (p^N / 2) */
        fmpz_pow_ui(B, p, N);
        fmpz_fdiv_q_ui(l, B, 2);
        fmpz_neg(l, l);
        fmpz_add(u, l, B);

        for (i = 0; i <= (b + 2) / 2; i++)
        {
            if (padic_val((padic_struct *) _vec_entry(f, i, ctx)) < 0)
            {
                printf("ERROR (diagfrob_revcharpoly).\n");
                printf("  Entry f[%ld] = ", i);
                ctx->print(ctx, _vec_entry(f, i, ctx));
                printf(" is not an integer.\n");
                printf("  Valuation ord_p(f[%ld]) = %ld.\n", 
                    i, padic_val((padic_struct *) _vec_entry(f, i, ctx)));
                abort();
            }

            padic_get_fmpz(a + i, 
                (padic_struct *) _vec_entry(f, i, ctx), ctx->pctx);
            fmpz_mod(a + i, a + i, B);
            if (fmpz_cmp(a + i, u) >= 0)
                fmpz_sub(a + i, a + i, B);
        }

        fmpz_clear(B);
        fmpz_clear(l);
        fmpz_clear(u);
    }

    /* Compute the top coefficient */
    {
        fmpz_t B, l, u, pow;

        fmpz_init(B);
        fmpz_init(l);
        fmpz_init(u);
        fmpz_init(pow);

        fmpz_pow_ui(B, p, (r * (n - 1) * b) / 2 + 2);
        fmpz_fdiv_q_ui(l, B, 2);
        fmpz_neg(l, l);
        fmpz_add(u, l, B);
        fmpz_pow_ui(pow, p, (r * (n - 1) * b) / 2);

        if (padic_val((padic_struct *) _vec_entry(f, b, ctx)) < 0)
        {
            printf("ERROR (diagfrob_revcharpoly).\n");
            printf("  Entry f[%ld] = ", b);
            ctx->print(ctx, _vec_entry(f, b, ctx));
            printf(" is not an integer.\n");
            printf("  Valuation ord_p(f[%ld]) = %ld.\n", 
                b, padic_val((padic_struct *) _vec_entry(f, b, ctx)));
            abort();
        }

        {
            padic_ctx_t ctx2;

            padic_ctx_init(ctx2, p, (r * (n - 1) * b) / 2 + 2, PADIC_SERIES);
            padic_get_fmpz(a + b, (void *) _vec_entry(f, b, ctx), ctx2);
            padic_ctx_clear(ctx2);
        }

        fmpz_mod(a + b, a + b, B);
        if (fmpz_cmp(a + b, u) >= 0)
            fmpz_sub(a + b, a + b, B);

        if (fmpz_cmpabs(a + b, pow) != 0)
        {
            printf("ERROR (diagfrob_revcharpoly).  Top coefficient of wrong size.\n");
            abort();
        }

        fmpz_clear(B);
        fmpz_clear(l);
        fmpz_clear(u);
        fmpz_clear(pow);
    }

    /* Use a[i] = sgn * a[b-i] * q^{(n-1) i - (n-1) b / 2} */
    {
        int sgn = fmpz_sgn(a + b);
        long i;
        fmpz_t pow, q;

        fmpz_init(pow);
        fmpz_init(q);
        fmpz_pow_ui(q, p, r);

        for (i = (b + 2) / 2 + 1; i < b; i++)
        {
            fmpz_pow_ui(pow, q, (n - 1) * i - ((n - 1) * b) / 2);
            fmpz_mul(a + i, a + (b - i), pow);
            if (sgn == -1)
                fmpz_neg(a + i, a + i);
        }

        fmpz_clear(pow);
        fmpz_clear(q);
    }
}

