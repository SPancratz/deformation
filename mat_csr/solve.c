#include <assert.h>

#include "mat.h"
#include "vec.h"

#include "mat_csr.h"

#define DEBUG  0

void mat_csr_solve(char *x, const mat_csr_solve_t s, const char *b, 
                   const ctx_t ctx)
{
    char *c, *t;
    long i, k, m;

    assert(s->m == s->n);

    #if (DEBUG > 0)
    printf("mat_csr_solve()\n"), fflush(stdout);
    #endif

    m = s->m;

    c = _vec_init(m, ctx);
    t = _vec_init(1, ctx);

    _vec_permute(c, b, m, s->pi, ctx);

    for (k = 0; k < s->nb; k++)
    {
        const long i1  = s->B[k];
        const long i2  = s->B[k + 1];
        const long len = i2 - i1;

        #if (DEBUG > 0)
        printf("  Block %ld out of %ld\n", k, s->nb), fflush(stdout);
        #endif

        _mat_lup_solve(x + i1 * ctx->size, s->LU + i1, len, len, 
                       s->P + i1, c + i1 * ctx->size, ctx);

        /* Update c.  Subtract A{bk} y{k} from c{b} for b > k */
        for (i = i2; i < m; i++)
        {
            long q;
            const long p1 = s->p[i];
            const long p2 = s->p[i] + s->lenr[i];

            for (q = p1; q < p2; q++)
            {
                const long j = s->j[q];

                if (i1 <= j && j < i2)
                {
                    ctx->mul(ctx, t, s->x + q * ctx->size, x + j * ctx->size);
                    ctx->sub(ctx, c + i * ctx->size, c + i * ctx->size, t);
                }
            }
        }
    }

    /* Find Q^{-1} x */

    _vec_permute(x, x, m, s->qi, ctx);

    _vec_clear(c, m, ctx);
    _vec_clear(t, 1, ctx);
}

