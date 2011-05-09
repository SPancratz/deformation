#include "mpoly.h"

void mpoly_randhom(mpoly_t rop, flint_rand_t state, long d, long N, 
                   const mat_ctx_t ctx)
{
    long i, n = rop->n;

    mpoly_zero(rop, ctx);

    for (i = 0; i < N; i++)
    {
        int ins;
        mon_t m1, m2;
        char *c1;
        void *c2;

        mon_init(m1);
        {
            long j;
            ulong e, sum = 0;

            for (j = 0; j < n - 1; j++)
            {
                e = n_randint(state, d - sum + 1);
                mon_set_exp(m1, j, e);
                sum += e;
            }

            mon_set_exp(m1, j, d - sum);
        }

        c1 = malloc(ctx->size);
        ctx->init(c1);
        ctx->randtest_not_zero(c1, state);

        ins = RBTREE_INSERT(mpoly, &m2, &c2, rop->dict, m1, c1, &mon_cmp);

        if (ins)
        {
            mon_clear(m2);
            ctx->clear(c2);
            free(c2);
        }
    }
}
