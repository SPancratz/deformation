#include "mpoly.h"

void mpoly_randtest(mpoly_t rop, flint_rand_t state, long d, long N, 
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
        mon_randtest(m1, state, n, d);

        c1 = malloc(ctx->size);
        ctx->init(ctx, c1);
        ctx->randtest_not_zero(ctx, c1, state);

        ins = RBTREE_INSERT(mpoly, &m2, &c2, rop->dict, m1, c1, &mon_cmp);

        if (ins)
        {
            mon_clear(m2);
            ctx->clear(ctx, c2);
            free(c2);
        }
    }
}

void mpoly_randtest_not_zero(mpoly_t rop, flint_rand_t state, long d, long N, 
                             const mat_ctx_t ctx)
{
    long i, n = rop->n;

    if (!(d >= 0 && N > 0))
    {
        printf("ERROR (mpoly_randtest_not_zero).\n");
        abort();
    }

    do 
        mpoly_randtest(rop, state, d, N, ctx);
    while (mpoly_is_zero(rop, ctx));
}
