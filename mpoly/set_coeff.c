#include "mpoly.h"

void mpoly_set_coeff(mpoly_t rop, const mon_t m, const void *c,
                     const ctx_t ctx)
{
    if (ctx->is_zero(ctx, c))
    {
        int del;
        mon_t m2;
        void *c2;

        del = RBTREE_DELETE(mpoly, &m2, &c2, rop->dict, m, &mon_cmp);
        if (del)
        {
            mon_clear(m2);
            ctx->clear(ctx, c2);
        }
    }
    else
    {
        int ins;
        mon_t m2;
        char *c1;
        void *c2;

        c1 = malloc(ctx->size);
        ctx->init(ctx, c1);
        ctx->set(ctx, c1, c);

        ins = RBTREE_INSERT(mpoly, &m2, &c2, rop->dict, m, c1, &mon_cmp);
        if (ins)
        {
            mon_clear(m2);
            ctx->clear(ctx, c2);
        }
    }
}

