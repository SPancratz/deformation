#include "mpoly.h"

void mpoly_set(mpoly_t rop, const mpoly_t op, const mat_ctx_t ctx)
{
    if (rop != op)
    {
        mpoly_iter_t iter;
        mpoly_term t;

        mpoly_clear(rop, ctx);
        mpoly_init(rop, op->n, ctx);

        mpoly_iter_init(iter, op);
        while ((t = mpoly_iter_next(iter)))
        {
            mon_t m2;
            char *c1;
            void *c2;

            c1 = malloc(ctx->size);
            ctx->init(ctx, c1);
            ctx->set(ctx, c1, t->val);

            RBTREE_INSERT(mpoly, &m2, &c2, rop->dict, t->key, c1, &mon_cmp);
        }
        mpoly_iter_clear(iter);
    }
}
