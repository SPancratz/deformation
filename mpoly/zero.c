#include "mpoly.h"

void mpoly_zero(mpoly_t rop, const mat_ctx_t ctx)
{
    mpoly_iter_t iter;
    mpoly_term t;

    mpoly_iter_init(iter, rop);

    while ((t = mpoly_iter_next(iter)))
    {
        ctx->clear(ctx, t->val);
        free(t->val);
    }

    mpoly_iter_clear(iter);

    RBTREE_CLEAR(mpoly, rop->dict, NULL, NULL);
}

