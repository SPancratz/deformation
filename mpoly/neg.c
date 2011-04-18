#include "mpoly.h"

void mpoly_neg(mpoly_t rop, const mpoly_t op, const mat_ctx_t ctx)
{
    mpoly_iter_t iter;
    mpoly_term t;

    mpoly_set(rop, op, ctx);

    if (mpoly_is_zero(rop, ctx))
        return;

    mpoly_iter_init(iter, rop);
    while ((t = mpoly_iter_next(iter)))
    {
        ctx->neg(ctx, t->val, t->val);
    }
    mpoly_iter_clear(iter);
}

