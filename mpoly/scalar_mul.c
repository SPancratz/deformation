#include "mpoly.h"

void mpoly_scalar_mul(mpoly_t rop, const mpoly_t op, const void *x, 
                      const ctx_t ctx)
{
    mpoly_iter_t iter;
    mpoly_term t;

    if (ctx->is_zero(x))
    {
        mpoly_zero(rop, ctx);
        rop->n = op->n;
        return;
    }

    if (rop != op)
        mpoly_set(rop, op, ctx);

    mpoly_iter_init(iter, rop);
    while ((t = mpoly_iter_next(iter)))
    {
        ctx->mul(t->val, t->val, x);
    }
    mpoly_iter_clear(iter);
}

