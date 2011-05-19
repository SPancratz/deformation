#include "mpoly.h"


void mpoly_scalar_div(mpoly_t rop, const mpoly_t op, const void *x, 
                      const ctx_t ctx)
{
    mpoly_iter_t iter;
    mpoly_term t;

    if (ctx->is_zero(ctx, x))
    {
        printf("ERROR (mpoly_scalar_div).  Division by zero.\n");
        abort();
    }

    if (rop != op)
        mpoly_set(rop, op, ctx);
    
    mpoly_iter_init(iter, rop);
    while ((t = mpoly_iter_next(iter)))
    {
        ctx->div(ctx, t->val, t->val, x);
    }
    mpoly_iter_clear(iter);
}

