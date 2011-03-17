#include "mpoly.h"

void mpoly_mul_mon(mpoly_t rop, const mpoly_t op, const mon_t m,
                   const mat_ctx_t ctx)
{
    mpoly_iter_t iter;
    mpoly_term t;
    
    if (rop != op)
        mpoly_set(rop, op, ctx);

    mpoly_iter_init(iter, rop);
    while ((t = mpoly_iter_next(iter)))
    {
        mon_mul(t->key, t->key, m);
    }
    mpoly_iter_clear(iter);
}

