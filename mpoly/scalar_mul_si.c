#include "mpoly.h"

void mpoly_scalar_mul_si(mpoly_t rop, const mpoly_t op, long x, 
                         const mat_ctx_t ctx)
{
    char *c;

    c = malloc(ctx->size);
    ctx->init(c);
    ctx->set_si(c, x);

    mpoly_scalar_mul(rop, op, c, ctx);

    ctx->clear(c);
    free(c);
}

