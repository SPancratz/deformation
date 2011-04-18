#include "vec.h"

void _vec_scalar_mul(char *res, const char *vec, long n, const char *x, 
                     const mat_ctx_t ctx)
{
    long i;

    for (i = 0; i < n; i++)
        ctx->mul(ctx, res + i * ctx->size, vec + i * ctx->size, x);
}

