#include "vec.h"

void _vec_scalar_div(char *res, const char *vec, long n, const char *x, 
                     const mat_ctx_t ctx)
{
    long i;

    for (i = 0; i < n; i++)
        ctx->div(ctx, res + i * ctx->size, vec + i * ctx->size, x);
}

