#include "vec.h"

void _vec_zero(char *vec, long n, const mat_ctx_t ctx)
{
    long i;

    for (i = 0; i < n; i++)
        ctx->zero(ctx, vec + i * ctx->size);
}

