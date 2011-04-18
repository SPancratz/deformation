#include "vec.h"

void _vec_sub(char *res, const char *vec1, const char *vec2, long n, 
              const mat_ctx_t ctx)
{
    long i;

    for (i = 0; i < n; i++)
        ctx->sub(ctx, res + i * ctx->size, vec1 + i * ctx->size, vec2 + i * ctx->size);
}

