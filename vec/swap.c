#include "vec.h"

void _vec_swap(char *vec1, char *vec2, long n, const ctx_t ctx)
{
    long i;

    for (i = 0; i < n; i++)
        ctx->swap(ctx, vec1 + i * ctx->size, vec2 + i * ctx->size);
}

