#include "vec.h"

void _vec_set(char *vec1, const char *vec2, long n, const ctx_t ctx)
{
    long i;

    for (i = 0; i < n; i++)
        ctx->set(vec1 + i * ctx->size, vec2 + i * ctx->size);
}

