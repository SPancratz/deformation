#include "vec.h"

int _vec_is_zero(const char *vec, long n, const mat_ctx_t ctx)
{
    long i;

    for (i = 0; i < n; i++)
        if (!ctx->is_zero(ctx, vec + i * ctx->size))
            return 0;
    return 1;
}

