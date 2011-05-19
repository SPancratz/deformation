#include "vec.h"

void _vec_scalar_div(char *res, const char *vec, long n, const char *x, 
                     const ctx_t ctx)
{
    long i;

    for (i = 0; i < n; i++)
        ctx->div(res + i * ctx->size, vec + i * ctx->size, x);
}

