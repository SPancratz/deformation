#include <string.h>

#include "vec.h"

void 
_vec_permute(char *res, const char *vec, long n, long *pi, const mat_ctx_t ctx)
{
    long i;

    if (res == vec)
    {
        char *t = malloc(n * ctx->size);

        memcpy(t, vec, n * ctx->size);

        for (i = 0; i < n; i++)
            memcpy(res + i * ctx->size, t + pi[i] * ctx->size, ctx->size);

        free(t);
    }
    else
    {
        for (i = 0; i < n; i++)
            ctx->set(res + i * ctx->size, vec + pi[i] * ctx->size);
    }
}

