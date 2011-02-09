#include <stdlib.h>
#include <assert.h>

#include "vec.h"

char * _vec_init(long n, const mat_ctx_t ctx)
{
    char *vec;
    long i;

    assert(n > 0);

    vec = malloc(n * ctx->size);

    if (!vec)
    {
        printf("ERROR (_vec_init).\n\n");
        abort();
    }

    for (i = 0; i < n; i++)
        ctx->init(vec + i * ctx->size);

    return vec;
}

