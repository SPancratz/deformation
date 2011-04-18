#include <stdlib.h>
#include <assert.h>

#include "vec.h"

void _vec_clear(char *vec, long n, const mat_ctx_t ctx)
{
    long i;

    assert(n > 0);

    for (i = 0; i < n; i++)
        ctx->clear(ctx, vec + i * ctx->size);

    free(vec);
}

