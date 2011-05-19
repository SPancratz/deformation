#include "mat.h"

void mat_randtest(mat_t mat, flint_rand_t state, 
                        const ctx_t ctx)
{
    long i;

    for (i = 0; i < mat->m * mat->n; i++)
        ctx->randtest(ctx, mat->entries + i * ctx->size, state);
}

