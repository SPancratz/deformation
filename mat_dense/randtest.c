#include "mat_dense.h"

void mat_dense_randtest(mat_dense_t mat, flint_rand_t state, 
                        const mat_ctx_t ctx)
{
    long i;

    for (i = 0; i < mat->m * mat->n; i++)
        ctx->randtest(mat->entries + i * ctx->size, state);
}

