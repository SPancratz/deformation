#include "mat_dense.h"

void mat_dense_clear(mat_dense_t mat, const mat_ctx_t ctx)
{
    long i;

    for (i = 0; i < mat->m * mat->n; i++)
        ctx->clear(mat->entries + i * ctx->size);

    free(mat->entries);
    free(mat->rows);
}

