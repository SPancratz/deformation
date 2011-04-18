#include "mat.h"

void mat_clear(mat_t mat, const mat_ctx_t ctx)
{
    long i;

    for (i = 0; i < mat->m * mat->n; i++)
        ctx->clear(ctx, mat->entries + i * ctx->size);

    free(mat->entries);
    free(mat->rows);
}

