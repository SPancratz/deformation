#include "mat_dense.h"

void mat_dense_zero(mat_dense_t mat, const mat_ctx_t ctx)
{
    long i;

    for (i = 0; i < mat->m * mat->n; i++)
        ctx->zero(mat->entries + i * ctx->size);
}
