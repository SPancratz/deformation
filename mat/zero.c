#include "mat.h"

void mat_zero(mat_t mat, const ctx_t ctx)
{
    long i;

    for (i = 0; i < mat->m * mat->n; i++)
        ctx->zero(mat->entries + i * ctx->size);
}
