#include "mat.h"

int mat_is_zero(const mat_t mat, const mat_ctx_t ctx)
{
    long i;

    for (i = 0; i < mat->m * mat->n; i++)
        if (!ctx->is_zero(mat->entries + i * ctx->size))
            return 0;
    return 1;
}

