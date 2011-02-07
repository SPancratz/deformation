#include "mat_dense.h"

int mat_dense_is_one(const mat_dense_t mat, const mat_ctx_t ctx)
{
    long i, j;

    if (mat->m != mat->n)
        return 0;

    for (i = 0; i < mat->m; i++)
    {
        for (j = 0; j < i; j++)
            if (!ctx->is_zero(mat_dense_entry(mat, i, j, ctx)))
                return 0;
        if (!ctx->is_one(mat_dense_entry(mat, i, j, ctx)))
            return 0;
        for (j++; j < mat->n; j++)
            if (!ctx->is_zero(mat_dense_entry(mat, i, j, ctx)))
                return 0;
    }
    return 1;
}

