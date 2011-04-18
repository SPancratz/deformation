#include "mat.h"

int mat_is_one(const mat_t mat, const mat_ctx_t ctx)
{
    long i, j;

    if (mat->m != mat->n)
        return 0;

    for (i = 0; i < mat->m; i++)
    {
        for (j = 0; j < i; j++)
            if (!ctx->is_zero(ctx, mat_entry(mat, i, j, ctx)))
                return 0;
        if (!ctx->is_one(ctx, mat_entry(mat, i, j, ctx)))
            return 0;
        for (j++; j < mat->n; j++)
            if (!ctx->is_zero(ctx, mat_entry(mat, i, j, ctx)))
                return 0;
    }
    return 1;
}

