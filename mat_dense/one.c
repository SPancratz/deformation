#include <assert.h>

#include "mat_dense.h"

void mat_dense_one(mat_dense_t mat, const mat_ctx_t ctx)
{
    long i, j;

    assert(mat->m == mat->n);

    for (i = 0; i < mat->m; i++)
    {
        for (j = 0; j < i; j++)
        {
            ctx->zero(mat_dense_entry(mat, i, j, ctx));
            ctx->zero(mat_dense_entry(mat, j, i, ctx));
        }
        ctx->one(mat_dense_entry(mat, i, i, ctx));
    }
}

