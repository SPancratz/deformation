#include <assert.h>

#include "mat.h"

void mat_one(mat_t mat, const ctx_t ctx)
{
    long i, j;

    assert(mat->m == mat->n);

    for (i = 0; i < mat->m; i++)
    {
        for (j = 0; j < i; j++)
        {
            ctx->zero(mat_entry(mat, i, j, ctx));
            ctx->zero(mat_entry(mat, j, i, ctx));
        }
        ctx->one(mat_entry(mat, i, i, ctx));
    }
}

