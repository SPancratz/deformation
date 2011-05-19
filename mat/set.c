#include <assert.h>

#include "mat.h"

void mat_set(mat_t mat1, const mat_t mat2, const ctx_t ctx)
{
    assert(mat1->m == mat2->m && mat1->n == mat2->n);

    if (mat1 != mat2)
    {
        long i, j;

        for (i = 0; i < mat1->m; i++)
            for (j = 0; j < mat1->n; j++)
                ctx->set(ctx, mat_entry(mat1, i, j, ctx), 
                         mat_entry(mat2, i, j, ctx));
    }
}

