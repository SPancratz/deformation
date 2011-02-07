#include <assert.h>

#include "mat_dense.h"

int mat_dense_equal(const mat_dense_t mat1, const mat_dense_t mat2, const mat_ctx_t ctx)
{
    assert(mat1->m == mat2->m && mat1->n == mat2->n);

    if (mat1 != mat2)
    {
        long i, j;

        for (i = 0; i < mat1->m; i++)
            for (j = 0; j < mat1->n; j++)
                if (!ctx->equal(mat_dense_entry(mat1, i, j, ctx), 
                                mat_dense_entry(mat2, i, j, ctx)))
                    return 0;
    }

    return 1;
}

