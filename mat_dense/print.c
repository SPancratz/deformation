#include <stdlib.h>
#include <stdio.h>

#include "mat_dense.h"

int mat_dense_print(const mat_dense_t mat, const mat_ctx_t ctx)
{
    long i, j;

    for (i = 0; i < mat->m; i++)
    {
        printf("[ ");
        for (j = 0; j < mat->n; j++)
        {
            ctx->print(mat_dense_entry(mat, i, j, ctx));
            printf(" ");
        }
        printf("]");
        if (i != mat->m - 1)
            printf("\n");
    }

    return 1;
}

