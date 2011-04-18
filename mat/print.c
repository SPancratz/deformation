#include <stdlib.h>
#include <stdio.h>

#include "mat.h"

int _mat_print(char ** const rows, long m, long n, const mat_ctx_t ctx)
{
    long i, j;

    for (i = 0; i < m; i++)
    {
        printf("[ ");
        for (j = 0; j < n; j++)
        {
            ctx->print(ctx, rows[i] + j * ctx->size);
            printf(" ");
        }
        printf("]");
        if (i != m - 1)
            printf("\n");
    }

    return 1;
}

int mat_print(const mat_t mat, const mat_ctx_t ctx)
{
    return _mat_print(mat->rows, mat->m, mat->n, ctx);
}

