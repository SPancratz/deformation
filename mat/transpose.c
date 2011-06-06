#include "mat.h"

void _mat_transpose(char **rowsB, char ** const rowsA, long m, long n, 
                    const ctx_t ctx)
{
    long i, j;

    if (rowsB == rowsA)
    {
        for (i = 0; i < m; i++)
            for (j = 0; j < i; j++)
                ctx->swap(ctx, 
                          rowsB[j] + i * ctx->size, rowsB[i] + j * ctx->size);
    }
    else
    {
        for (i = 0; i < m; i++)
            for (j = 0; j < n; j++)
                ctx->set(ctx, 
                         rowsB[j] + i * ctx->size, rowsA[i] + j * ctx->size);
    }
}

void mat_transpose(mat_t B, const mat_t A, const ctx_t ctx)
{
    if (!(A->m == B->n && A->n == B->m))
    {
        printf("ERROR (mat_transpose).\n\n");
        abort();
    }

    _mat_transpose(B->rows, A->rows, A->m, A->n, ctx);
}

