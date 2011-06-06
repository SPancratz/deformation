#include "mat.h"

void _mat_neg(char **rowsB, char ** const rowsA, long m, long n, 
              const ctx_t ctx)
{
    long i, j;

    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
            ctx->neg(ctx, rowsB[i] + j * ctx->size, 
                          rowsA[i] + j * ctx->size);
}

void mat_neg(mat_t B, const mat_t A, const ctx_t ctx)
{
    if (!(A->m == B->m && A->n == B->n))
    {
        printf("ERROR (mat_neg).\n\n");
        abort();
    }

    _mat_neg(B->rows, A->rows, A->m, A->n, ctx);
}

