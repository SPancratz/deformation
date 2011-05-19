#include "mat.h"

void _mat_add(char **rowsC, char ** const rowsA, char ** const rowsB, 
                    long m, long n, const ctx_t ctx)
{
    long i, j;

    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
            ctx->add(ctx, rowsC[i] + j * ctx->size, 
                     rowsA[i] + j * ctx->size, 
                     rowsB[i] + j * ctx->size);
}

void mat_add(mat_t C, const mat_t A, const mat_t B, 
                   const ctx_t ctx)
{
    if (!(A->m == B->m && A->m == C->m && A->n == B->n && A->n == C->n))
    {
        printf("ERROR (mat_add).\n\n");
        abort();
    }

    _mat_add(C->rows, A->rows, B->rows, A->m, A->n, ctx);
}

