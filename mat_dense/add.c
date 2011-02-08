#include "mat_dense.h"

void _mat_dense_add(char **rowsC, char ** const rowsA, char ** const rowsB, 
                    long m, long n, const mat_ctx_t ctx)
{
    long i, j;

    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
            ctx->add(rowsC[i] + j * ctx->size, 
                     rowsA[i] + j * ctx->size, 
                     rowsB[i] + j * ctx->size);
}

void mat_dense_add(mat_dense_t C, const mat_dense_t A, const mat_dense_t B, 
                   const mat_ctx_t ctx)
{
    if (!(A->m == B->m && A->m == C->m && A->n == B->n && A->n == C->n))
    {
        printf("ERROR (mat_dense_add).\n\n");
        abort();
    }

    _mat_dense_add(C->rows, A->rows, B->rows, A->m, A->n, ctx);
}

