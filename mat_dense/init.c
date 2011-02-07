#include <assert.h>

#include "mat_dense.h"

void mat_dense_init(mat_dense_t mat, long m, long n, const mat_ctx_t ctx)
{
    long i;

    assert(m > 0 && n > 0);

    mat->entries = malloc(m * n * ctx->size);
    mat->rows    = malloc(m * sizeof(char *));

    for (i = 0; i < m * n; i++)
        ctx->init(mat->entries + i * ctx->size);

    for (i = 0; i < m; i++)
        mat->rows[i] = mat->entries + i * n * ctx->size;

    mat->m = m;
    mat->n = n;
}

