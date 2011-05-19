#include <assert.h>

#include "mat.h"

void mat_init(mat_t mat, long m, long n, const ctx_t ctx)
{
    long i;

    assert(m > 0 && n > 0);

    mat->entries = malloc(m * n * ctx->size);
    mat->rows    = malloc(m * sizeof(char *));

    for (i = 0; i < m * n; i++)
        ctx->init(ctx, mat->entries + i * ctx->size);

    for (i = 0; i < m; i++)
        mat->rows[i] = mat->entries + i * n * ctx->size;

    mat->m = m;
    mat->n = n;
}

