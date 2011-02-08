#include <assert.h>

#include "mat_dense.h"

void _mat_dense_lup_solve(char *x, char ** const rows, long m, long n, 
                          const long *pi, const char *b, const mat_ctx_t ctx)
{
    long i, j;
    char *t;

    assert(m == n);

    t = malloc(ctx->size);

    if (!t)
    {
        printf("ERROR (_mat_dense_lup_solve).\n\n");
        abort();
    }

    ctx->init(t);

    /* Solve the lower unit-triangular system L y = P b */

    for (i = 0; i < m; i++)
    {
        ctx->set(x + i * ctx->size, b + pi[i] * ctx->size);
        for (j = 0; j < i; j++)
        {
            ctx->mul(t, rows[i] + j * ctx->size, x + j * ctx->size);
            ctx->sub(x + i * ctx->size, x + i * ctx->size, t);
        }
    }

    /* Solve the upper triangular system U x = y */

    for (j = m - 1; j >= 0; j--)
    {
        for (i = 0; i < j; i++)
        {
            ctx->mul(t, rows[i] + j * ctx->size, x + j * ctx->size);
            ctx->sub(x + i * ctx->size, x + i * ctx->size, t);
        }
        ctx->div(x + j * ctx->size, x + j * ctx->size, rows[j] + j * ctx->size);
    }

    ctx->clear(t);
    free(t);
}

void mat_dense_lup_solve(char *x, const mat_dense_t mat, const long *pi, 
                         const char *b, const mat_ctx_t ctx)
{
    _mat_dense_lup_solve(x, mat->rows, mat->m, mat->n, pi, b, ctx);
}

