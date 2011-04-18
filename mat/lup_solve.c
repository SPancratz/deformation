#include <assert.h>

#include "mat.h"
#include "vec.h"

void _mat_lup_solve(char *x, char ** const rows, long m, long n, 
                    const long *pi, const char *b, const mat_ctx_t ctx)
{
    long i, j;
    char *t;

    assert(m == n);
    assert(x != b);

    t = _vec_init(1, ctx);

    /* Solve the lower unit-triangular system L y = P b */

    for (i = 0; i < m; i++)
    {
        ctx->set(ctx, x + i * ctx->size, b + pi[i] * ctx->size);
        for (j = 0; j < i; j++)
        {
            if (!ctx->is_zero(ctx, rows[i] + j * ctx->size))
            {
                ctx->mul(ctx, t, rows[i] + j * ctx->size, x + j * ctx->size);
                ctx->sub(ctx, x + i * ctx->size, x + i * ctx->size, t);
            }
        }
    }

    /* Solve the upper triangular system U x = y */

    for (i = m - 1; i >= 0; i--)
    {
        for (j = i + 1; j < m; j++)
        {
            if (!ctx->is_zero(ctx, rows[i] + j * ctx->size))
            {
                ctx->mul(ctx, t, rows[i] + j * ctx->size, x + j * ctx->size);
                ctx->sub(ctx, x + i * ctx->size, x + i * ctx->size, t);
            }
        }
        ctx->div(ctx, x + i * ctx->size, x + i * ctx->size, rows[i] + i * ctx->size);
    }

    _vec_clear(t, 1, ctx);
}

void mat_lup_solve(char *x, const mat_t mat, const long *pi, 
                         const char *b, const mat_ctx_t ctx)
{
    if (mat->m == 1 && mat->n == 1)
        ctx->div(ctx, x, b, (mat->rows)[0]);
    else
        _mat_lup_solve(x, mat->rows, mat->m, mat->n, pi, b, ctx);
}

