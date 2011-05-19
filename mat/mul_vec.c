#include "mat.h"

void _mat_mul_vec(char *y, char ** const rows, long m, long n, 
                        const char *x, const ctx_t ctx)
{
    long i, j;
    char *t;

    t = malloc(ctx->size);

    if (!t)
    {
        printf("ERROR (_mat_mul_vec).\n\n");
        abort();
    }

    ctx->init(ctx, t);

    for (i = 0; i < m; i++)
    {
        ctx->zero(ctx, y + i * ctx->size);
        for (j = 0; j < n; j++)
        {
            ctx->mul(ctx, t, rows[i] + j * ctx->size, x + j * ctx->size);
            ctx->add(ctx, y + i * ctx->size, y + i * ctx->size, t);
        }
    }

    ctx->clear(ctx, t);
    free(t);
}

void mat_mul_vec(char *y, const mat_t A, const char *x, 
                       const ctx_t ctx)
{
    _mat_mul_vec(y, A->rows, A->m, A->n, x, ctx);
}

