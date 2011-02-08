#include "mat_dense.h"

int 
_mat_dense_lup_decompose(long *pi, char **rows, long m, const mat_ctx_t ctx)
{
    long i, j, k, t;
    char *p, *x;

    x = malloc(ctx->size);
    ctx->init(x);

    for (i = 0; i < m; i++)
        pi[i] = i;

    for (k = 0; k < m - 1; k++)
    {
        for (i = k; i < m; i++)
            if (!ctx->is_zero(rows[i] + k * ctx->size))
                break;
        if (i == m)
            return 1;
        t     = pi[k];
        pi[k] = pi[i];
        pi[i] = t;
        p       = rows[k];
        rows[k] = rows[i];
        rows[i] = p;
        for (i = k + 1; i < m; i++)
        {
            ctx->div(rows[i] + k * ctx->size, 
                     rows[i] + k * ctx->size, 
                     rows[k] + k * ctx->size);
            for (j = k + 1; j < m; j++)
            {
                ctx->mul(x, rows[i] + k * ctx->size, 
                            rows[k] + j * ctx->size);
                ctx->sub(rows[i] + j * ctx->size, 
                         rows[i] + j * ctx->size, x);
            }
        }
    }

    ctx->clear(x);
    free(x);

    return 0;
}

int mat_dense_lup_decompose(mat_dense_t out, long *pi, mat_dense_t mat, 
                            const mat_ctx_t ctx)
{
    if (!(mat->m == out->m && mat->n && out->n))
    {
        printf("ERROR (mat_dense_lup_decompose).\n\n");
        abort();
    }

    mat_dense_set(out, mat, ctx);

    return _mat_dense_lup_decompose(pi, mat->rows, mat->m, ctx);
}

