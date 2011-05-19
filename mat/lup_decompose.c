#include "mat.h"
#include "vec.h"

int 
_mat_lup_decompose(long *pi, char **rows, long m, const ctx_t ctx)
{
    long i, j, k, t;
    char *p, *x;

    if (m == 1)
    {
        pi[0] = 0;
        return 0;
    }

    x = _vec_init(1, ctx);

    for (i = 0; i < m; i++)
        pi[i] = i;

    for (k = 0; k < m - 1; k++)
    {
        /*
            Find a pivot.  We choose the non-zero one in column k 
            in the row with the largest number of zero elements
         */
        {
            long nz, nz_best = -1, i_best = -1;

            for (i = k; i < m; i++)
            {
                if (ctx->is_zero(ctx, rows[i] + k * ctx->size))
                    continue;

                nz = 0;
                for (j = k; j < m; j++)
                    if (ctx->is_zero(ctx, rows[i] + j * ctx->size))
                        nz++;
                if (nz > nz_best)
                {
                    nz_best = nz;
                    i_best  = i;
                }
            }
            i = i_best;
        }

        if (i < 0)
            return 1;

        t     = pi[k];
        pi[k] = pi[i];
        pi[i] = t;
        p       = rows[k];
        rows[k] = rows[i];
        rows[i] = p;

        for (i = k + 1; i < m; i++)
        {
            if (ctx->is_zero(ctx, rows[i] + k * ctx->size))
                continue;

            ctx->div(ctx, rows[i] + k * ctx->size, 
                     rows[i] + k * ctx->size, 
                     rows[k] + k * ctx->size);
            for (j = k + 1; j < m; j++)
            {
                if (ctx->is_zero(ctx, rows[k] + j * ctx->size))
                    continue;

                ctx->mul(ctx, x, rows[i] + k * ctx->size, 
                            rows[k] + j * ctx->size);
                ctx->sub(ctx, rows[i] + j * ctx->size, 
                         rows[i] + j * ctx->size, x);
            }
        }
    }

    _vec_clear(x, 1, ctx);

    return 0;
}

int mat_lup_decompose(mat_t out, long *pi, const mat_t mat, 
                            const ctx_t ctx)
{
    if (!(mat->m == out->m && mat->n && out->n && mat->m == mat->n))
    {
        printf("ERROR (mat_lup_decompose).\n\n");
        abort();
    }

    mat_set(out, mat, ctx);

    return _mat_lup_decompose(pi, out->rows, out->m, ctx);
}

