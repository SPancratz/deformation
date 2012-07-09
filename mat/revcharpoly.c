#include <assert.h>

#include "mat.h"
#include "vec.h"

void mat_revcharpoly(char *poly, mat_t mat, const ctx_t ctx)
{
    char *a, *A;
    char *r, *s;
    long i, j, k, t;
    long n = mat->n;

    assert(mat->m == mat->n);

    if (n == 0)
    {
        ctx->one(ctx, _vec_entry(poly, 0, ctx));
        return;
    }
    if (n == 1)
    {
        ctx->one(ctx, _vec_entry(poly, 0, ctx));
        ctx->neg(ctx, _vec_entry(poly, 1, ctx), mat_entry(mat, 0, 0, ctx));
        return;
    }

    a = _vec_init(n * n, ctx);
    A = _vec_entry(a, (n - 1) * n, ctx);

    r = _vec_init(1, ctx);
    s = _vec_init(1, ctx);

    ctx->neg(ctx, _vec_entry(poly, 0, ctx), mat_entry(mat, 0, 0, ctx));

    for (t = 1; t < n; t++)
    {
        for (i = 0; i <= t; i++)
        {
            ctx->set(ctx, _vec_entry(a, 0 * n + i, ctx), 
                          mat_entry(mat, i, t, ctx));
        }

        ctx->set(ctx, _vec_entry(A, 0, ctx), mat_entry(mat, t, t, ctx));

        for (k = 1; k < t; k++)
        {
            for (i = 0; i <= t; i++)
            {
                ctx->zero(ctx, s);
                for (j = 0; j <= t; j++)
                {
                    ctx->mul(ctx, r, mat_entry(mat, i, j, ctx), 
                                     _vec_entry(a, (k - 1) * n + j, ctx));
                    ctx->add(ctx, s, s, r);
                }
                ctx->set(ctx, _vec_entry(a, k * n + i, ctx), s);
            }

            ctx->set(ctx, _vec_entry(A, k, ctx), _vec_entry(a, k * n + t, ctx));
        }

        ctx->zero(ctx, s);
        for (j = 0; j <= t; j++)
        {
            ctx->mul(ctx, r, mat_entry(mat, t, j, ctx), 
                             _vec_entry(a, (t - 1) * n + j, ctx));
            ctx->add(ctx, s, s, r);
        }
        ctx->set(ctx, _vec_entry(A, t, ctx), s);

        for (k = 0; k <= t; k++)
        {
            ctx->set(ctx, s, _vec_entry(poly, k, ctx));
            for (j = 0; j < k; j++)
            {
                ctx->mul(ctx, r, _vec_entry(A, j, ctx), _vec_entry(poly, k - j - 1, ctx));
                ctx->sub(ctx, s, s, r);
            }
            ctx->sub(ctx, _vec_entry(poly, k, ctx), s, _vec_entry(A, k, ctx));
        }
    }

    /* Shift all coefficients up by one */

    for (i = n; i > 0; i--)
    {
        ctx->swap(ctx, _vec_entry(poly, i, ctx), _vec_entry(poly, i - 1, ctx));
    }
    ctx->one(ctx, _vec_entry(poly, 0, ctx));

    _vec_clear(r, 1, ctx);
    _vec_clear(s, 1, ctx);
    _vec_clear(a, n * n, ctx);
}

