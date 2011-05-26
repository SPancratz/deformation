#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "vec.h"
#include "mat.h"

long mat_inv(mat_t B, const mat_t A, const ctx_t ctx)
{
    long ans;

    mat_t LU, Linv, Uinv;
    long *pi;
    long i, j, k, n = A->m;
    char *t;

    assert(B->m == B->n && A->m == A->n && A->m == B->m);

    mat_init(LU, n, n, ctx);
    mat_init(Linv, n, n, ctx);
    mat_init(Uinv, n, n, ctx);
    t = _vec_init(1, ctx);

    pi = malloc(n * sizeof(long));

    ans = mat_lup_decompose(LU, pi, A, ctx);
    if (ans)
    {
        goto exit;   
    }

    for (j = 0; j < n; j++)
        for (i = 0; i < n; i++)
        {
            if (i == j)
            {
                ctx->one(ctx, mat_entry(Linv, i, j, ctx));
            }
            else
            {
                ctx->zero(ctx, mat_entry(Linv, i, j, ctx));
            }
            for (k = 0; k < i; k++)
            {
                ctx->mul(ctx, t, mat_entry(LU, i, k, ctx), 
                                 mat_entry(Linv, k, j, ctx));
                ctx->sub(ctx, mat_entry(Linv, i, j, ctx), 
                              mat_entry(Linv, i, j, ctx), t);
            }
        }

    for (j = 0; j < n; j++)
        for (i = n - 1; i >= 0; i--)
        {
            if (i == j)
            {
                ctx->one(ctx, mat_entry(Uinv, i, j, ctx));
            }
            else
            {
                ctx->zero(ctx, mat_entry(Uinv, i, j, ctx));
            }
            for (k = i + 1; k < n; k++)
            {
                ctx->mul(ctx, t, mat_entry(LU, i, k, ctx), 
                                 mat_entry(Uinv, k, j, ctx));
                ctx->sub(ctx, mat_entry(Uinv, i, j, ctx), 
                              mat_entry(Uinv, i, j, ctx), t);
            }
            ctx->div(ctx, mat_entry(Uinv, i, j, ctx), 
                          mat_entry(Uinv, i, j, ctx), 
                          mat_entry(LU, i, i, ctx));
        }

    mat_mul(B, Uinv, Linv, ctx);

    /* B[i,pi[j]] = B[i,j] */
    {
        char *w = malloc(n * ctx->size);

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
                memcpy(w + j * ctx->size, mat_entry(B, i, j, ctx), ctx->size);

            for (j = 0; j < n; j++)
                memcpy(mat_entry(B, i, pi[j], ctx), w + j * ctx->size, ctx->size);
        }

        free(w);
    }

  exit: 

    mat_clear(LU, ctx);
    mat_clear(Linv, ctx);
    mat_clear(Uinv, ctx);
    _vec_clear(t, 1, ctx);
    free(pi);

    return ans;
}

