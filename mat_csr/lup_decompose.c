#include <assert.h>

#include "mat_csr.h"

void 
mat_csr_lupdecompose(mat_csr_t L, mat_csr_t U, long *pi, const mat_csr_t A, 
                     const mat_ctx_t ctx)
{
    char *mem, *x;
    long i, j, k, n, q;

    assert(A->m == A->n);
    assert(L->m == A->m && L->n && A->n);
    assert(U->m == A->m && U->n && A->n);

    n = A->m;
    mem = malloc(n * n * ctx->size);

    for (k = 0; k < n * n; k++)
        ctx->zero(mem + k * ctx->size);

    for (i = 0; i < n; i++)
        for (q = A->p[i]; q < A->p[i] + A->lenr[i]; q++)
        {
            j = A->j[q];
            x = A->x + (q * ctx->size);
            ctx->set(mem + (i * n + j) * ctx->size, x);
        }

    /* Dense LUP decomposition */
    x = malloc(ctx->size);
    ctx->init(x);
    for (i = 0; i < n; i++)
        pi[i] = i;
    for (k = 0; k < n - 1; k++)
    {
        for (i = k; i < n; i++)
            if (!ctx->is_zero(mem + (i * n + k) * ctx->size))
                break;
        assert(i < n);
        if (i != k)
        {
            long t;

            t     = pi[i];
            pi[i] = pi[k];
            pi[k] = t;

            for (j = 0; j < n; j++)
                ctx->swap(mem + (k * n + j) * ctx->size, 
                          mem + (i * n + j) * ctx->size);
        }
        for (i = k + 1; i < n; i++)
        {
            ctx->div(mem + (i * n + k) * ctx->size, 
                     mem + (i * n + k) * ctx->size, 
                     mem + (k * n + k) * ctx->size);
            for (j = k + 1; j < n; j++)
            {
                ctx->mul(x, mem + (i * n + k) * ctx->size, 
                            mem + (k * n + j) * ctx->size);
                ctx->sub(mem + (i * n + j) * ctx->size, 
                         mem + (i * n + j) * ctx->size, x);
            }
        }
    }
    ctx->clear(x);
    free(x);



    /* No `clearing' of individual coefficients */
    free(mem);
}

