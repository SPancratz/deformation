#include "mat_csr.h"

int mat_csr_debug(const mat_csr_t A, const mat_ctx_t ctx)
{
    long i, k;

    printf("m     = %ld\n", A->m);
    printf("n     = %ld\n", A->n);
    printf("alloc = %ld\n", A->alloc);

    printf("x     = [ ");
    for (k = 0; k < A->alloc; k++)
    {
        ctx->print(ctx, A->x + k * (ctx->size));
        printf(" ");
    }
    printf("]\n");

    printf("j     = [ ");
    for (k = 0; k < A->alloc; k++)
    {
        printf("%ld ", A->j[k]);
    }
    printf("]\n");

    printf("p     = [ ");
    for (i = 0; i < A->m; i++)
    {
        printf("%ld ", A->p[i]);
    }
    printf("]\n");

    printf("lenr  = [ ");
    for (i = 0; i < A->m; i++)
    {
        printf("%ld ", A->j[i]);
    }
    printf("]\n");

    return 1;
}

