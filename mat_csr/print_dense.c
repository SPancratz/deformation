#include "mat_csr.h"

int 
_mat_csr_print_dense(long m, long n, const char *x, const long *j, 
                                     const long *p, const long *lenr, 
                                     const mat_ctx_t ctx)
{
    long i, c, k;

    for (i = 0; i < m; i++)
    {
        printf("[");
        for (c = 0; c < n; c++)
        {
            for (k = p[i]; k < p[i] + lenr[i]; k++)
                if (j[k] == c)
                    break;
            printf(" ");
            if (k < p[i] + lenr[i])
                ctx->print(ctx, x + k * ctx->size);
            else
                printf("0");
        }
        if (c != n - 1)
            printf(" ");
        printf("]");
        if (i != m - 1)
            printf("\n");
    }

    return 1;
}

int mat_csr_print_dense(const mat_csr_t A, const mat_ctx_t ctx)
{
    return _mat_csr_print_dense(A->m, A->n, A->x, A->j, A->p, A->lenr, ctx);
}

