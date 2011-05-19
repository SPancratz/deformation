#include "mat_coo.h"

int 
_mat_coo_print_dense(long m, long n, const char *list, long len, 
                     const ctx_t ctx)
{
    long i, j, k, u = 2 * sizeof(long) + ctx->size;
    char *x;

    for (i = 0; i < m; i++)
    {
        printf("[");
        for (j = 0; j < n; j++)
        {
            x = NULL;
            for (k = 0; k < len; k++)
                if (   i == *(long *) (list + k * u) 
                    && j == *(long *) (list + k * u + sizeof(long)))
                {
                    x = (char *) list + k * u + 2 * sizeof(long);
                }
            printf(" ");
            if (x)
                ctx->print(x);
            else
                printf("0");
        }
        printf(" ]");
        if (i != m - 1)
            printf("\n");
    }

    return 1;
}

int 
mat_coo_print_dense(const mat_coo_t A, const ctx_t ctx)
{
    return _mat_coo_print_dense(A->m, A->n, A->list, A->length, ctx);
}
