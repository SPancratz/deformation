#include <stdlib.h>
#include <string.h>

#include "mat_csr.h"

static 
void _mat_csr_sort_row(char *x, long *j, long len, const mat_ctx_t ctx)
{
    long i, k, key;
    char *y;

    if (len == 1)
        return;

    y = malloc(ctx->size);
    if (!y)
    {
        printf("ERROR (_mat_csr_sort_row).\n\n");
        abort();
    }

    for (i = 1; i < len; i++)
    {
        key = j[i];
        memcpy(y, x + i * ctx->size, ctx->size);

        for (k = i - 1; (k >= 0) && (j[k] > key); k--)
        {
            j[k + 1] = j[k];
            memcpy(x + (k + 1) * ctx->size, x + k * ctx->size, ctx->size);
        }
        j[k + 1] = key;
        memcpy(x + (k + 1) * ctx->size, y, ctx->size);
    }

    free(y);
}

void _mat_csr_sort_rows(long m, char *x, long *j, long *p, long *lenr, 
                        const mat_ctx_t ctx)
{
    long i;

    for (i = 0; i < m; i++)
        _mat_csr_sort_row(x + p[i] * ctx->size, j + p[i], lenr[i], ctx);
}

void mat_csr_sort_rows(mat_csr_t A, const mat_ctx_t ctx)
{
    _mat_csr_sort_rows(A->m, A->x, A->j, A->p, A->lenr, ctx);
}

