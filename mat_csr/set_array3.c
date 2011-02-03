#include <string.h>

#include "mat_csr.h"

void mat_csr_set_array3(mat_csr_t A, char *mem, long len, int copy, const mat_ctx_t ctx)
{
    long i, j, k, u;
    long *lenr;

    /* Also sets A->lenr[i] to zero */
    mat_csr_zero(A, ctx);

    if (len == 0)
        return;

    u = 2 * sizeof(long) + ctx->size;

    lenr = calloc(A->m, sizeof(long));

    if (!lenr)
    {
        printf("ERROR (mat_csr_set_array3).\n\n");
        abort();
    }

    mat_csr_fit_length(A, len, ctx);

    for (k = 0; k < len; k++)
    {
        if (!ctx->is_zero(mem + k * u + 2 * sizeof(long)))
        {
            i = *(long *) (mem + k * u);
            lenr[i] ++;
        }
    }

    A->p[0] = 0;
    for (i = 1; i < A->m; i++)
        A->p[i] = A->p[i - 1] + lenr[i - 1];

    for (k = 0; k < len; k++)
    {
        char *x, *y = mem + k * u + 2 * sizeof(long);

        if (!ctx->is_zero(y))
        {
            i = *(long *) (mem + k * u);
            j = *(long *) (mem + k * u + sizeof(long));
            x = A->x + (A->p[i] + A->lenr[i]) * ctx->size;

            if (copy)
            {
                ctx->set(x, y);
            }
            else
            {
                ctx->clear(x);
                memcpy(x, y, ctx->size);
            }
            A->j[A->p[i] + A->lenr[i]] = j;
            A->lenr[i] ++;
        }
        else
            if (!copy)
                ctx->clear(y);
    }

    free(lenr);
}
