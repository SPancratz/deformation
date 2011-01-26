#include <string.h>

#include "mat_csr.h"

/*
    Sets the matrix A to matrix \code{(mem, len)} given in triplet form. 
    Each entry in the matrix has width \code{2 * sizeof(long) + ctx->size}.

    If \code{copy} evaluates to true, the entry values are deep-copied 
    using the context's \code{set()} function.  Otherwise, they are only 
    shallow copies and it is safe for the caller to free \code{mem}.

    Assumes the co-ordinates in \code{(mem, len)} are distinct.
 */

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
        if (!ctx->is_zero(mem + k * u + 2 * sizeof(long)))
        {
            i = *(long *) (mem + k * u);
            j = *(long *) (mem + k * u + sizeof(long));

            if (copy)
            {
                ctx->set(A->x + (A->p[i] + A->lenr[i]) * ctx->size, 
                         mem + k * u + 2 * sizeof(long));
            }
            else
            {
                memcpy(A->x + (A->p[i] + A->lenr[i]) * ctx->size, 
                       mem + k * u + 2 * sizeof(long), ctx->size);
            }
            A->j[A->p[i] + A->lenr[i]] = j;
            A->lenr[i] ++;
        }
    }

    free(lenr);
}
