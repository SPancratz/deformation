#include <assert.h>

#include "mat_coo.h"

void 
mat_coo_set_entry(mat_coo_t A, long i, long j, const void *x, 
                  const mat_ctx_t ctx)
{
    long k, u = 2 * sizeof(long) + ctx->size;

    assert(!ctx->is_zero(x));

    for (k = 0; k < A->length; k++)
        if (   i == *(long *) (A->list + k * u) 
            && j == *(long *) (A->list + k * u + sizeof(long)))
            break;

    if (k < A->length)
        ctx->set(A->list + k * u + 2 * sizeof(long), x);
    else
    {
        mat_coo_fit_length(A, A->length + 1, ctx);
        A->length += 1;
        *(long *) (A->list + k * u) = i;
        *(long *) (A->list + k * u + sizeof(long)) = j;
        ctx->set(A->list + k * u + 2 * sizeof(long), x);
    }
}
