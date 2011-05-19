#include "mat_coo.h"

void mat_coo_clear(mat_coo_t A, int clear, const ctx_t ctx)
{
    if (A->alloc)
    {
        if (clear)
        {
            long k, u = 2 * sizeof(long) + ctx->size;

            for (k = 0; k < A->length; k++)
                ctx->clear(ctx, A->list + k * u + 2 * sizeof(long));
        }

        A->alloc  = 0;
        A->length = 0;
        free(A->list);
    }
}

