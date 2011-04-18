#include "mat_coo.h"

void mat_coo_zero(mat_coo_t A, const mat_ctx_t ctx)
{
    long k, u = 2 * sizeof(long) + ctx->size;

    for (k = 0; k < A->length; k++)
        ctx->clear(ctx, A->list + k * u + 2 * sizeof(long));

    A->length = 0;
}

