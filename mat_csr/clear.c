#include "mat_csr.h"

void mat_csr_clear(mat_csr_t A, const mat_ctx_t ctx)
{
    free(A->p);
    free(A->lenr);

    if (A->alloc)
    {
        long k;

        for (k = 0; k < A->alloc; k++)
            ctx->clear(A->x + k * (ctx->size));

        A->alloc = 0;
        free(A->x);
        free(A->j);
    }
}
