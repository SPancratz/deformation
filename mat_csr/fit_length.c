#include "mat_csr.h"

void mat_csr_fit_length(mat_csr_t A, long len, const mat_ctx_t ctx)
{
    if (len > A->alloc)
    {
        if (len < 2 * A->alloc)
            len = 2 * A->alloc;
        mat_csr_realloc(A, len, ctx);
    }
}

