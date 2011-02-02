#include "mat_coo.h"

void mat_coo_fit_length(mat_coo_t A, long len, const mat_ctx_t ctx)
{
    if (len > A->alloc)
    {
        if (len < 2 * A->alloc)
            len = 2 * A->alloc;
        mat_coo_realloc(A, len, ctx);
    }
}

