#include "mat_csr.h"

int mat_csr_is_zero(const mat_csr_t A, const ctx_t ctx)
{
    long i;

    for (i = 0; i < A->m; i++)
        if (A->lenr[i])
            return 0;
    return 1;
}

