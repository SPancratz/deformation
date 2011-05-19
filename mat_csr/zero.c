#include "mat_csr.h"

void mat_csr_zero(mat_csr_t A, const ctx_t ctx)
{
    long i;

    for (i = 0; i < A->m; i++)
        A->lenr[i] = 0;
}

