#include <stdlib.h>

#include "vec.h"

#include "mat_csr.h"

void mat_csr_solve_clear(mat_csr_solve_t s, const ctx_t ctx)
{
    long k, len, sum = 0;

    for (k = 0; k < s->nb; k++)
    {
        len  = s->B[k + 1] - s->B[k];
        sum += len * len;
    }

    _vec_clear(s->entries, sum, ctx);

    free(s->j);
    free(s->LU);
}

