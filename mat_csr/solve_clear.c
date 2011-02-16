#include <stdlib.h>

#include "vec.h"

#include "mat_csr.h"

void mat_csr_solve_clear(mat_csr_solve_t s, const mat_ctx_t ctx)
{
    free(s->pi);
    free(s->LU);

    _vec_clear(s->entries, s->alloc, ctx);
}

