#include "mat.h"

void mat_randrank(mat_t mat, flint_rand_t state, long rank, 
                        const ctx_t ctx)
{
    long i;

    mat_zero(mat, ctx);

    for (i = 0; i < rank; i++)
        ctx->randtest_not_zero(mat_entry(mat, i, i, ctx), state);
}
