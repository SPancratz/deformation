#include "mat_dense.h"

void mat_dense_randrank(mat_dense_t mat, flint_rand_t state, long rank, 
                        const mat_ctx_t ctx)
{
    long i;

    mat_dense_zero(mat, ctx);

    for (i = 0; i < rank; i++)
        ctx->randtest_not_zero(mat_dense_entry(mat, i, i, ctx), state);
}
