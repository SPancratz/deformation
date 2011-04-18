#include "mat_coo.h"

#define N  1000

void mat_coo_randtest(mat_coo_t A, flint_rand_t state, double d, const mat_ctx_t ctx)
{
    long i, j;
    long f;
    long k = 0, u = 2 * sizeof(long) + ctx->size, z;

    d = FLINT_MAX(d, 0.0);
    d = FLINT_MIN(d, 1.0);
    f = N * d;
    z = d * A->m * A->n;

    mat_coo_zero(A, ctx);

    if (z == 0)
        return;

    mat_coo_fit_length(A, z, ctx);

    for (i = 0; z && i < A->m; i++)
        for (j = 0; z && j < A->n; j++)
            if (n_randint(state, N) <= f)
            {
                char *off = A->list + k * u;

                *(long *) (off) = i;
                *(long *) (off + sizeof(long)) = j;

                ctx->init(ctx, off + 2 * sizeof(long));
                ctx->randtest_not_zero(ctx, off + 2 * sizeof(long), state);

                k++;
                z--;
            }

    A->length = k;
}
