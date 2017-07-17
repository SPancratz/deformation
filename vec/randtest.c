#include <assert.h>

#include "vec.h"

#include "flint/flint.h"
#include "flint/ulong_extras.h"

void _vec_randtest(char *vec, long n, flint_rand_t state, const ctx_t ctx)
{
    long i;

    for (i = 0; i < n; i++)
        ctx->randtest(ctx, vec + i * ctx->size, state);
}

void _vec_randtest_not_zero(char *vec, long n, flint_rand_t state, const ctx_t ctx)
{
    long i;

    assert(n > 0);

    for (i = 0; i < n; i++)
        ctx->randtest(ctx, vec + i * ctx->size, state);

    i = n_randint(state, n);
    ctx->randtest_not_zero(ctx, vec + i * ctx->size, state);
}

