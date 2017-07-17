#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>

#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/ulong_extras.h"

#include "mpoly.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;
    ctx_t ctx;

    printf("set/equal... ");
    fflush(stdout);

    _randinit(state);

    ctx_init_mpq(ctx);

    /* Equal polynomials */
    for (i = 0; i < 1000; i++)
    {
        mpoly_t a, b;
        long n, d, N;

        n = n_randint(state, MON_MAX_VARS) + 1;
        d = n_randint(state, 50) + 1;
        N = n_randint(state, 50) + 1;

        mpoly_init(a, n, ctx);
        mpoly_init(b, n, ctx);
        mpoly_randtest(a, state, d, N, ctx);

        mpoly_set(b, a, ctx);

        result = mpoly_equal(a, b, ctx);
        if (!result)
        {
            printf("FAIL:\n");
            mpoly_print(a, ctx); printf("\n");
            mpoly_print(b, ctx); printf("\n");
            abort();
        }

        mpoly_clear(a, ctx);
        mpoly_clear(b, ctx);
    }

    /* Non-equal polynomials */
    for (i = 0; i < 1000; i++)
    {
        mpoly_t a, b;
        long n, d, N;
        mon_t m;
        char *x;

        n = n_randint(state, MON_MAX_VARS) + 1;
        d = n_randint(state, 50) + 1;
        N = n_randint(state, 50) + 1;

        mon_init(m);
        mon_randtest(m, state, n, d);
        x = malloc(ctx->size);
        ctx->init(ctx, x);
        ctx->randtest_not_zero(ctx, x, state);

        mpoly_init(a, n, ctx);
        mpoly_init(b, n, ctx);
        mpoly_randtest(a, state, d, N, ctx);

        mpoly_set(b, a, ctx);
        mpoly_add_coeff(b, m, x, ctx);

        result = (!mpoly_equal(a, b, ctx));
        if (!result)
        {
            printf("FAIL:\n");
            mpoly_print(a, ctx); printf("\n");
            mpoly_print(b, ctx); printf("\n");
            mon_print(m, n); printf("\n");
            ctx->print(ctx, x); printf("\n");
            abort();
        }

        mpoly_clear(a, ctx);
        mpoly_clear(b, ctx);
        mon_clear(m);
        ctx->clear(ctx, x);
        free(x);
    }

    ctx_clear(ctx);

    _randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
