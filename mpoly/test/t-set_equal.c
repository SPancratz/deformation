#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>

#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

#include "mpoly.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;
    ctx_t ctx;

    printf("set/equal... ");
    fflush(stdout);

    flint_randinit(state);

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
        ctx->init(x);
        ctx->randtest_not_zero(x, state);

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
            ctx->print(x); printf("\n");
            abort();
        }

        mpoly_clear(a, ctx);
        mpoly_clear(b, ctx);
        mon_clear(m);
        ctx->clear(x);
        free(x);
    }

    ctx_clear(ctx);

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
