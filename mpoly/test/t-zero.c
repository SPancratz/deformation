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

    printf("is_zero/ zero... ");
    fflush(stdout);

    flint_randinit(state);

    ctx_init_mpq(ctx);

    for (i = 0; i < 1000; i++)
    {
        mpoly_t a;
        long n, d, N;

        n = n_randint(state, MON_MAX_VARS) + 1;
        d = n_randint(state, 50) + 1;
        N = n_randint(state, 50) + 1;

        mpoly_init(a, n, ctx);
        mpoly_randtest(a, state, d, N, ctx);

        mpoly_zero(a, ctx);

        result = (mpoly_is_zero(a, ctx));
        if (!result)
        {
            printf("FAIL:\n");
            mpoly_print(a, ctx); printf("\n");
            abort();
        }

        mpoly_clear(a, ctx);
    }

    for (i = 0; i < 1000; i++)
    {
        mpoly_t a;
        long n, d, N;

        n = n_randint(state, MON_MAX_VARS) + 1;
        d = n_randint(state, 50) + 1;
        N = n_randint(state, 50) + 1;

        mpoly_init(a, n, ctx);
        mpoly_randtest_not_zero(a, state, d, N, ctx);

        result = (!mpoly_is_zero(a, ctx));
        if (!result)
        {
            printf("FAIL:\n");
            mpoly_print(a, ctx); printf("\n");
            abort();
        }

        mpoly_clear(a, ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
