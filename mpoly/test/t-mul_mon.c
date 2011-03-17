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
    mat_ctx_t ctx;

    printf("mul_mon... ");
    fflush(stdout);

    flint_randinit(state);

    mat_ctx_init_mpq(ctx);

    /* Check aliasing of a and b */
    for (i = 0; i < 1000; i++)
    {
        mpoly_t a, b;
        mon_t m;
        long n, d, N;

        n = n_randint(state, MON_MAX_VARS) + 1;
        d = n_randint(state, 50) + 1;
        N = n_randint(state, 50) + 1;

        mon_init(m);
        mon_randtest(m, state, n, d);

        mpoly_init(a, n, ctx);
        mpoly_init(b, n, ctx);
        mpoly_randtest(a, state, d, N, ctx);

        mpoly_mul_mon(b, a, m, ctx);
        mpoly_mul_mon(a, a, m, ctx);

        result = (mpoly_equal(a, b, ctx));
        if (!result)
        {
            printf("FAIL:\n");
            mpoly_print(a, ctx); printf("\n");
            mpoly_print(b, ctx); printf("\n");
            mon_print(m, n); printf("\n");
            abort();
        }

        mpoly_clear(a, ctx);
        mpoly_clear(b, ctx);
        mon_clear(m);
    }

    /* Check (b*c) + (b*d) = b*(c+d) */
    for (i = 0; i < 1000; i++)
    {
        mpoly_t  a1, a2, c1, c2;
        mon_t b;
        long n, d, N;

        n = n_randint(state, MON_MAX_VARS) + 1;
        d = n_randint(state, 50) + 1;
        N = n_randint(state, 50) + 1;

        mon_init(b);
        mon_randtest(b, state, n, d);

        mpoly_init(a1, n, ctx);
        mpoly_init(a2, n, ctx);
        mpoly_init(c1, n, ctx);
        mpoly_init(c2, n, ctx);
        mpoly_randtest(c1, state, d, N, ctx);
        mpoly_randtest(c2, state, d, N, ctx);

        mpoly_mul_mon(a1, c1, b, ctx);
        mpoly_mul_mon(a2, c2, b, ctx);
        mpoly_add(a1, a1, a2, ctx);

        mpoly_add(c1, c1, c2, ctx);
        mpoly_mul_mon(a2, c1, b, ctx);

        result = (mpoly_equal(a1, a2, ctx));
        if (!result)
        {
            printf("FAIL:\n");
            mpoly_print(a1, ctx); printf("\n");
            mpoly_print(a2, ctx); printf("\n");
            mon_print(b, n); printf("\n");
            mpoly_print(c1, ctx); printf("\n");
            mpoly_print(c2, ctx); printf("\n");
            abort();
        }

        mpoly_clear(a1, ctx);
        mpoly_clear(a2, ctx);
        mpoly_clear(c1, ctx);
        mpoly_clear(c2, ctx);
        mon_clear(b);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
