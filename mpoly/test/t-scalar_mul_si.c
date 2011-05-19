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

    printf("scalar_mul_si... ");
    fflush(stdout);

    flint_randinit(state);

    ctx_init_mpq(ctx);

    /* Check aliasing of a and b */
    for (i = 0; i < 1000; i++)
    {
        mpoly_t a, b;
        long n, d, N;
        long x;

        n = n_randint(state, MON_MAX_VARS) + 1;
        d = n_randint(state, 50) + 1;
        N = n_randint(state, 50) + 1;

        mpoly_init(a, n, ctx);
        mpoly_init(b, n, ctx);
        mpoly_randtest(b, state, d, N, ctx);
        x = z_randtest(state);

        mpoly_scalar_mul_si(a, b, x, ctx);
        mpoly_scalar_mul_si(b, b, x, ctx);

        result = (mpoly_equal(a, b, ctx));
        if (!result)
        {
            printf("FAIL:\n");
            mpoly_print(a, ctx); printf("\n");
            mpoly_print(b, ctx); printf("\n");
            printf("%ld\n", x);
            abort();
        }

        mpoly_clear(a, ctx);
        mpoly_clear(b, ctx);
    }

    /* Check that x (a + b) == x * a + x * b */
    for (i = 0; i < 1000; i++)
    {
        mpoly_t a, b, c1, c2;
        long n, d, N;
        long x;

        n = n_randint(state, MON_MAX_VARS) + 1;
        d = n_randint(state, 50) + 1;
        N = n_randint(state, 50) + 1;

        mpoly_init(a, n, ctx);
        mpoly_init(b, n, ctx);
        mpoly_init(c1, n, ctx);
        mpoly_init(c2, n, ctx);
        mpoly_randtest(a, state, d, N, ctx);
        mpoly_randtest(b, state, d, N, ctx);
        x = z_randtest(state);

        mpoly_scalar_mul_si(c1, a, x, ctx);
        mpoly_scalar_mul_si(c2, b, x, ctx);
        mpoly_add(c2, c1, c2, ctx);

        mpoly_add(c1, a, b, ctx);
        mpoly_scalar_mul_si(c1, c1, x, ctx);

        result = (mpoly_equal(c1, c2, ctx));
        if (!result)
        {
            printf("FAIL:\n");
            mpoly_print(a, ctx); printf("\n");
            mpoly_print(b, ctx); printf("\n");
            mpoly_print(c1, ctx); printf("\n");
            mpoly_print(c2, ctx); printf("\n");
            printf("%ld\n", x);
            abort();
        }

        mpoly_clear(a, ctx);
        mpoly_clear(b, ctx);
        mpoly_clear(c1, ctx);
        mpoly_clear(c2, ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
