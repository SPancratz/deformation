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

    printf("derivative... ");
    fflush(stdout);

    _randinit(state);

    ctx_init_mpq(ctx);

    /* Check aliasing of a and b */
    for (i = 0; i < 1000; i++)
    {
        mpoly_t a, b;
        long n, d, N;
        int var;

        n = n_randint(state, MON_MAX_VARS) + 1;
        d = n_randint(state, 50) + 1;
        N = n_randint(state, 50) + 1;

        var = n_randint(state, n);

        mpoly_init(a, n, ctx);
        mpoly_init(b, n, ctx);
        mpoly_randtest(a, state, d, N, ctx);

        mpoly_derivative(b, a, var, ctx);
        mpoly_derivative(a, a, var, ctx);

        result = (mpoly_equal(a, b, ctx));
        if (!result)
        {
            printf("FAIL:\n");
            printf("a = "), mpoly_print(a, ctx), printf("\n");
            printf("b = "), mpoly_print(b, ctx), printf("\n");
            printf("var = %d\n", var);
            abort();
        }

        mpoly_clear(a, ctx);
        mpoly_clear(b, ctx);
    }

    /* Check d(a + b) == da + db */
    for (i = 0; i < 1000; i++)
    {
        mpoly_t a, b, c1, c2;
        long n, d, N;
        int var;

        n = n_randint(state, MON_MAX_VARS) + 1;
        d = n_randint(state, 50) + 1;
        N = n_randint(state, 50) + 1;

        var = n_randint(state, n);

        mpoly_init(a, n, ctx);
        mpoly_init(b, n, ctx);
        mpoly_init(c1, n, ctx);
        mpoly_init(c2, n, ctx);
        mpoly_randtest(a, state, d, N, ctx);
        mpoly_randtest(b, state, d, N, ctx);

        mpoly_derivative(c1, a, var, ctx);
        mpoly_derivative(c2, b, var, ctx);
        mpoly_add(c2, c1, c2, ctx);

        mpoly_add(c1, a, b, ctx);
        mpoly_derivative(c1, c1, var, ctx);

        result = (mpoly_equal(c1, c2, ctx));
        if (!result)
        {
            printf("FAIL:\n");
            printf("a = "), mpoly_print(a, ctx), printf("\n");
            printf("b = "), mpoly_print(b, ctx), printf("\n");
            printf("c = "), mpoly_print(c1, ctx), printf("\n");
            printf("d = "), mpoly_print(c2, ctx), printf("\n");
            printf("var = %d\n", var);
            abort();
        }

        mpoly_clear(a, ctx);
        mpoly_clear(b, ctx);
        mpoly_clear(c1, ctx);
        mpoly_clear(c2, ctx);
    }

    _randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
