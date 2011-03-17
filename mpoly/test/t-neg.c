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

    printf("neg... ");
    fflush(stdout);

    flint_randinit(state);

    mat_ctx_init_mpq(ctx);

    /* Check aliasing of a and b */
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

        mpoly_neg(b, a, ctx);
        mpoly_neg(b, b, ctx);

        result = (mpoly_equal(a, b, ctx));
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

    /* Check a + (-a) == 0 */
    for (i = 0; i < 1000; i++)
    {
        mpoly_t a, b, c;
        long n, d, N;

        n = n_randint(state, MON_MAX_VARS) + 1;
        d = n_randint(state, 50) + 1;
        N = n_randint(state, 50) + 1;

        mpoly_init(a, n, ctx);
        mpoly_init(b, n, ctx);
        mpoly_init(c, n, ctx);
        mpoly_randtest(a, state, d, N, ctx);

        mpoly_neg(b, a, ctx);
        mpoly_add(c, a, b, ctx);

        result = (mpoly_is_zero(c, ctx));
        if (!result)
        {
            printf("FAIL:\n");
            mpoly_print(a, ctx); printf("\n");
            mpoly_print(b, ctx); printf("\n");
            mpoly_print(c, ctx); printf("\n");
            abort();
        }

        mpoly_clear(a, ctx);
        mpoly_clear(b, ctx);
        mpoly_clear(c, ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
