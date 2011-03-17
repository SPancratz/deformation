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

    printf("degree... ");
    fflush(stdout);

    flint_randinit(state);

    mat_ctx_init_mpq(ctx);

    /* Check deg(a) + deg(b) == deg(ab) for a, b != 0 */
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
        mpoly_randtest_not_zero(a, state, d, N, ctx);
        mpoly_randtest_not_zero(b, state, d, N, ctx);

        mpoly_mul(c, a, b, ctx);

        result = (mpoly_degree(c, -1, ctx) == mpoly_degree(a, -1, ctx) + mpoly_degree(b, -1, ctx));
        if (!result)
        {
            printf("FAIL:\n");
            printf("n d N = %ld %ld %ld\n", n, d, N);
            printf("a = "), mpoly_print(a, ctx), printf("\n");
            printf("b = "), mpoly_print(b, ctx), printf("\n");
            printf("c = "), mpoly_print(c, ctx), printf("\n");
            printf("deg(a) = %ld\n", mpoly_degree(a, -1, ctx));
            printf("deg(b) = %ld\n", mpoly_degree(b, -1, ctx));
            printf("deg(c) = %ld\n", mpoly_degree(c, -1, ctx));
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
