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

    printf("sub_coeff... ");
    fflush(stdout);

    flint_randinit(state);

    ctx_init_mpq(ctx);

    for (i = 0; i < 1000; i++)
    {
        mpoly_t a;
        long n, d, N;
        mon_t m;
        char *x, *y, *z;

        n = n_randint(state, MON_MAX_VARS) + 1;
        d = n_randint(state, 50) + 1;
        N = n_randint(state, 50) + 1;

        mon_init(m);
        mon_randtest(m, state, n, d);
        x = malloc(ctx->size);
        y = malloc(ctx->size);
        z = malloc(ctx->size);
        ctx->init(x);
        ctx->init(y);
        ctx->init(z);
        ctx->randtest(x, state);

        mpoly_init(a, n, ctx);
        mpoly_randtest(a, state, d, N, ctx);

        mpoly_get_coeff(y, a, m, ctx);
        ctx->sub(y, y, x);

        mpoly_sub_coeff(a, m, x, ctx);
        mpoly_get_coeff(z, a, m, ctx);

        result = (ctx->equal(y, z));
        if (!result)
        {
            printf("FAIL:\n");
            mpoly_print(a, ctx); printf("\n");
            ctx->print(x); printf("\n");
            ctx->print(y); printf("\n");
            ctx->print(z); printf("\n");
            abort();
        }

        mpoly_clear(a, ctx);
        mon_clear(m);
        ctx->clear(x);
        ctx->clear(y);
        ctx->clear(z);
        free(x);
        free(y);
        free(z);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
