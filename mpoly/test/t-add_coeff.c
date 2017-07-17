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

    printf("add_coeff... ");
    fflush(stdout);

    _randinit(state);

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
        ctx->init(ctx, x);
        ctx->init(ctx, y);
        ctx->init(ctx, z);
        ctx->randtest(ctx, x, state);

        mpoly_init(a, n, ctx);
        mpoly_randtest(a, state, d, N, ctx);

        mpoly_get_coeff(y, a, m, ctx);
        ctx->add(ctx, y, y, x);

        mpoly_add_coeff(a, m, x, ctx);
        mpoly_get_coeff(z, a, m, ctx);

        result = (ctx->equal(ctx, y, z));
        if (!result)
        {
            printf("FAIL:\n");
            mpoly_print(a, ctx); printf("\n");
            ctx->print(ctx, x); printf("\n");
            ctx->print(ctx, y); printf("\n");
            ctx->print(ctx, z); printf("\n");
            abort();
        }

        mpoly_clear(a, ctx);
        mon_clear(m);
        ctx->clear(ctx, x);
        ctx->clear(ctx, y);
        ctx->clear(ctx, z);
        free(x);
        free(y);
        free(z);
    }

    _randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
