/******************************************************************************

    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <mpir.h>

#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

#include "mpoly.h"

int
main(void)
{
    int i;
    flint_rand_t state;
    mat_ctx_t ctx;

    printf("init/clear... ");
    fflush(stdout);

    flint_randinit(state);

    mat_ctx_init_mpq(ctx);
    for (i = 0; i < 1000; i++)
    {
        mpoly_t a;
        long n, d, N;

        n = n_randint(state, MON_MAX_VARS) + 1;
        d = n_randint(state, 50) + 1;
        N = n_randint(state, 50) + 1;

        mpoly_init(a, n, ctx);
        mpoly_randtest(a, state, d, N, ctx);
        mpoly_clear(a, ctx);
    }
    mat_ctx_clear(ctx);

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
