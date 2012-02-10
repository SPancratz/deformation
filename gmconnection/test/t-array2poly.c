#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>

#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

#include "mpoly.h"
#include "mat_csr.h"
#include "gmconnection.h"

#define RUNS 100

int
main(void)
{
    int i, result;
    flint_rand_t state;
    ctx_t ctx;

    printf("array2poly... ");
    fflush(stdout);

    _randinit(state);

    ctx_init_fmpz_poly_q(ctx);

    printf("[");
    for (i = 0; i < RUNS; i++)
    {
        mpoly_t a, b, c;
        long n, d, N;
        mon_t *m;
        long len;
        char *v;

        n = n_randint(state, MON_MAX_VARS) + 1;
        d = n_randint(state, 10) + 1;
        N = n_randint(state, 20) + 1;

        mpoly_init(a, n, ctx);
        mpoly_init(b, n, ctx);
        mpoly_init(c, n, ctx);
        mpoly_randhom(a, state, d, N, ctx);

        m = mon_generate_by_degree(&len, n, d);

        v = _vec_init(len, ctx);

        gmc_poly2array(v, a, m, len, ctx);
        gmc_array2poly(b, v, m, len, ctx);

        gmc_poly2array(v, b, m, len, ctx);
        gmc_array2poly(c, v, m, len, ctx);

        if (!mpoly_is_zero(b, ctx))
            printf("."), fflush(stdout);

        result = (mpoly_equal(b, c, ctx));
        if (!result)
        {
            long k;

            printf("FAIL:\n\n");
            printf("n d N = %ld %ld %ld\n", n, d, N);
            printf("m = { ");
            for (k = 0; k < len; k++)
                mon_print(m[k], n), printf(" ");
            printf("}, len = %ld\n", len);
            printf("a = "), mpoly_print(a, ctx), printf("\n");
            printf("b = "), mpoly_print(b, ctx), printf("\n");
            printf("c = "), mpoly_print(c, ctx), printf("\n");
            abort();
        }

        mpoly_clear(a, ctx);
        mpoly_clear(b, ctx);
        mpoly_clear(c, ctx);
        _vec_clear(v, len, ctx);
        free(m);
    }
    printf("] ");

    ctx_clear(ctx);

    _randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

