#include "mat.h"

#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

int
main(void)
{
    int i;
    flint_rand_t state;

    printf("init/ clear... ");
    fflush(stdout);

    _randinit(state);

    /* Unmanaged element type (long) */
    for (i = 0; i < 100; i++)
    {
        long m, n;
        ctx_t ctx;
        mat_t A;

        m = n_randint(state, 100) + 1;
        n = n_randint(state, 100) + 1;

        ctx_init_long(ctx);
        mat_init(A, m, n, ctx);

        mat_randtest(A, state, ctx);

        mat_clear(A, ctx);
        ctx_clear(ctx);
    }

    /* Managed element type (mpq_t) */
    for (i = 0; i < 100; i++)
    {
        long m, n;
        ctx_t ctx;
        mat_t A;

        m = n_randint(state, 100) + 1;
        n = n_randint(state, 100) + 1;

        ctx_init_mpq(ctx);
        mat_init(A, m, n, ctx);

        mat_randtest(A, state, ctx);

        mat_clear(A, ctx);
        ctx_clear(ctx);
    }

    _randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
