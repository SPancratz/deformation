#include "mat_csr.h"

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
        double d;
        ctx_t ctx;
        mat_csr_t A;

        m = n_randint(state, 100) + 1;
        n = n_randint(state, 100) + 1;
        d = (double) n_randint(state, 101) / (double) 100;

        ctx_init_long(ctx);
        mat_csr_init(A, m, n, ctx);

        mat_csr_randtest(A, state, d, ctx);

        mat_csr_clear(A, ctx);
        ctx_clear(ctx);
    }

    /* Managed element type (mpq_t) */
    for (i = 0; i < 1; i++)
    {
        long m, n;
        double d;
        ctx_t ctx;
        mat_csr_t A;

        m = n_randint(state, 100) + 1;
        n = n_randint(state, 100) + 1;
        d = (double) n_randint(state, 101) / (double) 100;

        ctx_init_mpq(ctx);
        mat_csr_init(A, m, n, ctx);

        mat_csr_randtest(A, state, d, ctx);

        mat_csr_clear(A, ctx);
        ctx_clear(ctx);
    }

    _randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
