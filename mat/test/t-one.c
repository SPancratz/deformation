#include "mat.h"

#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("one... ");
    fflush(stdout);

    flint_randinit(state);

    /* Unmanaged element type (long) */
    for (i = 0; i < 100; i++)
    {
        long m, n;
        ctx_t ctx;
        mat_t A;

        m = n_randint(state, 100) + 1;
        n = m;

        ctx_init_long(ctx);
        mat_init(A, m, n, ctx);

        mat_randtest(A, state, ctx);
        mat_one(A, ctx);

        result = mat_is_one(A, ctx);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("Matrix A\n");
            mat_print(A, ctx);
            printf("\n");
            abort();
        }

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
        n = m;

        ctx_init_mpq(ctx);
        mat_init(A, m, n, ctx);

        mat_randtest(A, state, ctx);
        mat_one(A, ctx);

        result = mat_is_one(A, ctx);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("Matrix A\n");
            mat_print(A, ctx);
            printf("\n");
            abort();
        }

        mat_clear(A, ctx);
        ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
