#include "mat_dense.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("zero\n");
    printf("----\n");
    fflush(stdout);

    flint_randinit(state);

    /* Unmanaged element type (long) */
    for (i = 0; i < 100; i++)
    {
        long m, n;
        mat_ctx_t ctx;
        mat_dense_t A;

        m = n_randint(state, 100) + 1;
        n = n_randint(state, 100) + 1;

        mat_ctx_init_long(ctx);
        mat_dense_init(A, m, n, ctx);

        mat_dense_randtest(A, state, ctx);
        mat_dense_zero(A, ctx);

        result = mat_dense_is_zero(A, ctx);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("Matrix A\n");
            mat_dense_print(A, ctx);
            printf("\n");
            abort();
        }

        mat_dense_clear(A, ctx);
        mat_ctx_clear(ctx);
    }

    /* Managed element type (mpq_t) */
    for (i = 0; i < 100; i++)
    {
        long m, n;
        mat_ctx_t ctx;
        mat_dense_t A;

        m = n_randint(state, 100) + 1;
        n = n_randint(state, 100) + 1;

        mat_ctx_init_mpq(ctx);
        mat_dense_init(A, m, n, ctx);

        mat_dense_randtest(A, state, ctx);
        mat_dense_zero(A, ctx);

        result = mat_dense_is_zero(A, ctx);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("Matrix A\n");
            mat_dense_print(A, ctx);
            printf("\n");
            abort();
        }

        mat_dense_clear(A, ctx);
        mat_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
