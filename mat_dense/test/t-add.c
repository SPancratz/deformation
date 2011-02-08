#include "mat_dense.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("add\n");
    printf("---\n");
    fflush(stdout);

    flint_randinit(state);

    /* Check that addition is abelian */

    /* Unmanaged element type (long) */
    for (i = 0; i < 100; i++)
    {
        long m, n;
        mat_ctx_t ctx;
        mat_dense_t A, B, C, D;

        m = n_randint(state, 50) + 1;
        n = n_randint(state, 50) + 1;

        mat_ctx_init_long(ctx);
        mat_dense_init(A, m, n, ctx);
        mat_dense_init(B, m, n, ctx);
        mat_dense_init(C, m, n, ctx);
        mat_dense_init(D, m, n, ctx);

        mat_dense_randtest(A, state, ctx);
        mat_dense_randtest(B, state, ctx);

        mat_dense_add(C, A, B, ctx);
        mat_dense_add(D, B, A, ctx);

        result = mat_dense_equal(C, D, ctx);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("Matrix A\n");
            mat_dense_print(A, ctx);
            printf("\n");
            printf("Matrix A\n");
            mat_dense_print(B, ctx);
            printf("\n");
        }

        mat_dense_clear(A, ctx);
        mat_dense_clear(B, ctx);
        mat_dense_clear(C, ctx);
        mat_dense_clear(D, ctx);
        mat_ctx_clear(ctx);
    }

    /* Managed element type (mpq_t) */
    for (i = 0; i < 100; i++)
    {
        long m, n;
        mat_ctx_t ctx;
        mat_dense_t A, B, C, D;

        m = n_randint(state, 50) + 1;
        n = n_randint(state, 50) + 1;

        mat_ctx_init_mpq(ctx);
        mat_dense_init(A, m, n, ctx);
        mat_dense_init(B, m, n, ctx);
        mat_dense_init(C, m, n, ctx);
        mat_dense_init(D, m, n, ctx);

        mat_dense_randtest(A, state, ctx);
        mat_dense_randtest(B, state, ctx);

        mat_dense_add(C, A, B, ctx);
        mat_dense_add(D, B, A, ctx);

        result = mat_dense_equal(C, D, ctx);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("Matrix A\n");
            mat_dense_print(A, ctx);
            printf("\n");
            printf("Matrix A\n");
            mat_dense_print(B, ctx);
            printf("\n");
        }

        mat_dense_clear(A, ctx);
        mat_dense_clear(B, ctx);
        mat_dense_clear(C, ctx);
        mat_dense_clear(D, ctx);
        mat_ctx_clear(ctx);
    }

    /* Check that the zero matrix does what it's supposed to do */

    /* Unmanaged element type (long) */
    for (i = 0; i < 100; i++)
    {
        long m, n;
        mat_ctx_t ctx;
        mat_dense_t A, B, C, D;

        m = n_randint(state, 50) + 1;
        n = n_randint(state, 50) + 1;

        mat_ctx_init_long(ctx);
        mat_dense_init(A, m, n, ctx);
        mat_dense_init(B, m, n, ctx);
        mat_dense_init(C, m, n, ctx);
        mat_dense_init(D, m, n, ctx);

        mat_dense_randtest(A, state, ctx);
        mat_dense_zero(B, ctx);

        mat_dense_add(C, A, B, ctx);
        mat_dense_add(D, B, A, ctx);

        result = mat_dense_equal(A, C, ctx) && mat_dense_equal(C, D, ctx);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("Matrix A\n");
            mat_dense_print(A, ctx);
            printf("\n");
            printf("Matrix A\n");
            mat_dense_print(B, ctx);
            printf("\n");
        }

        mat_dense_clear(A, ctx);
        mat_dense_clear(B, ctx);
        mat_dense_clear(C, ctx);
        mat_dense_clear(D, ctx);
        mat_ctx_clear(ctx);
    }

    /* Managed element type (mpq_t) */
    for (i = 0; i < 100; i++)
    {
        long m, n;
        mat_ctx_t ctx;
        mat_dense_t A, B, C, D;

        m = n_randint(state, 50) + 1;
        n = n_randint(state, 50) + 1;

        mat_ctx_init_mpq(ctx);
        mat_dense_init(A, m, n, ctx);
        mat_dense_init(B, m, n, ctx);
        mat_dense_init(C, m, n, ctx);
        mat_dense_init(D, m, n, ctx);

        mat_dense_randtest(A, state, ctx);
        mat_dense_zero(B, ctx);

        mat_dense_add(C, A, B, ctx);
        mat_dense_add(D, B, A, ctx);

        result = mat_dense_equal(A, C, ctx) && mat_dense_equal(C, D, ctx);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("Matrix A\n");
            mat_dense_print(A, ctx);
            printf("\n");
            printf("Matrix A\n");
            mat_dense_print(B, ctx);
            printf("\n");
        }

        mat_dense_clear(A, ctx);
        mat_dense_clear(B, ctx);
        mat_dense_clear(C, ctx);
        mat_dense_clear(D, ctx);
        mat_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
