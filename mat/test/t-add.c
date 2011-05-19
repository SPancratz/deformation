#include "mat.h"

#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("add... ");
    fflush(stdout);

    flint_randinit(state);

    /* Check that addition is abelian */

    /* Unmanaged element type (long) */
    for (i = 0; i < 100; i++)
    {
        long m, n;
        ctx_t ctx;
        mat_t A, B, C, D;

        m = n_randint(state, 50) + 1;
        n = n_randint(state, 50) + 1;

        ctx_init_long(ctx);
        mat_init(A, m, n, ctx);
        mat_init(B, m, n, ctx);
        mat_init(C, m, n, ctx);
        mat_init(D, m, n, ctx);

        mat_randtest(A, state, ctx);
        mat_randtest(B, state, ctx);

        mat_add(C, A, B, ctx);
        mat_add(D, B, A, ctx);

        result = mat_equal(C, D, ctx);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("Matrix A\n");
            mat_print(A, ctx);
            printf("\n");
            printf("Matrix A\n");
            mat_print(B, ctx);
            printf("\n");
        }

        mat_clear(A, ctx);
        mat_clear(B, ctx);
        mat_clear(C, ctx);
        mat_clear(D, ctx);
        ctx_clear(ctx);
    }

    /* Managed element type (mpq_t) */
    for (i = 0; i < 100; i++)
    {
        long m, n;
        ctx_t ctx;
        mat_t A, B, C, D;

        m = n_randint(state, 50) + 1;
        n = n_randint(state, 50) + 1;

        ctx_init_mpq(ctx);
        mat_init(A, m, n, ctx);
        mat_init(B, m, n, ctx);
        mat_init(C, m, n, ctx);
        mat_init(D, m, n, ctx);

        mat_randtest(A, state, ctx);
        mat_randtest(B, state, ctx);

        mat_add(C, A, B, ctx);
        mat_add(D, B, A, ctx);

        result = mat_equal(C, D, ctx);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("Matrix A\n");
            mat_print(A, ctx);
            printf("\n");
            printf("Matrix A\n");
            mat_print(B, ctx);
            printf("\n");
        }

        mat_clear(A, ctx);
        mat_clear(B, ctx);
        mat_clear(C, ctx);
        mat_clear(D, ctx);
        ctx_clear(ctx);
    }

    /* Check that the zero matrix does what it's supposed to do */

    /* Unmanaged element type (long) */
    for (i = 0; i < 100; i++)
    {
        long m, n;
        ctx_t ctx;
        mat_t A, B, C, D;

        m = n_randint(state, 50) + 1;
        n = n_randint(state, 50) + 1;

        ctx_init_long(ctx);
        mat_init(A, m, n, ctx);
        mat_init(B, m, n, ctx);
        mat_init(C, m, n, ctx);
        mat_init(D, m, n, ctx);

        mat_randtest(A, state, ctx);
        mat_zero(B, ctx);

        mat_add(C, A, B, ctx);
        mat_add(D, B, A, ctx);

        result = mat_equal(A, C, ctx) && mat_equal(C, D, ctx);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("Matrix A\n");
            mat_print(A, ctx);
            printf("\n");
            printf("Matrix A\n");
            mat_print(B, ctx);
            printf("\n");
        }

        mat_clear(A, ctx);
        mat_clear(B, ctx);
        mat_clear(C, ctx);
        mat_clear(D, ctx);
        ctx_clear(ctx);
    }

    /* Managed element type (mpq_t) */
    for (i = 0; i < 100; i++)
    {
        long m, n;
        ctx_t ctx;
        mat_t A, B, C, D;

        m = n_randint(state, 50) + 1;
        n = n_randint(state, 50) + 1;

        ctx_init_mpq(ctx);
        mat_init(A, m, n, ctx);
        mat_init(B, m, n, ctx);
        mat_init(C, m, n, ctx);
        mat_init(D, m, n, ctx);

        mat_randtest(A, state, ctx);
        mat_zero(B, ctx);

        mat_add(C, A, B, ctx);
        mat_add(D, B, A, ctx);

        result = mat_equal(A, C, ctx) && mat_equal(C, D, ctx);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("Matrix A\n");
            mat_print(A, ctx);
            printf("\n");
            printf("Matrix A\n");
            mat_print(B, ctx);
            printf("\n");
        }

        mat_clear(A, ctx);
        mat_clear(B, ctx);
        mat_clear(C, ctx);
        mat_clear(D, ctx);
        ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
