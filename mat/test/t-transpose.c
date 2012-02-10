#include "mat.h"

#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("transpose... ");
    fflush(stdout);

    _randinit(state);

    /* Check aliasing */

    /* Unmanaged element type (long) */
    for (i = 0; i < 100; i++)
    {
        long m, n;
        ctx_t ctx;
        mat_t A, B, C;

        m = n_randint(state, 50) + 1;
        n = m;

        ctx_init_long(ctx);
        mat_init(A, m, n, ctx);
        mat_init(B, m, n, ctx);
        mat_init(C, m, n, ctx);

        mat_randtest(A, state, ctx);

        mat_set(B, A, ctx);
        mat_transpose(C, B, ctx);
        mat_transpose(B, B, ctx);

        result = mat_equal(B, C, ctx);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("Matrix A\n");
            mat_print(A, ctx);
            printf("\n");
            printf("Matrix B\n");
            mat_print(B, ctx);
            printf("\n");
            printf("Matrix C\n");
            mat_print(C, ctx);
            printf("\n");
        }

        mat_clear(A, ctx);
        mat_clear(B, ctx);
        mat_clear(C, ctx);
        ctx_clear(ctx);
    }

    /* Check that (A^t)^t == A */

    /* Unmanaged element type (long) */
    for (i = 0; i < 100; i++)
    {
        long m, n;
        ctx_t ctx;
        mat_t A, B, C;

        m = n_randint(state, 50) + 1;
        n = n_randint(state, 50) + 1;

        ctx_init_long(ctx);
        mat_init(A, m, n, ctx);
        mat_init(B, n, m, ctx);
        mat_init(C, m, n, ctx);

        mat_randtest(A, state, ctx);

        mat_transpose(B, A, ctx);
        mat_transpose(C, B, ctx);

        result = mat_equal(A, C, ctx);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("Matrix A\n");
            mat_print(A, ctx);
            printf("\n");
            printf("Matrix B\n");
            mat_print(B, ctx);
            printf("\n");
            printf("Matrix C\n");
            mat_print(C, ctx);
            printf("\n");
        }

        mat_clear(A, ctx);
        mat_clear(B, ctx);
        mat_clear(C, ctx);
        ctx_clear(ctx);
    }

    _randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

