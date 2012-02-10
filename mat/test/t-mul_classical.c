#include "mat.h"

#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("mul_classical... ");
    fflush(stdout);

    _randinit(state);

    /* Check that the identity matrix does what it's supposed to do */

    /* Unmanaged element type (long) */
    for (i = 0; i < 100; i++)
    {
        long m, n;
        ctx_t ctx;
        mat_t A, B, C, D;

        m = n_randint(state, 50) + 1;
        n = m;

        ctx_init_long(ctx);
        mat_init(A, m, n, ctx);
        mat_init(B, m, n, ctx);
        mat_init(C, m, n, ctx);
        mat_init(D, m, n, ctx);

        mat_randtest(A, state, ctx);
        mat_one(B, ctx);

        mat_mul_classical(C, A, B, ctx);
        mat_mul_classical(D, B, A, ctx);

        result = mat_equal(A, C, ctx) && mat_equal(C, D, ctx);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("Matrix A\n");
            mat_print(A, ctx);
            printf("\n");
            printf("Matrix B\n");
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
        n = m;

        ctx_init_mpq(ctx);
        mat_init(A, m, n, ctx);
        mat_init(B, m, n, ctx);
        mat_init(C, m, n, ctx);
        mat_init(D, m, n, ctx);

        mat_randtest(A, state, ctx);
        mat_one(B, ctx);

        mat_mul_classical(C, A, B, ctx);
        mat_mul_classical(D, B, A, ctx);

        result = mat_equal(A, C, ctx) && mat_equal(C, D, ctx);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("Matrix A\n");
            mat_print(A, ctx);
            printf("\n");
            printf("Matrix B\n");
            mat_print(B, ctx);
            printf("\n");
        }

        mat_clear(A, ctx);
        mat_clear(B, ctx);
        mat_clear(C, ctx);
        mat_clear(D, ctx);
        ctx_clear(ctx);
    }

    /* Check that A * (B + C) == A * B + A * C */

    /* Unmanaged element type (long) */
    for (i = 0; i < 100; i++)
    {
        long ell, m, n;
        ctx_t ctx;
        mat_t A, B, C, D, E, T1, T2;

        ell = n_randint(state, 50) + 1;
        m   = n_randint(state, 50) + 1;
        n   = n_randint(state, 50) + 1;

        ctx_init_long(ctx);
        mat_init(A, ell, m, ctx);
        mat_init(B, m, n, ctx);
        mat_init(C, m, n, ctx);
        mat_init(D, ell, n, ctx);
        mat_init(E, ell, n, ctx);
        mat_init(T1, m, n, ctx);
        mat_init(T2, ell, n, ctx);

        mat_randtest(A, state, ctx);
        mat_randtest(B, state, ctx);
        mat_randtest(C, state, ctx);

        mat_add(T1, B, C, ctx);
        mat_mul_classical(D, A, T1, ctx);

        mat_mul_classical(E, A, B, ctx);
        mat_mul_classical(T2, A, C, ctx);
        mat_add(E, E, T2, ctx);

        result = mat_equal(D, E, ctx);
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
        mat_clear(D, ctx);
        mat_clear(E, ctx);
        mat_clear(T1, ctx);
        mat_clear(T2, ctx);
        ctx_clear(ctx);
    }

    /* Managed element type (mpq) */
    for (i = 0; i < 100; i++)
    {
        long ell, m, n;
        ctx_t ctx;
        mat_t A, B, C, D, E, T1, T2;

        ell = n_randint(state, 50) + 1;
        m   = n_randint(state, 50) + 1;
        n   = n_randint(state, 50) + 1;

        ctx_init_mpq(ctx);
        mat_init(A, ell, m, ctx);
        mat_init(B, m, n, ctx);
        mat_init(C, m, n, ctx);
        mat_init(D, ell, n, ctx);
        mat_init(E, ell, n, ctx);
        mat_init(T1, m, n, ctx);
        mat_init(T2, ell, n, ctx);

        mat_randtest(A, state, ctx);
        mat_randtest(B, state, ctx);
        mat_randtest(C, state, ctx);

        mat_add(T1, B, C, ctx);
        mat_mul_classical(D, A, T1, ctx);

        mat_mul_classical(E, A, B, ctx);
        mat_mul_classical(T2, A, C, ctx);
        mat_add(E, E, T2, ctx);

        result = mat_equal(D, E, ctx);
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
        mat_clear(D, ctx);
        mat_clear(E, ctx);
        mat_clear(T1, ctx);
        mat_clear(T2, ctx);
        ctx_clear(ctx);
    }

    _randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
