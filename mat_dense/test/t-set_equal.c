#include "mat_dense.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("set/ equal\n");
    printf("----------\n");
    fflush(stdout);

    flint_randinit(state);

    /* Unmanaged element type (long), equal */
    for (i = 0; i < 100; i++)
    {
        long m, n;
        mat_ctx_t ctx;
        mat_dense_t A, B;

        m = n_randint(state, 100) + 1;
        n = n_randint(state, 100) + 1;

        mat_ctx_init_long(ctx);
        mat_dense_init(A, m, n, ctx);
        mat_dense_init(B, m, n, ctx);

        mat_dense_randtest(A, state, ctx);
        mat_dense_set(B, A, ctx);

        result = mat_dense_equal(A, B, ctx);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("Matrix A\n");
            mat_dense_print(A, ctx);
            printf("\n");
            printf("Matrix B\n");
            mat_dense_print(B, ctx);
            printf("\n");
            abort();
        }

        mat_dense_clear(A, ctx);
        mat_dense_clear(B, ctx);
        mat_ctx_clear(ctx);
    }

    /* Unmanaged element type (long), unequal */
    for (i = 0; i < 100; i++)
    {
        long m, n;
        long r, c;
        mat_ctx_t ctx;
        mat_dense_t A, B;

        m = n_randint(state, 100) + 1;
        n = n_randint(state, 100) + 1;

        mat_ctx_init_long(ctx);
        mat_dense_init(A, m, n, ctx);
        mat_dense_init(B, m, n, ctx);

        mat_dense_randtest(A, state, ctx);
        mat_dense_set(B, A, ctx);

        r = n_randint(state, m);
        c = n_randint(state, n);
        if (ctx->is_zero(mat_dense_entry(B, r, c, ctx)))
            ctx->randtest_not_zero(mat_dense_entry(B, r, c, ctx), state);
        else
            ctx->zero(mat_dense_entry(B, r, c, ctx));

        result = !mat_dense_equal(A, B, ctx);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("Matrix A\n");
            mat_dense_print(A, ctx);
            printf("\n");
            printf("Matrix B\n");
            mat_dense_print(B, ctx);
            printf("\n");
            abort();
        }

        mat_dense_clear(A, ctx);
        mat_dense_clear(B, ctx);
        mat_ctx_clear(ctx);
    }

    /* Managed element type (mpq_t), equal */
    for (i = 0; i < 100; i++)
    {
        long m, n;
        mat_ctx_t ctx;
        mat_dense_t A, B;

        m = n_randint(state, 100) + 1;
        n = n_randint(state, 100) + 1;

        mat_ctx_init_mpq(ctx);
        mat_dense_init(A, m, n, ctx);
        mat_dense_init(B, m, n, ctx);

        mat_dense_randtest(A, state, ctx);
        mat_dense_set(B, A, ctx);

        result = mat_dense_equal(A, B, ctx);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("Matrix A\n");
            mat_dense_print(A, ctx);
            printf("\n");
            printf("Matrix B\n");
            mat_dense_print(B, ctx);
            printf("\n");
            abort();
        }

        mat_dense_clear(A, ctx);
        mat_dense_clear(B, ctx);
        mat_ctx_clear(ctx);
    }

    /* Managed element type (mpq_t), unequal */
    for (i = 0; i < 100; i++)
    {
        long m, n;
        long r, c;
        mat_ctx_t ctx;
        mat_dense_t A, B;

        m = n_randint(state, 100) + 1;
        n = n_randint(state, 100) + 1;

        mat_ctx_init_mpq(ctx);
        mat_dense_init(A, m, n, ctx);
        mat_dense_init(B, m, n, ctx);

        mat_dense_randtest(A, state, ctx);
        mat_dense_set(B, A, ctx);

        r = n_randint(state, m);
        c = n_randint(state, n);
        if (ctx->is_zero(mat_dense_entry(B, r, c, ctx)))
            ctx->randtest_not_zero(mat_dense_entry(B, r, c, ctx), state);
        else
            ctx->zero(mat_dense_entry(B, r, c, ctx));

        result = !mat_dense_equal(A, B, ctx);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("Matrix A\n");
            mat_dense_print(A, ctx);
            printf("\n");
            printf("Matrix B\n");
            mat_dense_print(B, ctx);
            printf("\n");
            abort();
        }

        mat_dense_clear(A, ctx);
        mat_dense_clear(B, ctx);
        mat_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
