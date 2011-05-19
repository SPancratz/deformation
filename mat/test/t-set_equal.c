#include "mat.h"

#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("set/ equal... ");
    fflush(stdout);

    flint_randinit(state);

    /* Unmanaged element type (long), equal */
    for (i = 0; i < 100; i++)
    {
        long m, n;
        ctx_t ctx;
        mat_t A, B;

        m = n_randint(state, 100) + 1;
        n = n_randint(state, 100) + 1;

        ctx_init_long(ctx);
        mat_init(A, m, n, ctx);
        mat_init(B, m, n, ctx);

        mat_randtest(A, state, ctx);
        mat_set(B, A, ctx);

        result = mat_equal(A, B, ctx);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("Matrix A\n");
            mat_print(A, ctx);
            printf("\n");
            printf("Matrix B\n");
            mat_print(B, ctx);
            printf("\n");
            abort();
        }

        mat_clear(A, ctx);
        mat_clear(B, ctx);
        ctx_clear(ctx);
    }

    /* Unmanaged element type (long), unequal */
    for (i = 0; i < 100; i++)
    {
        long m, n;
        long r, c;
        ctx_t ctx;
        mat_t A, B;

        m = n_randint(state, 100) + 1;
        n = n_randint(state, 100) + 1;

        ctx_init_long(ctx);
        mat_init(A, m, n, ctx);
        mat_init(B, m, n, ctx);

        mat_randtest(A, state, ctx);
        mat_set(B, A, ctx);

        r = n_randint(state, m);
        c = n_randint(state, n);
        if (ctx->is_zero(ctx, mat_entry(B, r, c, ctx)))
            ctx->randtest_not_zero(ctx, mat_entry(B, r, c, ctx), state);
        else
            ctx->zero(ctx, mat_entry(B, r, c, ctx));

        result = !mat_equal(A, B, ctx);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("Matrix A\n");
            mat_print(A, ctx);
            printf("\n");
            printf("Matrix B\n");
            mat_print(B, ctx);
            printf("\n");
            abort();
        }

        mat_clear(A, ctx);
        mat_clear(B, ctx);
        ctx_clear(ctx);
    }

    /* Managed element type (mpq_t), equal */
    for (i = 0; i < 100; i++)
    {
        long m, n;
        ctx_t ctx;
        mat_t A, B;

        m = n_randint(state, 100) + 1;
        n = n_randint(state, 100) + 1;

        ctx_init_mpq(ctx);
        mat_init(A, m, n, ctx);
        mat_init(B, m, n, ctx);

        mat_randtest(A, state, ctx);
        mat_set(B, A, ctx);

        result = mat_equal(A, B, ctx);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("Matrix A\n");
            mat_print(A, ctx);
            printf("\n");
            printf("Matrix B\n");
            mat_print(B, ctx);
            printf("\n");
            abort();
        }

        mat_clear(A, ctx);
        mat_clear(B, ctx);
        ctx_clear(ctx);
    }

    /* Managed element type (mpq_t), unequal */
    for (i = 0; i < 100; i++)
    {
        long m, n;
        long r, c;
        ctx_t ctx;
        mat_t A, B;

        m = n_randint(state, 100) + 1;
        n = n_randint(state, 100) + 1;

        ctx_init_mpq(ctx);
        mat_init(A, m, n, ctx);
        mat_init(B, m, n, ctx);

        mat_randtest(A, state, ctx);
        mat_set(B, A, ctx);

        r = n_randint(state, m);
        c = n_randint(state, n);
        if (ctx->is_zero(ctx, mat_entry(B, r, c, ctx)))
            ctx->randtest_not_zero(ctx, mat_entry(B, r, c, ctx), state);
        else
            ctx->zero(ctx, mat_entry(B, r, c, ctx));

        result = !mat_equal(A, B, ctx);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("Matrix A\n");
            mat_print(A, ctx);
            printf("\n");
            printf("Matrix B\n");
            mat_print(B, ctx);
            printf("\n");
            abort();
        }

        mat_clear(A, ctx);
        mat_clear(B, ctx);
        ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
