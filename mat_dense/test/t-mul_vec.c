#include "mat_dense.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("mul_vec... ");
    fflush(stdout);

    flint_randinit(state);

    /* Check that the identity matrix does what it's supposed to do */

    /* Unmanaged element type (long) */
    for (i = 0; i < 100; i++)
    {
        long m;
        mat_ctx_t ctx;
        mat_dense_t A;
        long k;
        char *x, *y;

        m = n_randint(state, 50) + 1;

        mat_ctx_init_long(ctx);
        mat_dense_init(A, m, m, ctx);
        mat_dense_one(A, ctx);

        x = malloc(m * ctx->size);
        for (k = 0; k < m; k++)
            ctx->init(x + k * ctx->size);
        y = malloc(m * ctx->size);
        for (k = 0; k < m; k++)
            ctx->init(y + k * ctx->size);
        for (k = 0; k < m; k++)
            ctx->randtest(x + k * ctx->size, state);

        mat_dense_mul_vec(y, A, x, ctx);

        for (k = 0; k < m; k++)
            if (!ctx->equal(x + k * ctx->size, y + k * ctx->size))
                break;

        result = (k == m);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("Matrix A\n");
            mat_dense_print(A, ctx);
            printf("\n");
        }

        for (k = 0; k < m; k++)
            ctx->clear(x + k * ctx->size);
        free(x);
        for (k = 0; k < m; k++)
            ctx->clear(y + k * ctx->size);
        free(y);

        mat_dense_clear(A, ctx);
        mat_ctx_clear(ctx);
    }

    /* Check that A * (x + y) == A * x + A * y */

    /* Unmanaged element type (long) */
    for (i = 0; i < 100; i++)
    {
        long m, n;
        mat_ctx_t ctx;
        mat_dense_t A;
        long k;
        char *x, *y, *z1, *z2;

        m = n_randint(state, 50) + 1;
        n = n_randint(state, 50) + 1;

        mat_ctx_init_long(ctx);
        mat_dense_init(A, m, n, ctx);
        mat_dense_randtest(A, state, ctx);

        x = malloc(n * ctx->size);
        for (k = 0; k < n; k++)
            ctx->init(x + k * ctx->size);
        y = malloc(n * ctx->size);
        for (k = 0; k < n; k++)
            ctx->init(y + k * ctx->size);
        for (k = 0; k < n; k++)
            ctx->randtest(x + k * ctx->size, state);
        for (k = 0; k < n; k++)
            ctx->randtest(y + k * ctx->size, state);

        z1 = malloc(m * ctx->size);
        for (k = 0; k < m; k++)
            ctx->init(z1 + k * ctx->size);
        z2 = malloc(m * ctx->size);
        for (k = 0; k < m; k++)
            ctx->init(z2 + k * ctx->size);

        mat_dense_mul_vec(z1, A, x, ctx);
        mat_dense_mul_vec(z2, A, y, ctx);
        for (k = 0; k < m; k++)
            ctx->add(z2 + k * ctx->size, 
                     z1 + k * ctx->size, z2 + k * ctx->size);

        for (k = 0; k < n; k++)
            ctx->add(x + k * ctx->size, x + k * ctx->size, y + k * ctx->size);
        mat_dense_mul_vec(z1, A, x, ctx);

        for (k = 0; k < m; k++)
            if (!ctx->equal(z1 + k * ctx->size, z2 + k * ctx->size))
                break;

        result = (k == m);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("Matrix A\n");
            mat_dense_print(A, ctx);
            printf("\n");
        }

        for (k = 0; k < n; k++)
            ctx->clear(x + k * ctx->size);
        free(x);
        for (k = 0; k < n; k++)
            ctx->clear(y + k * ctx->size);
        free(y);
        for (k = 0; k < m; k++)
            ctx->clear(z1 + k * ctx->size);
        free(z1);
        for (k = 0; k < m; k++)
            ctx->clear(z2 + k * ctx->size);
        free(z2);

        mat_dense_clear(A, ctx);
        mat_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
