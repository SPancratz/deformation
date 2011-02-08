#include "mat_dense.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("lup_solve... ");
    fflush(stdout);

    flint_randinit(state);

    /* Check that P * A == L * U */

    /* Managed element type (mpq_t) */
    for (i = 0; i < 100; i++)
    {
        long m;
        mat_ctx_t ctx;
        mat_dense_t A, LU;
        long *pi;
        char *x, *b, *c;
        long k;

        m = n_randint(state, 5) + 1;

        pi = malloc(m * sizeof(long));

        mat_ctx_init_mpq(ctx);
        mat_dense_init(A, m, m, ctx);
        mat_dense_init(LU, m, m, ctx);

        mat_dense_randrank(A, state, m, ctx);
        mat_dense_randops(A, state, 1.5 * m, ctx);

        x = malloc(m * ctx->size);
        for (k = 0; k < m; k++)
            ctx->init(x + k * ctx->size);
        b = malloc(m * ctx->size);
        for (k = 0; k < m; k++)
            ctx->init(b + k * ctx->size);
        c = malloc(m * ctx->size);
        for (k = 0; k < m; k++)
            ctx->init(c + k * ctx->size);

        for (k = 0; k < m; k++)
            ctx->randtest(b + k * ctx->size, state);

        mat_dense_lup_decompose(LU, pi, A, ctx);
        mat_dense_lup_solve(x, LU, pi, b, ctx);

        mat_dense_mul_vec(c, A, x, ctx);

        for (k = 0; k < m; k++)
            if (!ctx->equal(b + k * ctx->size, c + k * ctx->size))
                break;

        result = (k == m);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("Matrix A\n");
            mat_dense_print(A, ctx);
            printf("\n");
            abort();
        }

        free(pi);

        for (k = 0; k < m; k++)
            ctx->clear(x + k * ctx->size);
        for (k = 0; k < m; k++)
            ctx->clear(b + k * ctx->size);
        for (k = 0; k < m; k++)
            ctx->clear(c + k * ctx->size);

        mat_dense_clear(A, ctx);
        mat_dense_clear(LU, ctx);
        mat_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
