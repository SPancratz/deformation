#include "mat_dense.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("lup_decompose... ");
    fflush(stdout);

    flint_randinit(state);

    /* Check that P * A == L * U */

    /* Managed element type (mpq_t) */
    for (i = 0; i < 100; i++)
    {
        long m, n;
        mat_ctx_t ctx;
        mat_dense_t A, L, U, T;
        long i, j;
        long *pi;

        m = n_randint(state, 50) + 1;
        n = m;

        pi = malloc(m * sizeof(long));

        mat_ctx_init_mpq(ctx);
        mat_dense_init(A, m, n, ctx);
        mat_dense_init(L, m, n, ctx);
        mat_dense_init(U, m, n, ctx);
        mat_dense_init(T, m, n, ctx);

        mat_dense_randrank(A, state, m, ctx);
        mat_dense_randops(A, state, 1.5 * m, ctx);

        mat_dense_lup_decompose(L, pi, A, ctx);

        mat_dense_permute_rows(A, pi, ctx);

        for (i = 0; i < m; i++)
        {
            for (j = i; j < n; j++)
                ctx->swap(mat_dense_entry(L, i, j, ctx), 
                          mat_dense_entry(U, i, j, ctx));
            ctx->one(mat_dense_entry(L, i, i, ctx));
        }

        mat_dense_mul_classical(T, L, U, ctx);

        result = mat_dense_equal(A, T, ctx);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("Matrix A\n");
            mat_dense_print(A, ctx);
            printf("\n");
            printf("Matrix L\n");
            mat_dense_print(L, ctx);
            printf("Matrix U\n");
            mat_dense_print(U, ctx);
            printf("Matrix T\n");
            mat_dense_print(T, ctx);
            printf("\n");
        }

        free(pi);

        mat_dense_clear(A, ctx);
        mat_dense_clear(L, ctx);
        mat_dense_clear(U, ctx);
        mat_dense_clear(T, ctx);
        mat_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
