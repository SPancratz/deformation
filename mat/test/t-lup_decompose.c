#include "mat.h"

#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("lup_decompose... ");
    fflush(stdout);

    _randinit(state);

    /* Check that P * A == L * U */

    /* Managed element type (mpq_t) */
    for (i = 0; i < 100; i++)
    {
        long m, n;
        ctx_t ctx;
        mat_t A, L, U, T;
        long i, j;
        long *pi;

        m = n_randint(state, 50) + 1;
        n = m;

        pi = malloc(m * sizeof(long));

        ctx_init_mpq(ctx);
        mat_init(A, m, n, ctx);
        mat_init(L, m, n, ctx);
        mat_init(U, m, n, ctx);
        mat_init(T, m, n, ctx);

        mat_randrank(A, state, m, ctx);
        mat_randops(A, state, 1.5 * m, ctx);

        mat_lup_decompose(L, pi, A, ctx);

        mat_permute_rows(A, pi, ctx);

        for (i = 0; i < m; i++)
        {
            for (j = i; j < n; j++)
                ctx->swap(ctx, mat_entry(L, i, j, ctx), 
                          mat_entry(U, i, j, ctx));
            ctx->one(ctx, mat_entry(L, i, i, ctx));
        }

        mat_mul_classical(T, L, U, ctx);

        result = mat_equal(A, T, ctx);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("Matrix A\n");
            mat_print(A, ctx);
            printf("\n");
            printf("Matrix L\n");
            mat_print(L, ctx);
            printf("Matrix U\n");
            mat_print(U, ctx);
            printf("Matrix T\n");
            mat_print(T, ctx);
            printf("\n");
        }

        free(pi);

        mat_clear(A, ctx);
        mat_clear(L, ctx);
        mat_clear(U, ctx);
        mat_clear(T, ctx);
        ctx_clear(ctx);
    }

    _randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
