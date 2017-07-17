#include "mat.h"
#include "vec.h"

#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("lup_solve... ");
    fflush(stdout);

    _randinit(state);

    /* Check that P * A == L * U */

    /* Managed element type (mpq_t) */
    for (i = 0; i < 100; i++)
    {
        long m;
        ctx_t ctx;
        mat_t A, LU;
        long *pi;
        char *x, *b, *c;
        long k;

        m = n_randint(state, 5) + 1;

        pi = malloc(m * sizeof(long));

        ctx_init_mpq(ctx);
        mat_init(A, m, m, ctx);
        mat_init(LU, m, m, ctx);

        mat_randrank(A, state, m, ctx);
        mat_randops(A, state, 1.5 * m, ctx);

        x = _vec_init(m, ctx);
        b = _vec_init(m, ctx);
        c = _vec_init(m, ctx);

        _vec_randtest(b, m, state, ctx);

        mat_lup_decompose(LU, pi, A, ctx);
        mat_lup_solve(x, LU, pi, b, ctx);

        mat_mul_vec(c, A, x, ctx);

        result = _vec_equal(b, c, m, ctx);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("Matrix A\n");
            mat_print(A, ctx);
            printf("\n");
            printf("Matrix LU\n");
            mat_print(LU, ctx);
            printf("\n");
            printf("pi = {%ld ", m);
            for (k = 0; k < m; k++)
                printf(" %ld", pi[k]);
            printf("}\n");
            printf("b = {");
            _vec_print(b, m, ctx);
            printf("}\n");
            printf("c = {");
            _vec_print(c, m, ctx);
            printf("}\n");
            printf("x = {");
            _vec_print(x, m, ctx);
            printf("}\n");
            abort();
        }

        free(pi);

        _vec_clear(x, m, ctx);
        _vec_clear(b, m, ctx);
        _vec_clear(c, m, ctx);

        mat_clear(A, ctx);
        mat_clear(LU, ctx);
        ctx_clear(ctx);
    }

    _randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
