#include "mat_csr.h"
#include "mat.h"
#include "vec.h"

#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("solve... ");
    fflush(stdout);

    flint_randinit(state);

    /* Check that A x == b */

    /* Managed element type (mpq_t) */
    for (i = 0; i < 100; i++)
    {
        long m;
        mat_ctx_t ctx;
        mat_csr_t A;
        mat_csr_solve_t S;
        mat_t B;
        char *x, *b, *c;

        m = n_randint(state, 100) + 1;

        mat_ctx_init_mpq(ctx);
        mat_init(B, m, m, ctx);

        mat_randrank(B, state, m, ctx);
        mat_randops(B, state, 1.5 * m, ctx);

        mat_csr_init(A, m, m, ctx);
        mat_csr_set_mat(A, B, ctx);

        x = _vec_init(m, ctx);
        b = _vec_init(m, ctx);
        c = _vec_init(m, ctx);

        _vec_randtest(b, m, state, ctx);

        mat_csr_solve_init(S, A, ctx);
        mat_csr_solve(x, S, b, ctx);

        mat_csr_mul_vec(c, A, x, ctx);

        result = _vec_equal(b, c, m, ctx);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("Matrix A:\n"), mat_csr_print_dense(A, ctx), printf("\n");
            printf("Vector b = {"), _vec_print(b, m, ctx), printf("}\n");
            printf("Vector x = {"), _vec_print(x, m, ctx), printf("}\n");
            printf("Vector c = {"), _vec_print(c, m, ctx), printf("}\n");
            abort();
        }

        _vec_clear(x, m, ctx);
        _vec_clear(b, m, ctx);
        _vec_clear(c, m, ctx);

        mat_csr_clear(A, ctx);
        mat_csr_solve_clear(S, ctx);
        mat_clear(B, ctx);
        mat_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
