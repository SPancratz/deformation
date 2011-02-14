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

    printf("mul_vec... ");
    fflush(stdout);

    flint_randinit(state);

    /* Check that the identity matrix does what it's supposed to do */

    /* Unmanaged element type (long) */
    for (i = 0; i < 100; i++)
    {
        long m;
        mat_ctx_t ctx;
        mat_t A;
        char *x, *y;

        m = n_randint(state, 50) + 1;

        mat_ctx_init_long(ctx);
        mat_init(A, m, m, ctx);
        mat_one(A, ctx);

        x = _vec_init(m, ctx);
        y = _vec_init(m, ctx);
        _vec_randtest(x, m, state, ctx);

        mat_mul_vec(y, A, x, ctx);

        result = _vec_equal(x, y, m, ctx);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("Matrix A\n");
            mat_print(A, ctx);
            printf("\n");
        }

        _vec_clear(x, m, ctx);
        _vec_clear(y, m, ctx);

        mat_clear(A, ctx);
        mat_ctx_clear(ctx);
    }

    /* Check that A * (x + y) == A * x + A * y */

    /* Unmanaged element type (long) */
    for (i = 0; i < 100; i++)
    {
        long m, n;
        mat_ctx_t ctx;
        mat_t A;
        char *x, *y, *z1, *z2;

        m = n_randint(state, 50) + 1;
        n = n_randint(state, 50) + 1;

        mat_ctx_init_long(ctx);
        mat_init(A, m, n, ctx);
        mat_randtest(A, state, ctx);

        x = _vec_init(n, ctx);
        y = _vec_init(n, ctx);
        _vec_randtest(x, n, state, ctx);
        _vec_randtest(y, n, state, ctx);

        z1 = _vec_init(m, ctx);
        z2 = _vec_init(m, ctx);

        mat_mul_vec(z1, A, x, ctx);
        mat_mul_vec(z2, A, y, ctx);
        _vec_add(z2, z1, z2, m, ctx);

        _vec_add(x, x, y, n, ctx);
        mat_mul_vec(z1, A, x, ctx);

        result = _vec_equal(z1, z2, m, ctx);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("Matrix A\n");
            mat_print(A, ctx);
            printf("\n");
        }

        _vec_clear(x, n, ctx);
        _vec_clear(y, n, ctx);
        _vec_clear(z1, m, ctx);
        _vec_clear(z2, m, ctx);

        mat_clear(A, ctx);
        mat_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
