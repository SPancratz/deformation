#include "vec.h"

#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    int c;
    __ctx_struct *ctx;
    ctx_t *CTX;

    printf("scalar_div... ");
    fflush(stdout);

    flint_randinit(state);

    CTX = malloc(1 * sizeof(ctx_t));
    ctx_init_mpq(CTX[0]);

    for (c = 0; c < 1; c++)
    {
        ctx = CTX[c];

        /* (A + B) / x ==  A / x + B / x */
        for (i = 0; i < 100; i++)
        {
            char *A, *B, *C, *D, *x;
            long n;

            n = n_randint(state, 100) + 1;

            A = _vec_init(n, ctx);
            B = _vec_init(n, ctx);
            C = _vec_init(n, ctx);
            D = _vec_init(n, ctx);
            x = _vec_init(1, ctx);
            _vec_randtest(A, n, state, ctx);
            _vec_randtest(B, n, state, ctx);
            _vec_randtest(x, 1, state, ctx);

            _vec_scalar_div(C, A, n, x, ctx);
            _vec_scalar_div(D, B, n, x, ctx);
            _vec_add(D, C, D, n, ctx);

            _vec_add(C, A, B, n, ctx);
            _vec_scalar_div(C, C, n, x, ctx);

            result = _vec_equal(C, D, n, ctx);
            if (!result)
            {
                printf("FAIL:\n\n");
                abort();
            }

            _vec_clear(A, n, ctx);
            _vec_clear(B, n, ctx);
            _vec_clear(C, n, ctx);
            _vec_clear(D, n, ctx);
            _vec_clear(x, 1, ctx);
        }

        /* A / (x y) == (A / x) / y */
        for (i = 0; i < 100; i++)
        {
            char *A, *B, *C, *D, *x, *y, *z;
            long n;

            n = n_randint(state, 100) + 1;

            A = _vec_init(n, ctx);
            B = _vec_init(n, ctx);
            C = _vec_init(n, ctx);
            D = _vec_init(n, ctx);
            x = _vec_init(1, ctx);
            y = _vec_init(1, ctx);
            z = _vec_init(1, ctx);
            _vec_randtest(A, n, state, ctx);
            _vec_randtest_not_zero(x, 1, state, ctx);
            _vec_randtest_not_zero(y, 1, state, ctx);

            ctx->mul(ctx, z, x, y);
            _vec_scalar_div(C, A, n, z, ctx);

            _vec_scalar_div(B, A, n, x, ctx);
            _vec_scalar_div(D, B, n, y, ctx);

            result = _vec_equal(C, D, n, ctx);
            if (!result)
            {
                printf("FAIL:\n\n");
                abort();
            }

            _vec_clear(A, n, ctx);
            _vec_clear(B, n, ctx);
            _vec_clear(C, n, ctx);
            _vec_clear(D, n, ctx);
            _vec_clear(x, 1, ctx);
            _vec_clear(y, 1, ctx);
            _vec_clear(z, 1, ctx);
        }

    }

    ctx_clear(CTX[0]);
    free(CTX);

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
