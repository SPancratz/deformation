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

    printf("scalar_mul... ");
    fflush(stdout);

    flint_randinit(state);

    CTX = malloc(2 * sizeof(ctx_t));
    ctx_init_long(CTX[0]);
    ctx_init_mpq(CTX[1]);

    for (c = 0; c < 2; c++)
    {
        ctx = CTX[c];

        /* x (A + B) ==  x A + x B */
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

            _vec_scalar_mul(C, A, n, x, ctx);
            _vec_scalar_mul(D, B, n, x, ctx);
            _vec_add(D, C, D, n, ctx);

            _vec_add(C, A, B, n, ctx);
            _vec_scalar_mul(C, C, n, x, ctx);

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

        /* (x y) A == x (y A) */
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
            _vec_randtest(x, 1, state, ctx);
            _vec_randtest(y, 1, state, ctx);

            ctx->mul(z, x, y);
            _vec_scalar_mul(C, A, n, z, ctx);

            _vec_scalar_mul(B, A, n, y, ctx);
            _vec_scalar_mul(D, B, n, x, ctx);

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
    ctx_clear(CTX[1]);
    free(CTX);

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
