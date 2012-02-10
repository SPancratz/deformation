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

    printf("swap... ");
    fflush(stdout);

    _randinit(state);

    CTX = malloc(2 * sizeof(ctx_t));
    ctx_init_long(CTX[0]);
    ctx_init_mpq(CTX[1]);

    for (c = 0; c < 2; c++)
    {
        ctx = CTX[c];

        for (i = 0; i < 100; i++)
        {
            char *A, *B, *C, *D;
            long n;

            n = n_randint(state, 100) + 1;

            A = _vec_init(n, ctx);
            B = _vec_init(n, ctx);
            C = _vec_init(n, ctx);
            D = _vec_init(n, ctx);
            _vec_randtest(A, n, state, ctx);
            _vec_randtest(B, n, state, ctx);

            _vec_set(C, A, n, ctx);
            _vec_set(D, B, n, ctx);
            _vec_swap(A, B, n, ctx);

            result = _vec_equal(A, D, n, ctx) && _vec_equal(B, C, n, ctx);
            if (!result)
            {
                printf("FAIL:\n\n");
                abort();
            }

            _vec_clear(A, n, ctx);
            _vec_clear(B, n, ctx);
            _vec_clear(C, n, ctx);
            _vec_clear(D, n, ctx);
        }
    }

    ctx_clear(CTX[0]);
    ctx_clear(CTX[1]);
    free(CTX);

    _randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
