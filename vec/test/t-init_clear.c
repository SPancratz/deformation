#include "vec.h"

#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

int
main(void)
{
    int i;
    flint_rand_t state;

    int c;
    __ctx_struct *ctx;
    ctx_t *CTX;

    printf("init/ clear... ");
    fflush(stdout);

    flint_randinit(state);

    CTX = malloc(2 * sizeof(ctx_t));
    ctx_init_long(CTX[0]);
    ctx_init_mpq(CTX[1]);

    for (c = 0; c < 2; c++)
    {
        ctx = CTX[c];

        for (i = 0; i < 100; i++)
        {
            char *A;
            long n;

            n = n_randint(state, 100) + 1;

            A = _vec_init(n, ctx);
            _vec_randtest(A, n, state, ctx);
            _vec_clear(A, n, ctx);
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
