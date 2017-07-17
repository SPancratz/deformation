#include "mat_csr.h"

#include "mat.h"

#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/ulong_extras.h"

int
main(void)
{
    int i;
    flint_rand_t state;

    int c;
    __ctx_struct *ctx;
    ctx_t *CTX;

    printf("set_mat... ");
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
            long m, n;
            double d;
            mat_csr_t A;
            mat_t B;

            m = n_randint(state, 100) + 1;
            n = n_randint(state, 100) + 1;
            d = (double) n_randint(state, 101) / (double) 100;

            mat_csr_init(A, m, n, ctx);
            mat_init(B, m, n, ctx);

            mat_randtest(B, state, ctx);
            mat_csr_set_mat(A, B, ctx);

            mat_clear(B, ctx);
            mat_csr_clear(A, ctx);
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
