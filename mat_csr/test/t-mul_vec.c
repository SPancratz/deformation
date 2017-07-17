#include "mat_csr.h"
#include "vec.h"

#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    int c;
    __ctx_struct *ctx;
    ctx_t *CTX;

    printf("mul_vec... ");
    fflush(stdout);

    _randinit(state);

    CTX = malloc(2 * sizeof(ctx_t));
    ctx_init_long(CTX[0]);
    ctx_init_mpq(CTX[1]);

    for (c = 0; c < 2; c++)
    {
        ctx = CTX[c];

        /* A (x1 + x2) == A x1 + A x2 */
        for (i = 0; i < 100; i++)
        {
            mat_csr_t A;
            long m, n;
            char *x1, *x2, *y, *z;

            m = n_randint(state, 100) + 1;
            n = n_randint(state, 100) + 1;

            mat_csr_init(A, m, n, ctx);
            x1 = _vec_init(n, ctx);
            x2 = _vec_init(n, ctx);
            y = _vec_init(m, ctx);
            z = _vec_init(m, ctx);

            mat_csr_randtest(A, state, 0.3, ctx);

            mat_csr_mul_vec(y, A, x1, ctx);
            mat_csr_mul_vec(z, A, x2, ctx);
            _vec_add(z, y, z, m, ctx);

            _vec_add(x1, x1, x2, n, ctx);
            mat_csr_mul_vec(y, A, x1, ctx);

            result = _vec_equal(y, z, m, ctx);
            if (!result)
            {
                printf("FAIL:\n\n");
                abort();
            }

            mat_csr_clear(A, ctx);
            _vec_clear(x1, n, ctx);
            _vec_clear(x2, n, ctx);
            _vec_clear(y, m, ctx);
            _vec_clear(z, m, ctx);
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
