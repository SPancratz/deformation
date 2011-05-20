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

    printf("revcharpoly... ");
    fflush(stdout);

    flint_randinit(state);

    /* Compute a single example over QQ */
    {
        long m, n;
        ctx_t ctx;
        mat_t A;
        char *f;

        m = 2;
        n = 2;

        ctx_init_mpq(ctx);
        mat_init(A, m, n, ctx);
        f = _vec_init(n + 1, ctx);

        ctx->set_si(ctx, mat_entry(A, 0, 0, ctx), 1);
        ctx->set_si(ctx, mat_entry(A, 0, 1, ctx), 2);
        ctx->set_si(ctx, mat_entry(A, 1, 0, ctx), 3);
        ctx->set_si(ctx, mat_entry(A, 1, 1, ctx), 4);

        mat_revcharpoly(f, A, ctx);

        printf("Matrix = \n");
        mat_print(A, ctx);
        printf("\n");
        printf("Reverse charpoly = \n");
        _vec_print(f, n + 1, ctx);
        printf("\n");

        mat_clear(A, ctx);
        _vec_clear(f, n + 1, ctx);
        ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

