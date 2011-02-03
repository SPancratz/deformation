#include "mat_coo.h"

int
main(void)
{
    int i;
    flint_rand_t state;

    printf("set_entry\n");
    printf("---------\n");
    fflush(stdout);

    flint_randinit(state);

    /*
        Run a single example
     */
    {
        mat_ctx_t ctx;
        mat_coo_t A;
        long x;

        mat_ctx_init_long(ctx);

        /*
            [  0  -4     1 ]
            [  0   0  1  0 ]
            [  2   0  0 -3 ]
            [  0   1  0  0 ]
         */

        mat_coo_init(A, 4, 4, ctx);

        x = -4;
        mat_coo_set_entry(A, 0, 1, &x, ctx);
        x = 1;
        mat_coo_set_entry(A, 0, 3, &x, ctx);
        x = 1;
        mat_coo_set_entry(A, 1, 2, &x, ctx);
        x = 2;
        mat_coo_set_entry(A, 2, 0, &x, ctx);
        x = -3;
        mat_coo_set_entry(A, 2, 3, &x, ctx);
        x = 1;
        mat_coo_set_entry(A, 3, 1, &x, ctx);

        printf("Matrix A, debug:\n");
        mat_coo_debug(A, ctx);
        printf("\n");

        printf("Matrix A:\n");
        mat_coo_print_dense(A, ctx);
        printf("\n");

        mat_coo_clear(A, 1, ctx);
        mat_ctx_clear(ctx);
    }

    for (i = 0; i < 1000; i++)
    {
        long m, n, z;
        double d;
        mat_ctx_t ctx;
        mat_coo_t A;

        m = n_randint(state, 100) + 1;
        n = n_randint(state, 100) + 1;
        d = (double) n_randint(state, 101) / (double) 100;

        mat_ctx_init_long(ctx);
        mat_coo_init(A, m, n, ctx);

        mat_coo_randtest(A, state, d, ctx);

        for (z = 0; z < 10; z++)
        {
            long row, col;
            long x;

            row = n_randint(state, m);
            col = n_randint(state, n);
            x   = z_randtest_not_zero(state);
            mat_coo_set_entry(A, row, col, &x, ctx);
        }

        mat_coo_clear(A, 1, ctx);
        mat_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    return EXIT_SUCCESS;
}
