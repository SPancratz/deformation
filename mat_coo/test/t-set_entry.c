#include "mat_coo.h"

int
main(void)
{
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

    flint_randclear(state);
    _fmpz_cleanup();
    return EXIT_SUCCESS;
}
