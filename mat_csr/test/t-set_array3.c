#include "mat_csr.h"

#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/ulong_extras.h"

int
main(void)
{
    flint_rand_t state;

    printf("set_array3\n");
    printf("----------\n");
    fflush(stdout);

    _randinit(state);

    /*
        Run a single example
     */
    {
        ctx_t ctx;
        mat_csr_t A;

        char *mem;
        long k, len, u;

        ctx_init_long(ctx);

        u   = 2 * sizeof(long) + ctx->size;
        len = 6;
        mem = malloc(len * u);

        /*
            [  0  -4     1 ]
            [  0   0  1  0 ]
            [  2   0  0 -3 ]
            [  0   1  0  0 ]
         */

        k = 0;
        *(long *) (mem + k * u)                    =  0;
        *(long *) (mem + k * u + sizeof(long))     =  1;
        *(long *) (mem + k * u + 2 * sizeof(long)) = -4;
        k ++;
        *(long *) (mem + k * u)                    =  0;
        *(long *) (mem + k * u + sizeof(long))     =  3;
        *(long *) (mem + k * u + 2 * sizeof(long)) =  1;
        k ++;
        *(long *) (mem + k * u)                    =  1;
        *(long *) (mem + k * u + sizeof(long))     =  2;
        *(long *) (mem + k * u + 2 * sizeof(long)) =  1;
        k ++;
        *(long *) (mem + k * u)                    =  2;
        *(long *) (mem + k * u + sizeof(long))     =  0;
        *(long *) (mem + k * u + 2 * sizeof(long)) =  2;
        k ++;
        *(long *) (mem + k * u)                    =  2;
        *(long *) (mem + k * u + sizeof(long))     =  3;
        *(long *) (mem + k * u + 2 * sizeof(long)) = -3;
        k ++;
        *(long *) (mem + k * u)                    =  3;
        *(long *) (mem + k * u + sizeof(long))     =  1;
        *(long *) (mem + k * u + 2 * sizeof(long)) =  1;
        k ++;

        mat_csr_init(A, 4, 4, ctx);

        mat_csr_set_array3(A, mem, len, 0, ctx);

        mat_csr_debug(A, ctx);

        free(mem);
        mat_csr_clear(A, ctx);
        ctx_clear(ctx);
    }
    printf("... ");

    _randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
