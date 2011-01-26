#include "mat_csr.h"

int
main(void)
{
    flint_rand_t state;

    printf("zfdiagonal\n");
    printf("----------\n");
    fflush(stdout);

    flint_randinit(state);

    /*
        Run a single example
     */
    {
        long m, n;
        mat_csr_t A;
        mat_ctx_t ctx;

        char *mem;
        long k, len, u;

        long *pi, numnz;

        mat_ctx_init_long(ctx);
        m = n = 4;

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

        mat_csr_init(A, m, n, ctx);

        mat_csr_set_array3(A, mem, len, 0, ctx);

        mat_csr_print_dense(A, ctx);
        printf("\n");

        pi = malloc(m * sizeof(long));

        numnz = mat_csr_zfdiagonal(pi, A);

        printf("pi = {");
        _ulong_vec_print(pi, m);
        printf("}\n");
        printf("numnz = %ld\n", numnz);

        free(mem);
        free(pi);

        mat_csr_clear(A, ctx);
        mat_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    return EXIT_SUCCESS;
}
