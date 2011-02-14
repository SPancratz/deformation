#include "mat_csr.h"

#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("zfdiagonal\n");
    printf("----------\n");
    fflush(stdout);

    flint_randinit(state);

    /* Run a single example, with output(?) */
    {
        long m;
        mat_csr_t A;
        mat_ctx_t ctx;

        char *mem;
        long k, len, u;

        long *pi, numnz;

        mat_ctx_init_long(ctx);
        m = 4;

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

        mat_csr_init(A, m, m, ctx);

        mat_csr_set_array3(A, mem, len, 0, ctx);

        printf("Matrix A:\n");
        mat_csr_print_dense(A, ctx);
        printf("\n");

        pi = malloc(m * sizeof(long));

        numnz = mat_csr_zfdiagonal(pi, A);

        printf("pi = {");
        _perm_print(pi, m);
        printf("}\n");
        printf("numnz = %ld\n", numnz);

        printf("Matrix A:\n");
        mat_csr_permute_rows(A, pi, ctx);
        mat_csr_print_dense(A, ctx);
        printf("\n");

        free(mem);
        free(pi);

        mat_csr_clear(A, ctx);
        mat_ctx_clear(ctx);
    }

    /* Permutation matrices */
    for (i = 0; i < 1000; i++)
    {
        long m;
        mat_csr_t A;
        mat_ctx_t ctx;

        char *mem;
        long k, len, u;

        long *pi, numnz;

        mat_ctx_init_long(ctx);
        m = n_randint(state, 100) + 1;

        u   = 2 * sizeof(long) + ctx->size;
        len = m;
        mem = malloc(len * u);

        for (k = 0; k < len; k++)
        {
            *(long *) (mem + k * u)                    = k;
            *(long *) (mem + k * u + sizeof(long))     = k;
            *(long *) (mem + k * u + 2 * sizeof(long)) = 
                z_randtest_not_zero(state);
        }

        pi = malloc(m * sizeof(long));

        _perm_randtest(pi, m, state);

        mat_csr_init(A, m, m, ctx);
        mat_csr_set_array3(A, mem, len, 0, ctx);
        mat_csr_permute_rows(A, pi, ctx);

        numnz = mat_csr_zfdiagonal(pi, A);
        result = (numnz == m);

        if (!result)
        {
            printf("FAIL:\n\n");
            printf("Matrix A:\n");
            mat_csr_debug(A, ctx);
            printf("numnz = %ld\n", numnz);
            abort();
        }

        free(mem);
        free(pi);

        mat_csr_clear(A, ctx);
        mat_ctx_clear(ctx);
    }

    /* Matrices with at most m entries */
    for (i = 0; i < 1000; i++)
    {
        long m;
        mat_csr_t A;
        mat_ctx_t ctx;

        char *mem;
        long k, len, u;

        long *pi, numnz, count;

        mat_ctx_init_long(ctx);
        m = n_randint(state, 100) + 1;

        u   = 2 * sizeof(long) + ctx->size;
        len = m;
        mem = malloc(len * u);

        count = 0;
        for (k = 0; k < len; k++)
        {
            long o;

            o = z_randtest(state);
            count += (o != 0L);

            *(long *) (mem + k * u)                    = k;
            *(long *) (mem + k * u + sizeof(long))     = k;
            *(long *) (mem + k * u + 2 * sizeof(long)) = o;
        }

        pi = malloc(m * sizeof(long));

        _perm_randtest(pi, m, state);

        mat_csr_init(A, m, m, ctx);
        mat_csr_set_array3(A, mem, len, 0, ctx);
        mat_csr_permute_rows(A, pi, ctx);

        numnz = mat_csr_zfdiagonal(pi, A);
        result = (numnz == count);

        if (!result)
        {
            printf("FAIL:\n\n");
            printf("Matrix A:\n");
            mat_csr_debug(A, ctx);
            printf("numnz = %ld\n", numnz);
            abort();
        }

        free(mem);
        free(pi);

        mat_csr_clear(A, ctx);
        mat_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
