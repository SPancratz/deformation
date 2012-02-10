#include "mat_csr.h"

#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

static int 
_is_block_triangular(const mat_csr_t A, const long * B, long b)
{
    long i, i1, i2, k, p;

    for (k = 1; k < b; k++)
    {
        /*
            Ensure that the entry at (i, j) is zero 
            for all i in [i1, i2) and j in [i2, n). 
         */
        i1 = B[k - 1];
        i2 = B[k];
        for (i = i1; i < i2; i++)
            for (p = A->p[i]; p < A->p[i] + A->lenr[i]; p++)
                if (A->j[p] >= i2)
                    return 0;
    }
    return 1;
}

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("block_triangularise\n");
    printf("-------------------\n");
    fflush(stdout);

    _randinit(state);

    /* Run a single example, with output(?) */
    {
        long m, n;
        mat_csr_t A;
        ctx_t ctx;

        long *P, *Q, *B, nz, b;

        ctx_init_long(ctx);
        m = n = 10;

        P = malloc(n * sizeof(long));
        Q = malloc(n * sizeof(long));
        B = malloc(n * sizeof(long));

        mat_csr_init(A, m, n, ctx);
        mat_csr_randtest(A, state, 0.3, ctx);

        printf("Matrix A:\n");
        mat_csr_print_dense(A, ctx);
        printf("\n");

        nz = mat_csr_zfdiagonal(P, A);

        printf("Matrix A:\n");
        mat_csr_permute_rows(A, P, ctx);
        mat_csr_print_dense(A, ctx);
        printf("\n");

        b = mat_csr_block_triangularise(Q, B, A, ctx);

        printf("Matrix A:\n");
        mat_csr_permute_rows(A, Q, ctx);
        mat_csr_permute_cols(A, Q, ctx);
        mat_csr_print_dense(A, ctx);
        printf("\n");

        printf("Blocks B:");
        _long_vec_print(B, b);
        printf("\n");

        free(P);
        free(Q);
        free(B);

        mat_csr_clear(A, ctx);
        ctx_clear(ctx);
    }
    printf("... ");

    for (i = 0; i < 1000; i++)
    {
        long m, n;
        mat_csr_t A;
        ctx_t ctx;

        long *P, *Q, *B, nz, b;

        ctx_init_long(ctx);
        m = n = 100;

        P = malloc(n * sizeof(long));
        Q = malloc(n * sizeof(long));
        B = malloc(n * sizeof(long));

        mat_csr_init(A, m, n, ctx);
        mat_csr_randtest(A, state, 0.03, ctx);

        nz = mat_csr_zfdiagonal(P, A);

        mat_csr_permute_rows(A, P, ctx);

        b = mat_csr_block_triangularise(Q, B, A, ctx);

        mat_csr_permute_rows(A, Q, ctx);
        mat_csr_permute_cols(A, Q, ctx);

        result = _is_block_triangular(A, B, b);
        if (!result)
        {
            printf("FAIL:\n\n");
            abort();
        }

        free(P);
        free(Q);
        free(B);

        mat_csr_clear(A, ctx);
        ctx_clear(ctx);
    }

    _randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
