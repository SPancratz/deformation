#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>

#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

#include "mpoly.h"
#include "mat_csr.h"
#include "gmconnection.h"

#define RUNS 1000

int
main(void)
{
    int i, j, result;
    flint_rand_t state;
    mat_ctx_t ctx;

    printf("decompose_poly... ");
    fflush(stdout);

    flint_randinit(state);

    mat_ctx_init_mpq(ctx);

    {
        mpoly_t P;

        mon_t *B;
        long *iB, lenB, l, u, k;
        long n, d;

        printf("\n");
        fflush(stdout);

        mpoly_init(P, 3, ctx);
        mpoly_set_str(P, "3  [3 0 0] [0 3 0] [0 0 3] (2)[1 1 1]", ctx);

        n = P->n;
        d = mpoly_degree(P, -1, ctx);

        gmc_basis_sets(&B, &iB, &lenB, &l, &u, n, d);

        printf("P = "), mpoly_print(P, ctx), printf("\n");
        printf("n = %ld\n", n);
        printf("d = %ld\n", d);
        printf("l u = %ld %ld\n", l, u);
        printf("B = "), gmc_basis_print(B, iB, lenB, n, d), printf("\n");

        for (k = l + 1; k <= u + 1; k++)
        {
            mat_csr_t mat;
            mat_csr_solve_t s;
            mon_t *rows, *cols;
            long *p;

            mpoly_t *A, *D;

            p = malloc((n + 1) * sizeof(long));
            gmc_init_auxmatrix(mat, &rows, &cols, p, P, k, ctx);
            mat_csr_solve_init(s, mat, ctx);

            A = malloc(n * sizeof(mpoly_t));
            for (j = 0; j < n; j++)
                mpoly_init(A[j], n, ctx);

            D = malloc(n * sizeof(mpoly_t));
            for (j = 0; j < n; j++)
                mpoly_init(D[j], n, ctx);

            gmc_derivatives(D, P, ctx);

            printf("k = %ld\n", k);
            printf("[");
            for (i = 0; i < RUNS; i++)
            {
                mpoly_t poly1, poly2, poly3;
                char *zero;

                mpoly_init(poly1, n, ctx);
                mpoly_init(poly2, n, ctx);
                mpoly_init(poly3, n, ctx);

                zero = malloc(ctx->size);
                ctx->init(zero);
                ctx->zero(zero);

                mpoly_randhom(poly1, state, k * d - n, 20, ctx);
                for (j = iB[k]; j < iB[k + 1]; j++)
                    mpoly_set_coeff(poly1, B[j], zero, ctx);

                gmc_decompose_poly(A, poly1, s, rows, cols, p, ctx);

                for (j = 0; j < n; j++)
                    mpoly_addmul(poly2, A[j], D[j], ctx);
                for (j = iB[k]; j < iB[k + 1]; j++)
                    mpoly_set_coeff(poly2, B[j], zero, ctx);

                if (!mpoly_is_zero(poly1, ctx))
                    printf("."), fflush(stdout);

                result = (mpoly_equal(poly1, poly2, ctx));
                if (!result)
                {
                    printf("FAIL:\n\n");
                    printf("poly1 = "), mpoly_print(poly1, ctx), printf("\n");
                    printf("poly2 = "), mpoly_print(poly2, ctx), printf("\n");
                    for (j = 0; j < n; j++)
                        printf("D[%d] = ", j), mpoly_print(D[j], ctx), printf("\n");
                    for (j = 0; j < n; j++)
                        printf("A[%d] = ", j), mpoly_print(A[j], ctx), printf("\n");
                    abort();
                }

                mpoly_clear(poly1, ctx);
                mpoly_clear(poly2, ctx);
                mpoly_clear(poly3, ctx);
                ctx->clear(zero);
                free(zero);
            }
            printf("]\n");

            mat_csr_clear(mat, ctx);
            mat_csr_solve_clear(s, ctx);
            free(rows);
            free(cols);
            free(p);

            for (j = 0; j < n; j++)
                mpoly_clear(A[j], ctx);
            free(A);
            for (j = 0; j < n; j++)
                mpoly_clear(D[j], ctx);
            free(D);
        }

        mpoly_clear(P, ctx);
        free(B);
        free(iB);
    }

    mat_ctx_clear(ctx);

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

