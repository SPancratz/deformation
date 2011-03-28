#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>

#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

#include "mpoly.h"
#include "mat_csr.h"
#include "gmconnection.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;
    mat_ctx_t ctx;

    printf("all... ");
    fflush(stdout);

    flint_randinit(state);

    mat_ctx_init_fmpz_poly_q(ctx);

    /*
        Manually computing basis sets for 

            X^3 + Y^3 + Z^3 + t X Y Z.
     */
    {
        mpoly_t P;

        mon_t *B;
        long *iB, lenB, l, u, k;
        long n, d;

        printf("\n");
        fflush(stdout);

        mpoly_init(P, 3, ctx);
        mpoly_set_str(P, "3  [3 0 0] [0 3 0] [0 0 3] (2  0 1)[1 1 1]", ctx);

        n = P->n;
        d = mpoly_degree(P, -1, ctx);

        gmc_basis_sets(&B, &iB, &lenB, &l, &u, n, d);

        printf("P = "), mpoly_print(P, ctx), printf("\n");
        printf("n = %ld\n", n);
        printf("d = %ld\n", d);
        printf("l u = %ld %ld\n", l, u);

        for (k = l + 1; k <= u + 1; k++)
        {
            mat_csr_t mat;
            mon_t *rows, *cols;
            long *p;

            p = malloc((n + 1) * sizeof(long));
            gmc_init_auxmatrix(mat, &rows, &cols, p, P, k, ctx);

            printf("k = %ld\n", k);
            mat_csr_print_dense(mat, ctx);
            printf("\n");

            mat_csr_clear(mat, ctx);
            free(rows);
            free(cols);
            free(p);
        }

        mpoly_clear(P, ctx);
        free(B);
        free(iB);
    }

    /*
        Computing Gauss--Manin connection for 

            X^3 + Y^3 + Z^3 + t X Y Z.
     */
    {
        mpoly_t P;

        mat_t M;
        mon_t *rows, *cols;

        printf("\n");
        fflush(stdout);

        mpoly_init(P, 3, ctx);
        mpoly_set_str(P, "3  [3 0 0] [0 3 0] [0 0 3] (2  0 1)[1 1 1]", ctx);

        gmc_compute(M, &rows, &cols, P, ctx);

        printf("M = \n"), mat_print(M, ctx), printf("\n");

        mat_clear(M, ctx);
        free(rows);
        free(cols);

        mpoly_clear(P, ctx);
    }

    mat_ctx_clear(ctx);

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

