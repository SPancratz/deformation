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
    ctx_t ctx;

    printf("all... ");
    fflush(stdout);

    flint_randinit(state);

    ctx_init_fmpz_poly_q(ctx);

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

        n = P->n - 1;
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

            p = malloc((n + 2) * sizeof(long));
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

    /*
        Computing Gauss--Manin connection for 

            P(0) = W^3 + X^3 + Y^3 + Z^3
            P(1) = W^3 + X^3 + Y^3 + Z^3 
                 + (W + X)(W + 2Y)(W + 3Z) + 3XY(W + X + Z)
            Q(t) = A^3 - (1 - t)P(0) - t P(1)
     */
    if (FLINT_BITS == 64)
    {
        mpoly_t P;
        char *str;

        mat_t M;
        mon_t *rows, *cols;

        printf("\n");
        fflush(stdout);

        str = "5  [3 0 0 0 0] (2  -1 -1) [0 3 0 0 0] (1  -1) [0 0 3 0 0] (1  -1) [0 0 0 3 0] (1  -1) [0 0 0 0 3] (2  0 -1) [0 2 1 0 0] (2  0 -2) [0 2 0 1 0] (2  0 -5) [0 1 1 1 0] (2  0 -3) [0 0 2 1 0] (2  0 -3) [0 2 0 0 1] (2  0 -3) [0 1 1 0 1] (2  0 -6) [0 1 0 1 1] (2  0 -9) [0 0 1 1 1]";

        mpoly_init(P, 5, ctx);
        mpoly_set_str(P, str, ctx);

        gmc_compute(M, &rows, &cols, P, ctx);

        printf("M = \n"), mat_print(M, ctx), printf("\n");

        mat_clear(M, ctx);
        free(rows);
        free(cols);

        mpoly_clear(P, ctx);
    }

    ctx_clear(ctx);

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

