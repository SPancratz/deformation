#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>

#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

#include "mpoly.h"
#include "mat_csr.h"
#include "gmconnection.h"

int main(int argc, const char* argv[])
{
    ctx_t ctx;

    mpoly_t P;
    mat_t M;
    mon_t *rows, *cols;
    long b, n;

    if (argc != 2)
    {
        printf("Syntax: gmconnection <polynomial>\n");
        fflush(stdout);
        return EXIT_FAILURE;
    }

    ctx_init_fmpz_poly_q(ctx);

    n = atoi(argv[1]) - 1;
    mpoly_init(P, n, ctx);
    mpoly_set_str(P, argv[1], ctx);

    b = gmc_basis_size(n, mpoly_degree(P, -1, ctx));
    mat_init(M, b, b, ctx);
    gmc_compute(M, &rows, &cols, P, ctx);

    {
        long i, j;

        printf("Row index set:\n");
        for (i = 0; i < M->m; i++)
        {
            mon_print_pretty(rows[i], n + 1, NULL);
            if (i != M->m - 1)
                printf(", ");
        }
        printf("\n");
        printf("Column index set:\n");
        for (j = 0; j < M->n; j++)
        {
            mon_print_pretty(cols[j], n + 1, NULL);
            if (j != M->n - 1)
                printf(", ");
        }
        printf("\n");
    }

    printf("Connection matrix:\n");
    mat_print(M, ctx);
    printf("\n");

    mpoly_clear(P, ctx);
    mat_clear(M, ctx);
    free(rows);
    free(cols);

    ctx_clear(ctx);

    _fmpz_cleanup();
    return EXIT_SUCCESS;
}

