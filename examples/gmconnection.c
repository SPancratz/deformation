#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>

#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

#include "mpoly.h"
#include "mat_csr.h"
#include "gmconnection.h"

static 
int _gmc_mat_print(char ** const rows, long m, long n, const ctx_t ctx)
{
    long i, j;

    for (i = 0; i < m; i++)
    {
        printf("[");
        for (j = 0; j < n; j++)
        {
            ctx->print(ctx, rows[i] + j * ctx->size);
            if (j != n - 1)
                printf(" ");
        }
        if (i != m - 1)
            printf("],\n");
        else
            printf("]");
    }

    return 1;
}

static 
int gmc_mat_print(const mat_t mat, const ctx_t ctx)
{
    return _gmc_mat_print(mat->rows, mat->m, mat->n, ctx);
}


int 
main(int argc, const char* argv[])
{
    ctx_t ctx;

    mpoly_t P;
    mat_t M;
    mon_t *rows, *cols;

    if (argc != 2)
    {
        printf("Syntax: gmconnection <polynomial>\n");
        fflush(stdout);
        return EXIT_FAILURE;
    }

    ctx_init_fmpz_poly_q(ctx);

    mpoly_init(P, 1, ctx);
    mpoly_set_str(P, argv[1], ctx);

    gmc_compute(M, &rows, &cols, P, ctx);

    gmc_mat_print(M, ctx), printf("\n");

    mpoly_clear(P, ctx);
    mat_clear(M, ctx);
    free(rows);
    free(cols);

    ctx_clear(ctx);

    _fmpz_cleanup();
    return EXIT_SUCCESS;
}

