#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>

#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

#include "mpoly.h"
#include "mat_csr.h"
#include "gmconnection.h"
#include "fmpz_poly.h"

static int 
_mat_print_magma(char ** const rows, long m, long n, const ctx_t ctx)
{
    long i, j;

    printf("[");
    for (i = 0; i < m; i++)
    {
        printf("[ ");
        for (j = 0; j < n; j++)
        {
            fmpz_poly_q_print_pretty((fmpz_poly_q_struct *) (rows[i] + j * ctx->size), "t");
            if (j != n - 1)
                printf(", ");
        }
        printf("]");
        if (i != m - 1)
            printf(", ");
    }
    printf("]\n");

    return 1;
}

static int 
mat_print_magma(const mat_t mat, const ctx_t ctx)
{
    return _mat_print_magma(mat->rows, mat->m, mat->n, ctx);
}

static int
fmpz_poly_mat_print_latex(const fmpz_poly_mat_t mat)
{
    long i, j;

    for (i = 0; i < mat->r; i++)
    {
        for (j = 0; j < mat->c; j++)
        {
            fmpz_poly_print_pretty(fmpz_poly_mat_entry(mat, i, j), "t");
            if (j < mat->c - 1)
            {
                printf(" & ");
            }
        }
        if (i < mat->r - 1)
        {
            printf(" \\\\\n");
        }
    }

    return 1;
}

int main(int argc, const char* argv[])
{
    ctx_t ctx;

    mpoly_t P;
    mat_t M;
    mon_t *rows, *cols;
    long b, n;

    if (argc != 2)
    {
        printf("Syntax: gmconnection_pretty <polynomial>\n");
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
    fmpz_poly_mat_t numM;
    fmpz_poly_t denM;

    fmpz_poly_mat_init(numM, M->m, M->n);
    fmpz_poly_init(denM);

    gmc_convert(numM, denM, M, ctx);

    fmpz_poly_mat_print_latex(numM);

    fmpz_poly_mat_clear(numM);
    fmpz_poly_clear(denM);
}

    mat_print(M, ctx);
    printf("\n");

    mat_print_magma(M, ctx);
    printf("\n");

    mpoly_clear(P, ctx);
    mat_clear(M, ctx);
    free(rows);
    free(cols);

    ctx_clear(ctx);

    _fmpz_cleanup();
    return EXIT_SUCCESS;
}

