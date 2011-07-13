#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include <math.h>

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
    long b, i, n, N;

    time_t t0, t1;
    clock_t c0, c1;
    long double cputime, cputimepr;

    if (argc != 3)
    {
        printf("Syntax: gmconnection-p <polynomial> <integer>\n");
        fflush(stdout);
        return EXIT_FAILURE;
    }

    ctx_init_fmpz_poly_q(ctx);

    n = atoi(argv[1]) - 1;
    mpoly_init(P, n, ctx);
    mpoly_set_str(P, argv[1], ctx);
    b = gmc_basis_size(n, mpoly_degree(P, -1, ctx));

    N = atoi(argv[2]);

    t0 = time(NULL);
    c0 = clock();

    for (i = 0; i < N; i++)
    {
        mat_init(M, b, b, ctx);
        gmc_compute(M, &rows, &cols, P, ctx);

        mat_clear(M, ctx);
        free(rows);
        free(cols);
    }

    t1 = time(NULL);
    c1 = clock();
    cputime   = ((long double) (c1 - c0)) / ((long double) CLOCKS_PER_SEC);
    cputimepr = cputime / ((long double) N);

    printf ("Elapsed wall clock time:  %ld\n", (long) (t1 - t0));
    printf ("Elapsed CPU time:         %LG\n", cputime);
    printf ("Elapsed CPU time per run: %LG\n", cputimepr);

    mpoly_clear(P, ctx);
    ctx_clear(ctx);

    _fmpz_cleanup();
    return EXIT_SUCCESS;
}

