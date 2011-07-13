/******************************************************************************

    Copyright (C) 2010, 2011 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include <limits.h>
#include <math.h>

#include "generics.h"
#include "mat.h"
#include "gmconnection.h"
#include "gmde.h"

#include "flint.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"
#include "fmpz_poly_q.h"
#include "padic.h"
#include "padic_mat.h"

int main(void)
{

    char *str;  /* String for the input polynomial P */
    mpoly_t P;  /* Input polynomial P */
    long n;     /* Number of variables minus one */
    long N;     /* Required t-adic precision */
    long i, b;

    mat_t M;
    ctx_t ctxM;

    mon_t *rows, *cols;

    padic_mat_struct *C;
    fmpz_t p;
    padic_ctx_t pctx;

    /* Example 3-1-1 */
    /* str = "3  [3 0 0] [0 3 0] [0 0 3] (2  0 1)[1 1 1]"; */

    /* Example 4-4-2 */
    str = "4  [4 0 0 0] [0 4 0 0] [0 0 4 0] [0 0 0 4] (2  0 1)[3 1 0 0] (2  0 1)[1 0 1 2] (2  0 1)[0 1 0 3]";

    /* Example 3-3-6 */
    /* str = "3  [3 0 0] [0 3 0] [0 0 3] (2  0 314)[2 1 0] (2  0 42)[0 2 1] (2  0 271)[1 0 2] (2  0 -23)[1 1 1]"; */

    /* Example 3-3-2 */
    /* str = "3  [3 0 0] [0 3 0] [0 0 3] (2  0 1)[2 1 0] (2  0 1)[0 2 1] (2  0 1)[1 0 2]"; */

    n = atoi(str) - 1;
    N = 200;

    ctx_init_fmpz_poly_q(ctxM);

    mpoly_init(P, n + 1, ctxM);
    mpoly_set_str(P, str, ctxM);

    printf("P = "), mpoly_print(P, ctxM), printf("\n");

    b = gmc_basis_size(n, mpoly_degree(P, -1, ctxM));

    gmc_compute(M, &rows, &cols, P, ctxM);
    mat_print(M, ctxM);
    printf("\n");

    fmpz_init(p);
    fmpz_set_ui(p, 7);
    padic_ctx_init(pctx, p, 40, PADIC_VAL_UNIT);

    C = malloc(N * sizeof(padic_mat_struct));
    for(i = 0; i < N; i++)
        _padic_mat_init(C + i, b, b);

    gmde_solve(C, N, pctx, M, ctxM);

    printf("Valuations\n");

    for (i = 0; i < N; i++)
    {
        long v = padic_mat_val(C + i);

        if (v < LONG_MAX)
            printf("  i = %ld val = %ld val/log(i) = %f\n", i, v, 
                (i > 1) ? (double) v / log(i) : 0);
        else
            printf("  i = %ld val = +infty\n", i);
    }

    fmpz_clear(p);
    padic_ctx_clear(pctx);
    mpoly_clear(P, ctxM);
    mat_clear(M, ctxM);
    for (i = 0; i < N; i++)
        _padic_mat_clear(C + i);
    free(C);
    ctx_clear(ctxM);

    return EXIT_SUCCESS;
}

