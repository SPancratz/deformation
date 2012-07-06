/******************************************************************************

    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>

#include "generics.h"
#include "mat.h"
#include "gmconnection.h"
#include "gmde.h"

#include "flint.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"
#include "fmpz_poly_q.h"

static int 
__fmpz_poly_q_print_pretty(const struct __ctx_struct * ctx, const void *op)
{
    return fmpz_poly_q_print_pretty(op, "t");
}

static int 
__padic_poly_print_pretty(const struct __ctx_struct * ctx, const void *op)
{
    return padic_poly_print_pretty(op, "t", ctx->pctx);
}

int main(void)
{
    char *str;  /* String for the input polynomial P */
    mpoly_t P;  /* Input polynomial P */
    int n;      /* Number of variables minus one */
    int N;      /* Required t-adic precision */
    long b;     /* Matrix dimensions */
    long i, j, k;

    mat_t M;
    ctx_t ctxM;

    mon_t *rows, *cols;

    padic_mat_struct *C;
    padic_ctx_t pctx;
    fmpz_t p;

    mat_t B;
    ctx_t ctxZpt;

    printf("solve... \n");
    fflush(stdout);

    /* Example 3-1-1 */
    str = "3  [3 0 0] [0 3 0] [0 0 3] (2  0 1)[1 1 1]";

    /* Example 3-3-2 */
    /* str = "3  [3 0 0] [0 3 0] [0 0 3] (2  0 1)[2 1 0] (2  0 1)[0 2 1] (2  0 1)[1 0 2]"; */

    /* Example 4-4-2 */
    /* str = "4  [4 0 0 0] [0 4 0 0] [0 0 4 0] [0 0 0 4] (2  0 1)[3 1 0 0] (2  0 1)[1 0 1 2] (2  0 1)[0 1 0 3]"; */

    /* Example ... */
    /* str = "4  [4 0 0 0] [0 4 0 0] [0 0 4 0] [0 0 0 4] (2  0 1)[1 1 1 1]"; */

    n = atoi(str) - 1;
    N = 100;

    fmpz_init(p);
    fmpz_set_ui(p, 7);
    padic_ctx_init(pctx, p, 100, PADIC_VAL_UNIT);
    ctx_init_padic_poly(ctxZpt, pctx);
    ctx_init_fmpz_poly_q(ctxM);

    ctxZpt->print = &__padic_poly_print_pretty;
    ctxM->print   = &__fmpz_poly_q_print_pretty;

    mpoly_init(P, n + 1, ctxM);
    mpoly_set_str(P, str, ctxM);

    printf("P = "), mpoly_print(P, ctxM), printf("\n");

    b = gmc_basis_size(n, mpoly_degree(P, -1, ctxM));

    mat_init(M, b, b, ctxM);
    mat_init(B, b, b, ctxZpt);

    gmc_compute(M, &rows, &cols, P, ctxM);

    mat_print(M, ctxM);
    printf("\n");

    C = malloc(N * sizeof(padic_mat_struct));
    for(i = 0; i < N; i++)
        padic_mat_init(C + i, b, b);

    gmde_solve(C, N, pctx, M, ctxM);
    gmde_convert_soln(B, ctxZpt, C, N);

    printf("Solution to (d/dt + M) B = 0:\n");
    mat_print(B, ctxZpt);
    printf("\n");

    gmde_check_soln(B, ctxZpt, N, M, ctxM);

    mpoly_clear(P, ctxM);
    mat_clear(M, ctxM);
    free(rows);
    free(cols);

    mat_clear(B, ctxZpt);

    for (i = 0; i < N; i++)
        padic_mat_clear(C + i);
    free(C);

    ctx_clear(ctxM);
    ctx_clear(ctxZpt);
    padic_ctx_clear(pctx);
    fmpz_clear(p);

    _fmpz_cleanup();

    return EXIT_SUCCESS;
}

