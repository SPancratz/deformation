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

    mat_t B, LHS;
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
    mat_init(LHS, b, b, ctxZpt);

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

    for (i = 0; i < b; i++)
        for (j = 0; j < b; j++)
        {
            padic_poly_derivative((padic_poly_struct *) mat_entry(LHS, i, j, ctxZpt), 
                                  (padic_poly_struct *) mat_entry(B, i, j, ctxZpt), pctx);

            for (k = 0; k < b; k++)
            {
                fmpq_poly_t t1, t2;
                padic_poly_t t;

                fmpq_poly_init(t1);
                fmpq_poly_init(t2);
                padic_poly_init(t);

                fmpq_poly_set_fmpz_poly(t1, fmpz_poly_q_numref(
                    (fmpz_poly_q_struct *) mat_entry(M, i, k, ctxM)));
                fmpq_poly_set_fmpz_poly(t2, fmpz_poly_q_denref(
                    (fmpz_poly_q_struct *) mat_entry(M, i, k, ctxM)));

                fmpq_poly_inv_series(t2, t2, N);
                fmpq_poly_mul(t1, t1, t2);

                padic_poly_set_fmpq_poly(t, t1, pctx);
                padic_poly_mul(t, t, 
                    (padic_poly_struct *) mat_entry(B, k, j, ctxZpt), pctx);
                padic_poly_add((padic_poly_struct *) mat_entry(LHS, i, j, ctxZpt), 
                               (padic_poly_struct *) mat_entry(LHS, i, j, ctxZpt), 
                               t, pctx);

                padic_poly_truncate(
                    (padic_poly_struct *) mat_entry(LHS, i, j, ctxZpt), N - 1, 
                    pctx->p);

                fmpq_poly_clear(t1);
                fmpq_poly_clear(t2);
                padic_poly_clear(t);
            }
        }

    printf("(d/dt + M) * C (mod t^%d)\n", N - 1);
    mat_print(LHS, ctxZpt);
    printf("\n");

    mpoly_clear(P, ctxM);
    mat_clear(M, ctxM);
    free(rows);
    free(cols);

    mat_clear(B, ctxZpt);
    mat_clear(LHS, ctxZpt);

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

