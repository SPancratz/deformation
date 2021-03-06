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
__fmpq_poly_print_pretty(const struct __ctx_struct * ctx, const void *op)
{
    return fmpq_poly_print_pretty(op, "t");
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

    fmpq_mat_struct *C;

    mat_t B, LHS;
    ctx_t ctxB;

    printf("solve_fmpq... \n");
    fflush(stdout);

    /* Example 3-1-1 */
    str = "3  [3 0 0] [0 3 0] [0 0 3] (2  0 1)[1 1 1]";

    /* Example 3-3-2 */
    /* str = "3  [3 0 0] [0 3 0] [0 0 3] (2  0 1)[2 1 0] (2  0 1)[0 2 1] (2  0 1)[1 0 2]";*/

    /* Example ... */
    /* str = "4  [4 0 0 0] [0 4 0 0] [0 0 4 0] [0 0 0 4] (2  0 1)[1 1 1 1]"; */

    n = atoi(str) - 1;
    N = 100;

    ctx_init_fmpz_poly_q(ctxM);
    ctx_init_fmpq_poly(ctxB);

    ctxM->print = &__fmpz_poly_q_print_pretty;
    ctxB->print = &__fmpq_poly_print_pretty;

    mpoly_init(P, n + 1, ctxM);
    mpoly_set_str(P, str, ctxM);

    printf("P = "), mpoly_print(P, ctxM), printf("\n");

    b = gmc_basis_size(n, mpoly_degree(P, -1, ctxM));

    mat_init(M, b, b, ctxM);
    mat_init(B, b, b, ctxB);
    mat_init(LHS, b, b, ctxB);

    gmc_compute(M, &rows, &cols, P, ctxM);

    mat_print(M, ctxM);
    printf("\n");

    gmde_solve_fmpq(&C, N, M, ctxM);
    gmde_convert_soln_fmpq(B, ctxB, C, N);

    printf("B: \n");
    mat_print(B, ctxB);
    printf("\n");

    for (i = 0; i < b; i++)
        for (j = 0; j < b; j++)
        {
            fmpq_poly_derivative(
                (fmpq_poly_struct *) mat_entry(LHS, i, j, ctxB), 
                (fmpq_poly_struct *) mat_entry(B, i, j, ctxB));

            for (k = 0; k < n; k++)
            {
                fmpq_poly_t t1, t2;

                fmpq_poly_init(t1);
                fmpq_poly_init(t2);

                fmpq_poly_set_fmpz_poly(t1, fmpz_poly_q_numref(
                    (fmpz_poly_q_struct *) mat_entry(M, i, k, ctxM)));
                fmpq_poly_set_fmpz_poly(t2, fmpz_poly_q_denref(
                    (fmpz_poly_q_struct *) mat_entry(M, i, k, ctxM)));

                fmpq_poly_inv_series(t2, t2, N);
                fmpq_poly_mul(t1, t1, t2);

                fmpq_poly_addmul(
                    (fmpq_poly_struct *) mat_entry(LHS, i, j, ctxB), t1, 
                    (fmpq_poly_struct *) mat_entry(B, k, j, ctxB));
                fmpq_poly_truncate(
                    (fmpq_poly_struct *) mat_entry(LHS, i, j, ctxB), N - 1);

                fmpq_poly_clear(t1);
                fmpq_poly_clear(t2);
            }
        }

    printf("(d/dt + M) * C (mod t^%d)\n", N - 1);
    mat_print(LHS, ctxB);
    printf("\n");

    mpoly_clear(P, ctxM);
    mat_clear(M, ctxM);
    free(rows);
    free(cols);

    mat_clear(B, ctxB);
    mat_clear(LHS, ctxB);

    for (i = 0; i < N; i++)
        fmpq_mat_clear(C + i);
    free(C);

    ctx_clear(ctxM);
    ctx_clear(ctxB);

    return EXIT_SUCCESS;
}

