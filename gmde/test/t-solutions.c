/******************************************************************************

    Copyright (C) 2010, 2011 Sebastian Pancratz

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

int main(void)
{

    char *str;  /* String for the input polynomial P */
    mpoly_t P;  /* Input polynomial P */
    long n;     /* Number of variables minus one */
    long N;     /* Required t-adic precision */
    long b;     /* Matrix dimensions */
    long i, j, k;

    mat_t M, Mt;
    ctx_t ctxM;

    mon_t *rows, *cols;

    fmpq_mat_struct *C, *Cinv;

    mat_t A, B, Binv;
    ctx_t ctxB;

    /* Example 3-1-1 */
    /* str = "3  [3 0 0] [0 3 0] [0 0 3] (2  0 1)[1 1 1]"; */

    /* Example 3-3-2 */
    str = "3  [3 0 0] [0 3 0] [0 0 3] (2  0 1)[2 1 0] (2  0 1)[0 2 1] (2  0 1)[1 0 2]";

    n = atoi(str) - 1;
    N = 100;

    ctx_init_fmpz_poly_q(ctxM);
    ctx_init_fmpq_poly(ctxB);

    mpoly_init(P, n + 1, ctxM);
    mpoly_set_str(P, str, ctxM);

    printf("P = "), mpoly_print(P, ctxM), printf("\n");

    b = gmc_basis_size(n + 1, mpoly_degree(P, -1, ctxM));

    mat_init(M, b, b, ctxM);
    mat_init(Mt, b, b, ctxM);
    mat_init(A, b, b, ctxB);
    mat_init(B, b, b, ctxB);
    mat_init(Binv, b, b, ctxB);

    gmc_compute(M, &rows, &cols, P, ctxM);
    mat_print(M, ctxM);
    printf("\n");

    mat_transpose(Mt, M, ctxM);
    mat_neg(Mt, Mt, ctxM);

    C    = malloc(N * sizeof(fmpq_mat_struct));
    Cinv = malloc(N * sizeof(fmpq_mat_struct));
    for(i = 0; i < N; i++)
    {
        fmpq_mat_init(C + i, b, b);
        fmpq_mat_init(Cinv + i, b, b);
    }

    gmde_solve_fmpq(C, N, M, ctxM);
    gmde_solve_fmpq(Cinv, N, Mt, ctxM);
    gmde_convert_soln_fmpq(B, ctxB, C, N);
    gmde_convert_soln_fmpq(Binv, ctxB, Cinv, N);

    mat_transpose(Binv, Binv, ctxB);

    mat_mul(A, B, Binv, ctxB);

    printf("A = B * Binv: \n");
    mat_print(A, ctxB);
    printf("\n");

    mpoly_clear(P, ctxM);
    mat_clear(M, ctxM);
    mat_clear(Mt, ctxM);
    free(rows);
    free(cols);

    mat_clear(A, ctxB);
    mat_clear(B, ctxB);
    mat_clear(Binv, ctxB);

    for (i = 0; i < N; i++)
    {
        fmpq_mat_clear(C + i);
        fmpq_mat_clear(Cinv + i);
    }
    free(C);
    free(Cinv);

    ctx_clear(ctxM);
    ctx_clear(ctxB);

    return EXIT_SUCCESS;
}

