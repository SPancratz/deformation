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

int main(void)
{
    char *str;  /* String for the input polynomial P */
    mpoly_t P;  /* Input polynomial P */
    int n;      /* Number of variables minus one */
    long K;     /* Required t-adic precision */
    long N, Nw;
    long b;     /* Matrix dimensions */
    long i, j, k;

    mat_t M;
    ctx_t ctxM;

    mon_t *rows, *cols;

    padic_mat_struct *C;
    fmpz_t p;

    fmpz_poly_mat_t B;
    long vB;

    printf("solve... \n");
    fflush(stdout);

    /* Example 3-1-1 */
    /* str = "3  [3 0 0] [0 3 0] [0 0 3] (2  0 1)[1 1 1]"; */

    /* Example 3-3-2 */
    /* str = "3  [3 0 0] [0 3 0] [0 0 3] (2  0 1)[2 1 0] (2  0 1)[0 2 1] (2  0 1)[1 0 2]"; */

    /* Example 4-4-2 */
    /* str = "4  [4 0 0 0] [0 4 0 0] [0 0 4 0] [0 0 0 4] (2  0 1)[3 1 0 0] (2  0 1)[1 0 1 2] (2  0 1)[0 1 0 3]"; */

    /* Example ... */
    /* str = "4  [4 0 0 0] [0 4 0 0] [0 0 4 0] [0 0 0 4] (2  0 1)[1 1 1 1]"; */

    /* Example from AKR */
    str = "4  (1  3)[0 3 0 0] (2  0 3)[0 1 2 0] "
          "(2  0 -1)[1 1 1 0] (2  0 3)[1 1 0 1] "
          "(2  0 -1)[2 1 0 0] [0 0 3 0] (2  0 -1)[1 0 2 0] "
          "(1  2)[0 0 0 3] [3 0 0 0]";

    n  = atoi(str) - 1;
    K  = 616;
    N  = 10;
    Nw = 22;

    fmpz_init(p);
    fmpz_set_ui(p, 5);
    ctx_init_fmpz_poly_q(ctxM);

    ctxM->print   = &__fmpz_poly_q_print_pretty;

    mpoly_init(P, n + 1, ctxM);
    mpoly_set_str(P, str, ctxM);

    printf("P = "), mpoly_print(P, ctxM), printf("\n");

    b = gmc_basis_size(n, mpoly_degree(P, -1, ctxM));

    mat_init(M, b, b, ctxM);
    fmpz_poly_mat_init(B, b, b);
    vB = 0;

    gmc_compute(M, &rows, &cols, P, ctxM);

    mat_print(M, ctxM);
    printf("\n");

    gmde_solve(&C, K, p, N, Nw, M, ctxM);
    gmde_convert_soln(B, &vB, C, K, p);

    printf("Solution to (d/dt + M) C = 0:\n");
    fmpz_poly_mat_print(B, "t");
    printf("vB = %ld\n", vB);

    gmde_check_soln(B, vB, p, N, K, M, ctxM);

    mpoly_clear(P, ctxM);
    mat_clear(M, ctxM);
    free(rows);
    free(cols);
    fmpz_poly_mat_clear(B);
    ctx_clear(ctxM);
    fmpz_clear(p);

    for (i = 0; i < K; i++)
        padic_mat_clear(C + i);
    free(C);

    _fmpz_cleanup();

    return EXIT_SUCCESS;
}

