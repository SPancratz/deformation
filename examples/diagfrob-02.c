/******************************************************************************

    Copyright (C) 2011, 2012 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>

#include "flint/flint.h"
#include "flint/ulong_extras.h"
#include "gmconnection.h"
#include "diagfrob.h"

int main(void)
{
    /* Example 2 */
    long n = 3;
    long d = 5;
    long N = 370;
    fmpz_t p = {2L};
    fmpz a[4] = {1, 1, 1, 1};

    padic_ctx_t pctx;

    padic_mat_t F;

    long i, lenB = gmc_basis_size(n, d);

    padic_ctx_init(pctx, p, FLINT_MAX(0, N), N, PADIC_VAL_UNIT);

    padic_mat_init2(F, lenB, lenB, N);
    diagfrob(F, a, n, d, N, pctx, 1);

    padic_mat_print_pretty(F, pctx);
    printf("\n\n");

    /* Clean-up */
    fmpz_clear(p);
    padic_mat_clear(F);
    padic_ctx_clear(pctx);

    return EXIT_SUCCESS;
}

