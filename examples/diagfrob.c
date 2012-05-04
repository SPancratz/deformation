/******************************************************************************

    Copyright (C) 2011, 2012 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>

#include "flint.h"
#include "ulong_extras.h"
#include "gmconnection.h"
#include "diagfrob.h"

int main(void)
{
    /* Example 1 */
    /*
    long n = 3;
    long d = 4;
    long N = 374;
    fmpz_t p = {3L};
    fmpz a[4] = {1, 1, 1, 1};
     */

    /* Example 2 */
    long n = 3;
    long d = 5;
    long N = 355;
    fmpz_t p = {2L};
    fmpz a[4] = {1, 1, 1, 1};

    padic_ctx_t pctx;
    ctx_t ctx;

    mat_t F;

    long i, lenB = gmc_basis_size(n, d);

    padic_ctx_init(pctx, p, N, PADIC_VAL_UNIT);

    ctx_init_padic(ctx, pctx);

    mat_init(F, lenB, lenB, ctx);
    diagfrob(F, a, n, d, ctx, 1);

    mat_print(F, ctx);
    printf("\n\n");

    /* Clean-up */

    fmpz_clear(p);
    mat_clear(F, ctx);
    padic_ctx_clear(pctx);
    ctx_clear(ctx);

    return EXIT_SUCCESS;
}

