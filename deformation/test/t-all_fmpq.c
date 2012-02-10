#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>

#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

#include "mpoly.h"
#include "mat.h"
#include "gmconnection.h"
#include "deformation.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    mat_t F;
    mpoly_t P;
    ctx_t ctxF;
    ctx_t ctxFracQt;
    padic_ctx_t pctx;
    fmpz_t p;

    printf("all_fmpq... ");
    fflush(stdout);

    _randinit(state);

    /*
        Consider the example 

            X^3 + Y^3 + Z^3 + t X Y Z.
     */
    {
        long b, n = 2;

        fmpz_init(p);
        fmpz_set_ui(p, 7);

        padic_ctx_init(pctx, p, 20, PADIC_VAL_UNIT);
        ctx_init_fmpz_poly_q(ctxFracQt);
        ctx_init_padic_poly(ctxF, pctx);

        mpoly_init(P, n + 1, ctxFracQt);

        mpoly_set_str(P, "3  [3 0 0] [0 3 0] [0 0 3] (2  0 1)[1 1 1]", ctxFracQt);
        printf("P = "), mpoly_print(P, ctxFracQt), printf("\n");

        b = gmc_basis_size(n, mpoly_degree(P, -1, ctxFracQt));
        mat_init(F, b, b, ctxF);

        frob_with_precisions_fmpq(F, ctxF, P, ctxFracQt, 30, 2000);

        printf("Matrix F:\n");
        mat_print(F, ctxF);
        printf("\n");

        mat_clear(F, ctxF);
        mpoly_clear(P, ctxFracQt);

        padic_ctx_clear(pctx);
        ctx_clear(ctxF);
        ctx_clear(ctxFracQt);
    }

    _randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

