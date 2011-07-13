#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>

#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

#include "mpoly.h"
#include "mat.h"
#include "gmconnection.h"
#include "diagfrob.h"
#include "deformation.h"

#include "padic_mat.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    mpoly_t P;
    mat_t F0;
    ctx_t ctxFracQt;
    ctx_t ctxZp;
    padic_ctx_t pctx;
    fmpz_t p;

    printf("diagfrob... ");
    fflush(stdout);

    flint_randinit(state);

    /*
        Consider the example 

            W^4 + X^4 + Y^4 + Z^4 + t W X Y Z
     */
    {
        char *str = "4  [4 0 0 0] [0 4 0 0] [0 0 4 0] [0 0 0 4] (2  0 1)[1 1 1 1]";

        long b, d, n = atoi(str) - 1;

        fmpz_init(p);
        fmpz_set_ui(p, 3);

        padic_ctx_init(pctx, p, 374, PADIC_VAL_UNIT);
        ctx_init_padic(ctxZp, pctx);
        ctx_init_fmpz_poly_q(ctxFracQt);

        mpoly_init(P, n + 1, ctxFracQt);

        mpoly_set_str(P, str, ctxFracQt);
        printf("P = "), mpoly_print(P, ctxFracQt), printf("\n");

        d = mpoly_degree(P, -1, ctxFracQt);
        b = gmc_basis_size(n, d);

        mat_init(F0, b, b, ctxZp);

        /* Find F(0) */
        {
            fmpz * a = _fmpz_vec_init(n + 1);

            mpoly_diagonal_fibre(a, P, ctxFracQt);

            diagfrob(F0, a, n, d, ctxZp);

            printf("a = {"), _fmpz_vec_print(a, n + 1), printf("}\n");

            _fmpz_vec_clear(a, n + 1);
        }

        printf("Matrix F(0):\n");
        mat_print(F0, ctxZp);
        printf("\n");

        mat_clear(F0, ctxZp);
        mpoly_clear(P, ctxFracQt);
        padic_ctx_clear(pctx);
        ctx_clear(ctxZp);
        ctx_clear(ctxFracQt);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

