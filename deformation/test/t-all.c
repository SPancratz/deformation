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

#include "padic_mat.h"

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

    printf("all... ");
    fflush(stdout);

    _randinit(state);

    {
        const char *str = "4  [4 0 0 0] [0 4 0 0] [0 0 4 0] [0 0 0 4] (2  0 1)[1 1 1 1]";
        /* const char *str = "3  [3 0 0] [0 3 0] [0 0 3] (2  0 1)[1 1 1]"; */
        const long n = atoi(str) - 1;
        long b;

        fmpz_init(p);
        fmpz_set_ui(p, 7);

        padic_ctx_init(pctx, p, 20, PADIC_VAL_UNIT);
        ctx_init_fmpz_poly_q(ctxFracQt);
        ctx_init_padic_poly(ctxF, pctx);

        mpoly_init(P, n + 1, ctxFracQt);
        mpoly_set_str(P, str, ctxFracQt);
        printf("P = "), mpoly_print(P, ctxFracQt), printf("\n");

        b = gmc_basis_size(n, mpoly_degree(P, -1, ctxFracQt));
        mat_init(F, b, b, ctxF);

        frob_with_precisions(F, ctxF, P, ctxFracQt, 20, 70, 0);

        printf("Matrix r(t)^{3*21} F(t):\n");
        mat_print(F, ctxF);
        printf("\n");
        {
            long i, j;

            printf("Ad hoc prettier output:\n");
            for (i = 0; i < 2; i++)
                for (j = 0; j < 2; j++)
                {
                    printf("(%ld, %ld): ", i, j);
                    padic_poly_print_pretty(
                        (padic_poly_struct *) mat_entry(F, i, j, ctxF), "t", pctx);
                    printf("\n");
                }
        }

        {
            long i, j;
            padic_t one, x;
            fmpz_t y, z;

            padic_mat_t F1;

            padic_mat_init(F1, b, b);
            _padic_init(one);
            _padic_init(x);
            fmpz_set_si(padic_unit(one), -1);

            fmpz_init(z);
            fmpz_set_ui(z, 26);
            fmpz_init(y);
            fmpz_pow_ui(y, z, 3 * 21);

            for (i = 0; i < b; i++)
                for (j = 0; j < b; j++)
                {
                    padic_poly_evaluate_padic(x, (padic_poly_struct *) mat_entry(F, i, j, ctxF), one, ctxF->pctx);
                    padic_mat_set_entry_padic(F1, i, j, x, ctxF->pctx);
                }

            padic_mat_scalar_div_fmpz(F1, F1, y, ctxF->pctx);

            fmpz_clear(y);
            fmpz_clear(z);

            printf("Matrix F1:\n");
            padic_mat_print_pretty(F1, ctxF->pctx);
            printf("\n\n");

            _padic_clear(x);
            _padic_clear(one);
            padic_mat_clear(F1);
        }

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

