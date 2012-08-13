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

    _randinit(state);

    {
        const char *str = 
            "3  [5 0 0] [0 5 0] [0 0 5] (2  0 1)[1 1 3]";
        const long n = atoi(str) - 1;

        mpoly_t P;
        ctx_t ctxFracQt;
        qadic_ctx_t Qq;
        fmpz_t p  = {3L};
        long d    = 40;
        qadic_t t1;
        prec_t prec, prec_in;

/*
prec_in.N0   = 9;
prec_in.N1   = 9;
prec_in.N2   = 9;
prec_in.N3   = 13;
prec_in.N3i  = 14;
prec_in.N3w  = 23;
prec_in.N3iw = 22;
prec_in.N4   = 18;
prec_in.m    = 29;
prec_in.K    = 178;
prec_in.r    = 0;
prec_in.s    = 0;
 */

        ctx_init_fmpz_poly_q(ctxFracQt);
        qadic_ctx_init_conway(Qq, p, d, 1, "X", PADIC_SERIES);

        qadic_init(t1);
        qadic_gen(t1, Qq);

        mpoly_init(P, n + 1, ctxFracQt);
        mpoly_set_str(P, str, ctxFracQt);

        frob(P, ctxFracQt, t1, Qq, &prec, NULL, 1);

        qadic_clear(t1);
        mpoly_clear(P, ctxFracQt);
        ctx_clear(ctxFracQt);
        qadic_ctx_clear(Qq);
    }

    _randclear(state);
    _fmpz_cleanup();
    return EXIT_SUCCESS;
}

