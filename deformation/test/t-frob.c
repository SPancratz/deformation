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

    printf("frob... ");
    fflush(stdout);

    _randinit(state);

    {
        const char *str = 
            "4  [4 0 0 0] [0 4 0 0] [0 0 4 0] [0 0 0 4] (2  0 1)[1 1 1 1]";
        const long n = atoi(str) - 1;

        mpoly_t P;
        ctx_t ctxFracQt;
        fmpz_t p = {7L};

        ctx_init_fmpz_poly_q(ctxFracQt);

        mpoly_init(P, n + 1, ctxFracQt);
        mpoly_set_str(P, str, ctxFracQt);
        printf("P = "), mpoly_print(P, ctxFracQt), printf("\n");

        frob(P, ctxFracQt, p);

        mpoly_clear(P, ctxFracQt);
        ctx_clear(ctxFracQt);
    }

    _randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

