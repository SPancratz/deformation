/******************************************************************************

    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>

#include "diagfrob.h"

#include "flint.h"
#include "ulong_extras.h"

int main(void)
{
    int i, result;
    flint_rand_t state;
   
    printf("falling_fac_mpq... ");
    fflush(stdout);
   
    _randinit(state);

    /* Verify (u)_r = (u + r - 1)! / (u - 1)! for positive u */
    for (i = 0; i < 10000; i++)
    {
        mpq_t lhs, rhs;
        long u, d = 1, r;

        mpq_init(lhs);
        mpq_init(rhs);

        u = n_randint(state, 500) + 1;
        r = n_randint(state, 500);

        diagfrob_falling_fac_mpq(lhs, u, d, r);
        mpz_fac_ui(mpq_numref(rhs), u + r - 1);
        mpz_fac_ui(mpq_denref(rhs), u - 1);
        mpq_canonicalize(rhs);

        result = (mpq_cmp(lhs, rhs) == 0);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("u d r %ld %ld %ld\n", u, d, r);
            gmp_printf("lhs rhs %Qd %Qd\n", lhs, rhs);
        }

        mpq_clear(lhs);
        mpq_clear(rhs);
    }

    _randclear(state);

    printf("PASS\n");
    return EXIT_SUCCESS;
}

