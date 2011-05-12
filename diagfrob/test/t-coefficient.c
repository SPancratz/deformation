/******************************************************************************

    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>

#include "diagfrob.h"

#include "flint.h"
#include "ulong_extras.h"

/*
    Computes the coefficient $lambda_m$ in two ways, using rational 
    and $p$-adic arithmetic, and then compares the results.
 */

int main(void)
{
    int i, result;
    flint_rand_t state;
   
    printf("coefficient... ");
    fflush(stdout);
   
    flint_randinit(state);

    for (i = 0; i < 500; i++)
    {
        long N;
        long prime;
        fmpz_t p;
        padic_ctx_t ctx;

        long m;
        padic_t a, b;
        mpq_t c;

        fmpz_init(p);
        prime = n_randprime(state, 5, 1);
        fmpz_set_ui(p, prime);
        N = n_randint(state, 100) + 1;
        padic_ctx_init(ctx, p, N, PADIC_SERIES);

        m = n_randint(state, 1000);

        padic_init(a, ctx);
        padic_init(b, ctx);
        mpq_init(c);

        diagfrob_coefficient(a, m, ctx);
        diagfrob_coefficient_mpq(c, m, prime);
        padic_set_mpq(b, c, ctx);

        result = (padic_equal(a, b, ctx));
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("p = %ld\n", prime);
            printf("N = %ld\n", N);
            printf("m = %ld\n", m);
            printf("a = "), padic_print(a, ctx), printf("\n");
            printf("b = "), padic_print(b, ctx), printf("\n");
            gmp_printf("c = %Qd\n", c);
            abort();
        }

        fmpz_clear(p);
        padic_clear(a, ctx);
        padic_clear(b, ctx);
        padic_ctx_clear(ctx);
        mpq_clear(c);
    }

    flint_randclear(state);

    printf("PASS\n");
    return EXIT_SUCCESS;
}

