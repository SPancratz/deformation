/******************************************************************************

    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

#include "diagfrob.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("flog... ");
    fflush(stdout);

    _randinit(state);

    for (i = 0; i < 1000000; i++)
    {
        fmpz_t a, b, c, d;
        long e;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(d);

        /* a >= 1, b >= 2 */
        fmpz_randtest_not_zero(a, state, 200);
        fmpz_abs(a, a);
        fmpz_randtest_unsigned(b, state, 20);
        while (fmpz_cmp_ui(b, 2) < 0)
            fmpz_randtest_unsigned(b, state, 20);

        e = fmpz_flog(a, b);

        if (e >= 0)
        {
            fmpz_pow_ui(c, b, e);
            fmpz_pow_ui(d, b, e + 1);

            result = (fmpz_cmp(c, a) <= 0) && (fmpz_cmp(a, d) < 0);
        }
        else
        {
            result = 0;
        }

        if (!result)
        {
            printf("FAIL:\n");
            printf("a = "), fmpz_print(a), printf("\n");
            printf("b = "), fmpz_print(b), printf("\n");
            printf("e = %ld\n", e);
            if (e >= 0)
            {
                printf("c = "), fmpz_print(c), printf("\n");
                printf("d = "), fmpz_print(d), printf("\n");
            }
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(d);
    }

    _randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

