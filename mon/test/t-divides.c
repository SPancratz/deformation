#include <stdlib.h>
#include <stdio.h>

#include "mon.h"
#include "fmpz.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("div/ divides... ");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing of a and c */
    for (i = 0; i < 100; i++)
    {
        int n;
        mon_t a, b, c;

        n = n_randint(state, 4) + 1;

        mon_init(a);
        mon_init(b);
        mon_init(c);
        mon_randtest(a, state, n, 255 / 2);
        mon_randtest(b, state, n, 255 / 2);
        mon_mul(a, a, b);

        mon_div(c, a, b);
        mon_div(a, a, b);

        result = (mon_equal(a, c));
        if (!result)
        {
            printf("FAIL:\n");
            mon_print(a, n), printf("\n\n");
            mon_print(b, n), printf("\n\n");
            mon_print(c, n), printf("\n\n");
            abort();
        }

        mon_clear(a);
        mon_clear(b);
        mon_clear(c);
    }

    /* Check aliasing of b and c */
    for (i = 0; i < 100; i++)
    {
        int n;
        mon_t a, b, c;

        n = n_randint(state, 4) + 1;

        mon_init(a);
        mon_init(b);
        mon_init(c);
        mon_randtest(a, state, n, 255 / 2);
        mon_randtest(b, state, n, 255 / 2);
        mon_mul(a, a, b);

        mon_div(c, a, b);
        mon_div(b, a, b);

        result = (mon_equal(b, c));
        if (!result)
        {
            printf("FAIL:\n");
            mon_print(a, n), printf("\n\n");
            mon_print(b, n), printf("\n\n");
            mon_print(c, n), printf("\n\n");
            abort();
        }

        mon_clear(a);
        mon_clear(b);
        mon_clear(c);
    }

    /* Check divides */
    for (i = 0; i < 100; i++)
    {
        int k, n;
        mon_t a, b, c;

        n = n_randint(state, 4) + 1;
        k = n_randint(state, n);

        mon_init(a);
        mon_init(b);
        mon_init(c);
        mon_randtest(a, state, n, 255 / 2);

        mon_inc_exp(a, k, 1);
        mon_set(b, a);
        mon_set(c, a);
        mon_inc_exp(c, k, 1);

        result = (mon_divides(b, a) && !mon_divides(c, a));
        if (!result)
        {
            printf("FAIL:\n");
            mon_print(a, n), printf("\n\n");
            mon_print(b, n), printf("\n\n");
            mon_print(c, n), printf("\n\n");
            abort();
        }

        mon_clear(a);
        mon_clear(b);
        mon_clear(c);
    }

    /* Check (ab) / b == a */
    for (i = 0; i < 100; i++)
    {
        int k, n;
        mon_t a, b, c, d;

        n = n_randint(state, 4) + 1;
        k = n_randint(state, n);

        mon_init(a);
        mon_init(b);
        mon_init(c);
        mon_init(d);
        mon_randtest(a, state, n, 255 / 2);
        mon_randtest(b, state, n, 255 / 2);

        mon_mul(c, a, b);
        mon_div(d, c, b);

        result = (mon_divides(b, c) && mon_equal(d, a));
        if (!result)
        {
            printf("FAIL:\n");
            mon_print(a, n), printf("\n\n");
            mon_print(b, n), printf("\n\n");
            mon_print(c, n), printf("\n\n");
            mon_print(d, n), printf("\n\n");
            abort();
        }

        mon_clear(a);
        mon_clear(b);
        mon_clear(c);
        mon_clear(d);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
