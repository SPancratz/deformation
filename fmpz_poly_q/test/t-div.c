#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>

#include "fmpz_poly_q.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("div... ");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing of a and b */
    for (i = 0; i < 100; i++)
    {
        fmpz_poly_q_t a, b, c;

        fmpz_poly_q_init(a);
        fmpz_poly_q_init(b);
        fmpz_poly_q_init(c);
        fmpz_poly_q_randtest(b, state, n_randint(state, 50), 50, n_randint(state, 50), 50);
        fmpz_poly_q_randtest_not_zero(c, state, n_randint(state, 50), 50, n_randint(state, 50), 50);

        fmpz_poly_q_div(a, b, c);
        fmpz_poly_q_div(b, b, c);

        result = fmpz_poly_q_equal(a, b);
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_poly_q_print(a), printf("\n\n");
            fmpz_poly_q_print(b), printf("\n\n");
            abort();
        }

        fmpz_poly_q_clear(a);
        fmpz_poly_q_clear(b);
        fmpz_poly_q_clear(c);
    }

    /* Check aliasing of a and c */
    for (i = 0; i < 100; i++)
    {
        fmpz_poly_q_t a, b, c;

        fmpz_poly_q_init(a);
        fmpz_poly_q_init(b);
        fmpz_poly_q_init(c);
        fmpz_poly_q_randtest(b, state, n_randint(state, 50), 50, n_randint(state, 50), 50);
        fmpz_poly_q_randtest_not_zero(c, state, n_randint(state, 50), 50, n_randint(state, 50), 50);

        fmpz_poly_q_div(a, b, c);
        fmpz_poly_q_div(c, b, c);

        result = (fmpz_poly_q_equal(a, c));
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_poly_q_print(a), printf("\n\n");
            fmpz_poly_q_print(c), printf("\n\n");
            abort();
        }

        fmpz_poly_q_clear(a);
        fmpz_poly_q_clear(b);
        fmpz_poly_q_clear(c);
    }

    /* Check (c/b)+(d/b) = (c+d)/b */
    for (i = 0; i < 100; i++)
    {
        fmpz_poly_q_t a1, a2, b, c, d;

        fmpz_poly_q_init(a1);
        fmpz_poly_q_init(a2);
        fmpz_poly_q_init(b);
        fmpz_poly_q_init(c);
        fmpz_poly_q_init(d);
        fmpz_poly_q_randtest_not_zero(b, state, n_randint(state, 50), 50, n_randint(state, 50), 50);
        fmpz_poly_q_randtest(c, state, n_randint(state, 50), 50, n_randint(state, 50), 50);
        fmpz_poly_q_randtest(d, state, n_randint(state, 50), 50, n_randint(state, 50), 50);

        fmpz_poly_q_div(a1, c, b);
        fmpz_poly_q_div(a2, d, b);
        fmpz_poly_q_add(a1, a1, a2);

        fmpz_poly_q_add(c, c, d);
        fmpz_poly_q_div(a2, c, b);

        result = fmpz_poly_q_equal(a1, a2);
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_poly_q_print(a1), printf("\n\n");
            fmpz_poly_q_print(a2), printf("\n\n");
            abort();
        }

        fmpz_poly_q_clear(a1);
        fmpz_poly_q_clear(a2);
        fmpz_poly_q_clear(b);
        fmpz_poly_q_clear(c);
        fmpz_poly_q_clear(d);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
