#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>

#include "fmpz_poly_q.h"
#include "long_extras.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("scalar_mul_mpz... ");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing of a and b */
    for (i = 0; i < 100; i++)
    {
        fmpz_poly_q_t a, b;
        fmpz_t x;
        mpz_t y;

        fmpz_poly_q_init(a);
        fmpz_poly_q_init(b);
        fmpz_init(x);
        mpz_init(y);

        fmpz_poly_q_randtest(b, state, n_randint(state, 50), 50, n_randint(state, 50), 50);
        fmpz_randtest(x, state, 50);
        fmpz_get_mpz(y, x);

        fmpz_poly_q_scalar_mul_mpz(a, b, y);
        fmpz_poly_q_scalar_mul_mpz(b, b, y);

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
        fmpz_clear(x);
        mpz_clear(y);
    }

    /* Check that x (a + b) == x * a + x * b */
    for (i = 0; i < 100; i++)
    {
        fmpz_poly_q_t a, b, c, d;
        fmpz_t x;
        mpz_t y;

        fmpz_poly_q_init(a);
        fmpz_poly_q_init(b);
        fmpz_poly_q_init(c);
        fmpz_poly_q_init(d);
        fmpz_init(x);
        mpz_init(y);

        fmpz_poly_q_randtest(a, state, n_randint(state, 50), 50, n_randint(state, 50), 50);
        fmpz_poly_q_randtest(b, state, n_randint(state, 50), 50, n_randint(state, 50), 50);
        fmpz_randtest(x, state, 50);
        fmpz_get_mpz(y, x);

        fmpz_poly_q_scalar_mul_mpz(c, a, y);
        fmpz_poly_q_scalar_mul_mpz(d, b, y);
        fmpz_poly_q_add(d, c, d);

        fmpz_poly_q_add(c, a, b);
        fmpz_poly_q_scalar_mul_mpz(c, c, y);

        result = fmpz_poly_q_equal(c, d);
        if (!result)
        {
            printf("FAIL:\n");
            printf("a = "), fmpz_poly_q_print(a), printf("\n\n");
            printf("b = "), fmpz_poly_q_print(b), printf("\n\n");
            printf("c = "), fmpz_poly_q_print(c), printf("\n\n");
            printf("d = "), fmpz_poly_q_print(d), printf("\n\n");
            gmp_printf("y = %Zd\n\n", y);
            abort();
        }

        fmpz_poly_q_clear(a);
        fmpz_poly_q_clear(b);
        fmpz_poly_q_clear(c);
        fmpz_poly_q_clear(d);
        fmpz_clear(x);
        mpz_clear(y);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
