/******************************************************************************

    Copyright (C) 2013 Sebastian Pancratz

******************************************************************************/

#include "diagfrob.h"

int main(void)
{
    int i, result;
    flint_rand_t state;

    printf("nmod_mat_charpoly... ");
    fflush(stdout);

    flint_randinit(state);

    /* Verify the charpoly of the empty matrix */
    {
        mp_limb_t mod = n_randtest_prime(state, 0);
        nmod_mat_t A;
        nmod_poly_t f;

        nmod_mat_init(A, 0, 0, mod);
        nmod_poly_init(f, mod);
        nmod_mat_charpoly(f, A);

        result = (f->length == 1) && (f->coeffs[0] == 1);
        if (!result)
        {
            printf("FAIL:\n");
            printf("A:\n"), nmod_mat_print_pretty(A), printf("\n");
            printf("f   = "), nmod_poly_print(f), printf("\n");
            printf("mod = %lu\n", mod);
            abort();
        }
        nmod_mat_clear(A);
        nmod_poly_clear(f);
    }

    /* Verify the trace */
    for (i = 0; i < 1000; i++)
    {
        long n;
        mp_limb_t mod;
        nmod_mat_t A;
        nmod_poly_t f;
        mp_limb_t tr;

        n   = n_randint(state, 20) + 1;
        mod = n_randtest_prime(state, 0);

        nmod_mat_init(A, n, n, mod);
        nmod_mat_randtest(A, state);

        nmod_poly_init(f, mod);
        nmod_mat_charpoly(f, A);
        tr = nmod_mat_trace(A);
        tr = n_negmod(tr, mod);

        result = (f->length == n+1) && (f->coeffs[n-1] == tr);
        if (!result)
        {
            printf("FAIL:\n");
            printf("A:\n"), nmod_mat_print_pretty(A), printf("\n");
            printf("f   = "), nmod_poly_print(f), printf("\n");
            printf("tr  = %lu\n", tr);
            printf("mod = %lu\n", mod);
            abort();
        }
        nmod_mat_clear(A);
        nmod_poly_clear(f);
    }
    flint_randclear(state);

    printf("PASS\n");
    return EXIT_SUCCESS;
}

