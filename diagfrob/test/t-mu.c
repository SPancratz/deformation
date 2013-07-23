/******************************************************************************

    Copyright (C) 2013 Sebastian Pancratz

******************************************************************************/

#include "diagfrob.h"

int main(void)
{
    int i, result;
    flint_rand_t state;

    /* Common variables, to avoid re-allocation */
    long MMAX;
    fmpz *mu1, *mu2;
    long pbits, pMAX, *C;

    printf("mu... ");
    fflush(stdout);

    MMAX = 1000;
    mu1  = _fmpz_vec_init(MMAX + 1);
    mu2  = _fmpz_vec_init(MMAX + 1);

    pbits = 10;
    pMAX  = (1L << (pbits + 1)) - 1;
    C     = flint_malloc(pMAX * sizeof(long));

    for (i = 0; i < pMAX; i++)
        C[i] = i;

    flint_randinit(state);

    for (i = 0; i < 100; i++)
    {
        mp_limb_t p;
        long M, N;

        p = n_randprime(state, n_randint(state, pbits - 2) + 2, 0);
        M = n_randint(state, MMAX + 1);
        N = n_randint(state, 100) + 1;

        precompute_mu(  mu1, M, C, p, p, N);
        precompute_muex(mu2, M, C, p, p, N);

        result = _fmpz_vec_equal(mu1, mu2, p);
        if (!result)
        {
            printf("FAIL:\n");
            printf("p   = %lu\n", p);
            printf("M   = %ld\n", M);
            printf("N   = %ld\n", N);
            printf("mu1 = {"), _fmpz_vec_print(mu1, M + 1), printf("}\n");
            printf("mu2 = {"), _fmpz_vec_print(mu2, M + 1), printf("}\n");
            abort();
        }
    }

    flint_randclear(state);

    _fmpz_vec_clear(mu1, MMAX + 1);
    _fmpz_vec_clear(mu2, MMAX + 1);
    flint_free(C);

    printf("PASS\n");
    return EXIT_SUCCESS;
}

