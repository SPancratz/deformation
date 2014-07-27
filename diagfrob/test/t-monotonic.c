/******************************************************************************

    Copyright (C) 2013 Sebastian Pancratz

******************************************************************************/

#include "diagfrob.h"

int main(void)
{
    int i, result;
    flint_rand_t state;

    printf("monotonic... ");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 250; i++)
    {
        const long n = n_randint(state, 3) + 2;             /* n in [2,4] */
        const long d = (n % 2L) ? n_randint(state, 5) + 2   /* d in [2,6] */
                                : n_randint(state, 4) + 3;  /* d in [3,6] */
        const long a = 1;

        const long lenB = gmc_basis_size(n, d);
        long j, N1, N2;
        fmpz_t p;
        padic_ctx_t ctx1, ctx2;
        fmpz *A;
        padic_mat_t F1, F2, F3;

        fmpz_init(p);
        do 
        {
            ulong bits = n_randint(state, 5) + 2;
            *p = n_randprime(state, bits, 1);
        }
        while (d % *p == 0);

        A = _fmpz_vec_init(n + 1);
        for (j = 0; j <= n; j++)
            fmpz_set_ui(A + j, n_randint(state, *p - 1) + 1);

        /* 1 <= N1 <= N2 <= 100 */
        N2 = n_randint(state, 50) + 1;
        N1 = n_randint(state, N2) + 1;

        padic_ctx_init(ctx1, p, FLINT_MAX(0, N1-5), N1, PADIC_SERIES);
        padic_ctx_init(ctx2, p, FLINT_MAX(0, N2-5), N2, PADIC_SERIES);

        padic_mat_init2(F1, lenB, lenB, N1);
        padic_mat_init2(F2, lenB, lenB, N2);
        padic_mat_init2(F3, lenB, lenB, N1);

        diagfrob(F1, A, n, d, N1, ctx1, 0);
        diagfrob(F2, A, n, d, N2, ctx2, 0);

        padic_mat_set(F3, F2, ctx1);

        result = padic_mat_equal(F1, F3);
        if (!result)
        {
            printf("FAIL:\n");
            printf("A  = {"), _fmpz_vec_print(A, n + 1), printf("}\n");
            printf("d  = %ld\n", d);
            printf("n  = %ld\n", n);
            printf("p  = %ld\n", *p);
            printf("N1 = %ld\n", N1);
            printf("N2 = %ld\n", N2);
            printf("F1 = \n");
            padic_mat_print_pretty(F1, ctx1);
            printf("F2 = \n");
            padic_mat_print_pretty(F2, ctx2);
            printf("F3 = \n");
            padic_mat_print_pretty(F3, ctx1);
            abort();
        }

        /* Clean-up */
        _fmpz_vec_clear(A, n + 1);
        padic_mat_clear(F1);
        padic_mat_clear(F2);
        padic_mat_clear(F3);
        padic_ctx_clear(ctx1);
        padic_ctx_clear(ctx2);
        fmpz_clear(p);
    }

    flint_randclear(state);

    printf("PASS\n");
    return EXIT_SUCCESS;
}

