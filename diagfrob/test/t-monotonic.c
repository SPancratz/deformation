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

    for (i = 0; i < 200; i++)
    {
        const long d = 3;
        const long n = 2;
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

        padic_ctx_init(ctx1, p, N1, PADIC_SERIES);
        padic_ctx_init(ctx2, p, N2, PADIC_SERIES);

        padic_mat_init(F1, lenB, lenB);
        padic_mat_init(F2, lenB, lenB);
        padic_mat_init(F3, lenB, lenB);

        diagfrob(F1, A, n, d, ctx1, 0);
        diagfrob(F2, A, n, d, ctx2, 0);

        padic_mat_set(F3, F2);
        padic_mat_reduce(F3, ctx1);

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

