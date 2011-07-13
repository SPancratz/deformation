/******************************************************************************

    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>

#include "flint.h"
#include "ulong_extras.h"
#include "gmconnection.h"
#include "diagfrob.h"

#define DEBUG  0

/*
    Computes the matrix $\alpha$, which is closely related to the 
    matrix of Frobenius on the diagonal fibre, using both $p$-padic 
    and rational arithmetic and then compares the results.
 */

int main(void)
{
    int i, result;
    flint_rand_t state;

    /* Data for the dimension/ degree and basis */
    long n = 2, d = 3;
    mon_t *B;
    long *iB, lenB, l, u;
   
    printf("alpha... ");
    fflush(stdout);
   
    flint_randinit(state);

    gmc_basis_sets(&B, &iB, &lenB, &l, &u, n, d);

    for (i = 0; i < 10; i++)
    {
        fmpz_t p;
        long prime, N;
        padic_ctx_t ctx;
        long j, u, v;

        fmpz *a;

        mpq_t *F1;
        padic_t *F2, *F3;

        fmpz_init(p);
        do 
            prime = n_randprime(state, 5, 1);
        while (d % prime == 0);
        fmpz_set_ui(p, prime);
        N = n_randint(state, 50) + 1;

        padic_ctx_init(ctx, p, N, PADIC_SERIES);

        a = _fmpz_vec_init(n + 1);
        for (j = 0; j <= n; j++)
            fmpz_set_ui(a + j, n_randint(state, prime - 1) + 1);

        if (DEBUG)
        {
            printf("  Test case (p = %ld, d = %ld, n = %ld)\n", prime, d, n);
            printf("    a = {");
            _fmpz_vec_print(a, n + 1);
            printf("}\n");
        }

        F1 = (mpq_t *) malloc(lenB * lenB * sizeof(mpq_t));
        F2 = (padic_t *) malloc(lenB * lenB * sizeof(padic_t));
        F3 = (padic_t *) malloc(lenB * lenB * sizeof(padic_t));
        for (j = 0; j < lenB * lenB; j++)
        {
            mpq_init(F1[j]);
            padic_init(F2[j], ctx);
            padic_init(F3[j], ctx);
        }

        for (u = 0; u < lenB; u++)
        {
            for (v = 0; v < lenB; v++)
            {
                diagfrob_alpha_mpq(F1[u*lenB + v], a, n, d, B[u], B[v], ctx);
                diagfrob_alpha(F2[u*lenB + v], a, n, d, B[u], B[v], ctx);
                padic_set_mpq(F3[u*lenB + v], F1[u*lenB + v], ctx);

                result = (padic_equal(F2[u*lenB + v], F3[u*lenB + v], ctx));
                if (!result)
                {
                    printf("FAIL\n\n");
                    printf("u v = %ld %ld\n", u, v);
                    printf("F2(u,v) = "), padic_print(F2[u*lenB + v], ctx), printf("\n");
                    printf("F3(u,v) = "), padic_print(F3[u*lenB + v], ctx), printf("\n");
                    abort();
                }
            }
        }

        for (j = 0; j < lenB * lenB; j++)
        {
            mpq_clear(F1[j]);
            padic_clear(F2[j], ctx);
            padic_clear(F3[j], ctx);
        }
        free(F1);
        free(F2);
        free(F3);

        _fmpz_vec_clear(a, n + 1);
        padic_ctx_clear(ctx);
        fmpz_clear(p);
    }

    free(B);
    free(iB);

    flint_randclear(state);

    printf("PASS\n");
    return EXIT_SUCCESS;
}

