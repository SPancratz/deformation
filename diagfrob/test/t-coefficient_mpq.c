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
    Compare the first ten coefficients with the results obtained from 
    the following Sage code

    sage: p = 7
    sage: R.<t> = QQ[]
    sage: K.<pi> = NumberField(t**(p-1) + p)
    sage: s = sum([((pi * (t - t**p))**i / factorial(i)) for i in range(10)])
    sage: L = [s[i] * pi**(-(i % (p-1))) for i in range(10)]
    sage: L
    [1, 1, 1/2, 1/6, 1/24, 1/120, -7/720, -721/720, -5761/5760, -25921/51840]
 */

int main(void)
{
    int i, result;
    flint_rand_t state;
   
    printf("coefficient_mpq... ");
    fflush(stdout);
   
    flint_randinit(state);

    {
        long m, N = 10;
        long p = 7;
        mpq_t *c, *d;

        c = (mpq_t *) malloc(N * sizeof(mpq_t));
        d = (mpq_t *) malloc(N * sizeof(mpq_t));

        for (m = 0; m < N; m++)
        {
            mpq_init(c[m]);
            diagfrob_coefficient_mpq(c[m], m, p);
        }

        for (m = 0; m < N; m++)
            mpq_init(d[m]);
        
        gmp_sscanf("1", "%Qd", d[0]);
        gmp_sscanf("1", "%Qd", d[1]);
        gmp_sscanf("1/2", "%Qd", d[2]);
        gmp_sscanf("1/6", "%Qd", d[3]);
        gmp_sscanf("1/24", "%Qd", d[4]);
        gmp_sscanf("1/120", "%Qd", d[5]);
        gmp_sscanf("-7/720", "%Qd", d[6]);
        gmp_sscanf("-721/720", "%Qd", d[7]);
        gmp_sscanf("-5761/5760", "%Qd", d[8]);
        gmp_sscanf("-25921/51840", "%Qd", d[9]);

        result = 1;
        for (m = 0; m < N; m++)
            result = result && (mpq_cmp(c[m], d[m]) == 0);

        if (!result)
        {
            printf("\n\n");
            for (m = 0; m < N; m++)
                gmp_printf("%ld %Qd %Qd\n", m, c[m], d[m]);
            abort();
        }

        for (m = 0; m < N; m++)
        {
            mpq_clear(c[m]);
            mpq_clear(d[m]);
        }
        free(c);
        free(d);
    }

    flint_randclear(state);

    printf("PASS\n");
    return EXIT_SUCCESS;
}

