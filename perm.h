#ifndef PERM_H
#define PERM_H

#include <stdlib.h>
#include <stdio.h>

#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

/* Randomisation *************************************************************/

/*
    Generates a random permutation vector of length $n$.

    This function uses the Knuth shuffle to generate a uniformly 
    random permutation without retries.
 */

static void _long_vec_randperm(long *vec, long n, flint_rand_t state)
{
    long i, j, t;

    for (i = 0; i < n; i++)
        vec[i] = i;

    for (i = n - 1; i > 0; i--)
    {
        j = n_randint(state, i + 1);
        t = vec[i];
        vec[i] = vec[j];
        vec[j] = t;
    }
}

static int _long_vec_print(const long *vec, long len)
{
    long i;

    printf("%ld", len);
    if (len > 0)
    {
        printf(" ");
        for (i = 0; i < len; i++)
            printf(" %ld", vec[i]);
    }

    return 1;
}

#endif

