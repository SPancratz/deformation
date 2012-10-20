/*
    Copyright (C) 2010, 2011 Sebastian Pancratz
 */

#ifndef PERM_H
#define PERM_H

#include <stdlib.h>
#include <stdio.h>

#include "flint.h"

/* Memory management *********************************************************/

static long * _perm_init(long n)
{
    long i, *vec;

    vec = malloc(n * sizeof(long));

    if (!vec)
    {
        printf("ERROR (_perm_init).\n\n");
        abort();
    }

    for (i = 0; i < n; i++)
        vec[i] = i;

    return vec;
}

static void _perm_clear(long * vec)
{
    free(vec);
}

/* Assignment ****************************************************************/

static void _perm_set(long *res, const long *vec, long n)
{
    long i;

    for (i = 0; i < n; i++)
        res[i] = vec[i];
}

static void _perm_set_one(long *vec, long n)
{
    long i;

    for (i = 0; i < n; i++)
        vec[i] = i;
}

static void _perm_inv(long *res, const long *vec, long n)
{
    long i;

    if (res == vec)
    {
        long *t = malloc(n * sizeof(long));

        if (!t)
        {
            printf("ERROR (_perm_inv).\n\n");
            abort();
        }

        for (i = 0; i < n; i++)
            t[i] = vec[i];
        for (i = 0; i < n; i++)
            res[t[i]] = i;

        free(t);
    }
    else
    {
        for (i = 0; i < n; i++)
            res[vec[i]] = i;
    }
}

/* Composition ***************************************************************/

static void _perm_compose(long *res, const long *vec1, const long *vec2, long n)
{
    long i;

    for (i = 0; i < n; i++)
        res[i] = vec1[vec2[i]];
}

/* Randomisation *************************************************************/

void _perm_randtest(long * vec, long n, flint_rand_t state);

/* Input and output **********************************************************/

static int _long_vec_print(const long * vec, long len)
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

static int _perm_print(const long * vec, long n)
{
    long i;

    printf("%ld", n);
    if (n > 0)
    {
        printf(" ");
        for (i = 0; i < n; i++)
            printf(" %ld", vec[i]);
    }

    return 1;
}

#endif

