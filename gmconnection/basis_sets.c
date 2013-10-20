/******************************************************************************

    Copyright (C) 2013 Sebastian Pancratz
 
******************************************************************************/

#include "gmconnection.h"

void gmc_basis_sets(mon_t **B, long **iB, long *lenB, long *l, long *u, 
                    long n, long d)
{
    long j, k;

    *u = ((n + 1) * (d - 1)) / d;
    *l = (n + 1) - *u;

    *lenB = gmc_basis_size(n, d);
    *B    = malloc(*lenB * sizeof(mon_t));
    *iB   = malloc((n + 2) * sizeof(long));

    for (k = 0; k < *l; k++)
        (*iB)[k] = 0;

    j = 0;
    for (k = *l; k <= *u; k++)
    {
        mon_t *L;
        long i, len, var;

        (*iB)[k] = j;
        L = mon_generate_by_degree_invlex(&len, n + 1, k * d - (n + 1));
        for (i = 0; i < len; i++)
        {
            for (var = 0; var <= n; var++)
                if (mon_get_exp(L[i], var) >= d - 1)
                    break;
            if (var == n + 1)
            {
                mon_init((*B)[j]);
                mon_set((*B)[j], L[i]);
                j++;
            }
        }
        for (i = 0; i < len; i++)
            mon_clear(L[i]);
        free(L);
    }

    for (k = *u + 1; k <= n + 1; k++)
        (*iB)[k] = *lenB;
}

