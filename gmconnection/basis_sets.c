#include "gmconnection.h"

void gmc_basis_sets(mon_t **B, long **iB, long *lenB, long *l, long *u, 
                    long n, long d)
{
    long j, k;

    *u = (n * (d - 1)) / d;
    *l = n - *u;

    *lenB = gmc_basis_size(n, d);
    *B    = malloc(*lenB * sizeof(mon_t));

    *iB = malloc((*u + 2) * sizeof(long));
    for (k = 0; k < *l; k++)
        (*iB)[k] = 0;

    j = 0;
    for (k = *l; k <= *u; k++)
    {
        mon_t *L;
        long i, len, var;

        (*iB)[k] = j;
        L = mon_generate_by_degree_invlex(&len, n, k * d - n);
        for (i = 0; i < len; i++)
        {
            for (var = 0; var < n; var++)
                if (mon_get_exp(L[i], var) >= d - 1)
                    break;
            if (var == n)
            {
                mon_init((*B)[j]);
                mon_set((*B)[j], L[i]);
                j++;
            }
        }
        /* Clean up */
        for (i = 0; i < len; i++)
            mon_clear(L[i]);
        free(L);
    }
    (*iB)[*u + 1] = *lenB;
}

