/******************************************************************************

    Copyright (C) 2013 Sebastian Pancratz
 
******************************************************************************/

#include "gmconnection.h"

void gmc_basis_print(const mon_t *B, const long *iB, long lenB, long n, long d)
{
    const long u = ((n + 1) * (d - 1)) / d;
    const long l = (n + 1) - u;

    long i, k;

    printf("[");
    for (k = 1; k < l; k++)
        printf(" |");
    for ( ; k <= u; k++)
    {
        for (i = iB[k]; i < iB[k + 1]; i++)
        {
            printf(i == iB[k] ? " " : ", ");
            mon_print(B[i], n + 1);
        }
        printf(" |");
    }
    for ( ; k < n; k++)
        printf(" |");
    printf(" ]");
}

