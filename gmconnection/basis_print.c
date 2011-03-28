#include "gmconnection.h"

void gmc_basis_print(const mon_t *B, const long *iB, long lenB, long n, long d)
{
    long i, k;

    const long u = (n * (d - 1)) / d;
    const long l = n - u;

    printf("[");
    for (k = 1; k < l; k++)
        printf(" |");
    for ( ; k <= u; k++)
    {
        for (i = iB[k]; i < iB[k + 1]; i++)
        {
            printf(i == iB[k] ? " " : ", ");
            mon_print(B[i], n);
        }
        printf(" |");
    }
    for ( ; k < n - 1; k++)
        printf(" |");
    printf(" ]");
}

