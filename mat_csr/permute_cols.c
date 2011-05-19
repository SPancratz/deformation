#include "mat_csr.h"

void _mat_csr_permute_cols(long m, long n, long *j, long *p, long *lenr, const long *pi)
{
    long i, k, *pi_inv;

    pi_inv = malloc(n * sizeof(long));

    if (!pi_inv)
    {
        printf("ERROR (_mat_csr_permute_cols).\n\n");
        abort();
    }

    for (k = 0; k < n; k++)
        pi_inv[pi[k]] = k;

    for (i = 0; i < m; i++)
        for (k = p[i]; k < p[i] + lenr[i]; k++)
            j[k] = pi_inv[j[k]];

    free(pi_inv);
}

void mat_csr_permute_cols(mat_csr_t A, const long *pi, const ctx_t ctx)
{
    _mat_csr_permute_cols(A->m, A->n, A->j, A->p, A->lenr, pi);
}

