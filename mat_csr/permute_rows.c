#include "mat_csr.h"

void _mat_csr_permute_rows(long m, long *p, long *lenr, const long *pi)
{
    long i, *v, *w;

    v = malloc(2 * m * sizeof(long));
    w = v + m;

    if (!v)
    {
        printf("ERROR (_mat_csr_permute_rows).\n\n");
        abort();
    }

    for (i = 0; i < m; i++)
    {
        v[i] = p[i];
        w[i] = lenr[i];
    }
    for (i = 0; i < m; i++)
    {
        p[i]    = v[pi[i]];
        lenr[i] = w[pi[i]];
    }

    free(v);
}

void mat_csr_permute_rows(mat_csr_t A, const long *pi, const ctx_t ctx)
{
    _mat_csr_permute_rows(A->m, A->p, A->lenr, pi);
}
