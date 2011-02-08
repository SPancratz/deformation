#include "mat_dense.h"

void _mat_dense_permute_rows(char **rows, long m, const long *pi, char **w)
{
    long i;

    for (i = 0; i < m; i++)
        w[i] = rows[i];

    for (i = 0; i < m; i++)
        rows[i] = w[pi[i]];
}

void mat_dense_permute_rows(mat_dense_t mat, const long *pi, const mat_ctx_t ctx)
{
    char **w;

    w = malloc(mat->m * sizeof(char *));

    if (!w)
    {
        printf("ERROR (mat_dense_permute_rows).\n\n");
        abort();
    }

    _mat_dense_permute_rows(mat->rows, mat->m, pi, w);

    free(w);
}

