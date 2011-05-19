#include "mat.h"

void _mat_permute_rows(char **rows, long m, const long *pi, char **w)
{
    long i;

    for (i = 0; i < m; i++)
        w[i] = rows[i];

    for (i = 0; i < m; i++)
        rows[i] = w[pi[i]];
}

void mat_permute_rows(mat_t mat, const long *pi, const ctx_t ctx)
{
    char **w;

    w = malloc(mat->m * sizeof(char *));

    if (!w)
    {
        printf("ERROR (mat_permute_rows).\n\n");
        abort();
    }

    _mat_permute_rows(mat->rows, mat->m, pi, w);

    free(w);
}

