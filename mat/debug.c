#include <stdlib.h>
#include <stdio.h>

#include "mat.h"

int mat_debug(const mat_t mat, const mat_ctx_t ctx)
{
    long i;

    printf("m = %ld\n", mat->m);
    printf("n = %ld\n", mat->n);
    printf("entries = { ");
    for (i = 0; i < mat->m * mat->n; i++)
    {
        ctx->print(mat->entries + i * ctx->size);
        printf(" ");
    }
    printf("}\n");
    printf("rows = { ");
    for (i = 0; i < mat->m; i++)
        printf("%p ", mat->rows[i]);
    printf("}");

    return 1;
}
