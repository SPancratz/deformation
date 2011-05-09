#include "mat_coo.h"

int mat_coo_debug(const mat_coo_t A, const mat_ctx_t ctx)
{
    long k, u = 2 * sizeof(long) + ctx->size;

    printf("m      = %ld\n", A->m);
    printf("n      = %ld\n", A->n);
    printf("alloc  = %ld\n", A->alloc);
    printf("length = %ld\n", A->length);
    printf("{");
    for (k = 0; k < A->length; k++)
    {
        printf("(%ld ", (long) *(A->list + k * u));
        printf( "%ld ", (long) *(A->list + k * u + sizeof(long)));
        ctx->print(A->list + k * u + 2 * sizeof(long));
        printf(")");
    }
    printf("}");

    return 1;
}

