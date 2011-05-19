#include "mat_coo.h"

void mat_coo_realloc(mat_coo_t A, long alloc, const ctx_t ctx)
{
    long k, u = 2 * sizeof(long) + ctx->size;

    if (alloc == 0)
    {
        mat_coo_clear(A, 1, ctx);
        mat_coo_init(A, A->m, A->n, ctx);
        return;
    }

    if (A->alloc)
    {
        for (k = alloc; k < A->length; k++)
            ctx->clear(A->list + k * u + 2 * sizeof(long));

        A->list = realloc(A->list, alloc * u);

        if (!(A->list))
        {
            printf("ERROR (mat_coo_realloc).  Could not allocate memory.\n\n");
            abort();
        }

        if (alloc < A->length)
            A->length = alloc;
        A->alloc = alloc;
    }
    else
    {
        A->alloc  = alloc;
        A->length = 0;
        A->list   = malloc(alloc * u);

        if (!(A->list))
        {
            printf("ERROR (mat_coo_realloc).  Could not allocate memory.\n\n");
            abort();
        }
    }
}

