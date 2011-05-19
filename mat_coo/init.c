#include "mat_coo.h"

void mat_coo_init(mat_coo_t A, long m, long n, const ctx_t ctx)
{
    A->m = m;
    A->n = n;

    A->alloc  = 0;
    A->length = 0;
    A->list   = NULL;
}

void mat_coo_init2(mat_coo_t A, long m, long n, long alloc, const ctx_t ctx)
{
    A->m = m;
    A->n = n;

    if (alloc)
    {
        long k, u = 2 * sizeof(long) + ctx->size;

        A->alloc  = alloc;
        A->length = 0;
        A->list   = malloc(alloc * u);

        if (!(A->list))
        {
            printf("ERROR (mat_coo_init2).  Could not allocate memory.\n\n");
            abort();
        }
    }
    else
    {
        A->alloc  = 0;
        A->length = 0;
        A->list   = NULL;
    }
}

