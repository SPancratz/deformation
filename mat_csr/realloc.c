#include "mat_csr.h"

void mat_csr_realloc(mat_csr_t A, long alloc, const mat_ctx_t ctx)
{
    long k;

    if (alloc == 0)
    {
        mat_csr_clear(A, ctx);
        mat_csr_init(A, A->m, A->n, ctx);
        return;
    }

    if (A->alloc)
    {
        for (k = alloc; k < A->alloc; k++)
            ctx->clear(A->x + k * (ctx->size));

        A->x = realloc(A->x, alloc * (ctx->size));
        A->j = realloc(A->j, alloc * sizeof(long));

        if (!(A->x) || !(A->j))
        {
            printf("ERROR (mat_csr_realloc).\n\n");
            abort();
        }

        for (k = A->alloc; k < alloc; k++)
            ctx->init(A->x + k * (ctx->size));
        for (k = A->alloc; k < alloc; k++)
            A->j[k] = 0L;

        A->alloc = alloc;
    }
    else
    {
        A->alloc = alloc;
        A->x     = calloc(alloc, ctx->size);
        A->j     = calloc(alloc, sizeof(long));

        if (!(A->x) || !(A->j) || !(A->p) || !(A->lenr))
        {
            printf("ERROR (mat_csr_realloc).\n\n");
            abort();
        }

        for (k = 0; k < alloc; k++)
            ctx->init(A->x + k * (ctx->size));
    }
}

