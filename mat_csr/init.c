#include "mat_csr.h"

void mat_csr_init(mat_csr_t A, long m, long n, const ctx_t ctx)
{
    A->m = m;
    A->n = n;

    A->alloc = 0;
    A->x     = NULL;
    A->j     = NULL;

    A->p     = calloc(m, sizeof(long));
    A->lenr  = calloc(m, sizeof(long));

    if (!(A->p) || !(A->lenr))
    {
        printf("ERROR(mat_csr_init).\n\n");
        abort();
    }
}

void mat_csr_init2(mat_csr_t A, long m, long n, long alloc, const ctx_t ctx)
{
    A->m = m;
    A->n = n;

    A->p     = calloc(m, sizeof(long));
    A->lenr  = calloc(m, sizeof(long));

    if (!(A->p) || !(A->lenr))
    {
        printf("ERROR(mat_csr_init2).\n\n");
        abort();
    }

    if (alloc)
    {
        long k;

        A->alloc = alloc;
        A->x     = malloc(alloc * ctx->size);
        A->j     = calloc(alloc, sizeof(long));

        if (!(A->x) || !(A->j))
        {
            printf("ERROR(mat_csr_init2).\n\n");
            abort();
        }

        for (k = 0; k < alloc; k++)
            ctx->init(ctx, A->x + k * (ctx->size));
    }
    else
    {
        A->alloc = 0;
        A->x     = NULL;
        A->j     = NULL;
    }
}
