#include <assert.h>

#include "mat.h"

#include "mat_csr.h"

void mat_csr_set_mat(mat_csr_t A, const mat_t mat, const mat_ctx_t ctx)
{
    long i, j, *lenr, sum;

    assert(A->m == mat->m && A->n == mat->n);

    mat_csr_zero(A, ctx);

    lenr = calloc(mat->m, sizeof(long));

    if (!lenr)
    {
        printf("ERROR (mat_csr_set_mat).\n\n");
        abort();
    }

    for (i = 0; i < mat->m; i++)
        for (j = 0; j < mat->n; j++)
            if (!ctx->is_zero(ctx, mat_entry(mat, i, j, ctx)))
                lenr[i] ++;

    sum = 0;
    for (i = 0; i < mat->m; i++)
        sum += lenr[i];

    mat_csr_fit_length(A, sum, ctx);

    for (i = 0; i < mat->m; i++)
    {
        A->p[i] = (i > 0) ? A->p[i - 1] + A->lenr[i - 1] : 0;

        for (j = 0; j < mat->n; j++)
            if (!ctx->is_zero(ctx, mat_entry(mat, i, j, ctx)))
            {
                long q = A->p[i] + A->lenr[i];

                ctx->set(ctx, A->x + q * ctx->size, mat_entry(mat, i, j, ctx));
                A->j[q] = j;
                A->lenr[i] ++;
            }
    }

    free(lenr);
}

