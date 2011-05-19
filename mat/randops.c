/******************************************************************************

    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include "mat.h"

#include "flint.h"
#include "ulong_extras.h"

void mat_randops(mat_t mat, flint_rand_t state, long count, 
                       const ctx_t ctx)
{
    long c, i, j, k, m, n;

    m = mat->m;
    n = mat->n;

    for (c = 0; c < count; c++)
    {
        if (n_randint(state, 2))
        {
            if ((i = n_randint(state, m)) == (j = n_randint(state, m)))
                continue;
            if (n_randint(state, 2))
                for (k = 0; k < n; k++)
                    ctx->add(mat_entry(mat, j, k, ctx), 
                             mat_entry(mat, j, k, ctx), 
                             mat_entry(mat, i, k, ctx));
            else
                for (k = 0; k < n; k++)
                    ctx->sub(mat_entry(mat, j, k, ctx), 
                             mat_entry(mat, j, k, ctx), 
                             mat_entry(mat, i, k, ctx));
        }
        else
        {
            if ((i = n_randint(state, n)) == (j = n_randint(state, n)))
                continue;
            if (n_randint(state, 2))
                for (k = 0; k < m; k++)
                    ctx->add(mat_entry(mat, k, j, ctx), 
                             mat_entry(mat, k, j, ctx), 
                             mat_entry(mat, k, i, ctx));
            else
                for (k = 0; k < m; k++)
                    ctx->sub(mat_entry(mat, k, j, ctx), 
                             mat_entry(mat, k, j, ctx), 
                             mat_entry(mat, k, i, ctx));
        }
    }
}

