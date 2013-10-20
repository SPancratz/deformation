/******************************************************************************

    Copyright (C) 2013 Sebastian Pancratz
 
******************************************************************************/

#include <stdlib.h>
#include <assert.h>

#include "gmconnection.h"

/*
    Performs a binary search for the monomial $x$ in the array $a$ 
    between $l$ and $u$, both inclusive.

    Assumes that the array is sorted in ascending order according 
    to the inverse lexicographical order, which coincides with the 
    word order of the underlying data.
 */
static long gmc_bsearch(const mon_t *a, long l, long u, const mon_t x)
{
    long m;

    while (l <= u)
    {
        m = l + ((u - l) >> 1);
        if (x < a[m])
            u = m - 1;
        else if (x > a[m])
            l = m + 1;
        else
            return m;
    }
    return -1;
}

/* Expects uninitialised matrix M, and p of length n + 1. */

void gmc_init_auxmatrix(mat_csr_t M, 
                        mon_t **R, mon_t **C, long *p, 
                        const mpoly_t P, long k, 
                        const ctx_t ctx)
{
    const long d = mpoly_degree(P, -1, ctx);
    const long n = P->n;

    long i, j, var;
    mpoly_t *DP;
    mon_t *rows, *cols;
    long nrows, ncols;

    char *mem;
    long c;
    const long u = 2 * sizeof(long) + ctx->size;
    
    assert(n - 1 > 0);
    assert(d > 0);
    assert((1 <= k) && (k <= n));
    assert((k - 1) * d >= n - 1);
    
    DP = malloc(n * sizeof(mpoly_t));
    for (var = 0; var < n; var++)
    {
        mpoly_init(DP[var], n, ctx);
        mpoly_derivative(DP[var], P, var, ctx);
    }
    
    /* Step 1.  Compute the row and column index sets */
    {
        mon_t *rows0, *cols0;
        long nrows0, ncols0;

        /* Supersets of the rows */
        rows0 = mon_generate_by_degree_invlex(&nrows0, n, k * d - n);
        cols0 = mon_generate_by_degree_invlex(&ncols0, n, (k - 1) * d - (n - 1));
        
        /* Setting up the rows */
        rows = malloc(nrows0 * sizeof(mon_t));
        nrows = 0;
        for (i = 0; i < nrows0; i++)
        {
            for (var = 0; var < n; var++)
            {
                if (mon_get_exp(rows0[i], var) >= d - 1)
                {
                    mon_init(rows[nrows]);
                    mon_set(rows[nrows], rows0[i]);
                    nrows++;
                    break;
                }
            }
            mon_clear(rows0[i]);
        }
        
        /* Setting up the columns */
        cols  = malloc(nrows * sizeof(mon_t));
        ncols = 0;
        p[0]  = 0;
        for (j = 0; j < ncols0; j++)
        {
            mon_init(cols[ncols]);
            mon_set(cols[ncols++], cols0[j]);
            mon_clear(cols0[j]);
        }
        for (var = 1; var < n; var++)
        {
            p[var] = ncols;
            for (j = p[var - 1]; j < p[var]; j++)
            {
                if (mon_get_exp(cols[j], var - 1) < d - 1)
                {
                    mon_init(cols[ncols]);
                    mon_set(cols[ncols++], cols[j]);
                }
            }
        }
        p[n] = ncols;

        free(rows0);
        free(cols0);
    }

    /* Step 2.  Set-up the auxiliary matrix */
    {
        mem = malloc(u * nrows * nrows);
        c = 0;

        for (j = 0; j < n; j++)
        {
            long q;

            for (q = p[j]; q < p[j + 1]; q++)
            {
                mpoly_iter_t iter;
                mpoly_term mt;

                /* Look at the column (j, cols[q]) */
                mpoly_iter_init(iter, DP[j]);
                while ((mt = mpoly_iter_next(iter)))
                {
                    mon_t f;
                    int f_is_in_rows;

                    mon_mul(f, mt->key, cols[q]);
                    f_is_in_rows = 0;
                    while (f)
                    {
                        if ((f & MON_BITMASK_BLOCK) >= d - 1)
                        {
                            f_is_in_rows = 1;
                            break;
                        }
                        f >>= MON_BITS_PER_EXP;
                    }
                    if (f_is_in_rows)
                    {
                        mon_mul(f, mt->key, cols[q]);
                        i = gmc_bsearch(rows, 0, nrows - 1, f);

                        /* Add the entry mt->value in position (i, q) */
                        *(long *) (mem + c * u) = i;
                        *(long *) (mem + c * u + sizeof(long)) = q;
                        ctx->init(ctx, mem + c * u + 2 * sizeof(long));
                        ctx->set(ctx, mem + c * u + 2 * sizeof(long), mt->val);
                        c++;
                    }
                }
                mpoly_iter_clear(iter);
            }
        }

        mat_csr_init2(M, nrows, ncols, c, ctx);
        mat_csr_set_array3(M, mem, c, 0, ctx);
        free(mem);
    }

    /* Copy output */
    *R  = rows;
    *C  = cols;

    /* Cleaning up */
    for (i = 0; i < n; i++)
        mpoly_clear(DP[i], ctx);
    free(DP);
}

