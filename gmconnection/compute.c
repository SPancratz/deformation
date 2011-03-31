#include "gmconnection.h"

#define DEBUG  0

/*
    Sets the polynomial \code{rop} to the derivative of \code{op} 
    with respect to the indeterminate $t$ of the base field.
 */
static void mpoly_tderivative(mpoly_t rop, const mpoly_t op, const mat_ctx_t ctx)
{
    mpoly_t temp;
    
    mpoly_iter_t iter;
    mpoly_term t;
    
    if (mpoly_is_zero(op, ctx))
    {
        mpoly_zero(rop, ctx);
        rop->n = op->n;
        return;
    }
    
    mpoly_init(temp, op->n, ctx);
    
    mpoly_iter_init(iter, op);
    while ((t = mpoly_iter_next(iter)))
    {
        char *c;

        c = malloc(ctx->size);
        ctx->init(c);
        ctx->derivative(c, t->val);

        if (ctx->is_zero(c))
        {
            ctx->clear(c);
            free(c);
        }
        else
        {
            mon_t m2;
            void *c2;

            /* Note that no entry with this key is present in temp yet */
            RBTREE_INSERT(mpoly, &m2, &c2, temp->dict, t->key, c, &mon_cmp);
        }
    }
    mpoly_iter_clear(iter);
    
    mpoly_swap(rop, temp, ctx);
    mpoly_clear(temp, ctx);
}

void gmc_compute(mat_t M, mon_t **rows, mon_t **cols, 
                 const mpoly_t P, const mat_ctx_t ctx)
{
    long i, j, k, colk, rowk;
    mpoly_t *dP, dPdt;

    mon_t *B;
    long *iB, l, u, lenB;

    const long n = P->n;
    const long d = mpoly_degree(P, -1, ctx);
    
    /*
        Arrays for the auxiliary matrices;  the arrays are of length u + 1 
        but they are unused when k < l
     */
    mat_csr_t *aux;
    mon_t **aux_rows, **aux_cols;
    long **aux_p;
    mat_csr_solve_t *aux_s;
    
    #if (DEBUG > 0)
    printf("Input:\n");
    printf("  P = "), mpoly_print(P, ctx), printf("\n");
    #endif
    
    /*
        Compute all partial derivatives of P, and the derivative of P with 
        respect to t, the variable of the coefficient field
     */
    dP = malloc(n * sizeof(mpoly_t));
    for (i = 0; i < n; i++)
        mpoly_init(dP[i], n, ctx);
    gmc_derivatives(dP, P, ctx);
    mpoly_init(dPdt, n, ctx);
    mpoly_tderivative(dPdt, P, ctx);
    
    #if (DEBUG > 0)
    printf("Derivatives:\n");
    for (i = 0; i < n; i++)
        printf("  dPdX %ld = ", i), mpoly_print(dP[i], ctx), printf("\n");
    printf("  dPdt = "), mpoly_print(dPdt, ctx), printf("\n");
    #endif
    
    /* Construct the index set B */
    gmc_basis_sets(&B, &iB, &lenB, &l, &u, n, d);
    
    #if (DEBUG > 0)
    printf("Basis sets:\n  ");
    gmc_basis_print(B, iB, lenB, n, d), printf("\n");
    fflush(stdout);
    #endif
    
    /* Connection matrix init */
    mat_init(M, lenB, lenB, ctx);

    /*
        Construct the auxiliary matrices

        Note that a priori we might need them for k = 1, ..., n, n+1. 
        The rows and columns sets are only non-empty if (k-1)*d >= n.

        To ease indexing, we include the empty set in the case k = 0.
     */
    aux      = malloc((n + 1) * sizeof(mat_csr_t));
    aux_rows = malloc((n + 1) * sizeof(mon_t *));
    aux_cols = malloc((n + 1) * sizeof(mon_t *));
    aux_p    = malloc((n + 1) * sizeof(long *));
    aux_s    = malloc((n + 1) * sizeof(mat_csr_solve_t));

    #if (DEBUG > 0)
    printf("Computing auxiliary matrices..\n");
    #endif

    for (k = ((n - 1) + (d - 1)) / d + 1; k <= u + 1; k++)
    {
        
        #if (DEBUG > 0)
        printf("  k = %ld\n", k), fflush(stdout);
        #endif
        
        aux_p[k] = malloc((n + 1) * sizeof(long));

        gmc_init_auxmatrix(aux[k], aux_rows + k, aux_cols + k, aux_p[k], 
                           P, k, ctx);
        mat_csr_sort_rows(aux[k], ctx);
        mat_csr_solve_init(aux_s[k], aux[k], ctx);
    }
    
    /* Construct the Gauss--Manin connection matrix */

    #if (DEBUG > 0)
    printf("Computing columns..\n");
    #endif

    for (colk = l; colk <= u; colk++)
    {
        for (j = iB[colk]; j < iB[colk + 1]; j++)
        {
            mpoly_t Q, *R;

            mpoly_init(Q, n, ctx);
            R = malloc((n + 1) * sizeof(mpoly_t));
            for (i = 0; i <= n; i++)
                mpoly_init(R[i], n, ctx);

            #if (DEBUG > 0)
            printf("  j = %ld\n", j);
            #endif
            
            /* Set Q to -colk B[j] dPdt, and then reduce */
            mpoly_mul_mon(Q, dPdt, B[j], ctx);
            mpoly_scalar_mul_si(Q, Q, -colk, ctx);
            
            gmc_reduce(R, Q, colk + 1, d, dP, 
                       aux_s, aux_rows, aux_cols, aux_p, l, u, ctx);

            /* Extract the column vector */
            
            /* Iterate over all rows.. */
            
            for (rowk = l; rowk <= FLINT_MIN(u, colk + 1); rowk++)
            {
                for (i = iB[rowk]; i < iB[rowk + 1]; i++)
                {
                    mpoly_get_coeff(mat_entry(M, i, j, ctx), 
                                    R[rowk], B[i], ctx);
                }
            }

            mpoly_clear(Q, ctx);
            for (i = 0; i <= n; i++)
                mpoly_clear(R[i], ctx);
            free(R);
        }
    }

    /* Copy row and column index sets */

    *rows = malloc(lenB * sizeof(mon_t));
    *cols = malloc(lenB * sizeof(mon_t));

    for (i = 0; i < lenB; i++)
    {
        (*rows)[i] = B[i];
        (*cols)[i] = B[i];
    }

    /* Clean up */

    k = ((n - 1) + (d - 1)) / d + 1;  /* k = ceil(n / d) + 1 */

    for ( ; k <= u + 1; k++)
    {
        free(aux_rows[k]);
        free(aux_cols[k]);
        free(aux_p[k]);
        mat_csr_clear(aux[k], ctx);
        mat_csr_solve_clear(aux_s[k], ctx);
    }

    free(aux);
    free(aux_rows);
    free(aux_cols);
    free(aux_p);
    free(aux_s);

    free(B);
    free(iB);

    for (i = 0; i < n; i++)
        mpoly_clear(dP[i], ctx);
    free(dP);
    mpoly_clear(dPdt, ctx);
}

