#include "gmconnection.h"

#define DEBUG  1

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

void gmc_compute(mat_t M, mon_t *rows, mon_t *cols, 
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
        but they are NULL when k < l
     */
    mat_csr_t *aux;
    mon_t **aux_rows, **aux_cols;
    long **aux_p;
    mat_csr_solve_t *aux_s;
    
    #if DEBUG > 0
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
    
    #if DEBUG > 0
        printf("Derivatives:\n");
        for (i = 0; i < n; i++)
        {
            printf("  dPdX %ld = ", i), mpoly_print(dP[i], ctx), printf("\n");
        }
        printf("  dPdt = "), mpoly_print(dPdt, ctx), printf("\n");
    #endif
    
    /* Construct the index set B */
    gmc_basis_sets(&B, &iB, &lenB, &l, &u, n, d);
    
    #if DEBUG > 0
        printf("Basis sets:\n");
        printf("  [");
        for (k = 1; k < l; k++)
            printf(" |");
        for ( ; k <= u; k++)
        {
            for (i = iB[k]; i < iB[k + 1]; i++)
            {
                printf(i == iB[k] ? " " : ", ");
                mon_print(B[i], n);
                printf("\n");
            }
            printf(" |");
        }
        for ( ; k < n - 1; k++)
            printf(" |");
        printf(" ]\n");
        fflush(stdout);
    #endif
    
    /* Connection matrix init */
    mat_clear(M, ctx);
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

    /*
        While (k-1)*d < n OR kd < n+1, which is equivalent to saying 
        while (k-1)*d < n, since d > 1
     */
    for (k = 0; (k - 1) * d < (n - 1); k++) ;
    
    for ( ; k <= u + 1; k++)
    {
        
        #if (DEBUG > 0)
            printf("Computing auxiliary matrix for k = %ld.\n", k);
        #endif
        
        aux_p[k] = malloc((n + 1) * sizeof(long));
        gmc_init_auxmatrix(aux[k], aux_rows + k, aux_cols + k, aux_p[k], 
                           P, k, ctx);
        mat_csr_solve_init(aux_s[k], aux[k], ctx);
    }
    
    /* Construct the Gauss--Manin connection matrix */
    
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
            printf("Computing column j = %ld.\n", j);
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
    
    /* Clean up */
    for (i = 0; i < n; i++)
        mpoly_clear(dP[i], ctx);
    free(dP);
    mpoly_clear(dPdt, ctx);
    
    free(iB);
    
    k = ((n - 1) + (d - 1)) / d + 1;  /* k = ceil(n / d) + 1 */

    for ( ; k <= u + 1; k++)
    {
        #if (DEBUG > 0)
            printf("Computing auxiliary matrix for k = %ld.\n", k);
        #endif
        
        free(aux_p[k]);
        mat_csr_clear(aux[k], ctx);
        mat_csr_solve_clear(aux_s[k], ctx);
    }
}

