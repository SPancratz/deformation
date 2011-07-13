#include "gmconnection.h"
#include "diagfrob.h"

void diagfrob(mat_t F, const fmpz *a, long n, long d, const ctx_t ctx)
{
    mon_t *B;
    long *iB, lenB, l, u;

    /*
        Whenever $p > n - 1$ the monomial basis for $H_{rig}^n(U)$ is 
        crystalline and Frobenius is integral.  Thus, we can set the 
        bound on the valuations of the denominators to 0.
     */
    long bound = 0;

    long i, j;

    gmc_basis_sets(&B, &iB, &lenB, &l, &u, n, d);
    
    mat_clear(F, ctx);
    mat_init(F, lenB, lenB, ctx);
    
    for (i = 0; i < lenB; i++)
        for (j = 0; j < lenB; j++)
        {
            diagfrob_entry((void *) mat_entry(F, i, j, ctx), a, n, d, 
                           B[i], B[j], bound, ctx->pctx);
        }

    free(B);
    free(iB);
}

