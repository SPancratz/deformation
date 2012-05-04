#include <stdlib.h>
#include <mpir.h>

#include "gmconnection.h"
#include "diagfrob.h"
#include "gmde.h"

#include "flint.h"
#include "fmpq_mat.h"

#include "deformation.h"

#define DEBUG 1

/*
    Step 1.

    Computes the Gauss--Manin connection $M$ over $\mathbf{Q}(t)$ 
    with denominator $r$ over $\mathbf{Z}$.

    Step 2.

    Computes $F(0)$ to precision $N_3$ for $p^{-1} F_p$ 
    on the diagonal fibre.

    Step 3.

    Compute solution $C(t)$ over $\mathbf{Q}_p[[t]]$ modulo 
    $p^{N_2}$ and $t^{K}$.  Also compute $C^{-1}(t^p)$ to 
    the same precision.

    Step 4.

    Compute matrix $F_t$ for $p^{-1} F_p$ on the generic 
    fibre via 
    \begin{equation*}
    F(t) = C(t) F(0) C^{-1}(t^p)
    \end{equation*}
    modulo $p^{N_1}$ and $t^{K}$.

    Step 5.

    Compute matrix $G(t) = r(t)^m F(t)$ over $\mathbf{Q}_p[[t]]$ 
    modulo $p^{N_1}$ and $t^{K}$.

    Note that $m$ should be about $1.10 \times p N_1$.

    Step 6.

    Evaluate this at $\hat{z}_1$, the Teichmuller lift of $z_1$ 
    by computing $F(1) = r(\hat{z}_1)^{-m} G(\hat{z}_1)$ all 
    modulo $p^{N_1}$.

    Assumptions:

        - F is a matrix over $\mathbf{Q}_p[t]$
        - p is a word-sized prime
 */

void frob(const mpoly_t P, const ctx_t ctxFracQt, const fmpz_t p)
{
    const long n  = P->n - 1;
    const long d  = mpoly_degree(P, -1, ctxFracQt);
    const long b  = gmc_basis_size(n, d);

    long N0, N1, N2, N3, K, m;

    /* Diagonal fibre */
    padic_ctx_t pctx_F0;
    ctx_t ctxZp_F0;
    mat_t F0;

    /* Gauss--Manin Connection */
    mat_t M;
    mon_t *bR, *bC;
    fmpz_poly_t r;

    /* Local solution */
    padic_ctx_t pctx_C;
    ctx_t ctxZpt_C;
    mat_t C, Cinv;

    /* Frobenius */
    padic_ctx_t pctx_F;
    ctx_t ctxZpt_F;
    mat_t F;

    /* Step 1 {M, r} *********************************************************/

    mat_init(M, b, b, ctxFracQt);
    fmpz_poly_init(r);

    gmc_compute(M, &bR, &bC, P, ctxFracQt);

    {
        long i, j;
        fmpz_poly_t t;

        fmpz_poly_init(t);
        fmpz_poly_set_ui(r, 1);
        for (i = 0; i < M->m; i++)
            for (j = 0; j < M->n; j++)
            {
                fmpz_poly_lcm(t, r, fmpz_poly_q_denref(
                    (fmpz_poly_q_struct *) mat_entry(M, i, j, ctxFracQt)));
                fmpz_poly_swap(r, t);
            }
        fmpz_poly_clear(t);
    }

#if(DEBUG == 1)
printf("Gauss--Manin connection M:\n");
mat_print(M, ctxFracQt);
printf("\n");
printf("r = ");
fmpz_poly_print_pretty(r, "t");
printf("\n");
#endif

    /* Precisions ************************************************************/

    N0 = deformation_prec_zeta_function(p, 1, n, d);
    N1 = deformation_prec_frob(p, n, N0);
    m  = deformation_prec_pole_order(p, N1);
    K  = deformation_prec_tadic(p, 4, N1);
    N2 = deformation_prec_local_solution(p, n, N1, K);
    N3 = deformation_prec_diagfrob(p, n, N1, K);

    printf("Precisions:\n");
    printf("N0 = %ld\n", N0);
    printf("N1 = %ld\n", N1);
    printf("N2 = %ld\n", N2);
    printf("N3 = %ld\n", N3);
    printf("m  = %ld\n", m);
    printf("K  = %ld\n", K);

    /* Initialisation ********************************************************/

    padic_ctx_init(pctx_F0, p, N3, PADIC_VAL_UNIT);
    ctx_init_padic(ctxZp_F0, pctx_F0);
    mat_init(F0, b, b, ctxZp_F0);

    padic_ctx_init(pctx_C, p, N2, PADIC_VAL_UNIT);
    ctx_init_padic_poly(ctxZpt_C, pctx_C);

    mat_init(C, b, b, ctxZpt_C);
    mat_init(Cinv, b, b, ctxZpt_C);

    padic_ctx_init(pctx_F, p, N1, PADIC_VAL_UNIT);
    ctx_init_padic_poly(ctxZpt_F, pctx_F);
    mat_init(F, b, b, ctxZpt_F);

    /* Step 2 {F0} ***********************************************************/

    {
        fmpz *t = _fmpz_vec_init(n + 1);

        mpoly_diagonal_fibre(t, P, ctxFracQt);
        diagfrob(F0, t, n, d, ctxZp_F0, 0);

#if(DEBUG == 1)
printf("Diagonal fibre F(0):\n");
printf("{"), _fmpz_vec_print(t, n + 1), printf("}\n");
mat_print(F0, ctxZp_F0);
printf("\n");
#endif

        _fmpz_vec_clear(t, n + 1);
    }

    /* Step 3 {C, Cinv} ******************************************************/
    /*
        Compute C as a matrix over Z_p[[t]].  A is the same but as a series 
        of matrices over Z_p.  Mt is the matrix -M^t, and C^{-1} is the 
        local solution of the differential equation replacing M by Mt.
     */

    {
        long i;
        mat_t Mt;
        padic_mat_struct *A, *Ainv;

        A    = malloc(K * sizeof(padic_mat_struct));
        Ainv = malloc(((K + (*p) - 1) / (*p)) * sizeof(padic_mat_struct));
        for(i = 0; i < K; i++)
            padic_mat_init(A + i, b, b);
        for(i = 0; i < (K + (*p) - 1) / (*p); i++)
            padic_mat_init(Ainv + i, b, b);

        mat_init(Mt, b, b, ctxFracQt);
        mat_transpose(Mt, M, ctxFracQt);
        mat_neg(Mt, Mt, ctxFracQt);
        gmde_solve(A, K, pctx_C, M, ctxFracQt);
        gmde_solve(Ainv, (K + (*p) - 1) / (*p), pctx_C, Mt, ctxFracQt);
        gmde_convert_soln(C, ctxZpt_C, A, K);
        gmde_convert_soln(Cinv, ctxZpt_C, Ainv, (K + (*p) - 1) / (*p));
        mat_transpose(Cinv, Cinv, ctxZpt_C);

        for(i = 0; i < K; i++)
            padic_mat_clear(A + i);
        for(i = 0; i < (K + (*p) - 1) / (*p); i++)
            padic_mat_clear(Ainv + i);
        free(A);
        free(Ainv);
        mat_clear(Mt, ctxFracQt);
    }


    {
        long i, j;

        /* Replace t by t^p in C^{-1} */
        for (i = 0; i < b; i++)
            for (j = 0; j < b; j++)
            {    
                padic_poly_compose_pow(
                    (padic_poly_struct *) mat_entry(Cinv, i, j, ctxZpt_C), 
                    (padic_poly_struct *) mat_entry(Cinv, i, j, ctxZpt_C), fmpz_get_si(p), pctx_C);
            }
    }

    /* Step 4 {F} ************************************************************/

    /* Form the product C(t) F(0) C(t^p)^{-1} */
    {
        long i, j, k;
        mat_t RHS;

        mat_init(RHS, b, b, ctxZpt_C);

        for (i = 0; i < b; i++)
        {
            /* Find the unique k s.t. F0(i,k) is non-zero */
            for (k = 0; k < b; k++) 
                if (!_padic_is_zero((padic_struct *) mat_entry(F0, i, k, ctxZp_F0)))
                    break;

            if (k == b)
            {
                printf("Exception (frob). Bad F0.\n\n");
                abort();
            }

            for (j = 0; j < b; j++)
            {
                padic_poly_scalar_mul_padic(
                    (padic_poly_struct *) mat_entry(RHS, i, j, ctxZpt_C), 
                    (padic_poly_struct *) mat_entry(Cinv, k, j, ctxZpt_C), 
                    (padic_struct *) mat_entry(F0, i, k, ctxZp_F0), pctx_C);
            }
        }

        mat_mul(F, C, RHS, ctxZpt_C);

        for (i = 0; i < b; i++)
            for (j = 0; j < b; j++)
            {
                padic_poly_truncate(
                    (padic_poly_struct *) mat_entry(F, i, j, ctxZpt_F), K, p);
                padic_poly_reduce(
                    (padic_poly_struct *) mat_entry(F, i, j, ctxZpt_F), pctx_F);
            }

        mat_clear(RHS, ctxZpt_C);
    }

    /* Step 5 {G = r(t)^m F(t) ***********************************************/

    {
        long i, j;
        padic_poly_t t;

        padic_poly_init(t);

        padic_poly_set_fmpz_poly(t, r, pctx_F);
        padic_poly_pow(t, t, m, pctx_F);

        for (i = 0; i < b; i++)
            for (j = 0; j < b; j++)
            {
                padic_poly_mul(
                    (padic_poly_struct *) mat_entry(F, i, j, ctxZpt_F), 
                    (padic_poly_struct *) mat_entry(F, i, j, ctxZpt_F), 
                    t, pctx_F);

                padic_poly_truncate(
                    (padic_poly_struct *) mat_entry(F, i, j, ctxZpt_F), K, p);
            }

        padic_poly_clear(t);
    }

    /* Clean up **************************************************************/

    mat_clear(F0, ctxZp_F0);
    ctx_clear(ctxZp_F0);
    padic_ctx_clear(pctx_F0);

    mat_clear(M, ctxFracQt);
    free(bR);
    free(bC);
    fmpz_poly_clear(r);

    mat_clear(C, ctxZpt_C);
    mat_clear(Cinv, ctxZpt_C);
    ctx_clear(ctxZpt_C);
    padic_ctx_clear(pctx_C);

    mat_clear(F, ctxZpt_F);
    ctx_clear(ctxZpt_F);
    padic_ctx_clear(pctx_F);
}

