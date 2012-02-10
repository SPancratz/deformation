#include <stdlib.h>
#include <mpir.h>

#include "gmconnection.h"
#include "diagfrob.h"
#include "gmde.h"

#include "flint.h"
#include "fmpq_mat.h"

#include "deformation.h"

/*
    Computes $F(0)$ to precision \code{pF0}, $C(t)$ to 
    precisions \code{pC} and \code{tC}, and $C^{-1}(t)$ 
    to precisions \code{pCinv} and \code{tCinv}, finally 
    combining these to form the product 
    \begin{equation*}
    F(t) = C(t) F(0) C^{-1}(t^p), 
    \end{equation*}
    reducing the result modulo the precision found in 
    the context for $F$.

    Expects \code{ctxF} to be a context for objects in 
    the \code{padic_poly} module.
 */
void frob_with_precisions_fmpq(mat_t F, const ctx_t ctxF, 
                               const mpoly_t P, const ctx_t ctxFracQt, 
                               long NWork, long K)
{
    const long n = P->n - 1;
    const long d = mpoly_degree(P, -1, ctxFracQt);
    const long b = gmc_basis_size(n, mpoly_degree(P, -1, ctxFracQt));
    const long p = fmpz_get_si(ctxF->pctx->p);

    mon_t *bR, *bC;

    mat_t M, F0, C, Cinv;

    padic_ctx_t pctx;
    ctx_t ctxZp;
    ctx_t ctxZpt;

    long i, j, k;

    /* Init */
    padic_ctx_init(pctx, ctxF->pctx->p, NWork, PADIC_VAL_UNIT);
    ctx_init_padic(ctxZp, pctx);
    ctx_init_padic_poly(ctxZpt, pctx);

    mat_init(M, b, b, ctxFracQt);
    mat_init(F0, b, b, ctxZp);
    mat_init(C, b, b, ctxZpt);
    mat_init(Cinv, b, b, ctxZpt);

    /* Compute Gauss--Manin connection */
    gmc_compute(M, &bR, &bC, P, ctxFracQt);

    /* Find F(0) */
    {
        fmpz * a = _fmpz_vec_init(n + 1);

        mpoly_diagonal_fibre(a, P, ctxFracQt);

        diagfrob(F0, a, n, d, ctxZp);

        printf("Compute the diagonal fibre.\n");
        printf("a = {"), _fmpz_vec_print(a, n + 1), printf("}\n");

        _fmpz_vec_clear(a, n + 1);
    }

    /* Solve for C and Cinv */
    {
        ctx_t ctxQt;

        mat_t Mt;
        mat_t B, Binv;

        fmpq_mat_struct *A, *Ainv;

        A    = malloc(K * sizeof(fmpq_mat_struct));
        Ainv = malloc(((K + p - 1) / p) * sizeof(fmpq_mat_struct));
        for(i = 0; i < K; i++)
            fmpq_mat_init(A + i, b, b);
        for(i = 0; i < (K + p - 1) / p; i++)
            fmpq_mat_init(Ainv + i, b, b);

        ctx_init_fmpq_poly(ctxQt);

        mat_init(B, b, b, ctxQt);
        mat_init(Binv, b, b, ctxQt);

        /* Mt = - M^t */
        mat_init(Mt, b, b, ctxFracQt);
        mat_transpose(Mt, M, ctxFracQt);
        mat_neg(Mt, Mt, ctxFracQt);

        gmde_solve_fmpq(A, K, M, ctxFracQt);
        gmde_solve_fmpq(Ainv, (K + p - 1) / p, Mt, ctxFracQt);
        gmde_convert_soln_fmpq(B, ctxQt, A, K);
        gmde_convert_soln_fmpq(Binv, ctxQt, Ainv, (K + p - 1) / p);
        mat_transpose(Binv, Binv, ctxQt);

        /* Push B and Binv into Zp[t] */
        for (i = 0; i < b; i++)
            for (j = 0; j < b; j++)
            {
                padic_poly_set_fmpq_poly(
                    (padic_poly_struct *) mat_entry(C, i, j, ctxZpt), 
                    (fmpq_poly_struct *) mat_entry(B, i, j, ctxQt), 
                    pctx);
                padic_poly_set_fmpq_poly(
                    (padic_poly_struct *) mat_entry(Cinv, i, j, ctxZpt), 
                    (fmpq_poly_struct *) mat_entry(Binv, i, j, ctxQt), 
                    pctx);
            }

        for(i = 0; i < K; i++)
            fmpq_mat_clear(A + i);
        for(i = 0; i < (K + p - 1) / p; i++)
            fmpq_mat_clear(Ainv + i);
        free(A);
        free(Ainv);
        mat_clear(B, ctxQt);
        mat_clear(Binv, ctxQt);
        mat_clear(Mt, ctxFracQt);

        ctx_clear(ctxQt);
    }

    /* Replace t by t^p in C^{-1} */
    for (i = 0; i < b; i++)
        for (j = 0; j < b; j++)
        {    
            padic_poly_compose_pow(
                (padic_poly_struct *) mat_entry(Cinv, i, j, ctxZpt), 
                (padic_poly_struct *) mat_entry(Cinv, i, j, ctxZpt), p, pctx);
        }

    /* Form the product C(t) F(0) C(t^p)^{-1} */
    {
        mat_t RHS;

        mat_init(RHS, b, b, ctxZpt);

        for (i = 0; i < b; i++)
        {
            /* Find the unique k s.t. F0(i,k) is non-zero */
            for (k = 0; k < b; k++) 
                if (!_padic_is_zero((padic_struct *) mat_entry(F0, i, k, ctxZp)))
                    break;

            if (k == b)
            {
                printf("ERROR (frob_with_precisions). Bad F0.\n\n");
                abort();
            }

            for (j = 0; j < b; j++)
            {
                padic_poly_scalar_mul_padic(
                    (padic_poly_struct *) mat_entry(RHS, i, j, ctxZpt), 
                    (padic_poly_struct *) mat_entry(Cinv, k, j, ctxZpt), 
                    (padic_struct *) mat_entry(F0, i, k, ctxZp), pctx);
            }
        }

        mat_mul(F, C, RHS, ctxZpt);

        for (i = 0; i < b; i++)
            for (j = 0; j < b; j++)
            {
                padic_poly_truncate(
                    (padic_poly_struct *) mat_entry(F, i, j, ctxZpt), 
                    K, pctx->p);
            }

        mat_clear(RHS, ctxZpt);
    }

    /* Multiply F by r^{p prec(F) + 10%} */

    {
        fmpz_poly_t r, t;
        padic_poly_t s;

        fmpz_poly_init(r);
        fmpz_poly_init(t);
        padic_poly_init(s);

        fmpz_poly_set_ui(r, 1);
        for (i = 0; i < M->m; i++)
            for (j = 0; j < M->n; j++)
            {
                fmpz_poly_lcm(t, r, fmpz_poly_q_denref(
                    (fmpz_poly_q_struct *) mat_entry(M, i, j, ctxFracQt)));
                fmpz_poly_swap(r, t);
            }

        printf("r = ");
        fmpz_poly_print_pretty(r, "t");
        printf("\n");

        padic_poly_set_fmpz_poly(s, r, pctx);
        padic_poly_pow(s, s, ((p * ctxF->pctx->N) * 11) / 10, pctx);

        for (i = 0; i < b; i++)
            for (j = 0; j < b; j++)
            {
                padic_poly_mul(
                    (padic_poly_struct *) mat_entry(F, i, j, ctxF), 
                    (padic_poly_struct *) mat_entry(F, i, j, ctxF), 
                    s, ctxF->pctx);

                padic_poly_truncate(
                    (padic_poly_struct *) mat_entry(F, i, j, ctxF), K, pctx->p);
            }

        fmpz_poly_clear(r);
        fmpz_poly_clear(t);
        padic_poly_clear(s);
    }

    /* Clean up */
    free(bR);
    free(bC);

    mat_clear(M, ctxFracQt);
    mat_clear(F0, ctxZp);
    mat_clear(C, ctxZpt);
    mat_clear(Cinv, ctxZpt);

    padic_ctx_clear(pctx);
    ctx_clear(ctxZpt);
    ctx_clear(ctxZp);
}

