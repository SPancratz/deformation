#include <stdlib.h>
#include <mpir.h>
#include <limits.h>

#include "gmconnection.h"
#include "diagfrob.h"
#include "gmde.h"

#include "flint.h"
#include "fmpq_mat.h"
#include "padic_mat.h"
#include "fmpz_poly_mat.h"
#include "fmpz_mod_poly.h"

#include "deformation.h"

#define DEBUG 1

static void mat_print_sage(const mat_t mat, const ctx_t ctx)
{
    long i, j;

    printf("[");
    for (i = 0; i < mat->m; i++)
        for (j = 0; j < mat->n; j++)
        {
            ctx->print(ctx, mat_entry(mat, i, j, ctx));
            if (!(i == mat->m - 1 && j == mat->n - 1))
                printf(", ");
        }
    printf("]");
}

static 
void _fmpz_poly_compose_pow(fmpz *rop, const fmpz *op, long len, long k)
{
    if (k == 1)
    {
        if (rop != op)
        {
            _fmpz_vec_set(rop, op, len);
        }
    }
    else if (len == 1)
    {
        fmpz_set(rop, op);
    }
    else
    {
        long i, j, h;

        for (i = len - 1, j = (len - 1) * k ; i >= 0; i--, j -= k)
        {
            fmpz_set(rop + j, op + i);
            if (i != 0)
                for (h = 1; h < k; h++)
                    fmpz_zero(rop + (j - h));
        }
    }
}

static 
void fmpz_poly_compose_pow(fmpz_poly_t rop, const fmpz_poly_t op, long k)
{
    const long len  = op->length;
    const long lenr = (len - 1) * k + 1;

    if (len == 0)
    {
        fmpz_poly_zero(rop);
    }
    else
    {
        fmpz_poly_fit_length(rop, lenr);
        _fmpz_poly_compose_pow(rop->coeffs, op->coeffs, len, k);
        _fmpz_poly_set_length(rop, lenr);
    }
}

static 
void fmpz_poly_mat_compose_pow(fmpz_poly_mat_t B, const fmpz_poly_mat_t A, long k)
{
    long i, j;

    if (!(A->r == B->r && A->c == B->c))
    {
        printf("Exception (fmpz_poly_mat_compose_pow).  Incompatible dimensions.\n");
        abort();
    }

    for (i = 0; i < B->r; i++)
        for (j = 0; j < B->c; j++)
            fmpz_poly_compose_pow(fmpz_poly_mat_entry(B, i, j), 
                                  fmpz_poly_mat_entry(A, i, j), k);
}

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

    Evaluate this at $\hat{t}_1$, the Teichmuller lift of $t_1$ 
    by computing $F(1) = r(\hat{t}_1)^{-m} G(\hat{t}_1)$,  all 
    modulo $p^{N_1}$.

    Assumptions:

        - F is a matrix over $\mathbf{Q}_p[t]$
        - p is a word-sized prime
 */

void frob(const mpoly_t P, const fmpz_t t1, const ctx_t ctxFracQt, const fmpz_t p)
{
    const long n  = P->n - 1;
    const long d  = mpoly_degree(P, -1, ctxFracQt);
    const long b  = gmc_basis_size(n, d);
    const long a = 1;

    long i, j, k;

    prec_struct prec = {0};

    /* Diagonal fibre */
    padic_ctx_t pctx_F0;
    padic_mat_t F0;

    /* Gauss--Manin Connection */
    mat_t M;
    mon_t *bR, *bC;
    fmpz_poly_t r;

    /* Local solution */
    padic_ctx_t pctx_C;
    fmpz_poly_mat_t C, Cinv;
    long vC, vCinv;

    /* Frobenius */
    fmpz_poly_mat_t F;
    long vF;

    padic_mat_t F1;

    fmpz_poly_t cp;

#if(DEBUG == 1)
printf("Input:\n");
printf("P  = "), mpoly_print(P, ctxFracQt), printf("\n");
printf("p  = "), fmpz_print(p), printf("\n");
printf("t1 = "), fmpz_print(t1), printf("\n");
printf("\n");
#endif

    /* Step 1 {M, r} *********************************************************/

    mat_init(M, b, b, ctxFracQt);
    fmpz_poly_init(r);

    gmc_compute(M, &bR, &bC, P, ctxFracQt);

    {
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
printf("\n\n");
printf("Denominator r:\n");
fmpz_poly_print_pretty(r, "t");
printf("\n\n");
#endif

    /* Precisions ************************************************************/

    deformation_precisions(&prec, p, 1, n, d, fmpz_poly_degree(r));

#if(DEBUG == 1)
printf("Precisions:\n");
printf("N0   = %ld\n", prec.N0);
printf("N1   = %ld\n", prec.N1);
printf("N2   = %ld\n", prec.N2);
printf("N3   = %ld\n", prec.N3);
printf("N3i  = %ld\n", prec.N3i);
printf("N3w  = %ld\n", prec.N3w);
printf("N3iw = %ld\n", prec.N3iw);
printf("N4   = %ld\n", prec.N4);
printf("m    = %ld\n", prec.m);
printf("K    = %ld\n", prec.K);
printf("r    = %ld\n", prec.r);
printf("s    = %ld\n", prec.s);
printf("\n");
#endif

    /* Initialisation ********************************************************/

    padic_ctx_init(pctx_F0, p, prec.N4, PADIC_VAL_UNIT);
    padic_mat_init(F0, b, b);

    padic_ctx_init(pctx_C, p, prec.N3, PADIC_VAL_UNIT);

    fmpz_poly_mat_init(C, b, b);
    fmpz_poly_mat_init(Cinv, b, b);

    fmpz_poly_mat_init(F, b, b);
    vF = 0;

    padic_mat_init(F1, b, b);
    fmpz_poly_init(cp);

    /* Step 2 {F0} ***********************************************************/

    {
        fmpz *t = _fmpz_vec_init(n + 1);

        mpoly_diagonal_fibre(t, P, ctxFracQt);
        diagfrob(F0, t, n, d, pctx_F0, 0);

        padic_mat_transpose(F0, F0);

#if(DEBUG == 1)
printf("Diagonal fibre:\n");
printf("P(0) = {"), _fmpz_vec_print(t, n + 1), printf("}\n");
printf("Matrix F(0):\n");
padic_mat_print_pretty(F0, pctx_F0);
printf("\n\n");
#endif

        _fmpz_vec_clear(t, n + 1);
    }

    /* Step 3 {C, Cinv} ******************************************************/
    /*
        Compute C as a matrix over Z_p[[t]].  A is the same but as a series 
        of matrices over Z_p.  Mt is the matrix -M^t, and Cinv is C^{-1}^t, 
        the local solution of the differential equation replacing M by Mt.
     */

    {
        long K = prec.K;
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
        gmde_solve(A, K, p, prec.N3, prec.N3w, M, ctxFracQt);
        gmde_solve(Ainv, (K + (*p) - 1) / (*p), p, prec.N3i, prec.N3iw, Mt, ctxFracQt);
        gmde_convert_soln(C, &vC, A, K, p);
        gmde_convert_soln(Cinv, &vCinv, Ainv, (K + (*p) - 1) / (*p), p);

#if(DEBUG == 1)
printf("Local solution C(t):\n");
fmpz_poly_mat_print(C, "t");
printf("\n");
printf("Matrix C(t)^{-1}:\n");
fmpz_poly_mat_print(Cinv, "t");
printf("\n");
#endif

        fmpz_poly_mat_transpose(Cinv, Cinv);
        fmpz_poly_mat_compose_pow(Cinv, Cinv, *p);

        for(i = 0; i < K; i++)
            padic_mat_clear(A + i);
        for(i = 0; i < (K + (*p) - 1) / (*p); i++)
            padic_mat_clear(Ainv + i);
        free(A);
        free(Ainv);
        mat_clear(Mt, ctxFracQt);
    }

    /* Step 4 {F(t) := C(t) F(0) C(t^p)^{-1}} ********************************/

    /*
        Computes the product C(t) F(0) C(t^p)^{-1} modulo (p^{N_2}, t^K). 
        This is done by first computing the unit part of the product 
        exactly over the integers modulo t^K.
     */

    {
        fmpz_t pN;
        fmpz_poly_mat_t T;

        fmpz_init(pN);
        fmpz_poly_mat_init(T, b, b);

        for (i = 0; i < b; i++)
        {
            /* Find the unique k s.t. F0(i,k) is non-zero */
            for (k = 0; k < b; k++) 
                if (!fmpz_is_zero(padic_mat_unit(F0, i, k)))
                    break;
            if (k == b)
            {
                printf("Exception (frob). F0 is singular.\n\n");
                abort();
            }

            for (j = 0; j < b; j++)
            {
                fmpz_poly_scalar_mul_fmpz(fmpz_poly_mat_entry(T, i, j), 
                                          fmpz_poly_mat_entry(Cinv, k, j), 
                                          padic_mat_unit(F0, i, k));
            }
        }

        fmpz_poly_mat_mul(F, C, T);
        vF = vC + padic_mat_val(F0) + vCinv;

        fmpz_pow_ui(pN, p, prec.N2 - vF);

        for (i = 0; i < b; i++)
            for (j = 0; j < b; j++)
            {
                fmpz_poly_struct *poly = fmpz_poly_mat_entry(F, i, j);
                long len = poly->length;

                fmpz_poly_truncate(poly, prec.K);
                if (len != 0)
                {
                    _fmpz_vec_scalar_mod_fmpz(poly->coeffs, poly->coeffs, len, pN);
                    _fmpz_poly_normalise(poly);
                }
            }

        fmpz_clear(pN);
        fmpz_poly_mat_clear(T);
    }

#if(DEBUG == 1)
printf("Matrix Fp(t):\n");
fmpz_poly_mat_print(F, "t");
printf("and a factor p^%ld, not necessarily the valuation.\n\n", vF);
#endif

    /* Step 5 {G = r(t)^m F(t)} **********************************************/

    {
        fmpz_t pN;
        fmpz_poly_t t;

        fmpz_init(pN);
        fmpz_poly_init(t);

        fmpz_pow_ui(pN, p, prec.N2 - vF);

        /* TODO: Could reduce this mod p^{N2-vF} */
        fmpz_poly_pow(t, r, prec.m);

        fmpz_poly_mat_scalar_mul_fmpz_poly(F, F, t);

        for (i = 0; i < b; i++)
            for (j = 0; j < b; j++)
            {
                fmpz_poly_struct *poly = fmpz_poly_mat_entry(F, i, j);
                long len = poly->length;

                fmpz_poly_truncate(poly, prec.K);
                if (len != 0)
                {
                    _fmpz_vec_scalar_mod_fmpz(poly->coeffs, 
                                              poly->coeffs, len, pN);
                    _fmpz_poly_normalise(poly);
                }
            }

        fmpz_clear(pN);
        fmpz_poly_clear(t);
    }

    /* Step 6 {F(1) = r(t_1)^{-m} G(t_1)} ************************************/

    {
        long N;
        fmpz_t f, g, t, pN;

        fmpz_init(f);
        fmpz_init(g);
        fmpz_init(t);
        fmpz_init(pN);

        N = prec.N2 - vF;
        fmpz_pow_ui(pN, p, N);

        /* f := \hat{t_1}, g := r(\hat{t_1})^{-m} */
        _padic_teichmuller(f, t1, p, N);
        _fmpz_mod_poly_evaluate_fmpz(t, r->coeffs, r->length, f, pN);
        _padic_inv(t, t, p, N);
        fmpz_powm_ui(g, t, prec.m, pN);

        /* F1 := g G(\hat{t_1}) */
        for (i = 0; i < b; i++)
            for (j = 0; j < b; j++)
            {
                const fmpz_poly_struct *poly = fmpz_poly_mat_entry(F, i, j);
                const long len               = poly->length;

                if (len == 0)
                {
                    fmpz_zero(padic_mat_unit(F1, i, j));
                }
                else
                {
                    _fmpz_mod_poly_evaluate_fmpz(t, poly->coeffs, len, f, pN);
                    fmpz_mul(padic_mat_unit(F1, i, j), g, t);
                    fmpz_mod(padic_mat_unit(F1, i, j), padic_mat_unit(F1, i, j), pN);
                }
            }
        padic_mat_val(F1) = vF;
        _padic_mat_canonicalise(F1, pctx_F0);

#if(DEBUG == 1)
printf("Matrix Fp(1):\n");
padic_mat_print_pretty(F1, pctx_F0);
printf("\n\n");
#endif

        fmpz_clear(f);
        fmpz_clear(g);
        fmpz_clear(t);
        fmpz_clear(pN);
    }

    /* Step 7 {Norm} *********************************************************/

    /* Step 8 {Reverse characteristic polynomial} ****************************/

    deformation_revcharpoly(cp, F1, n, p, a, prec.N0, prec.r, prec.s);

    #if(DEBUG == 1)
    printf("Reverse characteristic polynomial:\n");
    fmpz_poly_print_pretty(cp, "T");
    printf("\n\n");
    #endif

    /* Clean up **************************************************************/

    padic_ctx_clear(pctx_F0);
    padic_mat_clear(F0);

    mat_clear(M, ctxFracQt);
    free(bR);
    free(bC);
    fmpz_poly_clear(r);

    padic_ctx_clear(pctx_C);

    fmpz_poly_mat_clear(C);
    fmpz_poly_mat_clear(Cinv);

    fmpz_poly_mat_clear(F);
    padic_mat_clear(F1);
    fmpz_poly_clear(cp);
}

