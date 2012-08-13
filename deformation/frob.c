#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include <mpir.h>

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

static 
void fmpz_poly_mat_get_fmpz_mat(fmpz_mat_t B, const fmpz_poly_mat_t A)
{
    long i, j;

    if (!(A->r == B->r && A->c == B->c))
    {
        printf("ERROR (fmpz_poly_mat_get_fmpz_mat).  Incompatible dimensions.\n");
        abort();
    }

    for (i = 0; i < A->r; i++)
        for (j = 0; j < A->c; j++)
        {
            const fmpz_poly_struct *poly = fmpz_poly_mat_entry(A, i, j);
            const long len               = poly->length;

            if (len == 0)
            {
                fmpz_zero(fmpz_mat_entry(B, i, j));
            }
            else if (len == 1)
            {
                fmpz_set(fmpz_mat_entry(B, i, j), poly->coeffs + 0);
            }
            else
            {
                printf("ERROR (fmpz_poly_mat_get_fmpz_mat).\n");
                printf("A contains a polynomial of length greater than 1.\n");
                abort();
            }
        }
}

static 
void fmpz_poly_mat_scalar_mod_fmpz(fmpz_poly_mat_t B, 
                                   const fmpz_poly_mat_t A, const fmpz_t p)
{
    long i, j;

    for (i = 0; i < A->r; i++)
        for (j = 0; j < A->c; j++)
            fmpz_poly_scalar_mod_fmpz(fmpz_poly_mat_entry(B, i, j), 
                                      fmpz_poly_mat_entry(A, i, j), p);
}

static 
void fmpz_poly_mat_scalar_divexact_fmpz(fmpz_poly_mat_t B, 
                                        const fmpz_poly_mat_t A, const fmpz_t f)
{
    long i, j;

    for (i = 0; i < A->r; i++)
        for (j = 0; j < A->c; j++)
            fmpz_poly_scalar_divexact_fmpz(fmpz_poly_mat_entry(B, i, j),
                                           fmpz_poly_mat_entry(A, i, j), f);
}

static 
long fmpz_poly_mat_ord_p(const fmpz_poly_mat_t A, const fmpz_t p)
{
    long i, j, v, w = LONG_MAX;

    for (i = 0; i < A->r; i++)
        for (j = 0; j < A->c; j++)
        {
            const fmpz_poly_struct *poly = fmpz_poly_mat_entry(A, i, j);

            if (!fmpz_poly_is_zero(poly))
            {
                v = _fmpz_vec_ord_p(poly->coeffs, poly->length, p);
                w = FLINT_MIN(v, w);
                if (w == 0)
                    return 0;
            }
        }
    return w;
}

static 
void fmpz_poly_mat_canonicalise(fmpz_poly_mat_t A, long *vA, const fmpz_t p)
{
    const long w = fmpz_poly_mat_ord_p(A, p);

    if (w == LONG_MAX)
    {
        *vA = 0;
    }
    else if (w > 0)
    {
        fmpz_t f;

        fmpz_init(f);
        fmpz_pow_ui(f, p, w);

        fmpz_poly_mat_scalar_divexact_fmpz(A, A, f);
        *vA += w;

        fmpz_clear(f);
    }
}

/* 
    Applies the operator $\sigma^e$ to all elements in the matrix \code{op}, 
    setting the corresponding elements in \code{rop} to the results.
 */
static 
void fmpz_poly_mat_frobenius(fmpz_poly_mat_t B, 
                             const fmpz_poly_mat_t A, long e, 
                             const fmpz_t p, long N, const qadic_ctx_t ctx)
{
    const long d = qadic_ctx_degree(ctx);

    e = e % d;
    if (e < 0)
        e += d;

    if (e == 0)
    {
        fmpz_t pN;

        fmpz_init(pN);
        fmpz_pow_ui(pN, p, N);

        fmpz_poly_mat_scalar_mod_fmpz(B, A, pN);

        fmpz_clear(pN);
    }
    else
    {
        long i, j;
        fmpz *t = _fmpz_vec_init(2 * d - 1);

        for (i = 0; i < B->r; i++)
            for (j = 0; j < B->c; j++)
            {
                const fmpz_poly_struct *a = fmpz_poly_mat_entry(A,  i, j);
                fmpz_poly_struct *b       = fmpz_poly_mat_entry(B, i, j);

                if (a->length == 0)
                {
                    fmpz_poly_zero(b);
                }
                else
                {
                    _qadic_frobenius(t, a->coeffs, a->length, e, 
                                     ctx->a, ctx->j, ctx->len, p, N);

                    fmpz_poly_fit_length(b, d);
                    _fmpz_vec_set(b->coeffs, t, d);
                    _fmpz_poly_set_length(b, d);
                    _fmpz_poly_normalise(b);
                }
            }

        _fmpz_vec_clear(t, 2 * d - 1);
    }
}

static 
void _qadic_mat_mul(fmpz_poly_mat_t C, 
                    const fmpz_poly_mat_t A, const fmpz_poly_mat_t B, 
                    const fmpz_t pN, const qadic_ctx_t ctx)
{
    long i, j;

    fmpz_poly_mat_mul(C, A, B);

    for (i = 0; i < C->r; i++)
        for (j = 0; j < C->c; j++)
        {
            fmpz_poly_struct *poly = fmpz_poly_mat_entry(C, i, j);

            _fmpz_mod_poly_reduce(poly->coeffs, poly->length, 
                                  ctx->a, ctx->j, ctx->len, pN);
        }
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
 */

void frob(const mpoly_t P, const ctx_t ctxFracQt, 
          const qadic_t t1, const qadic_ctx_t Qq, 
          prec_t *prec, const prec_t *prec_in, 
          int verbose)
{
    const padic_ctx_struct *Qp = &Qq->pctx;
    const fmpz *p = Qp->p;
    const long a  = qadic_ctx_degree(Qq);
    const long n  = P->n - 1;
    const long d  = mpoly_degree(P, -1, ctxFracQt);
    const long b  = gmc_basis_size(n, d);

    long i, j, k;

    /* Diagonal fibre */
    padic_mat_t F0;

    /* Gauss--Manin Connection */
    mat_t M;
    mon_t *bR, *bC;
    fmpz_poly_t r;

    /* Local solution */
    fmpz_poly_mat_t C, Cinv;
    long vC, vCinv;

    /* Frobenius */
    fmpz_poly_mat_t F;
    long vF;

    fmpz_poly_mat_t F1;
    long vF1;

    fmpz_poly_t cp;

    clock_t c0, c1;
    double c;

    if (verbose)
    {
        printf("Input:\n");
        printf("  P  = "), mpoly_print(P, ctxFracQt), printf("\n");
        printf("  p  = "), fmpz_print(p), printf("\n");
        printf("  t1 = "), qadic_print_pretty(t1, Qq), printf("\n");
        printf("\n");
        fflush(stdout);
    }

    /* Step 1 {M, r} *********************************************************/

    c0 = clock();

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

    c1 = clock();
    c  = (double) (c1 - c0) / CLOCKS_PER_SEC;

    if (verbose)
    {
        printf("Gauss-Manin connection:\n");
        printf("  r(t) = "), fmpz_poly_print_pretty(r, "t"), printf("\n");
        printf("  Time = %f\n", c);
        printf("\n");
        fflush(stdout);
    }

    /* Precisions ************************************************************/

    if (prec_in != NULL)
    {
        *prec = *prec_in;
    }
    else
    {
        deformation_precisions(prec, p, a, n, d, fmpz_poly_degree(r));
    }

    if (verbose)
    {
        printf("Precisions:\n");
        printf("  N0   = %ld\n", prec->N0);
        printf("  N1   = %ld\n", prec->N1);
        printf("  N2   = %ld\n", prec->N2);
        printf("  N3   = %ld\n", prec->N3);
        printf("  N3i  = %ld\n", prec->N3i);
        printf("  N3w  = %ld\n", prec->N3w);
        printf("  N3iw = %ld\n", prec->N3iw);
        printf("  N4   = %ld\n", prec->N4);
        printf("  m    = %ld\n", prec->m);
        printf("  K    = %ld\n", prec->K);
        printf("  r    = %ld\n", prec->r);
        printf("  s    = %ld\n", prec->s);
        printf("\n");
        fflush(stdout);
    }

    /* Initialisation ********************************************************/

    padic_mat_init(F0, b, b);

    fmpz_poly_mat_init(C, b, b);
    fmpz_poly_mat_init(Cinv, b, b);

    fmpz_poly_mat_init(F, b, b);
    vF = 0;

    fmpz_poly_mat_init(F1, b, b);
    vF1 = 0;

    fmpz_poly_init(cp);

    /* Step 2 {F0} ***********************************************************/

    {
        padic_ctx_t pctx_F0;
        fmpz *t;

        padic_ctx_init(pctx_F0, p, prec->N4, PADIC_VAL_UNIT);
        t = _fmpz_vec_init(n + 1);

        c0 = clock();

        mpoly_diagonal_fibre(t, P, ctxFracQt);
        diagfrob(F0, t, n, d, pctx_F0, 0);
        padic_mat_transpose(F0, F0);

        c1 = clock();
        c  = (double) (c1 - c0) / CLOCKS_PER_SEC;

        if (verbose)
        {
            printf("Diagonal fibre:\n");
            printf("  P(0) = {"), _fmpz_vec_print(t, n + 1), printf("}\n");
            printf("  Time = %f\n", c);
            printf("\n");
            fflush(stdout);
        }

        _fmpz_vec_clear(t, n + 1);
        padic_ctx_clear(pctx_F0);
    }

    /* Step 3 {C, Cinv} ******************************************************/
    /*
        Compute C as a matrix over Z_p[[t]].  A is the same but as a series 
        of matrices over Z_p.  Mt is the matrix -M^t, and Cinv is C^{-1}^t, 
        the local solution of the differential equation replacing M by Mt.
     */

    c0 = clock();
    {
        const long K = prec->K;
        padic_mat_struct *A;

        A    = malloc(K * sizeof(padic_mat_struct));
        for(i = 0; i < K; i++)
            padic_mat_init(A + i, b, b);

        gmde_solve(A, K, p, prec->N3, prec->N3w, M, ctxFracQt);
        gmde_convert_soln(C, &vC, A, K, p);

        for(i = 0; i < K; i++)
            padic_mat_clear(A + i);
        free(A);
    }
    c1 = clock();
    c  = (double) (c1 - c0) / CLOCKS_PER_SEC;
    if (verbose)
    {
        printf("Local solution:\n");
        printf("  Time for C      = %f\n", c);
        fflush(stdout);
    }

    c0 = clock();
    {
        const long K = (prec->K + (*p) - 1) / (*p);
        mat_t Mt;
        padic_mat_struct *Ainv;

        Ainv = malloc(K * sizeof(padic_mat_struct));
        for(i = 0; i < K; i++)
            padic_mat_init(Ainv + i, b, b);

        mat_init(Mt, b, b, ctxFracQt);
        mat_transpose(Mt, M, ctxFracQt);
        mat_neg(Mt, Mt, ctxFracQt);
        gmde_solve(Ainv, K, p, prec->N3i, prec->N3iw, Mt, ctxFracQt);
        gmde_convert_soln(Cinv, &vCinv, Ainv, K, p);

        fmpz_poly_mat_transpose(Cinv, Cinv);
        fmpz_poly_mat_compose_pow(Cinv, Cinv, *p);

        for(i = 0; i < K; i++)
            padic_mat_clear(Ainv + i);
        free(Ainv);
        mat_clear(Mt, ctxFracQt);
    }
    c1 = clock();
    c  = (double) (c1 - c0) / CLOCKS_PER_SEC;
    if (verbose)
    {
        printf("  Time for C^{-1} = %f\n", c);
        printf("\n");
        fflush(stdout);
    }

    /* Step 4 {F(t) := C(t) F(0) C(t^p)^{-1}} ********************************/
    /*
        Computes the product C(t) F(0) C(t^p)^{-1} modulo (p^{N_2}, t^K). 
        This is done by first computing the unit part of the product 
        exactly over the integers modulo t^K.
     */

    c0 = clock();
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
        fmpz_poly_mat_truncate(F, prec->K);
        vF = vC + padic_mat_val(F0) + vCinv;

        /* Canonicalise (F, vF) */
        {
            long v = fmpz_poly_mat_ord_p(F, p);

            if (v == LONG_MAX)
            {
                printf("ERROR (deformation_frob).  F(t) == 0.\n");
                abort();
            }
            else if (v > 0)
            {
                fmpz_pow_ui(pN, p, v);
                fmpz_poly_mat_scalar_divexact_fmpz(F, F, pN);
                vF = vF + v;
            }
        }

        /* Reduce (F, vF) modulo p^{N2} */
        fmpz_pow_ui(pN, p, prec->N2 - vF);
        fmpz_poly_mat_scalar_mod_fmpz(F, F, pN);

        fmpz_clear(pN);
        fmpz_poly_mat_clear(T);
    }
    c1 = clock();
    c  = (double) (c1 - c0) / CLOCKS_PER_SEC;
    if (verbose)
    {
        printf("Matrix for F(t):\n");
        printf("  Time = %f\n", c);
        printf("\n");
        fflush(stdout);
    }

    /* Step 5 {G = r(t)^m F(t)} **********************************************/

    c0 = clock();
    {
        fmpz_t pN;
        fmpz_poly_t t;

        fmpz_init(pN);
        fmpz_poly_init(t);

        fmpz_pow_ui(pN, p, prec->N2 - vF);

        /* TODO: Could reduce this mod p^{N2-vF} */
        {
            fmpz_mod_poly_t _t;

            fmpz_mod_poly_init(_t, pN);
            fmpz_mod_poly_set_fmpz_poly(_t, r);
            fmpz_mod_poly_pow(_t, _t, prec->m);
            fmpz_mod_poly_get_fmpz_poly(t, _t);
            fmpz_mod_poly_clear(_t);
        }

        fmpz_poly_mat_scalar_mul_fmpz_poly(F, F, t);
        fmpz_poly_mat_scalar_mod_fmpz(F, F, pN);

        /* TODO: This should not be necessary? */
        fmpz_poly_mat_truncate(F, prec->K);

        fmpz_clear(pN);
        fmpz_poly_clear(t);
    }
    c1 = clock();
    c  = (double) (c1 - c0) / CLOCKS_PER_SEC;
    if (verbose)
    {
        printf("Analytic continuation:\n");
        printf("  Time = %f\n", c);
        printf("\n");
        fflush(stdout);
    }

    /* Steps 6 and 7 *********************************************************/

    if (a == 1)
    {
        /* Step 6 {F(1) = r(t_1)^{-m} G(t_1)} ********************************/

        c0 = clock();
        {
            const long N = prec->N2 - vF;

            fmpz_t f, g, t, pN;

            fmpz_init(f);
            fmpz_init(g);
            fmpz_init(t);
            fmpz_init(pN);

            fmpz_pow_ui(pN, p, N);

            /* f := \hat{t_1}, g := r(\hat{t_1})^{-m} */
            _padic_teichmuller(f, t1->coeffs + 0, p, N);
            _fmpz_mod_poly_evaluate_fmpz(t, r->coeffs, r->length, f, pN);
            _padic_inv(t, t, p, N);
            fmpz_powm_ui(g, t, prec->m, pN);

            /* F1 := g G(\hat{t_1}) */
            for (i = 0; i < b; i++)
                for (j = 0; j < b; j++)
                {
                    const fmpz_poly_struct *poly = fmpz_poly_mat_entry(F, i, j);
                    const long len               = poly->length;

                    if (len == 0)
                    {
                        fmpz_poly_zero(fmpz_poly_mat_entry(F1, i, j));
                    }
                    else
                    {
                        fmpz_poly_fit_length(fmpz_poly_mat_entry(F1, i, j), 1);

                        _fmpz_mod_poly_evaluate_fmpz(t, poly->coeffs, len, f, pN);
                        fmpz_mul(fmpz_poly_mat_entry(F1, i, j)->coeffs + 0, g, t);
                        fmpz_mod(fmpz_poly_mat_entry(F1, i, j)->coeffs + 0, 
                                 fmpz_poly_mat_entry(F1, i, j)->coeffs + 0, pN);

                        _fmpz_poly_set_length(fmpz_poly_mat_entry(F1, i, j), 1);
                        _fmpz_poly_normalise(fmpz_poly_mat_entry(F1, i, j));
                    }
                }

            vF1 = vF;
            fmpz_poly_mat_canonicalise(F1, &vF1, p);

            fmpz_clear(f);
            fmpz_clear(g);
            fmpz_clear(t);
            fmpz_clear(pN);
        }
        c1 = clock();
        c  = (double) (c1 - c0) / CLOCKS_PER_SEC;
        if (verbose)
        {
            printf("Evaluation:\n");
            printf("  Time = %f\n", c);
            printf("\n");
            fflush(stdout);
        }
    }
    else
    {
        /* Step 6 {F(1) = r(t_1)^{-m} G(t_1)} ********************************/

        c0 = clock();
        {
            const long N = prec->N2 - vF;
            fmpz_t e, pN;
            fmpz *f, *g, *t;

            fmpz_init_set_ui(e, prec->m);
            fmpz_init(pN);

            f = _fmpz_vec_init(a);
            g = _fmpz_vec_init(2 * a - 1);
            t = _fmpz_vec_init(2 * a - 1);

            fmpz_pow_ui(pN, p, N);

            /* f := \hat{t_1}, g := r(\hat{t_1})^{-m} */
            _qadic_teichmuller(f, t1->coeffs, t1->length, Qq->a, Qq->j, Qq->len, p, N);
            _fmpz_mod_poly_compose_mod(g, r->coeffs, r->length, f, a, 
                                          Qq->a, Qq->j, Qq->len, pN);
            _qadic_inv(t, g, a, Qq->a, Qq->j, Qq->len, p, N);
            _qadic_pow(g, t, a, e, Qq->a, Qq->j, Qq->len, pN);

            /* F1 := g G(\hat{t_1}) */
            for (i = 0; i < b; i++)
                for (j = 0; j < b; j++)
                {
                    const fmpz_poly_struct *poly = fmpz_poly_mat_entry(F, i, j);
                    const long len               = poly->length;

                    fmpz_poly_struct *poly2 = fmpz_poly_mat_entry(F1, i, j);

                    if (len == 0)
                    {
                        fmpz_poly_zero(poly2);
                    }
                    else
                    {
                        _fmpz_mod_poly_compose_mod(t, poly->coeffs, len, f, a, 
                                                      Qq->a, Qq->j, Qq->len, pN);

                        fmpz_poly_fit_length(poly2, 2 * a - 1);
                        _fmpz_poly_mul(poly2->coeffs, g, a, t, a);
                        _fmpz_mod_poly_reduce(poly2->coeffs, 2 * a - 1, Qq->a, Qq->j, Qq->len, pN);
                        _fmpz_poly_set_length(poly2, a);
                        _fmpz_poly_normalise(poly2);
                    }
                }

            /* Now the matrix for p^{-1} F_p at t=t_1 is (F1, vF1). */
            vF1 = vF;
            fmpz_poly_mat_canonicalise(F1, &vF1, p);

            fmpz_clear(e);
            fmpz_clear(pN);
            _fmpz_vec_clear(f, a);
            _fmpz_vec_clear(g, 2 * a - 1);
            _fmpz_vec_clear(t, 2 * a - 1);
        }
        c1 = clock();
        c  = (double) (c1 - c0) / CLOCKS_PER_SEC;
        if (verbose)
        {
            printf("Evaluation:\n");
            printf("  Time = %f\n", c);
            printf("\n");
            fflush(stdout);
        }

        /* Step 7 {Norm} *****************************************************/
        /*
            Computes the matrix for $q^{-1} F_q$ at $t = t_1$ as the 
            product $F \sigma(F) \dotsm \sigma^{a-1}(F)$ up appropriate 
            transpositions because our convention of columns vs rows is 
            the opposite of that used by Gerkmann.

            Note that, in any case, transpositions do not affect 
            the characteristic polynomial.
         */

        c0 = clock();
        {
            const long N = prec->N1 - a * vF1;

            fmpz_t pN;
            fmpz_poly_mat_t T;

            fmpz_init(pN);
            fmpz_poly_mat_init(T, b, b);

            fmpz_pow_ui(pN, p, N);

            fmpz_poly_mat_frobenius(T, F1, 1, p, N, Qq);
            _qadic_mat_mul(F1, F1, T, pN, Qq);

            for (i = 2; i < a; i++)
            {
                fmpz_poly_mat_frobenius(T, T, 1, p, N, Qq);
                _qadic_mat_mul(F1, F1, T, pN, Qq);
            }

            vF1 = a * vF1;
            fmpz_poly_mat_canonicalise(F1, &vF1, p);

            fmpz_clear(pN);
            fmpz_poly_mat_clear(T);
        }
        c1 = clock();
        c  = (double) (c1 - c0) / CLOCKS_PER_SEC;
        if (verbose)
        {
            printf("Norm:\n");
            printf("  Time = %f\n", c);
            printf("\n");
            fflush(stdout);
        }
    }

    /* Step 8 {Reverse characteristic polynomial} ****************************/

    c0 = clock();

    deformation_revcharpoly(cp, F1, vF1, n, prec->N0, prec->r, prec->s, Qq);

    c1 = clock();
    c  = (double) (c1 - c0) / CLOCKS_PER_SEC;
    if (verbose)
    {
        printf("Reverse characteristic polynomial:\n");
        printf("  p(T) = "), fmpz_poly_print_pretty(cp, "T"), printf("\n");
        printf("  Time = %f\n", c);
        printf("\n");
        fflush(stdout);
    }

    /* Clean up **************************************************************/

    padic_mat_clear(F0);

    mat_clear(M, ctxFracQt);
    free(bR);
    free(bC);
    fmpz_poly_clear(r);

    fmpz_poly_mat_clear(C);
    fmpz_poly_mat_clear(Cinv);

    fmpz_poly_mat_clear(F);
    fmpz_poly_mat_clear(F1);
    fmpz_poly_clear(cp);
}

