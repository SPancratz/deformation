#include <stdlib.h>
#include <mpir.h>

#include "gmconnection.h"
#include "diagfrob.h"
#include "gmde.h"

#include "flint.h"
#include "fmpq_mat.h"

#include "deformation.h"

#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include <math.h>

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

void frob_with_precisions(mat_t F, const ctx_t ctxF, 
                          const mpoly_t P, const ctx_t ctxFracQt, 
                          long NWork, long Kfinite, long Kinfinite)
{
    const long n = P->n - 1;
    const long d = mpoly_degree(P, -1, ctxFracQt);
    const long b = gmc_basis_size(n, mpoly_degree(P, -1, ctxFracQt));
    const long p = fmpz_get_si(ctxF->pctx->p);

    long K;

    mon_t *bR, *bC;

    mat_t M, F0, C, Cinv;

    padic_ctx_t pctx;
    ctx_t ctxZp;
    ctx_t ctxZpt;

    long i, j, k;

time_t  t0, t1;  /* time_t is defined in <time.h> and <sys/types.h> as long */
clock_t c0, c1;  /* clock_t is defined in <time.h> and <sys/types.h> as int */
long double cputime;    /* CPU time, total */

    /* Init */
    padic_ctx_init(pctx, ctxF->pctx->p, NWork, PADIC_VAL_UNIT);
    ctx_init_padic(ctxZp, pctx);
    ctx_init_padic_poly(ctxZpt, pctx);

    mat_init(M, b, b, ctxFracQt);
    mat_init(F0, b, b, ctxZp);
    mat_init(C, b, b, ctxZpt);
    mat_init(Cinv, b, b, ctxZpt);

printf("Connection matrix:\n");
t0 = time(NULL);
c0 = clock();

    /* Compute Gauss--Manin connection */
    gmc_compute(M, &bR, &bC, P, ctxFracQt);

t1 = time(NULL);
c1 = clock();
cputime   = ((long double) (c1 - c0)) / ((long double) CLOCKS_PER_SEC);

printf ("Elapsed wall clock time:  %ld\n", (long) (t1 - t0));
printf ("Elapsed CPU time:         %LG\n", cputime);

    /* Find F(0) */
    {
        fmpz * a = _fmpz_vec_init(n + 1);

        mpoly_diagonal_fibre(a, P, ctxFracQt);

printf("Diagfrob:\n");
t0 = time(NULL);
c0 = clock();

        diagfrob(F0, a, n, d, ctxZp);

t1 = time(NULL);
c1 = clock();
cputime   = ((long double) (c1 - c0)) / ((long double) CLOCKS_PER_SEC);

printf ("Elapsed wall clock time:  %ld\n", (long) (t1 - t0));
printf ("Elapsed CPU time:         %LG\n", cputime);

        printf("Compute the diagonal fibre.\n");
        printf("a = {"), _fmpz_vec_print(a, n + 1), printf("}\n");

        _fmpz_vec_clear(a, n + 1);
    }

printf("Solve:\n");
t0 = time(NULL);
c0 = clock();

    /* Solve for C and Cinv *//*
    {
        mat_t Mt;

        padic_mat_struct *A, *Ainv;

        K = 4 * Kfinite + Kinfinite;

        A    = malloc(K * sizeof(padic_mat_struct));
        Ainv = malloc(((K + p - 1) / p) * sizeof(padic_mat_struct));
        for(i = 0; i < K; i++)
            _padic_mat_init(A + i, b, b);
        for(i = 0; i < (K + p - 1) / p; i++)
            _padic_mat_init(Ainv + i, b, b);

        /* Mt = - M^t *//*
        mat_init(Mt, b, b, ctxFracQt);
        mat_transpose(Mt, M, ctxFracQt);
        mat_neg(Mt, Mt, ctxFracQt);
        gmde_solve(A, K, pctx, M, ctxFracQt);
        gmde_solve(Ainv, (K + p - 1) / p, pctx, Mt, ctxFracQt);
        gmde_convert_soln(C, ctxZpt, A, K);
        gmde_convert_soln(Cinv, ctxZpt, Ainv, (K + p - 1) / p);
        mat_transpose(Cinv, Cinv, ctxZpt);

        for(i = 0; i < K; i++)
            _padic_mat_clear(A + i);
        for(i = 0; i < (K + p - 1) / p; i++)
            _padic_mat_clear(Ainv + i);
        free(A);
        free(Ainv);
        mat_clear(Mt, ctxFracQt);
    }
*/
    {
        ctx_t ctxQt;

        mat_t Mt;
        mat_t B, Binv;

        fmpq_mat_struct *A, *Ainv;

        K = 4 * Kfinite + Kinfinite;

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

t1 = time(NULL);
c1 = clock();
cputime   = ((long double) (c1 - c0)) / ((long double) CLOCKS_PER_SEC);

printf ("Elapsed wall clock time:  %ld\n", (long) (t1 - t0));
printf ("Elapsed CPU time:         %LG\n", cputime);

printf("C(t) F(0) C(t^p)^{-1}:\n");
t0 = time(NULL);
c0 = clock();

    /* Replace t by t^p in C^{-1} */
    for (i = 0; i < b; i++)
        for (j = 0; j < b; j++)
        {    
            _padic_poly_compose_pow(
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
                if (!_padic_is_zero((long *) mat_entry(F0, i, k, ctxZp)))
                    break;

            if (k == b)
            {
                printf("ERROR (frob_with_precisions). Bad F0.\n\n");
                abort();
            }

            for (j = 0; j < b; j++)
            {
                _padic_poly_scalar_mul_padic(
                    (padic_poly_struct *) mat_entry(RHS, i, j, ctxZpt), 
                    (padic_poly_struct *) mat_entry(Cinv, k, j, ctxZpt), 
                    (long *) mat_entry(F0, i, k, ctxZp));
            }
        }

        mat_mul(F, C, RHS, ctxZpt);

        for (i = 0; i < b; i++)
            for (j = 0; j < b; j++)
            {
                padic_poly_truncate(
                    (padic_poly_struct *) mat_entry(F, i, j, ctxZpt), 
                    K, ctxF->pctx);
            }

        mat_clear(RHS, ctxZpt);
    }

t1 = time(NULL);
c1 = clock();
cputime   = ((long double) (c1 - c0)) / ((long double) CLOCKS_PER_SEC);

printf ("Elapsed wall clock time:  %ld\n", (long) (t1 - t0));
printf ("Elapsed CPU time:         %LG\n", cputime);

printf("Rational function:\n");
t0 = time(NULL);
c0 = clock();

    /* Multiply F by r^{p prec(F) + 10%} */

    {
        fmpz_poly_t r, t;
        padic_poly_t s;

        fmpz_poly_init(r);
        fmpz_poly_init(t);
        _padic_poly_init(s);

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

        _padic_poly_set_fmpz_poly(s, r, pctx->p);
        padic_poly_pow(s, s, Kfinite, pctx);

        for (i = 0; i < b; i++)
            for (j = 0; j < b; j++)
            {
                padic_poly_mul(
                    (padic_poly_struct *) mat_entry(F, i, j, ctxF), 
                    (padic_poly_struct *) mat_entry(F, i, j, ctxF), 
                    s, ctxF->pctx);

                _padic_poly_truncate(
                    (padic_poly_struct *) mat_entry(F, i, j, ctxF), K, pctx->p);
            }

        fmpz_poly_clear(r);
        fmpz_poly_clear(t);
        _padic_poly_clear(s);
    }

t1 = time(NULL);
c1 = clock();
cputime   = ((long double) (c1 - c0)) / ((long double) CLOCKS_PER_SEC);

printf ("Elapsed wall clock time:  %ld\n", (long) (t1 - t0));
printf ("Elapsed CPU time:         %LG\n", cputime);

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

