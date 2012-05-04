/******************************************************************************

    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>

#include "flint.h"
#include "ulong_extras.h"
#include "gmconnection.h"
#include "diagfrob.h"

#define DEBUG  0

/*
    Returns the number of points on the elliptic curve 
    \begin{equation*}
    E \colon A X^3 + B Y^3 + C Z^3 = 0
    \end{equation*}
    in projective space over $\mathbf{F}_p$.

    Assumes that the number of points fits into a signed integer.
 */
static 
long points(const fmpz_t A, const fmpz_t B, const fmpz_t C, const fmpz_t p)
{
    fmpz_t a, b, Cinv, r, s, t, x, y;
    long N = 0;

    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(Cinv);
    fmpz_init(r);
    fmpz_init(s);
    fmpz_init(t);
    fmpz_init(x);
    fmpz_init(y);

    /* Deal with the points where Z = 0 */

    fmpz_invmod(b, B, p);
    fmpz_mul(a, A, b);

    for (fmpz_set_ui(x, 1); fmpz_cmp(x, p) < 0; fmpz_add_ui(x, x, 1))
    {
        fmpz_pow_ui(r, x, 3);
        fmpz_mul(r, a, r);
        fmpz_add_ui(r, r, 1);
        fmpz_mod(r, r, p);
        if (fmpz_is_zero(r))
            N++;
    }

    /* Deal with the points where Z != 0 */

    fmpz_invmod(Cinv, C, p);
    fmpz_mul(a, A, Cinv);
    fmpz_mul(b, B, Cinv);
    fmpz_mod(a, a, p);
    fmpz_mod(b, b, p);

    for (fmpz_set_ui(x, 0); fmpz_cmp(x, p) < 0; fmpz_add_ui(x, x, 1))
    {
        for (fmpz_set_ui(y, 0); fmpz_cmp(y, p) < 0; fmpz_add_ui(y, y, 1))
        {
            fmpz_pow_ui(r, x, 3);
            fmpz_mul(r, a, r);
            fmpz_pow_ui(s, y, 3);
            fmpz_mul(s, b, s);
            fmpz_add(t, r, s);
            fmpz_add_ui(t, t, 1);
            fmpz_mod(t, t, p);

            if (fmpz_is_zero(t))
                N++;
        }
    }

    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(Cinv);
    fmpz_clear(r);
    fmpz_clear(s);
    fmpz_clear(t);
    fmpz_clear(x);
    fmpz_clear(y);

    return N;
}

int main(void)
{
    int cases, result;
    flint_rand_t state;

    /*
       Consider the elliptic curve E in projective form 

           a X^3 + b Y^3 + c Z^3 = 0

       where, without loss of generality we choose -b = 1, with 
       model for the affine part 

           y^3 = a x^3 + c

       If 1 + p - a_p is the number of points of E over Fp, then 
       the zeta-function is given by 

           (1 - a_p T + p T^2) [(1 - T) (1 - pT)]^{-1}
     */

    printf("ell_curves... ");
    fflush(stdout);
    _randinit(state);

    for (cases = 0; cases < 100; cases++)
    {
        long d = 3;
        long n = 2;

        long N, r, s;
        long prime;
        fmpz_t p;
        padic_ctx_t pctx1, pctx2;

        ctx_t ctx;

        fmpz *a;

        mat_t F;
        char *f;
        fmpz *poly;

        long i, lenB = gmc_basis_size(n, d);

        long pts_naive;

        fmpz_init(p);
        do 
            prime = n_randprime(state, 5, 1);
        while (d % prime == 0);
        fmpz_set_ui(p, prime);

        a = _fmpz_vec_init(n + 1);
        for (i = 0; i <= n; i++)
            fmpz_set_ui(a + i, n_randint(state, prime - 1) + 1);

        /* Set precisions */

        N = diagfrob_charpoly_prec(n, d, p, 1);
        diagfrob_matrix_prec(&r, &s, n, p);

        padic_ctx_init(pctx2, p, N + r + s, PADIC_SERIES);
        padic_ctx_init(pctx1, p, N, PADIC_SERIES);

        ctx_init_padic(ctx, pctx2);

        mat_init(F, lenB, lenB, ctx);
        diagfrob(F, a, n, d, ctx, 0);

        f = _vec_init(lenB + 1, ctx);
        poly = _fmpz_vec_init(lenB + 1);

        ctx->pctx = pctx1;
        mat_revcharpoly(f, F, ctx);

        /*********************************************************************/
        #if DEBUG
        printf("p = %ld, ", prime);
        printf("a = {"), _fmpz_vec_print(a, n + 1), printf("}\n");

        printf("Precision: (N, r, s) = (%ld, %ld, %ld)\n", N, r, s);

        printf("F to precision N + r + s = %ld, \n", N + r + s);
        mat_print(F, ctx);
        printf("\n");

        printf("f = det(1 - t F) to precision N = %ld, \n", N);
        _vec_print(f, lenB + 1, ctx);
        printf("\n");

        fflush(stdout);
        #endif
        /*********************************************************************/

        diagfrob_revcharpoly(poly, f, n, lenB, ctx);

        /* Alternative computation */

        pts_naive = points(a + 0, a + 1, a + 2, p);

        result = (1 + prime - (- fmpz_get_si(poly + 1)) == pts_naive);
        if (!result)
        {
            printf("FAIL:\n");
            printf("a = {"), _fmpz_vec_print(a, n + 1), printf("}\n");
            printf("d = %ld\n", d);
            printf("n = %ld\n", n);
            printf("p = %ld\n", prime);
            printf("f = "), _fmpz_vec_print(poly, lenB + 1), printf("\n");
            abort();
        }

        /* Clean-up */

        fmpz_clear(p);
        _fmpz_vec_clear(a, n + 1);
        mat_clear(F, ctx);
        _vec_clear(f, lenB + 1, ctx);
        _fmpz_vec_clear(poly, lenB + 1);

        padic_ctx_clear(pctx1);
        padic_ctx_clear(pctx2);
        ctx_clear(ctx);
    }

    _randclear(state);

    printf("PASS\n");
    return EXIT_SUCCESS;
}

