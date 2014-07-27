/******************************************************************************

    Copyright (C) 2011, 2013 Sebastian Pancratz

******************************************************************************/

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

    /* Points (X : Y : Z) with Z = 0 */

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

    /* Points (X : Y : Z) with Z != 0 */

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
    int i, result;
    flint_rand_t state;

    printf("ell_curves... ");
    fflush(stdout);

    flint_randinit(state);

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

    for (i = 0; i < 5000; i++)
    {
        const long d = 3;
        const long n = 2;
        const long a = 1;

        const long lenB = gmc_basis_size(n, d);
        long j;

        fmpz_t p;

        long N;

        padic_ctx_t pctx;

        fmpz *A;

        long pts_naive;

        padic_mat_t F;
        fmpz_poly_t chi;

        fmpz_init(p);
        do
        {
            ulong bits = n_randint(state, 5) + 2;
            *p = n_randprime(state, bits, 1);
        }
        while (d % *p == 0);

        A = _fmpz_vec_init(n + 1);
        for (j = 0; j <= n; j++)
            fmpz_set_ui(A + j, n_randint(state, *p - 1) + 1);

        N = diagfrob_prec_phi(n, d, p, a);

        padic_ctx_init(pctx, p, FLINT_MAX(0, N-10), N+10, PADIC_SERIES);

        padic_mat_init2(F, lenB, lenB, N);
        diagfrob(F, A, n, d, padic_mat_prec(F), pctx, 0);

        fmpz_poly_init(chi);
        diagfrob_revcharpoly(chi, F, pctx);

        /*********************************************************************/
        #if DEBUG
        printf("i    = %ld\n", i);
        printf("n    = %ld\n", n);
        printf("d    = %ld\n", d);
        printf("lenB = %ld\n", lenB);
        printf("p    = %ld\n", *p);
        printf("A    = {"), _fmpz_vec_print(A, n + 1), printf("}\n");
        printf("N    = %ld\n", N);

        printf("F mod %ld^%ld:\n", *p, N);
        padic_mat_print_pretty(F, pctx);
        printf("\n");

        printf("chi:\n");
        fmpz_poly_print_pretty(chi, "T");
        printf("\n");

        fflush(stdout);
        #endif
        /*********************************************************************/

        diagfrob_zetafunction(chi, n, d, p, a);

        /*********************************************************************/
        #if DEBUG
        printf("Zeta function:\n");
        fmpz_poly_print_pretty(chi, "T");
        printf("\n");

        fflush(stdout);
        #endif
        /*********************************************************************/

        /* Alternative computation */
        pts_naive = points(A + 0, A + 1, A + 2, p);

        result = (1 + *p - (- fmpz_poly_get_coeff_si(chi, 1)) == pts_naive);
        if (!result)
        {
            printf("FAIL:\n");
            printf("A = {"), _fmpz_vec_print(A, n + 1), printf("}\n");
            printf("d = %ld\n", d);
            printf("n = %ld\n", n);
            printf("p = %ld\n", *p);
            printf("chi = "), fmpz_poly_print_pretty(chi, "T"), printf("\n");
            printf("Points = %ld\n", pts_naive);
            abort();
        }

        /* Clean-up */
        fmpz_clear(p);
        _fmpz_vec_clear(A, n + 1);
        padic_mat_clear(F);
        fmpz_poly_clear(chi);

        padic_ctx_clear(pctx);
    }

    flint_randclear(state);

    printf("PASS\n");
    return EXIT_SUCCESS;
}

