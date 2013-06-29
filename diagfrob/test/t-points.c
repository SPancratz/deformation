/******************************************************************************

    Copyright (C) 2013 Sebastian Pancratz

******************************************************************************/

#include "diagfrob.h"

#define DEBUG  1

static void _visit(mp_limb_t *a, long n, long d, mp_limb_t p, long i, int nz, 
                   mp_limb_t *sum, mp_limb_t *pow, long *count)
{
    long j;

    for (j = (i < n || nz) ? 0 : 1; j <= (nz ? p - 1 : 1); j++)
    {
        pow[i] = n_powmod(j, d, p);
        *sum   = n_addmod(*sum, pow[i], p);
        
        if (i < n)  /* Visit next co-ordinate */
        {
            _visit(a, n, d, p, i + 1, nz || (j > 0), sum, pow, count);
        }
        else  /* Check this point */
        {
            *count += (sum == 0) ? 1 : 0;
        }

        *sum = n_submod(*sum, pow[i], p);
    }
}

/*
    Returns the number of points of the diagonal hypersurface 
    in projective n-space over F_p.
 */
static long points(fmpz *a, long n, long d, mp_limb_t p)
{
    mp_limb_t sum = 0;
    mp_limb_t pow[MON_MAX_VARS] = {0};
    mp_limb_t b[MON_MAX_VARS];
    long count = 0;
    long i;

    /* We know 0 <= a[i] < p, and p is small. */
    for (i = 0; i <= n; i++)
        b[i] = a[i];

    _visit(b, n, d, p, 0, 0, &sum, pow, &count);

    return count;
}

int main(void)
{
    int i, result;
    flint_rand_t state;

    printf("points... ");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 100; i++)
    {
        const long d = 4; /* n_randint(state, 5) + 2;  /* d in [2,6] */
        const long n = 2; /* n_randint(state, 3) + 2;  /* n in [2,4] */
        const long a = 1;

        const long lenB = gmc_basis_size(n, d);
        long j;

        padic_ctx_t pctx;
        fmpz_t p;
        long N;

        fmpz *A;

        long pts1, pts2;

        padic_mat_t F;
        fmpz_poly_t chi;

        fmpz_init(p);
        do 
            *p = n_randprime(state, 5, 1);
        while (d % *p == 0);

        A = _fmpz_vec_init(n + 1);
        for (j = 0; j <= n; j++)
            fmpz_set_ui(A + j, n_randint(state, *p - 1) + 1);

        N = diagfrob_prec_phi(n, d, p, a);

        padic_ctx_init(pctx, p, N, PADIC_SERIES);

        padic_mat_init(F, lenB, lenB);
        diagfrob(F, A, n, d, pctx, 0);

        fmpz_poly_init(chi);
        diagfrob_revcharpoly(chi, F, pctx);

        /*********************************************************************/
        #if DEBUG
        printf("i    = %d\n", i);
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

        /* N_1 = (-1)^n (\chi'(T) / \chi(T))|_{T=0} + 1 + \dotsb + q^{n-1}   */
        pts1 = ((n % 2 == 0) ? 1 : -1) * (chi->coeffs[1] / chi->coeffs[0]) 
               + (n_pow(*p, n) - 1) / (*p - 1);
        pts2 = points(A, n, d, *p);

        result = ((chi->coeffs[1] % chi->coeffs[0] == 0) && (pts1 == pts2));
        if (!result)
        {
            printf("FAIL:\n");
            printf("A = {"), _fmpz_vec_print(A, n + 1), printf("}\n");
            printf("d = %ld\n", d);
            printf("n = %ld\n", n);
            printf("p = %ld\n", *p);
            printf("chi = "), fmpz_poly_print_pretty(chi, "T"), printf("\n");
            printf("pts1 = %ld\n", pts1);
            printf("pts2 = %ld\n", pts2);
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

