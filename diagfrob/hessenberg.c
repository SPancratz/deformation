#include "diagfrob.h"

#include <math.h>

#define M_LOG2E  1.44269504088896340736  /* log2(e) */

static __inline__ long double _log2(const long double x)
{
    return log(x) * M_LOG2E;
}

void nmod_mat_hessenberg(nmod_mat_t rop, const nmod_mat_t op)
{
    const long n = op->r;
    const nmod_t mod = op->mod;

    mp_limb_t **a;
    mp_limb_t inv, s, t;

    long i, j, m;

    nmod_mat_set(rop, op);

    a = rop->rows;

    for (m = 1; m < n - 1; m++)
    {
        for (i = m + 1; i < n && a[i][m - 1] == 0; i++) ;

        if (i < n)
        {
            if (a[m][m - 1] != 0)
                i = m;

            if (i > m)
            {
                /* Row swap in O(1) */
                mp_limb_t *r;
                r    = a[i];
                a[i] = a[m];
                a[m] = r;
                /* Column swap in O(n) */
                for (j = 0; j < n; j++)
                {
                    t       = a[j][i];
                    a[j][i] = a[j][m];
                    a[j][m] = t;
                }
            }

            inv = n_invmod(a[m][m - 1], mod.n);

            for (i = m + 1; i < n; i++)
            {
                if (a[i][m - 1] != 0)
                {
                    t = n_mulmod2_preinv(a[i][m - 1], inv, mod.n, mod.ninv);
                    s = n_negmod(t, mod.n);

                    for (j = 0; j < n; j++)
                        NMOD_ADDMUL(a[i][j], a[m][j], s, mod);

                    for (j = 0; j < n; j++)
                        NMOD_ADDMUL(a[j][m], a[j][i], t, mod);
                }
            }
        }
    }
}

void nmod_mat_charpoly(nmod_poly_t rop, const nmod_mat_t op)

{
    const long n = op->r;
    const nmod_t mod = op->mod;

    nmod_mat_t H, T;
    mp_limb_t **h, **t, x, y;

    long i, j, m;

    nmod_mat_init_set(H, op);
    nmod_mat_hessenberg(H, H);
    nmod_mat_init(T, n + 1, n + 1, op->mod.n);

    h = H->rows;
    t = T->rows;

    t[0][0] = 1;
    for (m = 1; m <= n; m++)
    {
        for (j = 1; j <= n; j++)
            t[m][j] = t[m-1][j-1];
        x = n_negmod(h[m-1][m-1], mod.n);
        for (j = 0; j <= n; j++)
            NMOD_ADDMUL(t[m][j], t[m-1][j], x, mod);
        y = 1;
        for (i = 1; i < m; i++)
        {
            y = n_mulmod2_preinv(y, h[m-i][m-i-1], mod.n, mod.ninv);
            x = n_mulmod2_preinv(y, h[m-i-1][m-1], mod.n, mod.ninv);
            x = n_negmod(x, mod.n);
            for (j = 0; j <= n; j++)
                NMOD_ADDMUL(t[m][j], t[m-i-1][j], x, mod);
        }
    }

    nmod_poly_fit_length(rop, n + 1);
    for (j = 0; j <= n; j++)
        rop->coeffs[j] = t[n][j];
    rop->length = n + 1;

    nmod_mat_clear(H);
    nmod_mat_clear(T);
}

void fmpz_mat_charpoly_modular(fmpz_poly_t rop, const fmpz_mat_t op)
{
    const long n = op->r;

    if (n < 4)  /* TODO: Write a small charpoly function */
    {
        fmpz_mat_charpoly(rop, op);
    }
    else
    {
        /*
            Lemma 4.1.

            If $A$ is an $n \times n$ matrix with $n \geq 4$ and 
            coefficients bounded in absolute value by $B > 1$ then 
            the coefficients of the characteristic polynomial have 
            less than $\ceil{n/2 (\log_2(n) + \log_2(B^2) + 1.6669)}$ 
            bits.
         */
        long bound;

        long pbits  = FLINT_BITS - 1;
        mp_limb_t p = (1UL << pbits);

        fmpz_t m;

        /* Determine the bound in bits */
        {
            long i, j;
            fmpz *ptr;
            double t;

            ptr = fmpz_mat_entry(op, 0, 0);
            for (i = 0; i < n; i++)
                for (j = 0; j < n; j++)
                    if (fmpz_cmpabs(ptr, fmpz_mat_entry(op, i, j)) < 0)
                        ptr = fmpz_mat_entry(op, i, j);

            if (fmpz_bits(ptr) == 0)  /* Zero matrix */
            {
                fmpz_poly_zero(rop);
                fmpz_poly_set_coeff_ui(rop, n, 1);
                return;
            }

            t = (fmpz_bits(ptr) <= FLINT_D_BITS) ? 
                _log2(FLINT_ABS(fmpz_get_d(ptr))) : fmpz_bits(ptr);

            bound = ceil( (n / 2.0) * (_log2(n) + 2.0 * t + 1.6669) );
        }

        fmpz_init_set_ui(m, 1);
        fmpz_poly_zero(rop);

        for ( ; fmpz_bits(m) < bound; )
        {
            nmod_mat_t mat;
            nmod_poly_t poly;

            p = n_nextprime(p, 0);

            nmod_mat_init(mat, n, n, p);
            nmod_poly_init(poly, p);

            fmpz_mat_get_nmod_mat(mat, op);
            nmod_mat_charpoly(poly, mat);

            fmpz_poly_CRT_ui(rop, rop, m, poly, 1);

            fmpz_mul_ui(m, m, p);

            nmod_mat_clear(mat);
            nmod_poly_clear(poly);
        }
        fmpz_clear(m);
    }
}

