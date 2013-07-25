/******************************************************************************

    Copyright (C) 2013 Sebastian Pancratz

******************************************************************************/

#include "diagfrob.h"

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

