/******************************************************************************

    Copyright (C) 2013 Sebastian Pancratz

******************************************************************************/

#include "diagfrob.h"

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

