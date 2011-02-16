#include <assert.h>
#include <stdlib.h>

#include "mat.h"
#include "vec.h"

#include "mat_csr.h"

void mat_csr_solve_init(mat_csr_solve_t s, const mat_csr_t mat, 
                        const mat_ctx_t ctx)
{
    long *w, *t;
    long *j, *p, *lenr;
    long i, k, m, nz;

    assert(mat->m == mat->n);

    m = mat->m;
    s->mat = mat;

    s->pi = malloc((4 * m + 1) * sizeof(long));
    s->LU = malloc(m * sizeof(char *));
    w     = malloc((2 * m + mat->alloc + 3 * m) * sizeof(long));

    if (!(s->pi) || !(s->LU) || !w)
    {
        printf("ERROR (mat_csr_solve_init).  Memory allocation.\n\n");
        abort();
    }

    s->qi = s->pi + 1 * m;
    s->P  = s->pi + 2 * m;
    s->B  = s->pi + 3 * m;

    nz = mat_csr_zfdiagonal(s->pi, mat);

    if (nz != m)
    {
        printf("ERROR (mat_csr_solve_init).  Singular matrix.\n\n");
        abort();
    }

    /* Copy the sparse structure of mat into {j, p, lenr} */

    lenr = w;
    p    = w + 1 * m;
    j    = w + 2 * m;
    t    = w + 2 * m + mat->alloc;

    for (i = 0; i < m; i++)
        lenr[i] = mat->lenr[i];
    for (i = 0; i < m; i++)
        p[i] = mat->p[i];
    for (i = 0; i < mat->alloc; i++)
        j[i] = mat->j[i];

    /* Set A := P A */

    _mat_csr_permute_rows(m, p, lenr, s->pi);

    /* Find Q s.t. Q P A Q^t is block triangular */

    s->nb = _mat_csr_block_triangularise(s->qi, s->B, m, j, p, lenr, t);
    s->B[s->nb] = m;

    _mat_csr_permute_rows(m, p, lenr, s->qi);
    _mat_csr_permute_cols(m, m, j, p, lenr, s->qi);

    {
        char *off, **rows;
        long len, sum;

        /* Allocate dense data blocks */

        sum = 0;
        for (k = 0; k < s->nb; k++)
        {
            len  = s->B[k + 1] - s->B[k];
            sum += len * len;
        }

        s->entries = _vec_init(sum, ctx);
        s->alloc   = sum;

        off = s->entries;
        for (k = 0; k < s->nb; k++)
        {
            long q, r, c;

            rows = s->LU + s->B[k];
            len  = s->B[k + 1] - s->B[k];

            for (r = 0; r < len; r++)
                rows[r] = off + r * len * ctx->size;

            for (r = 0; r < len; r++)
            {
                i = s->B[k] + r;

                for (q = p[i]; q < p[i] + lenr[i]; q++)
                {
                    c = j[q] - s->B[k];
                    ctx->set(rows[r] + c * ctx->size, mat->x + q * ctx->size);
                }
            }

            off += len * len * ctx->size;
        }

        /* Dense LUP decomposition */

        for (k = 0; k < s->nb; k++)
        {
            rows = s->LU + s->B[k];
            len  = s->B[k + 1] - s->B[k];

            _mat_lup_decompose(s->P + s->B[k], rows, len, ctx);
        }
    }

    /* Compose Q P, invert Q */

    _perm_compose(w, s->pi, s->qi, m);
    _perm_set(s->pi, w, m);

    _perm_inv(s->qi, s->qi, m);

    /* Clean-up temporary space */

    free(w);
}

