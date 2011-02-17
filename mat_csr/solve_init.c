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

#if !defined(NDEBUG)
    printf("mat_csr_solve_init(...)\n");
    printf("Matrix A:\n");
    mat_csr_print_dense(mat, ctx);
    printf("\n");
    fflush(stdout);
#endif

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

#if !defined(NDEBUG)
    printf("nz = %ld\n", nz);
    printf("pi = {"); _perm_print(s->pi, m); printf("}\n");
    fflush(stdout);
#endif

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

#if !defined(NDEBUG)
    printf("nb = %ld\n", s->nb);
    printf("B  = {"); _perm_print(s->B, s->nb); printf("}\n");
    printf("qi = {"); _perm_print(s->qi, m); printf("}\n");
    printf("Matrix Q P A Q^t:\n");
    _mat_csr_print_dense(m, m, s->mat->x, j, p, lenr, ctx);
    printf("\n");
    fflush(stdout);
#endif

    /* Allocate dense data blocks */
    {
        long len, sum = 0;

        for (k = 0; k < s->nb; k++)
        {
            len  = s->B[k + 1] - s->B[k];
            sum += len * len;
        }

        s->entries = _vec_init(sum, ctx);
        s->alloc   = sum;
    }

    /* Copy data of the square blocks */
    {
        char *off = s->entries;

        for (k = 0; k < s->nb; k++)
        {
            char **rows = s->LU + s->B[k];
            long len = s->B[k + 1] - s->B[k];
            long q, r;

            for (r = 0; r < len; r++)
                rows[r] = off + r * len * ctx->size;

            for (r = 0, i = s->B[k]; r < len; r++, i++)
                for (q = p[i]; q < p[i] + lenr[i]; q++)
                {
                    long c = j[q] - s->B[k];

                    if (0 <= c && c < len)
                        ctx->set(rows[r] + c * ctx->size, mat->x + q * ctx->size);
                }

            off += len * len * ctx->size;
        }
    }

    /* Dense LUP decomposition */

    for (k = 0; k < s->nb; k++)
    {
        char **rows = s->LU + s->B[k];
        long len = s->B[k + 1] - s->B[k];

        _mat_lup_decompose(s->P + s->B[k], rows, len, ctx);
    }

#if !defined(NDEBUG)
    for (k = 0; k < s->nb; k++)
    {
        long len = s->B[k + 1] - s->B[k];

        printf("Block %ld (of length %ld):\n", k, len);
        _mat_print(s->LU + s->B[k], len, len, ctx);
        printf("\n");
        fflush(stdout);
    }
#endif

    /* Compose Q P, invert Q */

    _perm_compose(w, s->pi, s->qi, m);
    _perm_set(s->pi, w, m);

    _perm_inv(s->qi, s->qi, m);

    /* Clean-up temporary space */

    free(w);
}

