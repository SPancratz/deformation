#include <assert.h>
#include <stdlib.h>

#include "mat.h"
#include "vec.h"

#include "mat_csr.h"

#define DEBUG  1

void mat_csr_solve_init(mat_csr_solve_t s, const mat_csr_t mat, 
                        const mat_ctx_t ctx)
{
    long *w;
    long *mem;
    long i, k, m, nz;

    assert(mat->m == mat->n);

    m = mat->m;

    #if (DEBUG > 0)
    printf("mat_csr_solve_init()\n");
    printf("Matrix A (%ld x %ld):\n", mat->m, mat->n);
    mat_csr_print_dense(mat, ctx);
    printf("\n");
    fflush(stdout);
    #endif

    mem   = malloc((mat->alloc + 6 * m + 1) * sizeof(long));
    s->LU = malloc(m * sizeof(char *));
    w     = malloc((3 * m) * sizeof(long));

    if (!mem || !(s->LU) || !w)
    {
        printf("ERROR (mat_csr_solve_init).  Memory allocation.\n\n");
        abort();
    }

    s->j    = mem;
    s->p    = mem + mat->alloc;
    s->lenr = mem + mat->alloc + m;
    s->pi   = mem + mat->alloc + 2 * m;
    s->qi   = mem + mat->alloc + 3 * m;
    s->P    = mem + mat->alloc + 4 * m;
    s->B    = mem + mat->alloc + 5 * m;

    #if (DEBUG > 0)
    printf("Calling mat_csr_zfdiagonal()\n");
    fflush(stdout);
    #endif

    nz = mat_csr_zfdiagonal(s->pi, mat);

    #if (DEBUG > 0)
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

    s->m = mat->m;
    s->n = mat->n;
    s->x = mat->x;

    for (i = 0; i < m; i++)
        s->lenr[i] = mat->lenr[i];
    for (i = 0; i < m; i++)
        s->p[i] = mat->p[i];
    for (i = 0; i < mat->alloc; i++)
        s->j[i] = mat->j[i];

    /* Set A := P A */

    #if (DEBUG > 0)
    printf("Calling _mat_csr_permute_rows()\n");
    fflush(stdout);
    #endif

    _mat_csr_permute_rows(m, s->p, s->lenr, s->pi);

    /* Find Q s.t. Q P A Q^t is block triangular */

    #if (DEBUG > 0)
    printf("Calling _mat_csr_block_triangularise()\n");
    fflush(stdout);
    #endif

    s->nb = _mat_csr_block_triangularise(s->qi, s->B, m, s->j, s->p, s->lenr, w);
    s->B[s->nb] = m;

    #if (DEBUG > 0)
    printf("Calling _mat_csr_permute_rows()\n");
    fflush(stdout);
    #endif

    _mat_csr_permute_rows(m, s->p, s->lenr, s->qi);

    #if (DEBUG > 0)
    printf("Calling _mat_csr_permute_cols()\n");
    fflush(stdout);
    #endif

    _mat_csr_permute_cols(m, m, s->j, s->p, s->lenr, s->qi);

    #if (DEBUG > 0)
    printf("nb = %ld\n", s->nb);
    printf("B  = {"); _perm_print(s->B, s->nb); printf("}\n");
    printf("qi = {"); _perm_print(s->qi, m); printf("}\n");
    printf("Matrix Q P A Q^t:\n");
    _mat_csr_print_dense(m, m, s->x, s->j, s->p, s->lenr, ctx);
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
                for (q = s->p[i]; q < s->p[i] + s->lenr[i]; q++)
                {
                    long c = s->j[q] - s->B[k];

                    if (0 <= c && c < len)
                        ctx->set(rows[r] + c * ctx->size, s->x + q * ctx->size);
                }

            off += len * len * ctx->size;
        }
    }

    /* Dense LUP decomposition */

    #if (DEBUG > 0)
    printf("LUP decomposition for dense kernels..\n");
    fflush(stdout);
    #endif

    for (k = 0; k < s->nb; k++)
    {
        char **rows = s->LU + s->B[k];
        long len = s->B[k + 1] - s->B[k];

        #if (DEBUG > 0)
        printf("  Block %ld out of %ld, of size %ld\n", k, s->nb, len);
        fflush(stdout);
        #endif

        _mat_lup_decompose(s->P + s->B[k], rows, len, ctx);
    }

    #if (DEBUG > 0)
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

