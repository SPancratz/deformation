#include "mat_csr.h"

#include "flint.h"
#include "ulong_extras.h"

void 
mat_csr_randtest(mat_csr_t A, 
                 flint_rand_t state, double d, const ctx_t ctx)
{
    char *mem;
    long f, i, j, k, len, u;

    d = FLINT_MAX(d, 0.0);
    d = FLINT_MIN(d, 1.0);
    f = 100 * d;

    mat_csr_zero(A, ctx);

    if (f == 0)
        return;

    u   = 2 * sizeof(long) + ctx->size;
    len = d * A->m * A->n;
    mem = malloc(u * len);

    if (!mem)
    {
        printf("ERROR (mat_csr_randtest).\n\n");
        abort();
    }

    k = 0;
    for (i = 0; (k < len) && (i < A->m); i++)
        for (j = 0; (k < len) && (j < A->n); j++)
        {
            if (n_randint(state, 100) <= f)
            {
                *(long *) (mem + k * u) = i;
                *(long *) (mem + k * u + sizeof(long)) = j;
                ctx->init(mem + k * u + 2 * sizeof(long));
                ctx->randtest_not_zero(mem + k * u + 2 * sizeof(long), state);
                k++;
            }
        }

    mat_csr_set_array3(A, mem, k, 0, ctx);

    free(mem);
}
