#include "gmconnection.h"

void gmc_decompose_poly(mpoly_t * A, const mpoly_t poly, 
                        const mat_csr_solve_t s, 
                        mon_t * const rows, 
                        mon_t * const cols, 
                        long * const p, 
                        const ctx_t ctx)
{
    char *x, *b;
    long j, var;

    x = _vec_init(s->m, ctx);
    b = _vec_init(s->n, ctx);

    gmc_poly2array(b, poly, rows, s->n, ctx);
    mat_csr_solve(x, s, b, ctx);

    for (var = 0; var < poly->n; var++)
    {
        mpoly_zero(A[var], ctx);
        for (j = p[var]; j < p[var + 1]; j++)
            mpoly_add_coeff(A[var], cols[j], x + j * ctx->size, ctx);
    }

    _vec_clear(x, s->m, ctx);
    _vec_clear(b, s->n, ctx);
}

