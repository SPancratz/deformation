#include "mat_csr.h"
#include "vec.h"

void 
_mat_csr_mul_vec(char *res, 
                 long m, long n, 
                 const char *x, const long *j, const long *p, const long *lenr, 
                 const char *vec, 
                 const mat_ctx_t ctx)
{
    char *t;
    long r, q;

    t = _vec_init(1, ctx);

    for (r = 0; r < m; r++)
    {
        ctx->zero(res + r * ctx->size);
        for (q = p[r]; q < p[r] + lenr[r]; q++)
        {
            ctx->mul(t, x + q * ctx->size, vec + j[q] * ctx->size);
            ctx->add(res + r * ctx->size, res + r * ctx->size, t);
        }
    }

    _vec_clear(t, 1, ctx);
}

void 
mat_csr_mul_vec(char *res, const mat_csr_t mat, const char *vec, 
                const mat_ctx_t ctx)
{
    char *in;

    in = (res == vec) ? _vec_init(mat->n, ctx) : (char *) vec;

    _mat_csr_mul_vec(res, mat->m, mat->n, mat->x, mat->j, mat->p, mat->lenr, 
                     vec, ctx);

    if (res == vec)
        _vec_clear(in, mat->n, ctx);

}

