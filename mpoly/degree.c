#include "mpoly.h"

long mpoly_degree(const mpoly_t op, long var, const mat_ctx_t ctx)
{
    long deg;
    mpoly_iter_t iter;
    mpoly_term t;

    if (!(var == -1 || (1 <= var && var <= 8)))
    {
        printf("ERROR (mpoly_degree).  var = %ld.\n", var);
        abort();
    }

    if (mpoly_is_zero(op, ctx))
        return -1;
    
    deg = -1;
    
    mpoly_iter_init(iter, op);
    while ((t = mpoly_iter_next(iter)))
    {
        long temp;

        temp = (var < 0) ? mon_degree(t->key) : mon_get_exp(t->key, var);
        if (deg < temp)
            deg = temp;
    }
    mpoly_iter_clear(iter);
    
    return deg;
}

