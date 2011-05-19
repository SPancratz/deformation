#include "mpoly.h"

void mpoly_get_coeff(void *rop, const mpoly_t op, const mon_t m, 
                     const ctx_t ctx)
{
    int find;
    mon_t m2;
    void *c2;

    find = RBTREE_FIND(mpoly, &m2, &c2, op->dict, m, &mon_cmp);

    if (find)
    {
        ctx->set(rop, c2);
    }
    else
    {
        ctx->zero(rop);
    }
}

