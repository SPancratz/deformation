#include "mpoly.h"


void mpoly_add_coeff(mpoly_t rop, const mon_t m, const void *x, 
                     const ctx_t ctx)
{
    int find;
    mon_t m2;
    void *x2;

    if (ctx->is_zero(x))
        return;

    find = RBTREE_FIND(mpoly, &m2, &x2, rop->dict, m, &mon_cmp);

    if (find)
    {
        ctx->add(x2, x2, x);
        if (ctx->is_zero(x2))
        {
            RBTREE_DELETE(mpoly, &m2, &x2, rop->dict, m, &mon_cmp);
            ctx->clear(x2);
            free(x2);
        }
    }
    else
    {
        char *x1;

        x1 = malloc(ctx->size);
        ctx->init(x1);
        ctx->set(x1, x);

        RBTREE_INSERT(mpoly, &m2, &x2, rop->dict, m, x1, &mon_cmp);
    }
}
