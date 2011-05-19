#include "mpoly.h"

int mpoly_is_one(const mpoly_t op, const ctx_t ctx)
{
    if (RBTREE_SIZE(mpoly, op->dict) != 1)
        return 0;
    else
    {
        mpoly_term t = RBTREE_ROOT(op->dict);

        return mon_is_one(t->key) && ctx->is_one(ctx, t->val);
    }
}
