#include "mpoly.h"

void mpoly_iter_init(mpoly_iter_t iter, const mpoly_t op)
{
    RBTREE_ITER_INIT(mpoly, iter, op->dict);
}

void mpoly_iter_clear(mpoly_iter_t iter)
{
    RBTREE_ITER_CLEAR(mpoly, iter);
}

mpoly_term mpoly_iter_next(mpoly_iter_t iter)
{
    return RBTREE_ITER_NEXT(mpoly, iter);
}

