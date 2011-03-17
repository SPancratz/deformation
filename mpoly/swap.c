#include "mpoly.h"

void mpoly_swap(mpoly_t op1, mpoly_t op2, const mat_ctx_t ctx)
{
    long tn;

    RBTREE_SWAP(mpoly, op1->dict, op2->dict);

    tn     = op1->n;
    op1->n = op2->n;
    op2->n = tn;
}

