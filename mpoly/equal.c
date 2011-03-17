#include "mpoly.h"

int mpoly_equal(const mpoly_t op1, const mpoly_t op2, const mat_ctx_t ctx)
{
    long len1, len2;
    mpoly_iter_t iter;
    mpoly_term t;

    if (op1->n != op2->n)
    {
        printf("ERROR (mpoly_equal).  op1->n = %ld, op2->n = %ld.\n", 
               op1->n, op2->n);
        abort();
    }

    len1 = RBTREE_SIZE(mpoly, op1->dict);
    len2 = RBTREE_SIZE(mpoly, op2->dict);

    if (len1 != len2)
        return 0;

    if (len1 == 0)
        return 1;

    mpoly_iter_init(iter, op1);
    while ((t = mpoly_iter_next(iter)))
    {
        int find, b;
        mon_t m;
        void *c;

        find = RBTREE_FIND(mpoly, &m, &c, op2->dict, t->key, &mon_cmp);

        if (!find || !ctx->equal(t->val, c))
        {
            mpoly_iter_clear(iter);
            return 0;
        }
    }
    mpoly_iter_clear(iter);
    return 1;
}

