#include "mpoly.h"

void _mpoly_sub_in_place(mpoly_t rop, const mpoly_t op, const mat_ctx_t ctx)
{
    mpoly_iter_t iter;
    mpoly_term t;

    if (rop->n != op->n)
    {
        printf("ERROR (_mpoly_sub_in_place).  rop->n = %ld, op->n = %ld.\n", 
                rop->n, op->n);
        abort();
    }

    if (mpoly_is_zero(op, ctx))
        return;

    if (rop == op)
    {
        mpoly_zero(rop, ctx);
        return;
    }

    mpoly_iter_init(iter, op);
    while ((t = mpoly_iter_next(iter)))
    {
        /* 
           For every term in op, we check whether its monomial is also 
           present in rop.  If it is, we subtract the coefficient in op 
           from that in rop.  Otherwise, we add a new term to rop.  In 
           either case, we are careful not to introduce a zero term.
         */

        int find;
        mon_t m;
        void *c;

        find = RBTREE_FIND(mpoly, &m, &c, rop->dict, t->key, &mon_cmp);

        if (find)
        {
            ctx->sub(c, c, t->val);
            if (ctx->is_zero(c))
            {
                RBTREE_DELETE(mpoly, &m, &c, rop->dict, t->key, &mon_cmp);
                ctx->clear(c);
                free(c);
            }
        }
        else
        {
            char *c2;

            c2 = malloc(ctx->size);
            ctx->init(c2);
            ctx->neg(c2, t->val);

            RBTREE_INSERT(mpoly, &m, &c, rop->dict, t->key, c2, &mon_cmp);
        }
    }
    mpoly_iter_clear(iter);
}

void mpoly_sub(mpoly_t rop, const mpoly_t op1, const mpoly_t op2, 
               const mat_ctx_t ctx)
{
    if (op1->n != op2->n)
    {
        printf("ERROR (mpoly_sub).  op1->n = %ld, op2->n = %ld.\n", 
               op1->n, op2->n);
        abort();
    }

    if (op1 == op2)
    {
        mpoly_zero(rop, ctx);
        return;
    }
    if (rop == op1)
    {
        _mpoly_sub_in_place(rop, op2, ctx);
        return;
    }
    if (rop == op2)
    {
        _mpoly_sub_in_place(rop, op1, ctx);
        mpoly_neg(rop, rop, ctx);
        return;
    }
    
    mpoly_set(rop, op1, ctx);
    _mpoly_sub_in_place(rop, op2, ctx);
}

