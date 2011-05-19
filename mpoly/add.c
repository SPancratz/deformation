#include "mpoly.h"

void _mpoly_add_in_place(mpoly_t rop, const mpoly_t op, const ctx_t ctx)
{
    mpoly_iter_t iter;
    mpoly_term t;

    if (rop->n != op->n)
    {
        printf("ERROR (_mpoly_add_in_place).\n");
        abort();
    }

    if (mpoly_is_zero(op, ctx))
        return;

    if (rop == op)
    {
        mpoly_scalar_mul_si(rop, rop, 2, ctx);
        return;
    }

    mpoly_iter_init(iter, op);
    while ((t = mpoly_iter_next(iter)))
    {
        /*
           For every term t in op, we check whether its monomial is also 
           present in rop.  If it is, we add the coefficient in op to that 
           in rop.  Otherwise, we add a new term to rop.  In either case, 
           we are careful not to introduce zero terms.
         */

        int find;
        mon_t m;
        void *c;

        find = RBTREE_FIND(mpoly, &m, &c, rop->dict, t->key, &mon_cmp);

        if (find)
        {
            ctx->add(c, c, t->val);
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
            ctx->set(c2, t->val);

            RBTREE_INSERT(mpoly, &m, &c, rop->dict, t->key, c2, &mon_cmp);
        }
    }
    mpoly_iter_clear(iter);
}

void mpoly_add(mpoly_t rop, const mpoly_t op1, const mpoly_t op2, 
               const ctx_t ctx)
{
    if (op1->n != op2->n)
    {
        printf("ERROR (mpoly_add).  op1->n = %ld, op2->n = %ld.\n", 
               op1->n, op2->n);
        abort();
    }
    
    if (op1 == op2)
    {
        mpoly_scalar_mul_si(rop, op1, 2, ctx);
        return;
    }
    if (rop == op1)
    {
        _mpoly_add_in_place(rop, op2, ctx);
        return;
    }
    if (rop == op2)
    {
        _mpoly_add_in_place(rop, op1, ctx);
        return;
    }
    
    mpoly_set(rop, op1, ctx);
    _mpoly_add_in_place(rop, op2, ctx);
}

