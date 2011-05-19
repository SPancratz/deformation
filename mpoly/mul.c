#include "mpoly.h"

void mpoly_mul(mpoly_t rop, const mpoly_t op1, const mpoly_t op2, 
               const ctx_t ctx)
{
    mpoly_t temp;
    
    mpoly_iter_t iter1, iter2;
    mpoly_term t1, t2;

    if (op1->n != op2->n)
    {
        printf("ERROR (mpoly_mul).  op1->n != op2->n.\n");
        abort();
    }

    if (mpoly_is_zero(op1, ctx) || mpoly_is_zero(op2, ctx))
    {
        mpoly_zero(rop, ctx);
        rop->n = op1->n;
        return;
    }

    mpoly_init(temp, op1->n, ctx);

    mpoly_iter_init(iter1, op1);
    while ((t1 = mpoly_iter_next(iter1)))
    {
        mpoly_iter_init(iter2, op2);
        while ((t2 = mpoly_iter_next(iter2)))
        {
            mon_t m;
            char *c;

            mon_init(m);
            c = malloc(ctx->size);
            ctx->init(c);

            mon_mul(m, t1->key, t2->key);
            ctx->mul(c, t1->val, t2->val);
            mpoly_add_coeff(temp, m, c, ctx);

            mon_clear(m);
            ctx->clear(c);
            free(c);
        }
        mpoly_iter_clear(iter2);
    }
    mpoly_iter_clear(iter1);
    
    mpoly_swap(rop, temp, ctx);
    mpoly_clear(temp, ctx);
}

