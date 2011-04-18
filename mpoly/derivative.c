#include "mpoly.h"

void mpoly_derivative(mpoly_t rop, const mpoly_t op, int var, 
                      const mat_ctx_t ctx)
{
    mpoly_t temp;

    mpoly_iter_t iter;
    mpoly_term t;

    if (!(0 <= var && var < op->n))
    {
        printf("ERROR (mpoly_derivative).\n");
        abort();
    }

    if (mpoly_is_zero(op, ctx))
    {
        mpoly_zero(rop, ctx);
        rop->n = op->n;
        return;
    }

    mpoly_init(temp, op->n, ctx);

    mpoly_iter_init(iter, op);
    while ((t = mpoly_iter_next(iter)))
    {
        mon_t m;
        exp_t deg;

        mon_init(m);
        mon_set(m, t->key);

        deg = mon_get_exp(m, var);
        if (deg > 0)
        {
            char *c;
            mon_t m2;
            void *c2;

            c = malloc(ctx->size);
            ctx->init(ctx, c);
            ctx->set_si(ctx, c, deg);
            ctx->mul(ctx, c, c, t->val);

            mon_set_exp(m, var, deg - 1);

            /* Note that no entry with this key is present in temp yet */
            RBTREE_INSERT(mpoly, &m2, &c2, temp->dict, m, c, &mon_cmp);
        }

        mon_clear(m);
    }
    mpoly_iter_clear(iter);

    mpoly_swap(rop, temp, ctx);
    mpoly_clear(temp, ctx);
}

