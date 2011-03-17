#include "mpoly.h"

void mpoly_submul(mpoly_t rop, const mpoly_t op1, const mpoly_t op2, 
                  const mat_ctx_t ctx)
{
    mpoly_t temp;

    if (rop->n != op1->n || rop->n != op2->n)
    {
        printf("ERROR (mpoly_submul).\n");
        abort();
    }

    mpoly_init(temp, rop->n, ctx);

    mpoly_mul(temp, op1, op2, ctx);
    mpoly_sub(rop, rop, temp, ctx);

    mpoly_clear(temp, ctx);
}

