#include "mpoly.h"

void mpoly_init(mpoly_t rop, long n, const mat_ctx_t ctx)
{
    if (!(1 <= n && n <= 8))
    {
        printf("ERROR (mpoly_init).  n == %ld.\n", n);
        abort();
    }

    RBTREE_INIT(mpoly, rop->dict);
    rop->n = n;
}
