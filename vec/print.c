#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "vec.h"

int _vec_print(const char *vec, long n, const ctx_t ctx)
{
    assert(n >= 0);

    printf("%ld", n);

    if (n)
    {
        long i;

        printf(" ");
        for (i = 0; i < n; i++)
        {
            printf(" ");
            ctx->print(ctx, vec + i * ctx->size);
        }
    }

    return 1;
}

