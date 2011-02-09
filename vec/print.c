#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "vec.h"

int _vec_print(const char *vec, long n, const mat_ctx_t ctx)
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
            ctx->print(vec + i * ctx->size);
        }
    }

    return 1;
}

