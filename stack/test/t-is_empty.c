/******************************************************************************

    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include <stdio.h>

#include "stack.h"
#include "generics.h"
#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/ulong_extras.h"

STACK_PROTOTYPE(ulong, ulong, static)

int
main(void)
{
    int i, j, result;
    flint_rand_t state;

    printf("is_empty... ");
    fflush(stdout);

    _randinit(state);

    for (i = 0; i < 10000; i++)
    {
        ulong_stack_t Q;
        long len;

        len = n_randint(state, 100) + 1;

        ulong_stack_init2(Q, len);

        result = (ulong_stack_is_empty(Q));
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("len = %ld\n", len);
            abort();
        }

        ulong_stack_clear(Q);
    }

    for (i = 0; i < 100; i++)
    {
        ulong_stack_t Q;
        long len;

        len = n_randint(state, 100) + 1;
        ulong_stack_init(Q);

        for (j = 0; j < len; j++)
            ulong_stack_push(Q, n_randtest(state));
        for (j = 0; j < len; j++)
            ulong_stack_pop(Q);

        result = (ulong_stack_is_empty(Q));
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("len = %ld\n", len);
            abort();
        }

        ulong_stack_clear(Q);
    }

    for (i = 0; i < 100; i++)
    {
        ulong_stack_t Q;
        long len;

        len = n_randint(state, 100) + 1;
        ulong_stack_init(Q);

        for (j = 0; j < len; j++)
            ulong_stack_push(Q, n_randtest(state));

        result = (!ulong_stack_is_empty(Q));
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("len = %ld\n", len);
            abort();
        }

        ulong_stack_clear(Q);
    }

    _randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
