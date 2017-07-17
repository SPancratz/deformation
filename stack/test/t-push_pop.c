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
    int i, j, k, result;
    flint_rand_t state;

    printf("push/ pop... ");
    fflush(stdout);

    _randinit(state);

    for (i = 0; i < 100; i++)
    {
        ulong_stack_t Q;

        ulong_stack_init(Q);

        for (j = 0; j < 100; j++)
            ulong_stack_push(Q, n_randtest(state));

        j = 0;
        while (!ulong_stack_is_empty(Q))
        {
            ulong_stack_pop(Q);
            j++;
        }

        result = (j == 100);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("j = %d\n", j);
            abort();
        }

        ulong_stack_clear(Q);
    }

    for (i = 0; i < 100; i++)
    {
        ulong_stack_t Q;

        ulong_stack_init(Q);

        for (j = 0; j < 100; j++)
            ulong_stack_push(Q, j);

        j = 100;
        while (!ulong_stack_is_empty(Q))
        {
            j--;

            k = ulong_stack_pop(Q);

            result = (k == j);
            if (!result)
            {
                printf("FAIL:\n\n");
                printf("Stack violates FILO behaviour.\n");
                abort();
            }
        }

        ulong_stack_clear(Q);
    }

    _randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
