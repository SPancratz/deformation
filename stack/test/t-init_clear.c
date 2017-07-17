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
    int i;
    flint_rand_t state;

    printf("init/ init2/ clear... ");
    fflush(stdout);

    _randinit(state);

    for (i = 0; i < 10000; i++)
    {
        ulong_stack_t Q;

        ulong_stack_init(Q);
        ulong_stack_clear(Q);
    }

    for (i = 0; i < 10000; i++)
    {
        ulong_stack_t Q;
        long len;

        len = n_randint(state, 100) + 1;

        ulong_stack_init2(Q, len);
        ulong_stack_clear(Q);
    }

    for (i = 0; i < 10000; i++)
    {
        ulong_stack_t Q;
        long len;

        len = n_randint(state, 100) + 1;

        ulong_stack_init(Q);
        ulong_stack_fit_size(Q, len);
        ulong_stack_clear(Q);
    }

    _randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
