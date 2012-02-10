/******************************************************************************

    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include <stdio.h>

#include "queue.h"
#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

QUEUE_PROTOTYPE(ulong, ulong, static)

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
        ulong_queue_t Q;

        ulong_queue_init(Q);
        ulong_queue_clear(Q);
    }

    for (i = 0; i < 10000; i++)
    {
        ulong_queue_t Q;
        long len;

        len = n_randint(state, 100) + 1;

        ulong_queue_init2(Q, len);
        ulong_queue_clear(Q);
    }

    for (i = 0; i < 10000; i++)
    {
        ulong_queue_t Q;
        long len;

        len = n_randint(state, 100) + 1;

        ulong_queue_init(Q);
        ulong_queue_fit_size(Q, len);
        ulong_queue_clear(Q);
    }

    _randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
