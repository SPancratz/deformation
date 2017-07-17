/******************************************************************************

    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include <stdio.h>

#include "queue.h"
#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/ulong_extras.h"

QUEUE_PROTOTYPE(ulong, ulong, static)

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
        ulong_queue_t Q;
        long len;

        len = n_randint(state, 100) + 1;

        ulong_queue_init2(Q, len);

        result = (ulong_queue_is_empty(Q));
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("len = %ld\n", len);
            abort();
        }

        ulong_queue_clear(Q);
    }

    for (i = 0; i < 100; i++)
    {
        ulong_queue_t Q;
        long len;

        len = n_randint(state, 100) + 1;
        ulong_queue_init(Q);

        for (j = 0; j < len; j++)
            ulong_queue_enqueue(Q, n_randtest(state));
        for (j = 0; j < len; j++)
            ulong_queue_dequeue(Q);

        result = (ulong_queue_is_empty(Q));
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("len = %ld\n", len);
            abort();
        }

        ulong_queue_clear(Q);
    }

    _randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
