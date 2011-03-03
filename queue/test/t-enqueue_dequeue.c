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
    int i, j, k, result;
    flint_rand_t state;

    printf("enqueue/ dequeue... ");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 100; i++)
    {
        ulong_queue_t Q;

        ulong_queue_init(Q);

        for (j = 0; j < 100; j++)
            ulong_queue_enqueue(Q, n_randtest(state));

        j = 0;
        while (!ulong_queue_is_empty(Q))
        {
            ulong_queue_dequeue(Q);
            j++;
        }

        result = (j == 100);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("j = %d\n", j);
            abort();
        }

        ulong_queue_clear(Q);
    }

    for (i = 0; i < 100; i++)
    {
        ulong_queue_t Q;

        ulong_queue_init(Q);

        for (j = 0; j < 100; j++)
            ulong_queue_enqueue(Q, j);

        j = 0;
        while (!ulong_queue_is_empty(Q))
        {
            k = ulong_queue_dequeue(Q);

            result = (k == j);
            if (!result)
            {
                printf("FAIL:\n\n");
                printf("Queue violates FIFO behaviour.\n");
                abort();
            }

            j++;
        }

        ulong_queue_clear(Q);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
