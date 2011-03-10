#include <stdlib.h>
#include <stdio.h>

#include "rbtree.h"

#include "mon.h"

#include "flint.h"
#include "fmpz.h"
#include "long_extras.h"
#include "ulong_extras.h"

static int mon_compare(const mon_t x, const mon_t y)
{
    return mon_cmp_invlex(x, y);
}

RBTREE_PROTOTYPE_H(mpoly, mon_t, int, static)
RBTREE_PROTOTYPE_C(mpoly, mon_t, int, static)

RBTREE_PROTOTYPE_DEBUG_H(mpoly, static)
RBTREE_PROTOTYPE_DEBUG_C(mpoly, static)

int main(void)
{
    long i, n = 3, N = 1000, result;

    flint_rand_t state;

    flint_randinit(state);

    printf("mpoly_rbtree... ");
    fflush(stdout);

    /* Initialise a tree, add 5N entries, remove 6N entries */
    {
        RBTREE_T(mpoly) T;

        RBTREE_INIT(mpoly, T);

        for (i = 0; i < 5 * N; i++)
        {
            int ins, c2, c4, c5;
            mon_t a, x;
            int b, y;

            mon_init(x);
            mon_randtest(x, state, n, 10);
            y = z_randtest(state);

            ins = RBTREE_INSERT(mpoly, &a, &b, T, x, y, &mon_compare);

            c2 = RBTREE_VERIFY2(mpoly, RBTREE_ROOT(T));
            c4 = RBTREE_VERIFY4(mpoly, RBTREE_ROOT(T));
            c5 = RBTREE_VERIFY5(mpoly, RBTREE_ROOT(T));

            if (!c2 || !c4 || !c5)
            {
                printf("FAIL:\n");
                printf("c2 c4 c5 = %d %d %d\n", c2, c4, c5);
                abort();
            }
        }

        {
            RBTREE_NODE(mpoly) *n;
            RBTREE_ITER_T(mpoly) iter;

            RBTREE_ITER_INIT(mpoly, iter, T);

            while ((n = RBTREE_ITER_NEXT(mpoly, iter))) ;

            RBTREE_ITER_CLEAR(mpoly, iter);
        }

        for (i = 0; i < 6 * N; i++)
        {
            int del, c2, c4, c5;

            mon_t a, x;
            int b;

            mon_init(x);
            mon_randtest(x, state, n, 10);

            del = RBTREE_DELETE(mpoly, &a, &b, T, x, &mon_compare);

            c2 = RBTREE_VERIFY2(mpoly, RBTREE_ROOT(T));
            c4 = RBTREE_VERIFY4(mpoly, RBTREE_ROOT(T));
            c5 = RBTREE_VERIFY5(mpoly, RBTREE_ROOT(T));

            if (!c2 || !c4 || !c5)
            {
                printf("FAIL:\n");
                printf("c2 c4 c5 = %d %d %d\n", c2, c4, c5);
                abort();
            }
        }

        RBTREE_CLEAR(mpoly, T, NULL, NULL);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

