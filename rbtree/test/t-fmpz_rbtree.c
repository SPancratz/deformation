#include <stdlib.h>
#include <stdio.h>

#include "rbtree.h"

#include "flint.h"
#include "fmpz.h"
#include "long_extras.h"
#include "ulong_extras.h"

static int fmpz_compare(const fmpz x, const fmpz y)
{
    return fmpz_cmp(&x, &y);
}

static void fmpz_delete(fmpz x)
{
    fmpz_clear(&x);
}

RBTREE_PROTOTYPE_H(fmpz, fmpz, fmpz, fmpz_compare, static)
RBTREE_PROTOTYPE_C(fmpz, fmpz, fmpz, fmpz_compare, static)

RBTREE_PROTOTYPE_DEBUG_H(fmpz, static)
RBTREE_PROTOTYPE_DEBUG_C(fmpz, static)

int main(void)
{
    long i, N = 1000, result;

    flint_rand_t state;

    flint_randinit(state);

    printf("fmpz_rbtree... ");
    fflush(stdout);

    /* Initialise a tree, add 5N entries, remove 6N entries */
    {
        fmpz_rbtree_t T;

        fmpz_rbtree_init(T);

        for (i = 0; i < 5 * N; i++)
        {
            int ins, c2, c4, c5;
            fmpz a, b, x, y;

            fmpz_init(&x);
            fmpz_init(&y);
            fmpz_randtest(&x, state, 100);
            fmpz_randtest(&y, state, 100);

            ins = fmpz_rbtree_insert(&a, &b, T, x, y);

            if (ins)
            {
                fmpz_clear(&a);
                fmpz_clear(&b);
            }

            c2 = fmpz_rbtree_verify2(RBTREE_ROOT(T));
            c4 = fmpz_rbtree_verify4(RBTREE_ROOT(T));
            c5 = fmpz_rbtree_verify5(RBTREE_ROOT(T));

            if (!c2 || !c4 || !c5)
            {
                printf("FAIL:\n");
                printf("c2 c4 c5 = %d %d %d\n", c2, c4, c5);
                abort();
            }
        }

        {
            fmpz_rbtree_node * n;
            fmpz_rbtree_iter_t iter;

            fmpz_rbtree_iter_init(iter, T);
            while ((n = fmpz_rbtree_iter_next(iter))) ;
            fmpz_rbtree_iter_clear(iter);
        }

        for (i = 0; i < 6 * N; i++)
        {
            int del, c2, c4, c5;
            fmpz a, b, x;

            fmpz_init(&x);
            fmpz_randtest(&x, state, 100);

            del = fmpz_rbtree_delete(&a, &b, T, x);

            if (del)
            {
                fmpz_clear(&a);
                fmpz_clear(&b);
            }

            fmpz_clear(&x);

            c2 = fmpz_rbtree_verify2(RBTREE_ROOT(T));
            c4 = fmpz_rbtree_verify4(RBTREE_ROOT(T));
            c5 = fmpz_rbtree_verify5(RBTREE_ROOT(T));

            if (!c2 || !c4 || !c5)
            {
                printf("FAIL:\n");
                printf("c2 c4 c5 = %d %d %d\n", c2, c4, c5);
                abort();
            }
        }

        fmpz_rbtree_clear(T, &fmpz_delete, &fmpz_delete);
    }

    /*
        Create a tree, add elements from 0, ..., L, check that 
        we can remove them one at a time, finally ending up with 
        an empty tree
     */
    {
        fmpz_rbtree_t T;
        int L = 1000;

        fmpz_rbtree_init(T);

        for (i = 0; i < L; i++)
        {
            int ins, c2, c4, c5;
            fmpz a, b, x, y;

            fmpz_init(&x);
            fmpz_init(&y);
            fmpz_set_si(&x, i);
            fmpz_randtest(&y, state, 100);

            ins = fmpz_rbtree_insert(&a, &b, T, x, y);

            if (ins)
            {
                fmpz_clear(&a);
                fmpz_clear(&b);
            }

            c2 = fmpz_rbtree_verify2(RBTREE_ROOT(T));
            c4 = fmpz_rbtree_verify4(RBTREE_ROOT(T));
            c5 = fmpz_rbtree_verify5(RBTREE_ROOT(T));

            if (!c2 || !c4 || !c5)
            {
                printf("FAIL:\n");
                printf("c2 c4 c5 = %d %d %d\n", c2, c4, c5);
                abort();
            }
        }

        for (i = L - 1; i >= 0; i--)
        {
            int del, c2, c4, c5;
            fmpz a, b, x;

            fmpz_init(&x);
            fmpz_set_si(&x, i);

            del = fmpz_rbtree_delete(&a, &b, T, x);

            if (del)
            {
                fmpz_clear(&a);
                fmpz_clear(&b);
            }

            result = (del == 1);
            if (!result)
            {
                printf("FAIL:\n\n");
                printf("Could not remove ");
                fmpz_print(&x);
                printf(" from the tree.\n");
                abort();
            }

            c2 = fmpz_rbtree_verify2(RBTREE_ROOT(T));
            c4 = fmpz_rbtree_verify4(RBTREE_ROOT(T));
            c5 = fmpz_rbtree_verify5(RBTREE_ROOT(T));

            if (!c2 || !c4 || !c5)
            {
                printf("FAIL:\n");
                printf("c2 c4 c5 = %d %d %d\n", c2, c4, c5);
                abort();
            }
        }

        result = fmpz_rbtree_is_empty(T);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("Tree is not empty.\n");
            abort();
        }

        fmpz_rbtree_clear(T, &fmpz_delete, &fmpz_delete);
    }


    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

