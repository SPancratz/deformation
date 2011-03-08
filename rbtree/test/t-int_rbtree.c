#include <stdlib.h>
#include <stdio.h>

#include "rbtree.h"

#include "flint.h"
#include "fmpz.h"
#include "long_extras.h"
#include "ulong_extras.h"

#define int_cmp(x, y)  ((x) - (y))
#define int_clear(x) 

RBTREE_PROTOTYPE_H(int, int, int, int_cmp, int_clear, int_clear, static)
RBTREE_PROTOTYPE_C(int, int, int, int_cmp, int_clear, int_clear, static)

RBTREE_PROTOTYPE_DEBUG_H(int, static)
RBTREE_PROTOTYPE_DEBUG_C(int, static)

int main(void)
{
    long i, N = 1000, result;

    flint_rand_t state;

    flint_randinit(state);

    printf("int_rbtree... ");
    fflush(stdout);

    /* Initialise a tree, add 5N entries, remove 6N entries */
    {
        int_rbtree_t T;

        int_rbtree_init(T);

        for (i = 0; i < 5 * N; i++)
        {
            int ins;
            int a, b, c2, c4, c5;
            int x = z_randtest(state);
            int y = z_randtest(state);

            ins = int_rbtree_insert(&a, &b, T, x, y);

            c2 = int_rbtree_verify2(RBTREE_ROOT(T));
            c4 = int_rbtree_verify4(RBTREE_ROOT(T));
            c5 = int_rbtree_verify5(RBTREE_ROOT(T));

            if (!c2 || !c4 || !c5)
            {
                printf("FAIL:\n");
                printf("c2 c4 c5 = %d %d %d\n", c2, c4, c5);
                abort();
            }
        }

        {
            int_rbtree_node_struct * n;
            int_rbtree_iter_t iter;

            int_rbtree_iter_init(iter, T);
            while ((n = int_rbtree_iter_next(iter))) ;
            int_rbtree_iter_clear(iter);
        }

        for (i = 0; i < 6 * N; i++)
        {
            int del;
            int a, b, c2, c4, c5;
            int x = z_randtest(state);

            del = int_rbtree_delete(&a, &b, T, x);

            c2 = int_rbtree_verify2(RBTREE_ROOT(T));
            c4 = int_rbtree_verify4(RBTREE_ROOT(T));
            c5 = int_rbtree_verify5(RBTREE_ROOT(T));

            if (!c2 || !c4 || !c5)
            {
                printf("FAIL:\n");
                printf("c2 c4 c5 = %d %d %d\n", c2, c4, c5);
                abort();
            }
        }

        int_rbtree_clear(T);
    }

    /*
        Create a tree, add elements from 0, ..., L, check that 
        we can remove them one at a time, finally ending up with 
        an empty tree
     */
    {
        int_rbtree_t T;
        int L = 1000;

        int_rbtree_init(T);

        for (i = 0; i < L; i++)
        {
            int ins;
            int a, b, c2, c4, c5;
            int x = i;
            int y = z_randtest(state);

            ins = int_rbtree_insert(&a, &b, T, x, y);

            c2 = int_rbtree_verify2(RBTREE_ROOT(T));
            c4 = int_rbtree_verify4(RBTREE_ROOT(T));
            c5 = int_rbtree_verify5(RBTREE_ROOT(T));

            if (!c2 || !c4 || !c5)
            {
                printf("FAIL:\n");
                printf("c2 c4 c5 = %d %d %d\n", c2, c4, c5);
                abort();
            }
        }

        for (i = L - 1; i >= 0; i--)
        {
            int del;
            int a, b, c2, c4, c5;
            int x = i;

            del = int_rbtree_delete(&a, &b, T, x);

            result = (del == 1);
            if (!result)
            {
                printf("FAIL:\n\n");
                printf("Could not remove %d from the tree.\n", x);
                abort();
            }

            c2 = int_rbtree_verify2(RBTREE_ROOT(T));
            c4 = int_rbtree_verify4(RBTREE_ROOT(T));
            c5 = int_rbtree_verify5(RBTREE_ROOT(T));

            if (!c2 || !c4 || !c5)
            {
                printf("FAIL:\n");
                printf("c2 c4 c5 = %d %d %d\n", c2, c4, c5);
                abort();
            }
        }

        result = int_rbtree_is_empty(T);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("Tree is not empty.\n");
            abort();
        }

        int_rbtree_clear(T);
    }


    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

