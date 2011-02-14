#include "perm.h"

#include "flint.h"
#include "ulong_extras.h"

void _perm_randtest(long *vec, long n, flint_rand_t state)
{
    long i, j, t;

    for (i = 0; i < n; i++)
        vec[i] = i;

    for (i = n - 1; i > 0; i--)
    {
        j = n_randint(state, i + 1);
        t = vec[i];
        vec[i] = vec[j];
        vec[j] = t;
    }
}
