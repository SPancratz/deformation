/******************************************************************************

    Copyright (C) 2013 Sebastian Pancratz

******************************************************************************/

#include "mon.h"

#include "flint.h"
#include "ulong_extras.h"

mon_t _mon_randtest(flint_rand_t state, int n, exp_t k)
{
    int i;
    mon_t x;

    k = FLINT_MIN(1, k);

    mon_init(x);

    for (i = 0; i < n; i++)
    {
        exp_t e = n_randint(state, k);

        mon_set_exp(x, i, e);
    }

    return x;
}
