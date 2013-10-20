/******************************************************************************

    Copyright (C) 2013 Sebastian Pancratz
 
******************************************************************************/

#include "gmconnection.h"

int gmc_basis_contains(const mpoly_t f, long d)
{
    mpoly_iter_t iter;
    mpoly_term mt;
    long var;
    
    mpoly_iter_init(iter, f);
    while ((mt = mpoly_iter_next(iter)))
        for (var = 0; var < f->n; var++)
            if (mon_get_exp(mt->key, var) >= d - 1)
            {
                mpoly_iter_clear(iter);
                return 0;
            }
    mpoly_iter_clear(iter);
    return 1;
}

