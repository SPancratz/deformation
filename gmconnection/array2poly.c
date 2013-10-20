/******************************************************************************

    Copyright (C) 2013 Sebastian Pancratz
 
******************************************************************************/

#include "gmconnection.h"

void gmc_array2poly(mpoly_t poly, const char *c, const mon_t *m, long len,
                    const ctx_t ctx)
{
    long i;

    mpoly_zero(poly, ctx);
    for (i = 0; i < len; i++)
        mpoly_add_coeff(poly, m[i], c + i * ctx->size, ctx);
}

