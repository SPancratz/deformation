/******************************************************************************

    Copyright (C) 2013 Sebastian Pancratz
 
******************************************************************************/

#include "gmconnection.h"

void gmc_poly2array(char *c, const mpoly_t poly, const mon_t *m, long len, 
                    const ctx_t ctx)
{
    long i;

    _vec_zero(c, len, ctx);

    for (i = 0; i < len; i++)
        mpoly_get_coeff(c + i * ctx->size, poly, m[i], ctx);
}

