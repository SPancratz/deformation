/******************************************************************************

    Copyright (C) 2013 Sebastian Pancratz
 
******************************************************************************/

#include "gmconnection.h"

void gmc_derivatives(mpoly_t *D, const mpoly_t P, const ctx_t ctx)
{
    long i;

    for (i = 0; i < P->n; i++)
        mpoly_derivative(D[i], P, i, ctx);
}

